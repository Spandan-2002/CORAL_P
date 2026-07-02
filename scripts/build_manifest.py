#!/usr/bin/env python3
"""Build the CoralShot pipeline manifest from the bcl2fastq run report + the FASTQ dir.

Emits (all under $BASE):
  manifest.tsv     one row per paired FASTQ library
  filelist.txt     two absolute fastq paths per library (drives the raw-QC array)
  asm_units.tsv    co-assembly units: unit_id kind member_sids
  barcode_audit.txt  i7/i5 cross-pool collision report (the index-hop risk)

hop_pair_flag=yes is set for any library whose i7 OR i5 sequence is REUSED across the
i7 and i5 pools (single-index-hop reachable on a patterned lane) — computed, not hard-coded.
"""
import csv, re, sys, glob, os

BASE = os.environ.get('BASE', '/scratch/sr7729/CORAL_P')
HTML = f'{BASE}/20260612_LH00284_0351_B22KCWMLT1.html'
SHOTGUN = f'{BASE}/Shotgun'
PART = {'C': 'core', 'S': 'shell', 'F': 'feces'}

rows = re.findall(r'<tr>(.*?)</tr>', open(HTML).read(), re.S)
recs = []
for r in rows:
    tds = re.findall(r'<td[^>]*>(.*?)</td>', r, re.S)
    if len(tds) != 12:
        continue
    lane, project, sample, barcode = (tds[0].strip(), tds[1].strip(), tds[2].strip(), tds[3].strip())
    if project != 'Ramadi' or sample == 'Undetermined':
        continue
    m = re.match(r'^([PC])(\d+)([MB])_([CSF])$', sample)
    if not m:
        sys.stderr.write(f'WARN unparsed: {sample}\n'); continue
    pc, num, arm, part = m.groups()
    i7, i5 = (barcode.split('+') + [''])[:2]
    recs.append(dict(sample_id=sample, subject=f'{pc}{num}', arm=arm, part=part,
                     type=PART[part], i7=i7, i5=i5))

by_sid = {r['sample_id']: r for r in recs}

observed_fastq_sids = sorted(
    re.sub(r'_S\d+_L002_R1_001\.fastq\.gz$', '', os.path.basename(p))
    for p in glob.glob(f'{SHOTGUN}/*_L002_R1_001.fastq.gz')
)
assert observed_fastq_sids, f'no R1 FASTQs found under {SHOTGUN}'

existing = {}
if os.path.exists(f'{BASE}/manifest.tsv'):
    with open(f'{BASE}/manifest.tsv', newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            existing[row['sample_id']] = row

if existing:
    # Preserve the active manifest order when regenerating in-place, then append
    # any newly discovered FASTQs deterministically.
    fastq_sids = [sid for sid in existing if sid in observed_fastq_sids]
    fastq_sids.extend(sid for sid in observed_fastq_sids if sid not in existing)
else:
    fastq_sids = observed_fastq_sids

missing_report = [sid for sid in fastq_sids if sid not in by_sid]
if missing_report:
    missing_existing = [sid for sid in missing_report if sid not in existing]
    if missing_existing:
        raise SystemExit(
            'FASTQs are present but missing from the bcl2fastq report and the existing manifest: '
            + ', '.join(missing_existing)
        )
    sys.stderr.write(
        'WARN: using existing manifest metadata for FASTQs absent from the bcl2fastq report: '
        + ', '.join(missing_report) + '\n'
    )
    for sid in missing_report:
        row = existing[sid]
        by_sid[sid] = {
            'sample_id': sid,
            'subject': row['subject'],
            'arm': row['arm'],
            'part': row['part'],
            'type': row['type'],
            'i7': row['i7'],
            'i5': row['i5'],
        }

extra_report = sorted(set(by_sid) - set(fastq_sids))
if extra_report:
    sys.stderr.write(
        'WARN: run-report samples without matching R1 FASTQ are excluded from manifest.tsv: '
        + ', '.join(extra_report) + '\n'
    )

recs = [by_sid[sid] for sid in fastq_sids]

# --- barcode cross-pool collision audit (index-hop risk) ---
i7set = {r['i7'] for r in recs}
i5set = {r['i5'] for r in recs}
shared = sorted(i7set & i5set)
audit = ['CoralShot barcode audit — i7/i5 cross-pool collisions', '=' * 55,
         f'{len(recs)} libraries; all i7+i5 combos unique: '
         f"{len({(r['i7'], r['i5']) for r in recs}) == len(recs)}",
         f'sequences reused across the i7 AND i5 pools: {len(shared)}', '']
flagged = set()
for seq in shared:
    as_i7 = [r['sample_id'] for r in recs if r['i7'] == seq]
    as_i5 = [r['sample_id'] for r in recs if r['i5'] == seq]
    flagged.update(as_i7 + as_i5)
    audit.append(f'  {seq}:  i7 of {as_i7}   |   i5 of {as_i5}')
audit.append('')
audit.append(f'hop_pair_flag=yes libraries ({len(flagged)}): {sorted(flagged)}')
audit.append('=> demux crosstalk matrix / 0-mismatch re-demux is a HARD GATE for these before any capture claim.')
open(f'{BASE}/barcode_audit.txt', 'w').write('\n'.join(audit) + '\n')

# --- locate fastqs ---
def fq(sid, rd):
    g = glob.glob(f'{SHOTGUN}/{sid}_S*_L002_{rd}_001.fastq.gz')
    assert len(g) == 1, f'{sid} {rd}: expected 1 file, found {g}'
    return g[0]

cols = ['sample_id', 'subject', 'arm', 'part', 'type', 'i7', 'i5', 'hop_pair_flag', 'humann', 'R1', 'R2']
with open(f'{BASE}/manifest.tsv', 'w') as fh, open(f'{BASE}/filelist.txt', 'w') as fl:
    fh.write('\t'.join(cols) + '\n')
    for r in recs:
        r['R1'], r['R2'] = fq(r['sample_id'], 'R1'), fq(r['sample_id'], 'R2')
        r['hop_pair_flag'] = 'yes' if r['sample_id'] in flagged else 'no'
        r['humann'] = existing.get(r['sample_id'], {}).get(
            'humann',
            'yes' if r['type'] == 'feces' else 'pending'   # pills decided post-S3
        )
        fh.write('\t'.join(str(r[c]) for c in cols) + '\n')
        fl.write(r['R1'] + '\n'); fl.write(r['R2'] + '\n')

# --- co-assembly units (per-participant: core+shell+feces; per-control x arm: core+shell) ---
units = {}
for r in recs:
    if r['subject'].startswith('P'):
        units.setdefault((r['subject'], 'participant'), []).append(r['sample_id'])
    else:
        units.setdefault((r['subject'] + r['arm'], 'control'), []).append(r['sample_id'])
with open(f'{BASE}/asm_units.tsv', 'w') as fh:
    fh.write('unit_id\tkind\tmember_sids\n')
    for (uid, kind), sids in sorted(units.items()):
        fh.write(f'{uid}\t{kind}\t{",".join(sorted(sids))}\n')

print(f'manifest.tsv: {len(recs)} samples | filelist.txt: {2 * len(recs)} fastqs | asm_units.tsv: {len(units)} co-assembly units')
print(f'barcode collisions: {len(shared)} shared seq(s) -> {len(flagged)} hop-flagged libraries: {sorted(flagged)}')
print(f'humann=yes: {sum(1 for r in recs if r["humann"]=="yes")}; humann!=yes: {sum(1 for r in recs if r["humann"]!="yes")}')
