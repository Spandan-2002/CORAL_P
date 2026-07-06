#!/usr/bin/env python3
"""Stacked genus relative-abundance composition. Controls in a SEPARATE block on the LEFT (each control
shown individually, cores then shells), a DASHED divider, then participants grouped STOOL | CORE | SHELL.
Run: conda run -p /scratch/sr7729/conda_envs/cs_viz python scripts/make_stacked_composition.py"""
import os, pandas as pd, numpy as np, re
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
B=os.environ.get("BASE","/scratch/sr7729/CORAL_P")   # portable: honor $BASE, else default to this repo
df=pd.read_csv(f"{B}/results/04_mpa/merged_sgb.tsv",sep="\t",skiprows=1)
df.columns=[c.replace(".profile","") for c in df.columns]; tc=df.columns[0]; cn=df[tc].astype(str)
g=df[cn.str.contains(r"\|g__",na=False)&~cn.str.contains(r"\|s__",na=False)].copy()
g.index=[x.split("g__")[-1].split("|")[0] for x in g[tc]]
cols=[c for c in df.columns if re.match(r'^[CP]\d+[MB]_[CSF]$',c)]
G=g[cols]/g[cols].sum(0).clip(lower=1e-9)
CLEAN={"P01":"M","P02":"M","P03":"M","P04":"M","P05":"M","P08":"B","P09":"B","P10":"B","P12":"B","P13":"M"}
def pcol(p,arm,c):
    s=f"{p}{arm}_{c}"; return G[s] if s in G.columns else None
ctrl_core=sorted([c for c in G.columns if c.startswith("C") and c.endswith("_C")])
ctrl_shell=sorted([c for c in G.columns if c.startswith("C") and c.endswith("_S")])
# CONTROLS block (left): cores then shells, individual
controls=[(f"{c[:-2]}·c", G[c], True) for c in ctrl_core] + [(f"{c[:-2]}·s", G[c], True) for c in ctrl_shell]
# participant compartment groups (no controls)
stool=[(p, pcol(p,a,"F"), False) for p,a in CLEAN.items() if pcol(p,a,"F") is not None]
core =[(p, pcol(p,a,"C"), False) for p,a in CLEAN.items() if pcol(p,a,"C") is not None]
shell=[(p, pcol(p,a,"S"), False) for p,a in CLEAN.items() if pcol(p,a,"S") is not None]
groups=[("CONTROLS",controls),("STOOL",stool),("CORE",core),("SHELL",shell)]
allser=[s for _,gb in groups for _,s,_ in gb]
top=pd.concat(allser,axis=1).mean(1).sort_values(ascending=False).index[:14].tolist()
cmap=plt.cm.tab20(np.linspace(0,1,len(top))); COL=dict(zip(top,cmap))
# layout (vertical columns, groups left→right)
recs=[]; xi=0.0; centers={}; gstart={}
for gname,gb in groups:
    gstart[gname]=xi; start=xi
    for lab,s,isc in gb: recs.append((xi,lab,s,isc)); xi+=1.0
    centers[gname]=(start+xi-1.0)/2; xi+=2.0
nb=len(recs); W=0.8
fig,ax=plt.subplots(figsize=(max(20,0.44*nb),8))
for xx,lab,s,isc in recs:
    topa=100.0                                    # stack TOP-DOWN: first genus anchored at the top
    for gen in top:
        h=s.get(gen,0.0)*100; ax.bar(xx,h,W,bottom=topa-h,color=COL[gen],edgecolor="white",lw=0.15); topa-=h
    ax.bar(xx,max(0,topa),W,bottom=0.0,color="#d9d9d9",edgecolor="white",lw=0.15)   # "Other" fills the bottom
    ax.text(xx,-1.5,lab,ha="center",va="top",rotation=90,fontsize=7,color="#c0392b" if isc else "#333")
# dashed divider between CONTROLS and STOOL
sep=gstart["STOOL"]-1.0
ax.axvline(sep,ls="--",color="#444",lw=1.6)
# group labels
for gname,cx in centers.items():
    ax.text(cx,-15,gname,ha="center",va="top",fontsize=15,fontweight="bold",
            color="#c0392b" if gname=="CONTROLS" else "#000")
ax.set_ylim(-17,100); ax.set_xlim(-1,xi-1.5); ax.set_xticks([])
ax.set_ylabel("genus relative abundance (%)")
ax.set_title("Stacked genus composition — controls (left, individual) vs participants (STOOL | CORE | SHELL)\n"
             "controls resemble a shared participant-gut blend (cross-talk), not a reagent/kit-ome profile",
             fontsize=13,fontweight="bold")
ax.legend(handles=[Patch(fc=COL[t],label=t) for t in top]+[Patch(fc="#d9d9d9",label="Other")],
          fontsize=9,ncol=8,loc="upper center",bbox_to_anchor=(0.5,-0.05),frameon=False,title="genus (top 14)")
for sp in ["top","right","bottom"]: ax.spines[sp].set_visible(False)
ax.spines["left"].set_bounds(0,100)   # y-axis line only spans 0-100, not the label area
fig.tight_layout(rect=[0,0.05,1,1]); fig.savefig(f"{B}/figures/stacked_composition.png",dpi=150,bbox_inches="tight")
print(f"wrote figures/stacked_composition.png ; CONTROLS {len(controls)} | STOOL {len(stool)} | CORE {len(core)} | SHELL {len(shell)}")
