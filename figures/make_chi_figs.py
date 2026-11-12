import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
"""Figures 10-11: uniform susceptibility and reduced-field scaling, from
results/susceptibility/summary.txt (branch partIX-ising-validation)."""
import re, sys
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":8.5,"axes.labelsize":9,"legend.fontsize":7.3,"xtick.labelsize":8,"ytick.labelsize":8,
  "font.family":"DejaVu Sans","mathtext.fontset":"dejavusans","axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":300})
BLUE,ORG,GRN,INK2,GRID="#2a78d6","#eb6834","#1baf7a","#52514e","#e4e3df"
SRC=sys.argv[1] if len(sys.argv)>1 else str(DATA/"susceptibility"/"summary.txt")
OUT=str(FIGS)+"/"
txt=open(SRC).read().splitlines()

def num(s):  # "0.0407(0.0001)[4]" -> (0.0407, 0.0001)
    m=re.match(r"\s*(-?[\d.]+)\((-?[\d.]+)\)",s); return float(m.group(1)),float(m.group(2))

# ---- T*chi table ----
rows=[l for l in txt if re.match(r"^\s+[\d.]+ \| ",l)]
T=np.array([float(l.split("|")[0]) for l in rows])
cols=[[num(c) for c in l.split("|")[1:]] for l in rows]
def col(i): return np.array([c[i][0] for c in cols]),np.array([c[i][1] for c in cols])
ann48,ann48e=col(1); ann10,ann10e=col(3); taf48,taf48e=col(5)
chi0=0.0407

def lab(ax,s,x=-0.15): ax.text(x,1.02,s,transform=ax.transAxes,fontsize=10,fontweight="bold",va="bottom")

fig,(a,b)=plt.subplots(1,2,figsize=(7.0,2.8))
TT=np.logspace(np.log10(0.16),1.05,200)
for y,e,c,mk,lb in [(taf48,taf48e,ORG,"s","fixed-bond TAFM ($L=48$)"),(ann48,ann48e,BLUE,"o","annealed, $\\mu=0$ ($L=48$)"),
                    (ann10,ann10e,GRN,"^","annealed, $\\mu=10$ ($L=24$)")]:
    a.errorbar(T,y/T,yerr=e/T,fmt=mk,ms=3.2,color=c,label=lb,zorder=3)
    b.errorbar(T,y,yerr=e,fmt=mk+"-",ms=3.2,lw=0.9,color=c,label=lb,zorder=3)
a.plot(TT,chi0/TT,color=INK2,ls="--",lw=0.9); a.text(2.1,0.0215,r"$\chi_0/T$",fontsize=8,color=INK2)
hi=TT>1.4
for th,c in [(3.0,ORG),(2.4,BLUE),(2.0,GRN)]:
    a.plot(TT[hi],0.25/(TT[hi]+th),color=c,ls=":",lw=1)
    b.plot(TT[hi],0.25*TT[hi]/(TT[hi]+th),color=c,ls=":",lw=1)
a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("$T$"); a.set_ylabel(r"$\chi=\mathrm{var}(n_1)/(N_{\rm site}T)$")
a.legend(loc="lower left",frameon=False); a.grid(color=GRID,lw=0.5); lab(a,"(a)")
b.axhline(chi0,color=INK2,ls="--",lw=0.8); b.text(3.0,chi0+0.006,r"$\chi_0=0.0407$",fontsize=7.5,color=INK2)
b.axhline(0.25,color=INK2,ls=":",lw=0.8); b.text(5.0,0.254,"free spins, 1/4",fontsize=7.5,color=INK2,va="bottom",ha="center")
b.text(4.2,0.105,"Curie–Weiss\n"+r"$\theta=3,\ 2.4,\ 2$",fontsize=7.2,color=INK2,ha="center",va="top")
b.set_xscale("log"); b.set_xlabel("$T$"); b.set_ylabel(r"$T\chi$"); b.set_ylim(0,0.27)
b.legend(loc="upper left",bbox_to_anchor=(0.0,0.88),frameon=False); b.grid(color=GRID,lw=0.5); lab(b,"(b)")
fig.tight_layout(w_pad=2.2); fig.savefig(OUT+"fig10_chi.pdf"); fig.savefig(OUT+"fig10_chi.png",dpi=200); plt.close(fig)

# ---- m(H) ----
mgs={}
for l in txt:
    m=re.match(r"#\s+H=([\d.]+)\s+L=24 ord: ([\d.]+)\(([\d.]+)\)\s+L=24 gs: \S+\s+L=48 ord: ([\d.]+)\(([\d.]+)\)",l)
    if m: mgs[float(m.group(1))]=(float(m.group(4)),float(m.group(5)))
dat={}
for l in txt:
    m=re.match(r"#\s+L=48 T=([\d.]+) H=([\d.]+)\s+tafm\s+([\d.]+)\(([\d.]+)\) d=([+-][\d.]+) \| annealed\s+([\d.]+)\(([\d.]+)\) d=([+-][\d.]+)",l)
    if m:
        t,h=float(m.group(1)),float(m.group(2))
        dat.setdefault(t,[]).append((h,float(m.group(3)),float(m.group(4)),float(m.group(5)),float(m.group(6)),float(m.group(7)),float(m.group(8))))
Hs=np.array(sorted(mgs)); mg=np.array([mgs[h][0] for h in Hs])
fig,(a,b)=plt.subplots(1,2,figsize=(7.0,2.9))
a.plot(Hs,mg,color="black",lw=1.2,label=r"$m_{\rm GS}(H)$",zorder=2)
for t,c in [(0.2,"#7b3fb5"),(0.5,"#b8322a")]:
    d=np.array(sorted(dat[t]))
    a.plot(d[:,0],d[:,1],"s",ms=3.6,mfc="none",mec=c,mew=1.0,label=f"TAFM, $T={t}$",zorder=3)
    a.plot(d[:,0],d[:,4],"o",ms=3.0,color=c,label=f"annealed, $T={t}$",zorder=4)
a.axvline(0.266,color=INK2,ls=":",lw=0.8); a.text(0.28,0.497,r"$H_c=0.266$",fontsize=7.2,color=INK2,va="top")
a.set_xlabel(r"$H=D/(2T)$"); a.set_ylabel(r"$m=\langle n_1\rangle/N_{\rm site}$  ($L=48$)")
a.legend(loc="lower left",frameon=False); a.grid(color=GRID,lw=0.5); lab(a,"(a)")
shades={0.3:0.45,0.4:0.72,0.5:1.0}
for t in (0.3,0.4,0.5):
    d=np.array(sorted(dat[t])); al=shades[t]
    b.plot(d[:,0],d[:,3]*1e3,"s-",ms=3,lw=1,color=ORG,alpha=al,label=f"TAFM, $T={t}$")
    b.plot(d[:,0],d[:,6]*1e3,"o--",ms=3,lw=1,color=BLUE,alpha=al,label=f"annealed, $T={t}$")
b.axhline(0,color=INK2,lw=0.6)
b.set_xlabel(r"$H=D/(2T)$"); b.set_ylabel(r"$(m-m_{\rm GS})\times10^{3}$  ($L=48$)")
b.legend(loc="upper left",frameon=False,ncol=2,columnspacing=1.0,handlelength=2.2); b.set_ylim(-1,22); b.grid(color=GRID,lw=0.5); lab(b,"(b)")
fig.tight_layout(w_pad=2.2); fig.savefig(OUT+"fig11_mH.pdf"); fig.savefig(OUT+"fig11_mH.png",dpi=200); plt.close(fig)
print("ok", len(T), len(Hs), sorted(dat))
