import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
exec(open(ROOT/"analysis"/"rg_line"/"rg_line.py").read().split("res={}")[0])   # reuse n_d mapping + anchor points
plt.rcParams.update({"font.size":8.5,"axes.labelsize":9,"legend.fontsize":7.2,"xtick.labelsize":8,"ytick.labelsize":8,
  "font.family":"DejaVu Sans","mathtext.fontset":"dejavusans","axes.spines.top":False,"axes.spines.right":False})
BLUE,ORG,INK2,GRID="#2a78d6","#eb6834","#52514e","#e4e3df"
Ds=np.r_[np.linspace(0.02,0.1,9),np.linspace(0.12,1.2,28)]
fig,(a,b)=plt.subplots(1,2,figsize=(7.0,2.8))
for HKT,ls in [(0.266,"-"),(0.305,"--")]:
    P=pts[:3]; dH=P[:,0]-HKT; A=np.c_[np.sqrt(dH),dH]; c,*_=np.linalg.lstsq(A,P[:,1],rcond=None)
    Tline=lambda H: c[0]*np.sqrt(max(H-HKT,0))+c[1]*(H-HKT); Hline=lambda Tp: brentq(lambda H: Tline(H)-Tp,HKT+1e-12,5)
    tt=np.array([brentq(lambda t: D/(2*t)-Hline(t),0.01,0.9) for D in Ds])
    ta=np.array([brentq(lambda t: D/(2*t)-Hline(Tprime(t)),0.012,0.75) for D in Ds])
    lb="" if ls=="--" else " (RG, $H_{\\rm KT}=0.266$)"
    a.plot(Ds,tt,color=ORG,ls=ls,lw=1.2,label="TAFM line"+(" (fit)" if ls=="-" else "") if ls=="-" else None)
    a.plot(Ds,ta,color=BLUE,ls=ls,lw=1.4,label="annealed, predicted" if ls=="-" else None)
    b.plot(Ds,ta/tt,color=INK2,ls=ls,lw=1.3,label=f"$H_{{\\rm KT}}={HKT}$")
    a.plot(Ds,Ds/(2*HKT),color=INK2,ls=":" if ls=="-" else (0,(1,3)),lw=0.8)
a.text(0.215,0.47,r"$D/(2H_{\rm KT})$",fontsize=7.5,color=INK2,rotation=63,ha="right")
# data
a.plot([0.1025,0.5,1.0],[0.125,0.3347,0.4713],"s",ms=4,mfc="white",mec=ORG,mew=1.1,label="TAFM, CTMRG / MC")
a.errorbar([1.0],[0.528],yerr=[0.002],fmt="o",ms=4.2,color=BLUE,label=r"annealed, MC ($\mu=0$)")
a.errorbar([0.5],[0.456],yerr=[0.013],fmt="o",ms=4.2,mfc="white",mec=BLUE,ecolor=BLUE,label=r"annealed, $L\leq64$ crossing")
a.set_xlim(0,1.2); a.set_ylim(0,0.7); a.set_xlabel("$D$"); a.set_ylabel(r"$T_c$")
a.legend(loc="lower right",frameon=False,fontsize=6.8); a.grid(color=GRID,lw=0.5)
a.text(-0.15,1.02,"(a)",transform=a.transAxes,fontsize=10,fontweight="bold",va="bottom")
b.errorbar([1.0],[0.528/0.4713],yerr=[0.005],fmt="o",ms=4.2,color=BLUE,label="MC, $D=1$")
b.axhline(1,color=INK2,lw=0.6); b.set_xlim(0,1.2); b.set_ylim(0.98,1.25)
b.set_xlabel("$D$"); b.set_ylabel(r"$T_c^{\rm ann}/T_c^{\rm TAFM}$"); b.legend(loc="upper right",frameon=False); b.grid(color=GRID,lw=0.5)
b.text(-0.15,1.02,"(b)",transform=b.transAxes,fontsize=10,fontweight="bold",va="bottom")
fig.tight_layout(w_pad=2.2)
for e in ("pdf","png"): fig.savefig(FIGS/f"fig12_rgline.{e}",dpi=200 if e=="png" else None)
