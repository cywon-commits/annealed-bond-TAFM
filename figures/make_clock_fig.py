"""Fig.: how the field D acts on the phase of the sqrt3 x sqrt3 order parameter.
(a) n1/N = 1/2 - cos(3 theta)/6 over the six manifold states; (b) a uniform XY field (p=1)
versus the threefold clock field (p=3); (c) relevance of exp(i p theta): x_p = p^2 eta/2."""
import sys as _sys, pathlib as _pl, os
_here=_pl.Path(__file__).resolve()
try:
    _sys.path.insert(0,str(_here.parents[1]/"src")); from repo_paths import FIGS
except Exception:
    FIGS=_pl.Path(os.environ.get("OUT_DIR","/home/claude/paper"))/"figs"
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.colors as mcol
from matplotlib.patches import Polygon
plt.rcParams.update({"font.size":8.5,"axes.labelsize":9,"legend.fontsize":7.2,"xtick.labelsize":8,"ytick.labelsize":8,
  "font.family":"DejaVu Sans","mathtext.fontset":"dejavusans","axes.spines.top":False,"axes.spines.right":False})
INK2,GRID="#52514e","#e4e3df"; SQ3=np.sqrt(3)
def hue(k,v=1.0): return mcol.hsv_to_rgb([(k/6)%1,0.75,v])
UPS={0:{0},1:{0,1},2:{1},3:{1,2},4:{2},5:{2,0}}
def patch_inset(ax,cx,cy,k,w=0.95,h=0.19):
    ia=ax.inset_axes([cx-w/2,cy-h/2,w,h],transform=ax.transData); ia.set_axis_off(); ia.set_aspect("equal")
    ia.set_xlim(-1.75,1.75); ia.set_ylim(-1.75,1.75); patch(ia,0,0,k,a=1.0)
def patch(ax,cx,cy,k,a):
    pts=[(i,j) for i in range(-1,2) for j in range(-1,2) if abs(i+j)<=1]
    P={q:np.array([cx+a*(q[0]+0.5*q[1]),cy+a*SQ3/2*q[1]]) for q in pts}
    hexv=[np.array([cx,cy])+a*1.55*np.array([np.cos(np.pi/6+m*np.pi/3),np.sin(np.pi/6+m*np.pi/3)]) for m in range(6)]
    ax.add_patch(Polygon(hexv,fc=hue(k),ec="none",alpha=0.55,zorder=1,transform=ax.transData,clip_on=False))
    for q in pts:
        for d in [(1,0),(0,1),(-1,1)]:
            r=(q[0]+d[0],q[1]+d[1])
            if r in P: ax.plot(*np.c_[P[q],P[r]],color="#9aa0a6",lw=0.5,zorder=2,clip_on=False)
    for q in pts:
        on=((q[0]-q[1])%3) in UPS[k]
        ax.scatter(*P[q],s=10,c="#1d2a44" if on else "white",edgecolors="#1d2a44",lw=0.7,zorder=4,clip_on=False)
def lab(ax,s): ax.text(-0.13,1.03,s,transform=ax.transAxes,fontsize=10,fontweight="bold",va="bottom")

fig=plt.figure(figsize=(7.0,2.5))
a=fig.add_axes([0.065,0.18,0.30,0.70]); b=fig.add_axes([0.445,0.18,0.24,0.70]); c=fig.add_axes([0.755,0.18,0.215,0.70])
# (a)
th=np.linspace(0,2*np.pi,400); a.plot(th,0.5-np.cos(3*th)/6,color=INK2,lw=1.2,zorder=3)
for k in range(6):
    t=k*np.pi/3; y=0.5-np.cos(3*t)/6
    a.scatter(t,y,s=34,color=hue(k),edgecolors="k",lw=0.6,zorder=5)
    patch_inset(a,t,y+(0.135 if k%2 else -0.135),k)
a.set_xlim(-0.45,2*np.pi+0.45); a.set_ylim(0.10,0.92)
# keep patches round: set aspect so x and y units comparable for the patch geometry
a.set_xticks([k*np.pi/3 for k in range(7)]); a.set_xticklabels(["0","60°","120°","180°","240°","300°","360°"])
a.set_yticks([1/3,1/2,2/3]); a.set_yticklabels(["1/3","1/2","2/3"])
a.set_xlabel(r"$\theta=\arg\psi_{\rm loc}$"); a.set_ylabel(r"$n_1/N_{\rm site}$")
a.text(0.56,1.03,r"$n_1/N_{\rm site}=\frac{1}{2}-\frac{1}{6}\cos3\theta$",transform=a.transAxes,ha="center",va="bottom",fontsize=8.4)
a.grid(color=GRID,lw=0.5,zorder=0); lab(a,"(a)")
# (b)
for p,ls,c_,lb in [(1,"--",INK2,"XY in uniform field, $-h\\cos\\theta$"),(3,"-","k","this model, $-\\frac{D}{6}\\cos3\\theta$")]:
    b.plot(th,-np.cos(p*th),ls=ls,color=c_,lw=1.2,label=lb)
for k in range(7):
    b.axvspan(k*np.pi/3-np.pi/6,k*np.pi/3+np.pi/6,color=hue(k),alpha=0.18,lw=0)
b.scatter([0,2*np.pi/3,4*np.pi/3],[-1,-1,-1],s=30,color=[hue(0),hue(2),hue(4)],edgecolors="k",lw=0.6,zorder=5)
b.scatter([0],[-1.0],s=60,facecolors="none",edgecolors=INK2,lw=1.0,ls="--",zorder=6)
b.set_xlim(-0.3,2*np.pi+0.3); b.set_ylim(-1.25,2.6); b.set_yticks([])
b.set_xticks([k*np.pi/3 for k in range(7)]); b.set_xticklabels(["0","","120°","","240°","","360°"])
b.set_xlabel(r"$\theta$"); b.set_ylabel("energy per site (arb.)")
b.legend(loc="upper left",frameon=True,framealpha=1,edgecolor="none",fontsize=6.4,handlelength=1.8,borderaxespad=0.15)
lab(b,"(b)")
# (c)
eta=np.linspace(0.3,0.6,50)
c.axvspan(0.3,4/9,color="#86ace0",alpha=0.22,lw=0)
c.plot(eta,eta/2,ls="--",color=INK2,lw=1.1,label=r"$p=1$: $\eta/2$")
c.plot(eta,9*eta/2,color="k",lw=1.3,label=r"$p=3$: $9\eta/2$")
c.axhline(2,color="#b8322a",lw=0.9); c.text(0.305,2.05,"marginal",fontsize=6.8,color="#b8322a",ha="left",va="bottom")
c.text(0.372,2.75,"locked",fontsize=7.2,color="#2a5aa0",ha="center"); c.text(0.545,2.75,"rough",fontsize=7.2,color=INK2,ha="center")
c.plot([0.5],[2.25],"o",color="k",ms=4,zorder=5); c.text(0.512,1.93,"TAFM, $H=0$",fontsize=6.4,va="top")
c.plot([4/9],[2.0],"s",color="#b8322a",ms=4,zorder=5); c.text(4/9-0.006,1.85,r"$H_{\rm KT}$",fontsize=6.8,ha="right",va="top",color="#b8322a")
c.annotate("",xy=(0.452,2.06),xytext=(0.494,2.25),arrowprops=dict(arrowstyle="-|>",color="#2a5aa0",lw=1.0))
c.text(0.47,2.42,"$D$ stiffens",fontsize=6.4,color="#2a5aa0",ha="center")
c.set_xlim(0.3,0.6); c.set_ylim(0,3.1)
c.set_xlabel(r"$\eta$  (stiffer to the left)"); c.set_ylabel(r"dimension of $e^{ip\theta}$")
c.legend(loc="center left",frameon=False,fontsize=6.6,bbox_to_anchor=(0.0,0.42)); c.grid(color=GRID,lw=0.5); lab(c,"(c)")
for e in ("pdf","png"): fig.savefig(FIGS/f"fig13_clock.{e}",dpi=200 if e=="png" else None)
print("ok")
