import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
SNAP = DATA / "snapshots"
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.colors as mc
from matplotlib.collections import PolyCollection
from rhombile_lattice import A_MAT, A_MAT_INV
plt.rcParams.update({"font.size":8.5,"font.family":"DejaVu Sans","mathtext.fontset":"dejavusans"})
G=np.load(SNAP/"geo48.npz"); pos,tri,nbr,sub=G["pos"],G["tri"],G["nbr"],G["sub"]; L=48; N=len(pos)
ref=A_MAT@np.array([L/2,L/2])
def wrap(p):
    d=(A_MAT_INV@(p-ref).T).T; d-=np.round(d/L)*L; return ref+(A_MAT@d.T).T
P=wrap(pos)
# coarse-graining clusters: graph distance <= R
def clusters(R):
    out=[]
    for v in range(N):
        seen={v}; fr=[v]
        for _ in range(R):
            fr=[u for x in fr for u in nbr[x] if u not in seen]; seen.update(fr)
        out.append(np.array(sorted(seen)))
    return out
CL=clusters(2); W=np.exp(2j*np.pi*sub/3)
def local_psi(s):
    a=(2*s-1)*W; return np.array([a[c].sum() for c in CL])/len(CL[0])
r=1/3; hexv=np.array([[r*np.cos(k*np.pi/3),r*np.sin(k*np.pi/3)] for k in range(6)])
def panel(ax,s,mate,title,RW=34,RH=26,show_fr=True,amp=0.6,dsz=9,fsz=1.6):
    psi=local_psi(s); hue=(np.angle(psi)/(2*np.pi))%1; val=np.clip(np.abs(psi)/amp,0,1)
    rgb=mc.hsv_to_rgb(np.c_[hue,0.75*np.ones(N),0.25+0.75*val])
    keep=(np.abs(P[:,0]-ref[0])<=RW/2+0.4)&(np.abs(P[:,1]-ref[1])<=RH/2+0.4)
    ax.add_collection(PolyCollection([P[i]+hexv for i in np.where(keep)[0]],facecolors=rgb[keep],edgecolors="face",linewidths=0.15))
    c=wrap(np.array([pos[t].mean(0) if np.ptp(pos[t][:,0])<2 and np.ptp(pos[t][:,1])<2 else pos[t[0]] for t in tri]))
    ss=s[tri]; dfe=(ss[:,0]==ss[:,1])&(ss[:,1]==ss[:,2]); fr=(mate<0)&~dfe
    inr=(np.abs(c[:,0]-ref[0])<=RW/2)&(np.abs(c[:,1]-ref[1])<=RH/2)
    if show_fr: ax.scatter(*c[fr&inr].T,s=fsz,c="white",edgecolors="none",zorder=3)
    ax.scatter(*c[dfe&inr].T,s=dsz,c="black",edgecolors="white",linewidths=0.3,zorder=4)
    ax.set_xlim(ref[0]-RW/2,ref[0]+RW/2); ax.set_ylim(ref[1]-RH/2,ref[1]+RH/2); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]); [sp.set_linewidth(0.5) for sp in ax.spines.values()]
    ax.set_title(title,fontsize=7.8)
    return dfe.mean(),fr.mean()
def wheel(ax,ticks=True):
    x=np.linspace(-1,1,201); X,Y=np.meshgrid(x,x); R=np.hypot(X,Y); H=(np.arctan2(Y,X)/(2*np.pi))%1
    rgb=mc.hsv_to_rgb(np.dstack([H,0.75*np.ones_like(H),0.25+0.75*np.clip(R,0,1)])); rgba=np.dstack([rgb,(R<=1).astype(float)])
    ax.imshow(rgba,extent=(-1,1,-1,1),origin="lower"); ax.set_axis_off(); ax.set_aspect("equal")
    if ticks:
        for k in range(6):
            a=k*np.pi/3; ax.plot([0.78*np.cos(a),1.0*np.cos(a)],[0.78*np.sin(a),np.sin(a)],color="k",lw=0.6)
