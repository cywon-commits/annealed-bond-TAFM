import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D
from annealed_thermo import build_arrays, total_energy, _seed, z3_order
from ising_sweep import sweep_ising
from rhombile_lattice import A_MAT, A_MAT_INV
from worm_monomer_walk import DimerState
plt.rcParams.update({"font.size":8.5,"font.family":"DejaVu Sans","savefig.dpi":300})
L=24; A=build_arrays(L); st=DimerState(L,L); pos,_=st.lat.site_positions()
tri=np.array([t["sites"] for t in st.tri]); ref0=A_MAT@np.array([L/2,L/2])
def near(p,ref):
    d=A_MAT_INV@(p-ref); d-=np.round(d/L)*L; return ref+A_MAT@d
def run(D2,T,n,seed,start="random"):
    _seed(seed); np.random.seed(seed); rng=np.random.default_rng(seed); dim,mate=A["dimer"].copy(),A["mate"].copy()
    s=np.where(A["sub"]==0,1,0).astype(np.int64) if start=="crystal" else rng.integers(2,size=A["ns"]).astype(np.int64)
    E=total_energy(s,dim,A["bi"],A["bj"],D2,0.0,0.0,mate)
    if n==0: return s,dim,mate,E
    Ts=[T] if (T>5 or D2>0) else list(np.linspace(2.0,T,10))
    for Tk in Ts:
        for _ in range(n//len(Ts)): E=sweep_ising(s,dim,mate,A["bi"],A["bj"],A["bta"],A["btb"],A["site_b"],A["hexb"],1/Tk,D2,0.0,E)
    return s,dim,mate,E
def stats(s,dim,mate):
    ss=s[tri]; df=((ss[:,0]==ss[:,1])&(ss[:,1]==ss[:,2])).mean()
    return (mate<0).mean(),df,z3_order(dim,A["site_b"],A["sub"])
X0,X1,Y0,Y1=6,17,6,15
RW,RH=10.0,7.6            # clean rectangular crop (lattice units), centred on ref0
DIMC="#86ace0"
def inrect(p,pad=0.0):
    return abs(p[0]-ref0[0])<=RW/2+pad and abs(p[1]-ref0[1])<=RH/2+pad
def tri_xy(t):
    sites=tri[t]; p0=pos[sites[0]]; p=np.array([near(pos[i],p0) for i in sites]); c=p.mean(0); cc=near(c,ref0)
    return p+(cc-c), cc
def draw(ax,s,dim,mate,rw=None,rh=None):
    global RW,RH
    RW0,RH0=RW,RH
    if rw: RW,RH=rw,rh
    for t in range(len(tri)):
        p,c=tri_xy(t)
        if not inrect(c,1.0): continue
        sv=s[list(tri[t])]
        if mate[t]<0: ax.add_patch(Polygon(p,fc="#f5b7b1",ec="none",zorder=1))
        if sv[0]==sv[1]==sv[2]: ax.add_patch(Polygon(p,fc="none",ec="#e0a100",lw=2.2,zorder=4))
    for k in range(A["nb"]):
        pi=pos[A["bi"][k]]; pj=near(pos[A["bj"][k]],pi); m=(pi+pj)/2; mm=near(m,ref0); pi,pj=pi+(mm-m),pj+(mm-m)
        if not inrect(mm,1.0): continue
        if dim[k]:                                   # dimer: segment joining the two paired triangles
            _,ca=tri_xy(A["bta"][k]); ca=near(ca,mm); _,cb=tri_xy(A["btb"][k]); cb=near(cb,mm)
            ax.plot([ca[0],cb[0]],[ca[1],cb[1]],color=DIMC,lw=1.1,zorder=2,solid_capstyle="round")
            continue
        eq=s[A["bi"][k]]==s[A["bj"][k]]
        ax.plot(*np.array([pi,pj]).T,color="#c0392b" if eq else "#6b7280",lw=2.2 if eq else 0.55,zorder=3 if eq else 2,solid_capstyle="round")
    for i in range(A["ns"]):
        p=near(pos[i],ref0)
        if inrect(p,0.5):
            ax.scatter(*p,s=8,c="#1d2a44" if s[i]==1 else "white",edgecolors="#1d2a44",lw=0.6,zorder=5)
    ax.set_xlim(ref0[0]-RW/2,ref0[0]+RW/2); ax.set_ylim(ref0[1]-RH/2,ref0[1]+RH/2); ax.set_aspect("equal"); ax.axis("off")
    RW,RH=RW0,RH0
H=[Line2D([],[],color=DIMC,lw=1.4,label="dimer: paired triangles (shared bond off)"),Patch(fc="#f5b7b1",label="frustrated (unpaired) triangle"),
   Patch(fc="none",ec="#e0a100",lw=2,label="defect triangle (three equal spins)"),Line2D([],[],color="#c0392b",lw=2,label="unsatisfied bond (active, equal spins)"),
   Line2D([],[],ls="",marker="o",ms=4,mfc="#1d2a44",mec="#1d2a44",label="spin 1"),Line2D([],[],ls="",marker="o",ms=4,mfc="white",mec="#1d2a44",label="spin 0")]
# Fig 2: ground states
fig,axs=plt.subplots(1,2,figsize=(7.0,2.7))
s,dim,mate,E=run(1.0,0,0,1,"crystal"); print("crystal E/Ntri",E/A["nt"],"expect",1/6, stats(s,dim,mate))
draw(axs[0],s,dim,mate); axs[0].set_title("(a)  $0<D<6$, $T=0$: unique rhombile crystal\n(hubs carry the $D$ state; $E/N_{\\rm tri}=D/6$)",fontsize=8)
def exact_ground_state(seed, nflip=400000):
    """Exact T=0 configuration in zero field: a random lozenge tiling (hexagon flips from the
    rhombile tiling, no frustrated triangles) carrying its zero-energy two-colouring."""
    rng=np.random.default_rng(seed); dim,mate=A["dimer"].copy(),A["mate"].copy(); hexb=A["hexb"]
    for v in rng.integers(A["ns"],size=nflip):
        d=dim[hexb[v]]
        if (d[0] and d[2] and d[4] and not (d[1] or d[3] or d[5])) or (d[1] and d[3] and d[5] and not (d[0] or d[2] or d[4])):
            for k in hexb[v]:
                dim[k]=not dim[k]
                if dim[k]: mate[A["bta"][k]]=k; mate[A["btb"][k]]=k
    s=-np.ones(A["ns"],np.int64); s[0]=int(rng.integers(2)); stack=[0]
    while stack:                                   # BFS 2-colouring of the active-bond graph
        v=stack.pop()
        for k in A["site_b"][v]:
            if dim[k]: continue
            j=A["bj"][k] if A["bi"][k]==v else A["bi"][k]
            if s[j]<0: s[j]=1-s[v]; stack.append(j)
            else: assert s[j]!=s[v], "active-bond graph not bipartite"
    E=total_energy(s,dim,A["bi"],A["bj"],0.0,0.0,0.0,mate)
    return s,dim,mate,E
s,dim,mate,E=exact_ground_state(5); fr,df,p2=stats(s,dim,mate); print("exact GS: E",E,"frustrated",fr,"defect",df,"psi2",p2,"lozenges on sublattice-hub sites", (A["hexb"].shape))
assert E==0 and fr==0
draw(axs[1],s,dim,mate); axs[1].set_title("(b)  $D=0$, $T=0$: a random lozenge tiling\nwith its zero-energy spins (a TAFM ground state)",fontsize=8)
fig.legend(handles=H,loc="lower center",ncol=3,fontsize=7,frameon=False,handlelength=1.6,columnspacing=1.2)
fig.subplots_adjust(left=0.01,right=0.99,top=0.84,bottom=0.16,wspace=0.05); fig.savefig(str(FIGS)+"/fig2_groundstates.pdf"); fig.savefig(str(FIGS)+"/fig2_groundstates.png",dpi=200); plt.close(fig)
# Fig 3: temperature series
cases=[(0.0,0.30,"(a)  $D=0$, $T=0.30$"),(0.0,0.45,"(b)  $D=0$, $T=0.45$"),(0.0,0.70,"(c)  $D=0$, $T=0.70$"),
       (1.0,0.48,"(d)  $D=1$, $T=0.48<T_c$"),(1.0,0.60,"(e)  $D=1$, $T=0.60>T_c$"),(0.0,1e6,"(f)  $T=\\infty$")]
fig,axs=plt.subplots(2,3,figsize=(7.0,3.75))
def heat_from_gs(T,n,seed):
    """Zero field, finite T: start from an exact ground state and equilibrate at T (heating),
    so that no defect pairs are left trapped by cooling."""
    s,dim,mate,E=exact_ground_state(100+seed)
    _seed(seed); np.random.seed(seed)
    for _ in range(n): E=sweep_ising(s,dim,mate,A["bi"],A["bj"],A["bta"],A["btb"],A["site_b"],A["hexb"],1/T,0.0,0.0,E)
    return s,dim,mate,E
for ax,(D2,T,tt),seed in zip(axs.ravel(),cases,[11,12,13,14,15,21]):
    s,dim,mate,E=heat_from_gs(T,60000,seed) if (D2==0 and T<5) else run(D2,T,60000,seed); fr,df,p2=stats(s,dim,mate); print(tt,fr,df,p2)
    draw(ax,s,dim,mate,7.6,5.8); ax.set_title(f"{tt}\nfrustr. {fr:.3f}, defect {df:.4f}, $|\\psi|^2$ {p2:.3f}",fontsize=7.8)
fig.legend(handles=H,loc="lower center",ncol=3,fontsize=7,frameon=False,handlelength=1.6,columnspacing=1.2)
fig.subplots_adjust(left=0.01,right=0.99,top=0.9,bottom=0.12,wspace=0.04,hspace=0.32); fig.savefig(str(FIGS)+"/fig3_snapshots.pdf"); fig.savefig(str(FIGS)+"/fig3_snapshots.png",dpi=200)
