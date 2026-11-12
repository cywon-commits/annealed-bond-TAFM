import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
import csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
D=str(FIGDATA)+"/"
plt.rcParams.update({"font.size":8.5,"axes.labelsize":9,"axes.titlesize":9,"legend.fontsize":7.5,
  "xtick.labelsize":8,"ytick.labelsize":8,"font.family":"DejaVu Sans","mathtext.fontset":"dejavusans",
  "axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":300})
BLUE,ORG,GRN,YEL="#2a78d6","#eb6834","#1baf7a","#eda100"; INK2="#52514e"; GRID="#e4e3df"
LCOL={24:YEL,32:BLUE,48:ORG,64:GRN,96:"#7b3fb5"}
def rd(f): return list(csv.DictReader(open(D+f)))
def save(fig,name):
    fig.savefig(FIGS/f"{name}.pdf"); fig.savefig(FIGS/f"{name}.png",dpi=200); plt.close(fig)
def lab(ax,s,x=-0.16,y=1.02): ax.text(x,y,s,transform=ax.transAxes,fontsize=10,fontweight="bold",va="bottom")

# ---------- Fig 1: mechanism schematic ----------
def rh(ax,cx,cy,mode,s=1.0):
    xl,xr,yt,yb=cx-0.5*s,cx+0.5*s,cy+0.866*s,cy-0.866*s
    red="#c0392b"; edge="#6b7280"
    if mode=="mono":
        for yy in (yt,yb): ax.add_patch(Polygon([(xl,cy),(xr,cy),(cx,yy)],fc="#f5b7b1",ec="none"))
    else:
        pass
    up=red if mode=="defect" else edge; w=2.4 if mode=="defect" else 1.0
    ax.plot([xl,cx,xr],[cy,yt,cy],color=up,lw=w,solid_capstyle="round",zorder=2)
    ax.plot([xl,cx,xr],[cy,yb,cy],color=edge,lw=1.0,zorder=2)
    if mode in("ann","defect"):
        ax.plot([xl,xr],[cy,cy],color=edge,lw=1.0,ls=(0,(2.5,2)),zorder=2)
        ax.plot([cx,cx],[cy+0.866*s/3,cy-0.866*s/3],color="#86ace0",lw=2.2,zorder=3,solid_capstyle="round")
    else: ax.plot([xl,xr],[cy,cy],color=red,lw=2.4,zorder=2)
    for (x,y,on) in [(xl,cy,1),(xr,cy,1),(cx,yt,mode=="defect"),(cx,yb,0)]:
        ax.scatter([x],[y],s=46,c="#1d2a44" if on else "white",edgecolors="#1d2a44",lw=0.9,zorder=4)
import matplotlib.colors as mcol
from matplotlib.patches import FancyArrowPatch, Circle
SQ3=np.sqrt(3)
def hue(k,v=1.0): return mcol.hsv_to_rgb([(k/6)%1,0.75,v])
def tri_xy(i,j,a=1.0): return np.array([a*(i+0.5*j),a*SQ3/2*j])
UPS={0:{0},1:{0,1},2:{1},3:{1,2},4:{2},5:{2,0}}          # sublattices in state 1 for clock state k
def site_state(i,j,k): return ((i-j)%3) in UPS[k]
def patch(ax,cx,cy,k,a,R=1,ms=16,bond=True,tint=None):
    pts=[(i,j) for i in range(-R,R+1) for j in range(-R,R+1) if abs(i)<=R and abs(j)<=R and abs(i+j)<=R]
    if tint is not None:
        hexv=[np.array([cx,cy])+a*(R+0.62)*np.array([np.cos(np.pi/6+m*np.pi/3),np.sin(np.pi/6+m*np.pi/3)]) for m in range(6)]
        ax.add_patch(Polygon(hexv,fc=hue(k),ec="none",alpha=tint,zorder=1))
    P={q:np.array([cx,cy])+tri_xy(*q,a) for q in pts}
    if bond:
        for q in pts:
            for d in [(1,0),(0,1),(-1,1)]:
                r=(q[0]+d[0],q[1]+d[1])
                if r in P: ax.plot(*np.c_[P[q],P[r]],color="#9aa0a6",lw=0.6,zorder=2)
    for q in pts:
        on=site_state(q[0],q[1],k)
        ax.scatter(*P[q],s=ms,c="#1d2a44" if on else "white",edgecolors="#1d2a44",lw=0.8,zorder=4)
    return P
def wheel_img(ax,ext,cores,lw_dark=True,rc=0.55):
    x=np.linspace(ext[0],ext[1],320); y=np.linspace(ext[2],ext[3],260); X,Y=np.meshgrid(x,y)
    th=np.zeros_like(X); amp=np.ones_like(X)
    for (x0,y0,q) in cores:
        th+=q*np.arctan2(Y-y0,X-x0); amp*=np.tanh(np.hypot(X-x0,Y-y0)/rc)
    rgb=mcol.hsv_to_rgb(np.dstack([(th/(2*np.pi))%1,0.75*np.ones_like(X),0.25+0.75*amp]))
    ax.imshow(rgb,extent=ext,origin="lower",zorder=0)
def defect_tri(ax,x,y,a=0.42,up=True):
    v=np.array([[x-a/2,y-a*SQ3/6],[x+a/2,y-a*SQ3/6],[x,y+a*SQ3/3]]) if up else np.array([[x-a/2,y+a*SQ3/6],[x+a/2,y+a*SQ3/6],[x,y-a*SQ3/3]])
    ax.add_patch(Polygon(v,fc="white",ec="#c0392b",lw=1.4,zorder=5,alpha=0.9))
    ax.scatter(*v.T,s=16,c="#1d2a44",edgecolors="#1d2a44",lw=0.6,zorder=6)
    ax.scatter([x],[y],s=10,c="black",edgecolors="white",lw=0.4,zorder=7)

def fig1():
    fig=plt.figure(figsize=(7.0,4.75))
    ttl=["fixed-bond TAFM","annealed: T = 0","frustrated-triangle pair","defect triangle"]
    sub=["equal-spin bond\nalways active: cost 1","bond switched off (dimer)\ncost 0","dimer released\ncost 1 + 2μ",
         "dimer hides only one bond;\ntwo equal-spin bonds exposed:\nisolated defect: cost 2 (TAFM 1)"]
    for n,(m,t,s_,l) in enumerate(zip(["fixed","ann","mono","defect"],ttl,sub,"abcd")):
        ax=fig.add_axes([0.01+n*0.25,0.53,0.23,0.44])
        rh(ax,0,0,m); ax.set_xlim(-1.25,1.25); ax.set_ylim(-1.95,1.25); ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(t,fontsize=8.5,pad=2); ax.text(0,-1.2,s_,ha="center",va="top",fontsize=7.3,color=INK2,linespacing=1.25)
        ax.text(-1.2,1.18,f"({l})",fontsize=9,fontweight="bold",va="top")
    fig.add_artist(plt.Line2D([0.03,0.97],[0.505,0.505],color=GRID,lw=0.8))
    B=[fig.add_axes([0.01+n*0.25,0.0,0.23,0.47]) for n in range(4)]
    for ax in B: ax.set_aspect("equal"); ax.axis("off")
    # (e) six manifold states = six colors
    ax=B[0]; ax.set_xlim(-2.3,2.3); ax.set_ylim(-3.05,2.45)
    ring=np.linspace(0,2*np.pi,200)
    for k in range(6):
        ang=k*np.pi/3; c=1.55*np.array([np.cos(ang),np.sin(ang)])
        patch(ax,*c,k,a=0.30,R=1,ms=9,tint=0.55)
        ax.plot([0.42*np.cos(ang),0.72*np.cos(ang)],[0.42*np.sin(ang),0.72*np.sin(ang)],color="k",lw=1.6 if k%2==0 else 0.6)
    for k in range(60):
        a0,a1=2*np.pi*k/60,2*np.pi*(k+1)/60
        ax.add_patch(Polygon([(0.38*np.cos(a0),0.38*np.sin(a0)),(0.38*np.cos(a1),0.38*np.sin(a1)),(0,0)],fc=hue(6*k/60),ec="none"))
    ax.text(0,-2.25,"six $\\sqrt{3}\\times\\sqrt{3}$ states\n= six colors of arg$\\,\\psi_{\\rm loc}$;\nthick: the three favored by $D$",ha="center",va="top",fontsize=7.1,color=INK2,linespacing=1.25)
    ax.set_title("states → color",fontsize=8.5,pad=2); ax.text(-2.25,2.4,"(e)",fontsize=9,fontweight="bold",va="top")
    # (f) a defect triangle is a vortex
    ax=B[1]; E=(-1.6,1.6,-1.25,1.25); wheel_img(ax,E,[(0,0,1)]); ax.set_xlim(-2.3,2.3); ax.set_ylim(-3.05,2.45)
    ax.add_patch(plt.Rectangle((E[0],E[2]),E[1]-E[0],E[3]-E[2],fc="none",ec="#9aa0a6",lw=0.5))
    defect_tri(ax,0,0,a=0.5)
    ax.add_patch(FancyArrowPatch((0.95,0.35),(-0.35,0.95),connectionstyle="arc3,rad=0.45",arrowstyle="-|>",mutation_scale=9,color="k",lw=0.9,zorder=6))
    ax.text(0,-1.45,"defect triangle (black dot) = vortex:\nthe six states follow in order,\nphase winds by ±2π; dark core\nwhere |ψ| is small",ha="center",va="top",fontsize=7.1,color=INK2,linespacing=1.25)
    ax.set_title("defect → vortex",fontsize=8.5,pad=2); ax.text(-2.25,2.4,"(f)",fontsize=9,fontweight="bold",va="top")
    # (g) frustrated pair: spins untouched, color unchanged
    ax=B[2]; ax.set_xlim(-2.3,2.3); ax.set_ylim(-3.05,2.45); a=0.52
    KG=4; P=patch(ax,0,0.05,KG,a=a,R=3,ms=13,tint=0.40)
    q=[(0,0),(1,-1),(1,0),(0,1)]
    ax.add_patch(Polygon([P[t] for t in q],fc="#f5b7b1",ec="#c0392b",lw=0.5,zorder=3))
    ax.plot(*np.c_[P[(0,0)],P[(1,0)]],color="#c0392b",lw=2.0,zorder=3)
    for t in q: ax.scatter(*P[t],s=13,c="#1d2a44" if site_state(*t,KG) else "white",edgecolors="#1d2a44",lw=0.8,zorder=5)
    for tt in ([(0,0),(1,0),(0,1)],[(0,0),(1,0),(1,-1)]):
        c=np.mean([P[t] for t in tt],axis=0); ax.scatter(*c,s=9,c="white",edgecolors="#555",lw=0.4,zorder=6)
    ax.text(0,-1.75,"frustrated pair (white dots):\nno spin changes, the color\nis unchanged (no phase charge)",ha="center",va="top",fontsize=7.1,color=INK2,linespacing=1.25)
    ax.set_title("frustrated pair → nothing",fontsize=8.5,pad=2); ax.text(-2.25,2.4,"(g)",fontsize=9,fontweight="bold",va="top")
    # (h) adjacent (bond-sharing) defect pair: opposite vortices cancel
    ax=B[3]; E=(-1.6,1.6,-1.25,1.25); aa=0.62; h=aa*SQ3/6
    wheel_img(ax,E,[(0,h,1),(0,-h,-1)],rc=0.28); ax.set_xlim(-2.3,2.3); ax.set_ylim(-3.05,2.45)
    ax.add_patch(plt.Rectangle((E[0],E[2]),E[1]-E[0],E[3]-E[2],fc="none",ec="#9aa0a6",lw=0.5))
    V=np.array([[-aa/2,0],[aa/2,0],[0,aa*SQ3/2],[0,-aa*SQ3/2]])
    ax.add_patch(Polygon(V[[0,2,1,3]],fc="white",ec="#c0392b",lw=1.3,zorder=5,alpha=0.9))
    ax.plot(*V[[0,1]].T,color="#c0392b",lw=1.3,zorder=5)
    ax.scatter(*V.T,s=16,c="#1d2a44",zorder=6)
    ax.scatter([0,0],[h,-h],s=10,c="black",edgecolors="white",lw=0.4,zorder=7)
    ax.text(0.40,h+0.02,"+",fontsize=9,fontweight="bold",va="center"); ax.text(0.40,-h-0.02,"−",fontsize=9,fontweight="bold",va="center")
    ax.text(0,-1.45,"adjacent pair (shared bond):\nwindings cancel; cost 3 at μ = 0\n(TAFM: 2). Only isolated\ndefects act as free vortices",ha="center",va="top",fontsize=7.1,color=INK2,linespacing=1.25)
    ax.set_title("adjacent pair → neutral",fontsize=8.5,pad=2); ax.text(-2.25,2.4,"(h)",fontsize=9,fontweight="bold",va="top")
    save(fig,"fig1_mechanism")

# ---------- Fig 4: xi and densities ----------
def fig4():
    x=rd("fig3a_xi.csv"); d=rd("fig3b_defects.csv")
    fig,(a,b,cc)=plt.subplots(1,3,figsize=(7.0,2.55),gridspec_kw={"width_ratios":[1,1,0.82]})
    for mdl,c,mk,lb in [("annealed",BLUE,"o","annealed (4 seeds)"),("tafm",ORG,"s","TAFM (1 seed)")]:
        r=[q for q in x if q["model"].lower().startswith(mdl[:4])]
        T=np.array([float(q["T"]) for q in r]); xi=np.array([float(q["xi"]) for q in r]); e=np.array([float(q["err"] or 0) for q in r])
        a.errorbar(1/T,xi,yerr=e,fmt=mk+"-",ms=3,lw=1,color=c,capsize=1.5,label=lb)
    u=np.linspace(1.6,3.0,10); u2=np.linspace(1.5,2.45,10)
    a.plot(u2,0.10*np.exp(2*u2),color=BLUE,ls=":",lw=1); a.text(2.05,0.10*np.exp(2*2.05)*1.9,r"$\propto e^{2/T}$",color=BLUE,fontsize=8,ha="right")
    a.plot(u,0.25*np.exp(u),color=ORG,ls=":",lw=1); a.text(2.75,0.25*np.exp(2.75)*0.52,r"$\propto e^{1/T}$",color=ORG,fontsize=8,va="top")
    a.axhline(16,color=INK2,lw=0.7,ls="--"); a.text(1.02,17.3,"L/2",fontsize=7.5,color=INK2)
    a.set_yscale("log"); a.set_xlabel("1/T"); a.set_ylabel(r"$\xi$  (L = 32, $D=0$)"); a.set_ylim(0.5,30)
    a.legend(loc="lower right",frameon=False); a.grid(color=GRID,lw=0.5); lab(a,"(a)")
    T=np.array([float(q["T"]) for q in d])
    g=rd("defect_classes_L32.csv")
    b.plot(1/T,[float(q["n_mon_annealed"]) for q in d],"^--",ms=3.3,lw=1,color=GRN,label="frustrated, annealed ($L=24$)")
    for mdl,c,mk in [("annealed",BLUE,"o"),("tafm",ORG,"s")]:
        q=sorted([z for z in g if z["model"]==mdl],key=lambda z:float(z["T"])); Tq=np.array([float(z["T"]) for z in q])
        ok=Tq>=0.40
        b.plot(1/Tq[ok],[float(z["non"]) for z,o in zip(q,ok) if o],mk+"-",ms=3.3,lw=1,color=c,label=("isolated defects, "+("annealed" if mdl=="annealed" else "TAFM")))
        b.plot(1/Tq[ok],[float(z["adj"]) for z,o in zip(q,ok) if o],mk+":",ms=3.0,mfc="white",lw=1,color=c,label=("adjacent pairs, "+("annealed" if mdl=="annealed" else "TAFM")))
    b.set_yscale("log"); b.set_xlabel("1/T"); b.set_ylabel("density per triangle")
    b.set_ylim(3e-7,2); b.legend(loc="lower left",frameon=False,fontsize=6.0,handlelength=2.2,borderaxespad=0.2); b.grid(color=GRID,lw=0.5); lab(b,"(b)")
    # (c) xi as a function of the density of isolated defects
    for mdl,c,mk in [("annealed",BLUE,"o"),("tafm",ORG,"s")]:
        r=sorted([q for q in x if q["model"].lower().startswith(mdl[:4])],key=lambda q:float(q["T"]))
        Tx=np.array([float(q["T"]) for q in r]); lx=np.log([float(q["xi"]) for q in r])
        q=sorted([z for z in g if z["model"]==mdl],key=lambda z:float(z["T"])); Tq=np.array([float(z["T"]) for z in q])
        nf=np.array([float(z["non"]) for z in q]); inr=(Tq>=Tx.min()-1e-9)&(Tq<=Tx.max()+1e-9)
        xi_i=np.exp(np.interp(1/Tq,(1/Tx)[::-1],lx[::-1]))
        ok=inr&((Tq>=0.45) if mdl=="annealed" else True)
        cc.plot(nf[ok]**-0.5,xi_i[ok],mk,ms=3.6,color=c,label=("annealed" if mdl=="annealed" else "TAFM"))
        if mdl=="annealed":
            o2=inr&~ok; cc.plot(nf[o2]**-0.5,xi_i[o2],mk,ms=3.6,mfc="none",mec=c)
    u=np.linspace(2,40,10); cc.plot(u,0.33*u,color=INK2,ls="--",lw=0.9); cc.text(16,0.33*16*0.6,r"$0.33\,n_{\rm iso}^{-1/2}$",fontsize=7.5,color=INK2,ha="left",va="top")
    cc.set_xscale("log"); cc.set_yscale("log"); cc.set_xlabel(r"$n_{\rm iso}^{-1/2}$"); cc.set_ylabel(r"$\xi$")
    cc.set_xlim(2.5,110); cc.set_ylim(0.5,30)
    from matplotlib.ticker import FixedLocator,NullLocator,FixedFormatter
    cc.xaxis.set_major_locator(FixedLocator([3,10,30,100])); cc.xaxis.set_major_formatter(FixedFormatter(["3","10","30","100"])); cc.xaxis.set_minor_locator(NullLocator())
    cc.legend(loc="upper left",frameon=False); cc.grid(color=GRID,lw=0.5); lab(cc,"(c)")
    fig.tight_layout(w_pad=1.4); save(fig,"fig4_xi_defects")

# ---------- Fig 5: Binder ----------
def fig5():
    r=rd("fig2_binder.csv")
    fig,axs=plt.subplots(1,2,figsize=(7.0,2.7),sharey=True)
    for ax,mdl,tc,tt,l in [(axs[0],"annealed",0.528,r"annealed bonds ($D=1$, $\mu=0$)","(a)"),(axs[1],"tafm",0.4713,r"fixed-bond TAFM ($h/J=1$)","(b)")]:
        rr=[q for q in r if q["model"].lower().startswith(mdl[:4])]
        for L in sorted({int(q["L"]) for q in rr}):
            q=[z for z in rr if int(z["L"])==L]; q.sort(key=lambda z:float(z["T"]))
            ax.plot([float(z["T"]) for z in q],[float(z["U4"]) for z in q],color=LCOL[L],lw=1.4,label=f"L = {L}")
        ax.axvline(tc,color=INK2,lw=0.7,ls=":"); ax.text(tc-0.0005*(1 if tc<0.5 else 2),0.03,f"$T_c$ = {tc} ",fontsize=7.5,color=INK2,ha="right",va="bottom")
        ax.set_title(tt,pad=3); ax.set_xlabel("T"); ax.grid(color=GRID,lw=0.5); ax.legend(loc="center left",frameon=False); lab(ax,l,-0.12)
    axs[0].set_ylabel(r"Binder cumulant $U_4$"); axs[0].set_ylim(0,0.53)
    fig.tight_layout(w_pad=1.5); save(fig,"fig5_binder")

# ---------- Fig 6: heat capacity ----------
def fig6():
    r=rd("fig4_heat.csv")
    fig,ax=plt.subplots(figsize=(3.4,2.6))
    a=[q for q in r if q["model"].lower().startswith("anne")]; t=[q for q in r if not q["model"].lower().startswith("anne")]
    ax.plot([float(q["T_over_Tc"]) for q in a],[float(q["C"]) for q in a],color=BLUE,lw=1.6,label=r"annealed $C_\mu$")
    ax.plot([float(q["T_over_Tc"]) for q in a],[float(q["C_n"]) for q in a],color=BLUE,lw=1.4,ls="--",label=r"annealed $C_n$ (fixed $n$)")
    ax.plot([float(q["T_over_Tc"]) for q in t],[float(q["C"]) for q in t],color=ORG,lw=1.6,label="fixed-bond TAFM")
    ax.set_xlabel(r"$T/T_c$"); ax.set_ylabel(r"$C/N_{\rm tri}$  (L = 64)"); ax.set_ylim(0.3,1.62)
    ax.legend(loc="upper left",frameon=False,ncol=1); ax.grid(color=GRID,lw=0.5)
    fig.tight_layout(); save(fig,"fig6_heat")

# ---------- Fig 7: Tc(mu) ----------
def fig7():
    r=rd("fig5_tc_mu.csv"); mu=[float(q["mu"]) for q in r]; tc=[float(q["T_c"]) for q in r]; er=[float(q["err"]) for q in r]
    fig,ax=plt.subplots(figsize=(3.4,2.7))
    ax.plot(mu,tc,"-",color=BLUE,lw=1.3,zorder=2)
    for m,t,e in zip(mu,tc,er):
        lo,hi,f=(e,1.16-t,False) if m==10 else (e,e,e<0.02)
        ax.errorbar([m],[t],yerr=[[lo],[hi]],fmt="o",ms=4.5,color=BLUE,mfc=BLUE if f else "white",mew=1.3,capsize=2.5,zorder=3)
    ax.axhline(1.2027,color=INK2,lw=1,ls=":"); ax.text(0.2,1.215,"dice-lattice Ising, 1.2027",fontsize=7,color=INK2,va="bottom")
    ax.axhline(0.4713,color=ORG,lw=1.2,ls="--"); ax.text(10.2,0.485,r"fixed-bond TAFM, $h/J=1$",fontsize=7,color=INK2,ha="right",va="bottom")
    ax.errorbar([],[],fmt="o",ms=4.5,color=BLUE,label="Binder crossing"); ax.errorbar([],[],fmt="o",ms=4.5,color=BLUE,mfc="white",mew=1.3,label="approximate")
    ax.legend(loc="center right",bbox_to_anchor=(1.0,0.36),frameon=False)
    ax.set_xlabel(r"$\mu$ (cost of a frustrated triangle)"); ax.set_ylabel(r"$T_c$  ($D=1$)"); ax.set_ylim(0.4,1.3); ax.set_xlim(-0.4,10.6)
    ax.grid(color=GRID,lw=0.5); fig.tight_layout(); save(fig,"fig7_tc_mu")

# ---------- Fig 8: eta_eff and cos 3theta ----------
def fig8():
    r=rd("fig6a_eta.csv"); c=rd("fig6b_cos3.csv")
    fig,(a,b)=plt.subplots(1,2,figsize=(7.0,2.8),gridspec_kw={"width_ratios":[1.35,1]})
    cols={1:"#86b6ef",2:"#5598e7",3:"#2a78d6",5:"#1c5cab",10:"#0d366b"}; mk={1:"o",2:"s",3:"^",5:"D",10:"v"}
    for m in [1,2,3,5,10]:
        q=[z for z in r if int(float(z["mu"]))==m]
        a.errorbar([float(z["T_over_Tc"]) for z in q],[float(z["eta_eff"]) for z in q],yerr=[float(z["err"]) for z in q],
                   fmt=mk[m]+"-",ms=3.5,lw=1,color=cols[m],capsize=1.5,label=f"μ = {m}")
    for v,s in [(4/15,"4/15 (Potts)"),(4/9,"4/9 (Z$_3$ locking)")]:
        a.axhline(v,color=INK2,lw=0.7,ls=":"); a.text(0.604,v+0.01,s,fontsize=6.6,color=INK2,va="bottom")
    a.set_xlim(0.6,1.06); a.set_ylim(-0.08,1.0); a.set_xlabel(r"$T/T_c$"); a.set_ylabel(r"effective $\eta$  ($\psi_s^2\propto L^{-\eta}$)")
    a.legend(loc="upper left",bbox_to_anchor=(0.0,1.0),frameon=False,ncol=3,columnspacing=1.0,handlelength=1.6); a.grid(color=GRID,lw=0.5); lab(a,"(a)",-0.12)
    tc=1.141
    for T,cl,mkk in [(0.7,BLUE,"o"),(1.04,ORG,"s"),(1.08,GRN,"^"),(1.14,YEL,"D")]:
        q=[z for z in c if abs(float(z["T"])-T)<1e-6]
        b.plot([int(z["L"]) for z in q],[float(z["cos3theta"]) for z in q],mkk+"-",ms=4,lw=1,color=cl,label=f"T/T$_c$ = {T/tc:.2f}")
    b.set_xscale("log"); b.set_xticks([32,48,64,96]); b.set_xticklabels(["32","48","64","96"]); b.minorticks_off()
    b.set_xlabel("L"); b.set_ylabel(r"$\langle\cos 3\theta\rangle$  (μ = 10)"); b.set_ylim(0,1.12)
    b.legend(loc="upper center",bbox_to_anchor=(0.5,0.84),frameon=False,ncol=2); b.grid(color=GRID,lw=0.5); lab(b,"(b)",-0.16)
    fig.tight_layout(w_pad=2.0); save(fig,"fig8_eta_cos3")

# ---------- Fig 9: plateau window: U4 vs T/Tc, mu=0 vs mu=10 ----------
def fig9():
    r=rd("fig2_binder.csv")
    fig,(a,b)=plt.subplots(1,2,figsize=(7.0,2.7),sharey=True)
    rr=[q for q in r if q["model"].lower().startswith("anne")]
    for L in sorted({int(q["L"]) for q in rr}):
        q=sorted([z for z in rr if int(z["L"])==L],key=lambda z:float(z["T"]))
        a.plot([float(z["T"])/0.528 for z in q],[float(z["U4"]) for z in q],color=LCOL[L],lw=1.4,label=f"L = {L}")
    a.set_title(r"$\mu=0$ (hub order parameter)",pad=3)
    T=[1.02,1.04,1.06,1.08,1.10,1.12,1.14,1.16,1.18,1.20]
    tab={24:[.464,.459,.459,.458,.447,.427,.384,.290,.147,.069],32:[.458,.460,.460,.456,.449,.435,.387,.237,.057,-.035],
         48:[.459,.461,.464,.460,.456,.450,.391,.137,-.012,-.077],64:[.460,.465,.460,.463,.461,.455,.3905,-.026,-.004,-.017],
         96:[None,.462,None,.455,None,.457,.418,None,None,None]}
    err={24:[2,1,1,1,3,4,8,11,14,21],32:[2,2,1,3,6,6,12,18,31,18],48:[3,4,3,2,4,4,6,20,31,27],64:[5,3,5,3,2,2,22,35,46,22],96:[0,3,0,8,0,6,9,0,0,0]}
    for L in [24,32,48,64,96]:
        xs=[t/1.141 for t,v in zip(T,tab[L]) if v is not None]; ys=[v for v in tab[L] if v is not None]
        es=[e/1000 for e,v in zip(err[L],tab[L]) if v is not None]
        b.errorbar(xs,ys,yerr=es,fmt="o-",ms=3,lw=1,color=LCOL[L],capsize=1.5,label=f"L = {L}")
    b.set_title(r"$\mu=10$ (spin order parameter)",pad=3)
    for ax,l in [(a,"(a)"),(b,"(b)")]:
        ax.axhline(0.5,color=INK2,lw=0.6,ls=":"); ax.axvspan(0.89,0.95,color="#dfe9f7",zorder=0)
        ax.set_xlim(0.83,1.06); ax.set_xlabel(r"$T/T_c$"); ax.grid(color=GRID,lw=0.5); lab(ax,l,-0.12)
        ax.legend(loc="lower left",frameon=False)
    a.set_ylabel(r"$U_4$"); a.set_ylim(-0.1,0.56)
    fig.tight_layout(w_pad=1.5); save(fig,"fig9_window")
for f in (fig1,fig4,fig5,fig6,fig7,fig8,fig9): f()
print("ok")
