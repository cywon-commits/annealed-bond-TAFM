import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
exec(open(_pl.Path(__file__).with_name("plot.py")).read())
CL1=clusters(1)
def vort(s,mate):
    a=(2*s-1)*W; psi=np.array([a[c].sum() for c in CL1]); ph=np.angle(psi)
    # elementary site triangles = the lattice triangles 'tri'
    dd=lambda x:(x+np.pi)%(2*np.pi)-np.pi
    wn=np.round((dd(ph[tri[:,1]]-ph[tri[:,0]])+dd(ph[tri[:,2]]-ph[tri[:,1]])+dd(ph[tri[:,0]]-ph[tri[:,2]]))/(2*np.pi)).astype(int)
    ss=s[tri]; dfe=(ss[:,0]==ss[:,1])&(ss[:,1]==ss[:,2])
    c=pos[tri].mean(1)
    # distance (torus) from each vortex triangle to nearest defect, in lattice-site spacings
    def tor(p,q):
        d=(A_MAT_INV@(p-q).T).T; d-=np.round(d/L)*L; return np.linalg.norm((A_MAT@d.T).T,axis=1)
    V=np.where(wn!=0)[0]; Dd=np.where(dfe)[0]
    near=[tor(c[Dd],np.repeat(c[[v]],len(Dd),0)).min()/0.577 if len(Dd) else np.nan for v in V]
    print(f"  vortex triangles {len(V)} (net {wn.sum()}), defects {len(Dd)}, frustrated {(mate<0).sum()}; vortices within 1.5 spacings of a defect: {np.mean(np.array(near)<1.5):.2f}")
    fr=np.where((mate<0)&~dfe)[0]
    if len(fr):
        nf=[tor(c[V],np.repeat(c[[f]],len(V),0)).min()/0.577 for f in fr[:300]]
        print(f"  frustrated triangles within 1.5 spacings of a vortex: {np.mean(np.array(nf)<1.5):.2f}")
for f in ("ann_T050.npz","tafm_T040.npz"):
    z=np.load(SNAP/f); print(f); vort(z["s"],z["mate"])
# random reference: fraction of random points within 1.5 of a vortex
z=np.load(SNAP/"ann_T050.npz"); s=z["s"]; a=(2*s-1)*W; ph=np.angle(np.array([a[c].sum() for c in CL1]))
dd=lambda x:(x+np.pi)%(2*np.pi)-np.pi
wn=np.round((dd(ph[tri[:,1]]-ph[tri[:,0]])+dd(ph[tri[:,2]]-ph[tri[:,1]])+dd(ph[tri[:,0]]-ph[tri[:,2]]))/(2*np.pi))
V=np.where(wn!=0)[0]; c=pos[tri].mean(1); rng=np.random.default_rng(0); samp=rng.choice(len(tri),400,replace=False)
def tor(p,q):
    d=(A_MAT_INV@(p-q).T).T; d-=np.round(d/L)*L; return np.linalg.norm((A_MAT@d.T).T,axis=1)
print("baseline random triangles within 1.5 of vortex:",np.mean([tor(c[V],np.repeat(c[[t]],len(V),0)).min()/0.577<1.5 for t in samp]))
