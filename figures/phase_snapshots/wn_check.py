import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
exec(open(_pl.Path(__file__).with_name("plot.py")).read())
CL1=clusters(1); dd=lambda x:(x+np.pi)%(2*np.pi)-np.pi
for f in ("ann_T050.npz","tafm_T040.npz"):
    s=np.load(SNAP/f)["s"]; a=(2*s-1)*W; ph=np.angle(np.array([a[c].sum() for c in CL1]))
    wn=np.round((dd(ph[tri[:,1]]-ph[tri[:,0]])+dd(ph[tri[:,2]]-ph[tri[:,1]])+dd(ph[tri[:,0]]-ph[tri[:,2]]))/(2*np.pi)).astype(int)
    print(f,np.unique(wn,return_counts=True))
