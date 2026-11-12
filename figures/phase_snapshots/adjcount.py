import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
exec(open(_pl.Path(__file__).with_name("plot.py")).read())
from emin import tri_arrays
from chi_basic import setup
A,_=setup(48); TB=tri_arrays(A); nt=A["nt"]
tadj=np.full((nt,3),-1)
for t in range(nt):
    for k,b in enumerate(TB[t]): tadj[t,k]=A["btb"][b] if A["bta"][b]==t else A["bta"][b]
for f in ("ann_T050.npz","tafm_T040.npz"):
    s=np.load(SNAP/f)["s"]; eq=(s[A["bi"]]==s[A["bj"]]); d=eq[TB].all(1); nb=d[tadj].any(1)&d
    print(f,"defects",d.sum(),"adjacent",nb.sum(),"non-adjacent",(d&~nb).sum())
