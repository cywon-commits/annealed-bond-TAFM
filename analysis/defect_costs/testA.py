import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
from emin import *
L=int(sys.argv[1]); A,nbr=setup(L); TB=tri_arrays(A)
s0,dim,mate=exact_gs(A,3)
print("GS: E_min(mu=0)=",emin(A,s0), " defects",len(defects(A,s0,TB)))
from collections import Counter
res=Counter()
for v in range(A["ns"]):
    s=s0.copy(); s[v]=1-s[v]; d=defects(A,s,TB)
    eq=(s[A["bi"]]==s[A["bj"]]); dE_tafm=int(eq.sum()-(s0[A["bi"]]==s0[A["bj"]]).sum())
    if len(d)==0: continue
    # adjacency of defects (share a bond)
    shared=len(set(TB[d[0]])&set(TB[d[1]])) if len(d)==2 else -1
    res[(len(d),shared,dE_tafm,emin(A,s,0.0),round(emin(A,s,0.5),2))]+=1
for k,v in sorted(res.items()): print("ndef=%d share=%d dE_TAFM=%d Emin(mu0)=%d Emin(mu0.5)=%s  count=%d"%(*k,v))
