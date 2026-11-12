import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
from emin import *
from collections import defaultdict
L=6; A,nbr=setup(L); TB=tri_arrays(A); s0,_,_=exact_gs(A,3)
res=defaultdict(set)
for v in range(A["ns"]):
    s=s0.copy(); s[v]=1-s[v]; d=defects(A,s,TB)
    if len(d)==0: continue
    key=(len(d), len(set(TB[d[0]])&set(TB[d[1]])) if len(d)==2 else -1)
    for mu in (0,0.25,0.5,1.0,3.0): res[key+(mu,)].add(round(emin(A,s,mu),3))
for k in sorted(res): print(k,sorted(res[k]))
