import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
"""Defect classes vs T: bond-sharing (adjacent) vs non-adjacent; L=32, given seeds."""
from emin import *
import time
from annealed_thermo import total_energy, _seed
from ising_sweep import sweep_ising
from chi_basic import tafm_sweeps
L=32; A,nbr=setup(L); TB=tri_arrays(A); nt=A["nt"]
tadj=np.full((nt,3),-1)
for t in range(nt):
    for k,b in enumerate(TB[t]): tadj[t,k]=A["btb"][b] if A["bta"][b]==t else A["bta"][b]
def classify(s):
    eq=(s[A["bi"]]==s[A["bj"]]); d=eq[TB].all(1); nb=d[tadj].any(1)&d
    return d.mean(), nb.mean(), (d&~nb).mean()
model=sys.argv[1]; seeds=[int(x) for x in sys.argv[2].split(",")]
Ts=[0.36,0.40,0.45,0.50,0.55,0.60,0.70,0.80,1.0] if model=="annealed" else [0.30,0.33,0.36,0.40,0.45,0.50,0.55,0.60,0.70,0.80,1.0]
for sd in seeds:
  for T in Ts:
    t0=time.time(); seed=1000*sd+int(1000*T); s,dim,mate=exact_gs(A,seed); acc=[]
    if model=="tafm":
        tafm_sweeps(s,nbr,1/T,10000,10001,seed)
        for k in range(1000): tafm_sweeps(s,nbr,1/T,20,21,seed+k+1); acc.append(classify(s))
    else:
        _seed(seed); np.random.seed(seed); E=total_energy(s,dim,A["bi"],A["bj"],0.0,0.0,0.0,mate)
        for _ in range(10000): E=sweep_ising(s,dim,mate,A["bi"],A["bj"],A["bta"],A["btb"],A["site_b"],A["hexb"],1/T,0.0,0.0,E)
        for k in range(20000):
            E=sweep_ising(s,dim,mate,A["bi"],A["bj"],A["bta"],A["btb"],A["site_b"],A["hexb"],1/T,0.0,0.0,E)
            if k%20==0: acc.append(classify(s))
    m=np.array(acc).mean(0)
    print(f"{model} seed={sd} T={T:.3f} n_d={m[0]:.5e} adj={m[1]:.5e} nonadj={m[2]:.5e} [{time.time()-t0:.0f}s]",flush=True)
