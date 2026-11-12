import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
SNAP = DATA / "snapshots"
import sys, time
import numpy as np
from annealed_thermo import build_arrays, total_energy, _seed
from ising_sweep import sweep_ising
from chi_basic import setup, exact_gs, tafm_sweeps
L=48; A,nbr=setup(L); sub=A["sub"]; w=np.exp(2j*np.pi/3); nsub=np.array([(sub==k).sum() for k in range(3)])
def psi_s(s):
    f=np.array([((s==1)&(sub==k)).sum() for k in range(3)])/nsub; return (f*w**np.arange(3)).sum()
def ann(T,D,mu,n,seed,start):
    _seed(seed); np.random.seed(seed); rng=np.random.default_rng(seed)
    if start=="gs": s,dim,mate=exact_gs(A,seed+100)
    else:
        dim,mate=A["dimer"].copy(),A["mate"].copy()
        s=np.where(sub==0,1,0).astype(np.int64) if start=="crystal" else rng.integers(2,size=A["ns"]).astype(np.int64)
    E=total_energy(s,dim,A["bi"],A["bj"],D,0.0,mu,mate); P=[]
    for k in range(n):
        E=sweep_ising(s,dim,mate,A["bi"],A["bj"],A["bta"],A["btb"],A["site_b"],A["hexb"],1/T,D,mu,E)
        if k>=n//2 and k%20==0: P.append(abs(psi_s(s))**2)
    return s,dim,mate,np.array(P)
job=sys.argv[1]; t0=time.time()
if job=="tafm":
    s,dim,mate=exact_gs(A,7); tafm_sweeps(s,nbr,1/0.40,20000,20001,7)
    np.savez(SNAP/"tafm_T040.npz",s=s,dim=np.zeros_like(A["dimer"]),mate=np.zeros(A["nt"],np.int64)); print("tafm",time.time()-t0)
elif job=="annnd":
    s,dim,mate,P=ann(0.50,0.0,0.0,20000,8,"gs"); np.savez(SNAP/"ann_T050.npz",s=s,dim=dim,mate=mate); print("ann",time.time()-t0)
else:
    T,D,mu,n,start=float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),int(sys.argv[5]),sys.argv[6]
    s,dim,mate,P=ann(T,D,mu,n,int(1000*T)+int(mu),start)
    np.savez(SNAP/f"win_mu{mu:g}_T{T:g}_{start}.npz",s=s,dim=dim,mate=mate,P=P)
    print(f"mu={mu} T={T} {start}: <psi^2>={P.mean():.4f} +- {P.std()/np.sqrt(len(P)/25):.4f}  last={P[-1]:.4f}  [{time.time()-t0:.0f}s]")
