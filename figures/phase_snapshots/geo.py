import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
import numpy as np, time
from worm_monomer_walk import DimerState
from chi_basic import setup
t=time.time(); L=48; st=DimerState(L,L); pos,_=st.lat.site_positions(); tri=np.array([t_["sites"] for t_ in st.tri])
A,nbr=setup(L)
d=np.linalg.norm(pos[nbr[0]]-pos[0],axis=1); print(pos.shape,tri.shape,A["nt"],d, time.time()-t)
print(pos.min(0),pos.max(0))
np.savez(DATA/"snapshots"/"geo48.npz",pos=pos,tri=tri,nbr=nbr,sub=A["sub"])
