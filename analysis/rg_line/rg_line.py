import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
"""Step 2: annealed critical line from the TAFM critical surface by matching vortex fugacities.
y ∝ sqrt(n_d) (xi ∝ 1/y, n_d ∝ xi^-2, same constant in both models: Fig. 4c).
Annealed (T, H) <-> TAFM (T', H) with n_d^TAFM(T') = n_d^ann(T)."""
import numpy as np, pandas as pd
from scipy.optimize import least_squares, brentq
# free-vortex proxy: density of defect triangles that share no bond with another defect (L=32, 4 seeds)
g=pd.read_csv(FIGDATA/"defect_classes_L32.csv")
ga=g[g.model=="annealed"]; gt=g[g.model=="tafm"]
sa=(ga["T"].values>=0.40)&(ga["T"].values<=0.80)   # T=0.36 not converged between seeds
qa=np.polyfit(ga["T"].values[sa],np.log(ga["non"].values[sa])+4/ga["T"].values[sa],2)
qt=np.polyfit(gt["T"].values,np.log(gt["non"].values)+2/gt["T"].values,2)
def ln_nd_ann(t): return np.polyval(qa,t)-4/t
def ln_nd_taf(t): return np.polyval(qt,t)-2/t
print("fit residuals ann",np.round(ln_nd_ann(ga["T"].values[sa])-np.log(ga["non"].values[sa]),3)," tafm",np.round(ln_nd_taf(gt["T"].values)-np.log(gt["non"].values),3))
def Tprime(t):  # TAFM temperature with the same defect density
    target=ln_nd_ann(t); return brentq(lambda u: ln_nd_taf(u)-target,0.003,1.0)
for t in (0.6,0.528,0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15): print(f"  T={t:.3f} -> T'={Tprime(t):.3f}  T'/T={Tprime(t)/t:.3f}")
# TAFM critical line (H, T'), model units T=T_J/2, H=D/(2T)
pts=np.array([[0.41,0.125],[0.747,0.3347],[1.057,0.4729],[2.233,0.6718]])  # Nyckees23, this work(D=.5), Qian04 h/J=1, this work D=3
res={}
for HKT in (0.266,0.305):
    # Nyckees Eq.28 to 2nd order: T' = c1 sqrt(dH) + c2 dH ; fit to the three low points
    P=pts[:3]; dH=P[:,0]-HKT
    A=np.c_[np.sqrt(dH),dH]; c,*_=np.linalg.lstsq(A,P[:,1],rcond=None)
    Tline=lambda H: c[0]*np.sqrt(max(H-HKT,0))+c[1]*(H-HKT)
    Hline=lambda Tp: brentq(lambda H: Tline(H)-Tp,HKT+1e-12,5)
    print(f"H_KT={HKT}: T'(H)= {c[0]:.4f} sqrt(dH) + {c[1]:.4f} dH ; residuals",np.round(A@c-P[:,1],4))
    rows=[]
    for D in (1.0,0.75,0.5,0.4,0.3,0.2,0.15,0.1,0.05):
        Ttaf=brentq(lambda t: D/(2*t)-Hline(t),0.02,0.9)
        Tann=brentq(lambda t: D/(2*t)-Hline(Tprime(t)),0.03,0.7)
        rows.append((D,Ttaf,Tann,Tann/Ttaf,D/(2*Ttaf),D/(2*Tann),D/(2*HKT),Tprime(Tann)))
    res[HKT]=rows
    print("    D   Tc_TAFM  Tc_ann  ratio  H_TAFM  H_ann  D/2HKT  T'(Tc_ann)")
    for r in rows: print("  %4.2f  %.4f  %.4f  %.3f  %.3f  %.3f  %.3f  %.3f"%r)
np.save(OUT/"rg_line_result.npy",res,allow_pickle=True)

# sensitivity: annealed fugacity off by a factor (1+eps) i.e. ln n_d shifted by 2 ln(1+eps)
HKT=0.266; P=pts[:3]; dH=P[:,0]-HKT; A=np.c_[np.sqrt(dH),dH]; c,*_=np.linalg.lstsq(A,P[:,1],rcond=None)
Tline=lambda H: c[0]*np.sqrt(max(H-HKT,0))+c[1]*(H-HKT); Hline=lambda Tp: brentq(lambda H: Tline(H)-Tp,HKT+1e-12,5)
def Tp_shift(t,s):
    target=ln_nd_ann(t)+s; return brentq(lambda u: ln_nd_taf(u)-target,0.003,1.0)
print("sensitivity (y_ann x 1.2 / x 1/1.2):")
for D in (1.0,0.5,0.2,0.1):
    out=[brentq(lambda t: D/(2*t)-Hline(Tp_shift(t,s)),0.03,0.7) for s in (2*np.log(1.2),0,-2*np.log(1.2))]
    print("  D=%.2f  Tc_ann = %.4f / %.4f / %.4f"%(D,*out))
