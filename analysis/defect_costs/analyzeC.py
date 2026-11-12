import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / "src"))
from repo_paths import ROOT, DATA, FIGDATA, FIGS, OUT
import numpy as np, re, pandas as pd, glob
rows=[]
for f in glob.glob(str(DATA/"defect_classes"/"C_*.log")):
    for l in open(f):
        m=re.match(r"(\w+) seed=(\d+) T=([\d.]+) n_d=([\d.e+-]+) adj=([\d.e+-]+) nonadj=([\d.e+-]+)",l)
        if m: rows.append([m.group(1),int(m.group(2))]+[float(x) for x in m.groups()[2:]])
df=pd.DataFrame(rows,columns=["model","seed","T","nd","adj","non"])
g=df.groupby(["model","T"]).agg(nd=("nd","mean"),nd_e=("nd","sem"),adj=("adj","mean"),adj_e=("adj","sem"),non=("non","mean"),non_e=("non","sem"),ns=("seed","count")).reset_index()
g.to_csv(FIGDATA/"defect_classes_L32.csv",index=False)
pd.set_option("display.width",200); print(g.to_string())
for m in ("annealed","tafm"):
    q=g[(g.model==m)&(g["T"]>=0.40)&(g["T"]<=0.60)]
    for k in ("nd","adj","non"):
        p=np.polyfit(1/q["T"],np.log(q[k]),1); print(m,k,"slope %.2f"%-p[0])
xi=pd.read_csv(FIGDATA/"fig3a_xi.csv")
X={}
for m in ("annealed","tafm"):
    x=xi[xi.model==m].sort_values("T"); q=g[g.model==m]
    X[m]=np.exp(np.interp(1/q["T"].values,(1/x["T"].values)[::-1],np.log(x["xi"].values)[::-1]))
a=g[g.model=="annealed"].reset_index(drop=True); t=g[g.model=="tafm"].reset_index(drop=True)
for k in ("nd","non"):
    o=np.argsort(np.log(t[k].values))
    xt=np.exp(np.interp(np.log(a[k].values),np.log(t[k].values)[o],np.log(X["tafm"])[o]))
    inr=(a[k].values>=t[k].min())
    print(k,"ratio xi_ann/xi_TAFM at equal n:",[f"T={T:.2f}:{r:.3f}{'' if ok else '(extrap)'}" for T,r,ok in zip(a["T"],X["annealed"]/xt,inr)])
    for m,q,Xm in (("annealed",a,X["annealed"]),("tafm",t,X["tafm"])):
        sel=q["T"]>=0.4; print("   xi*sqrt(n)",m,np.round(Xm[sel]*np.sqrt(q[k][sel]),3))
# amplitudes for RG text
for m,E,q in (("tafm",2,t),("annealed",4,a)):
    sel=(q["T"]>=0.4)&(q["T"]<=0.6); print(m,"n_f e^{%d/T}"%E,np.round(q["non"][sel]*np.exp(E/q["T"][sel]),2).values)
