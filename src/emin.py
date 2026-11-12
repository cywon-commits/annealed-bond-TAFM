import sys
import numpy as np, networkx as nx
from chi_basic import setup, exact_gs
def tri_arrays(A):
    nt=A["nt"]; tb=[[] for _ in range(nt)]
    for b in range(len(A["bi"])):
        tb[A["bta"][b]].append(b); tb[A["btb"][b]].append(b)
    return np.array(tb)
def emin(A,s,mu=0.0,return_M=False):
    """min over matchings of E = sum_active eq + mu*n_mon for spins s"""
    eq=(s[A["bi"]]==s[A["bj"]]).astype(int)
    G=nx.Graph()
    for b in range(len(eq)):
        w=eq[b]+2*mu
        if w>0: G.add_edge(int(A["bta"][b]),int(A["btb"][b]),weight=w,b=b)
    M=nx.max_weight_matching(G,maxcardinality=False)
    gain=sum(G[u][v]["weight"] for u,v in M)
    E=eq.sum()+mu*A["nt"]-gain
    return (E,M,eq) if return_M else E
def defects(A,s,TB):
    eq=(s[A["bi"]]==s[A["bj"]])
    return np.where(eq[TB].all(1))[0]
