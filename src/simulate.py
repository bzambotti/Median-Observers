import numpy as np

def Sn(k):
    S=np.zeros((k,k)); S[1:,: -1]=np.eye(k-1); return S

def sim_median(A,C,L,q,dseq,e0hist,T):
    N=2*q+1; e=list(e0hist); p = C.shape[0]
    d=[np.zeros(p) for _ in range(N)]; out=[]
    for k in range(T):
        d=[dseq[k]]+d[:-1]
        buf=np.array([[C[j]@e[i]+d[i][j] for i in range(N)] for j in range(p)])
        M=np.median(buf,axis=1)
        e=[A@e[0]-L@M]+e[:-1]; out.append(np.linalg.norm(e[0]))
    return np.array(out)

def sim_stubborn(A,C,L,nb,r,dseq,e0,T):
    S=Sn(nb); G=np.zeros((nb,C.shape[0])); G[0,:]=1.0
    e=e0.copy(); s=np.zeros(nb); out=[]
    for k in range(T):
        inn=C@e+dseq[k]; Delta=r*np.sum(s)*np.ones(C.shape[0])
        sat=np.clip(inn,-Delta,Delta)
        e=A@e-L@sat; s=S@s+G@np.abs(inn); out.append(np.linalg.norm(e))
    return np.array(out)

def burst(T,b,P,amp,rng,start=6,p=2):
    d=[np.zeros(p) for _ in range(T)]
    for k in range(start,T):
        if (k-start)%P < b: d[k]=amp*rng.standard_normal(p)
    return d
