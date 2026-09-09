import numpy as np, itertools, cvxpy as cp
np.set_printoptions(precision=4, suppress=True, linewidth=150)

# ---------------- model builders ----------------
def build(A, C, q):
    n = A.shape[0]; p = C.shape[0]; N = 2*q+1; Nn = N*n
    Abar = np.zeros((Nn, Nn)); Abar[:n, :n] = A
    for i in range(1, N): Abar[i*n:(i+1)*n, (i-1)*n:i*n] = np.eye(n)
    Bbar = np.zeros((Nn, n)); Bbar[:n, :] = np.eye(n)
    Cbar = np.zeros((p, Nn)); Cbar[:, :n] = C
    # Co  (2q p x Nn) : rows [C  0..  -C(at block i) .. 0]
    Co = np.zeros((2*q*p, Nn))
    for i in range(1, N):
        Co[(i-1)*p:i*p, :n] = C
        Co[(i-1)*p:i*p, i*n:(i+1)*n] = -C
    Qo = Co.T @ Co
    # modes
    Hs = []
    for idx in itertools.product(range(N), repeat=p):
        Cm = np.zeros((p, Nn))
        for j, i in enumerate(idx): Cm[j, i*n:(i+1)*n] = C[j, :]
        Hs.append(Cbar - Cm)          # phi_k = (Cbar - Cbar_m) E_k  <-- correct sign
    return dict(n=n, p=p, q=q, N=N, Nn=Nn, A=A, C=C, Abar=Abar, Bbar=Bbar,
                Cbar=Cbar, Co=Co, Qo=Qo, Hs=Hs)

def AoBo(S, L):
    return S['Abar'] - S['Bbar'] @ L @ S['Cbar'], S['Bbar'] @ L

SOLVE = dict(solver=cp.CLARABEL)

def feas(cons, eps_ok=True):
    try:
        prob = cp.Problem(cp.Minimize(0), cons); prob.solve(**SOLVE)
        return prob.status in ("optimal", "optimal_inaccurate")
    except Exception as e:
        return False

# ---------------- REDESIGN (given L) ----------------
def lure_redesign(S, L, eps=1e-6):
    Ao, Bo = AoBo(S, L); Nn, p = S['Nn'], S['p']
    P = cp.Variable((Nn, Nn), symmetric=True)
    Om = cp.bmat([[Ao.T@P@Ao - P + S['Qo'], Ao.T@P@Bo],
                  [Bo.T@P@Ao,               Bo.T@P@Bo - np.eye(p)]])   # tau = 1 (wlog)
    return feas([P >> eps*np.eye(Nn), Om << -eps*np.eye(Nn+p)])

def switched_redesign(S, L, eps=1e-6):
    Ao, Bo = AoBo(S, L); Nn = S['Nn']
    P = cp.Variable((Nn, Nn), symmetric=True)
    cons = [P >> eps*np.eye(Nn)]
    for H in S['Hs']:
        M = Ao + Bo@H
        cons.append(M.T@P@M - P << -eps*np.eye(Nn))
    return feas(cons)

# ---------------- CODESIGN ----------------
def G_struct(Nn, n):
    G1 = cp.Variable((n, n)); G2 = cp.Variable((n, Nn-n)); G3 = cp.Variable((Nn-n, Nn-n))
    G = cp.bmat([[G1, G2], [np.zeros((Nn-n, n)), G3]])
    return G, G1

def lure_codesign(S, eps=1e-6):
    Nn, n, p = S['Nn'], S['n'], S['p']
    P = cp.Variable((Nn, Nn), symmetric=True); G, G1 = G_struct(Nn, n)
    Y = cp.Variable((n, p)); Ab, Bb, Cb = S['Abar'], S['Bbar'], S['Cbar']
    Z = np.zeros
    M = cp.bmat([
        [-P + S['Qo'],       Z((Nn, p)),        Ab.T@G.T - Cb.T@Y.T@Bb.T],
        [Z((p, Nn)),         -np.eye(p),        Y.T@Bb.T],
        [G@Ab - Bb@Y@Cb,     Bb@Y,              P - G - G.T]])
    cons = [P >> eps*np.eye(Nn), M << -eps*np.eye(2*Nn+p)]
    try:
        prob = cp.Problem(cp.Minimize(0), cons); prob.solve(**SOLVE)
        if prob.status not in ("optimal", "optimal_inaccurate"): return False, None
        return True, np.linalg.solve(G1.value, Y.value)
    except Exception:
        return False, None

def switched_codesign(S, scheduled=True, eps=1e-6):
    """scheduled=True -> Y_m mode dependent ; False -> single static Y."""
    Nn, n, p = S['Nn'], S['n'], S['p']
    P = cp.Variable((Nn, Nn), symmetric=True); G, G1 = G_struct(Nn, n)
    Ab, Bb, Cb = S['Abar'], S['Bbar'], S['Cbar']
    cons = [P >> eps*np.eye(Nn)]
    Ys = []
    Y0 = cp.Variable((n, p))
    for H in S['Hs']:
        Y = cp.Variable((n, p)) if scheduled else Y0
        Ys.append(Y); Cm = Cb - H
        M = cp.bmat([[-P, Ab.T@G.T - Cm.T@Y.T@Bb.T],
                     [G@Ab - Bb@Y@Cm, P - G - G.T]])
        cons.append(M << -eps*np.eye(2*Nn))
    try:
        prob = cp.Problem(cp.Minimize(0), cons); prob.solve(**SOLVE)
        if prob.status not in ("optimal", "optimal_inaccurate"): return False, None
        Ls = [np.linalg.solve(G1.value, Y.value) for Y in Ys]
        return True, Ls
    except Exception:
        return False, None

def switched_codesign_pairs(S, eps=1e-6):
    """robustified scheduled design: LMI for every pair (gain mode m, actual output mode m')."""
    Nn, n, p = S['Nn'], S['n'], S['p']
    P = cp.Variable((Nn, Nn), symmetric=True); G, G1 = G_struct(Nn, n)
    Ab, Bb, Cb = S['Abar'], S['Bbar'], S['Cbar']
    cons = [P >> eps*np.eye(Nn)]; Ys = [cp.Variable((n, p)) for _ in S['Hs']]
    for Y, H in zip(Ys, S['Hs']):
        for H2 in S['Hs']:
            Cm = Cb - H2
            cons.append(cp.bmat([[-P, Ab.T@G.T - Cm.T@Y.T@Bb.T],
                                 [G@Ab - Bb@Y@Cm, P - G - G.T]]) << -eps*np.eye(2*Nn))
    try:
        prob = cp.Problem(cp.Minimize(0), cons); prob.solve(**SOLVE)
        ok = prob.status in ("optimal", "optimal_inaccurate")
        return ok, ([np.linalg.solve(G1.value, Y.value) for Y in Ys] if ok else None)
    except Exception:
        return False, None

# ---------------- Hinf of G(z) = Co (zI-Ao)^-1 Bo ----------------
def hinf(S, L, ngrid=4000):
    Ao, Bo = AoBo(S, L)
    if max(abs(np.linalg.eigvals(Ao))) >= 1: return np.inf
    th = np.linspace(0, np.pi, ngrid); best_val = 0
    I = np.eye(S['Nn'])
    for t in th:
        Gz = S['Co'] @ np.linalg.solve(np.exp(1j*t)*I - Ao, Bo)
        best_val = max(best_val, np.linalg.svd(Gz, compute_uv=False)[0])
    return best_val

# ---------------- RATE SYNTHESIS ----------------
def synth_rate(S,alpha,eps=1e-7):
    Nn,n,p=S['Nn'],S['n'],S['p']
    P=cp.Variable((Nn,Nn),symmetric=True); G,G1=G_struct(Nn,n); Y=cp.Variable((n,p))
    Ab,Bb,Cb=S['Abar'],S['Bbar'],S['Cbar']; cons=[P>>eps*np.eye(Nn)]
    for H in S['Hs']:
        Cm=Cb-H
        cons.append(cp.bmat([[-(alpha**2)*P, Ab.T@G.T-Cm.T@Y.T@Bb.T],
                             [G@Ab-Bb@Y@Cm, P-G-G.T]])<<-eps*np.eye(2*Nn))
    try:
        pr=cp.Problem(cp.Minimize(0),cons); pr.solve(solver=cp.CLARABEL)
        if pr.status not in ("optimal","optimal_inaccurate"): return None
        return np.linalg.solve(G1.value,Y.value)
    except Exception:
        return None

def best(S,lo=.05,hi=1.0,it=10):
    if synth_rate(S,hi) is None: return None,None
    for _ in range(it):
        m=.5*(lo+hi)
        if synth_rate(S,m) is not None: hi=m
        else: lo=m
    return hi, synth_rate(S,min(1.0,1.02*hi))

def synth(a,q,alpha,eps=1e-6):
    S=build(np.array([[a]]),np.array([[1.0]]),q); Nn=S['Nn']
    P=cp.Variable((Nn,Nn),symmetric=True); G,G1=G_struct(Nn,1); Y=cp.Variable((1,1))
    Ab,Bb,Cb=S['Abar'],S['Bbar'],S['Cbar']; cons=[P>>eps*np.eye(Nn)]
    for H in S['Hs']:
        Cm=Cb-H
        cons.append(cp.bmat([[-(alpha**2)*P,Ab.T@G.T-Cm.T@Y.T@Bb.T],
                             [G@Ab-Bb@Y@Cm,P-G-G.T]])<<-eps*np.eye(2*Nn))
    try:
        pr=cp.Problem(cp.Minimize(0),cons); pr.solve(solver=cp.CLARABEL)
        if pr.status not in ("optimal","optimal_inaccurate"): return None
        L=np.linalg.solve(G1.value,Y.value)
    except Exception: return None
    # verification: analysis LMI in P only, with the recovered gain
    Q=cp.Variable((Nn,Nn),symmetric=True); c2=[Q>>np.eye(Nn)]
    for H in S['Hs']:
        M=Ab-Bb@L@(Cb-H)
        if max(abs(np.linalg.eigvals(M)))>=alpha: return None
        c2.append(M.T@Q@M-(alpha**2)*Q<<-1e-8*np.eye(Nn))
    try:
        pr2=cp.Problem(cp.Minimize(0),c2); pr2.solve(solver=cp.CLARABEL)
        return L if pr2.status in ("optimal","optimal_inaccurate") else None
    except Exception: return None
