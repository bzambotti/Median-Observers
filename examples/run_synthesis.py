import numpy as np
import json
from src.config import A, C
from src.lmi import build, switched_codesign, best, synth

def main():
    print("---- Generic Synthesis ----")
    gains = {}
    for q in [1, 2]:
        S = build(A, C, q)
        ok, Ls = switched_codesign(S, scheduled=False)
        gains[q] = Ls[0]
        rho = max(abs(np.linalg.eigvals(A - Ls[0] @ C)))
        print(f"q={q} feasible={ok} rho(A-LC)={rho:.3f}", flush=True)

    np.save('gains.npy', np.array([gains[1], gains[2]]))
    print("L1=\n", gains[1], "\nL2=\n", gains[2])

    print("\n---- Optimal Rate Synthesis ----")
    for q in [1, 2]:
        S = build(A, C, q)
        a, L = best(S)
        rho = max(abs(np.linalg.eigvals(A - L @ C)))
        print(f"q={q}: alpha*={a:.4f} rho(A-LC)={rho:.4f}", flush=True)
        np.save(f'Lopt{q}.npy', L)

    print("\n---- Certified Rate Bound Arrays ----")
    qs = [1, 2, 3, 4, 5, 6]
    out = {"q": qs}

    for a in [0.5, 0.7, 0.9]:
        row = []
        for q in qs:
            if synth(a, q, 0.999) is None:
                row.append(None)
                continue
            lo, hi = 0.02, 0.999
            for _ in range(13):
                m = .5 * (lo + hi)
                if synth(a, q, m) is not None:
                    hi = m
                else:
                    lo = m
            row.append(hi)
        out[f"a{a}"] = row
        print(a, [None if r is None else round(r, 4) for r in row], flush=True)

    ceil = []
    for q in qs:
        lo, hi = 1.0, 2.0
        if synth(1.0, q, 0.999) is None:
            ceil.append(1.0)
            print("ceil", q, 1.0, flush=True)
            continue
        for _ in range(14):
            m = .5 * (lo + hi)
            if synth(m, q, 0.999) is not None:
                lo = m
            else:
                hi = m
        ceil.append(lo)
        print("ceil", q, round(lo, 4), flush=True)

    out["ceil"] = ceil
    with open('fig2c.json', 'w') as f:
        json.dump(out, f)
    print("done")

if __name__ == "__main__":
    main()
