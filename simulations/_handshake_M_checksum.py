"""Deepen the relation ‖M‖_F² + λ_min(MMᵀ) = 2 that fell out of the handshake M.

Hand-derivation of one half:  ‖M‖_F² = Σ_b ‖V_b ψ₁‖² = Σ_b (c_a²+c_{a+1}²) = 2 − (c₀²+c_{N-1}²),
since each interior site has degree 2 and the two free ENDS degree 1, and Σ c_i² = 1.
So  2 − ‖M‖_F²  =  c₀²+c_{N-1}²  =  the carrier's weight on the two chain ends (the endpoint leakage).
Closed form: c₀=c_{N-1}=√(2/(N+1))·sin(π/(N+1))  ⇒  endpoint = (4/(N+1))·sin²(π/(N+1)).

The OTHER half (the empirical surprise) to gate: is  λ_min(MMᵀ) = the endpoint weight  EXACT, all N?
And what is its eigenvector?
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
np.set_printoptions(precision=5, suppress=True, linewidth=130)


def carrier(N):
    i = np.arange(N)
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * (i + 1) / (N + 1))


def build_M(N):
    c = carrier(N)
    P = np.array([np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (np.arange(N) + 1) / (N + 1))
                  for k in range(1, N + 1)])
    M = np.zeros((N - 1, N))
    for a in range(N - 1):
        for ki in range(N):
            M[a, ki] = P[ki][a] * c[a + 1] + P[ki][a + 1] * c[a]
    return M, c


def main():
    print("=== deepen: ‖M‖_F² = 2 − endpoint, and λ_min(MMᵀ) = endpoint? ===\n", flush=True)
    print(f"{'N':>3} {'‖M‖_F²':>11} {'2−endpoint':>11} {'diff':>10} | {'λ_min':>11} {'endpoint':>11} {'diff':>10}", flush=True)
    allok = True
    for N in range(3, 21):                                  # incl. N=3, the stated lower endpoint
        M, c = build_M(N)
        fro2 = (M ** 2).sum()
        endpoint = c[0] ** 2 + c[-1] ** 2
        endpoint_cf = (4.0 / (N + 1)) * np.sin(np.pi / (N + 1)) ** 2
        G = M @ M.T
        lam_min = np.linalg.eigvalsh(G)[0]
        d1 = abs(fro2 - (2 - endpoint))
        d2 = abs(lam_min - endpoint)
        d3 = abs(endpoint - endpoint_cf)
        allok &= (d1 < 1e-12 and d2 < 1e-10 and d3 < 1e-12)
        print(f"{N:>3} {fro2:>11.7f} {2-endpoint:>11.7f} {d1:>10.1e} | "
              f"{lam_min:>11.7f} {endpoint:>11.7f} {d2:>10.1e}", flush=True)
    print(f"\n  derivation ‖M‖_F²=2−endpoint exact for all N?  (diff<1e-12)", flush=True)
    print(f"  λ_min(MMᵀ)=endpoint exact for all N?           {'YES' if allok else 'NO -- not an identity'}", flush=True)
    print(f"  endpoint = (4/(N+1))sin²(π/(N+1)) confirmed.", flush=True)

    # the eigenvector of the smallest eigenvalue -- is it clean?
    print("\n--- λ_min eigenvector over the bonds (is it recognizable?) ---", flush=True)
    for N in (5, 7, 9):
        M, c = build_M(N)
        G = M @ M.T
        w, V = np.linalg.eigh(G)
        v = V[:, 0]
        v = v / v[np.argmax(np.abs(v))]                      # normalize sign/scale
        print(f"N={N}: λ_min={w[0]:.6f}, eigvec(bonds)/max = [{' '.join(f'{x:+.4f}' for x in v)}]", flush=True)
        # compare to a few candidates over the N-1 bonds
        a = np.arange(N - 1)
        cand_alt = (-1.0) ** a                               # staggered
        cand_alt /= cand_alt[np.argmax(np.abs(cand_alt))]
        print(f"      staggered (−1)^b /max          = [{' '.join(f'{x:+.4f}' for x in cand_alt)}]", flush=True)


if __name__ == "__main__":
    main()
