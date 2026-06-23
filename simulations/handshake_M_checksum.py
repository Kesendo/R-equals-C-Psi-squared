"""Deepen the relation ‖M‖_F² + λ_min(MMᵀ) = 2 that fell out of the handshake M.

Hand-derivation of one half:  ‖M‖_F² = Σ_b ‖V_b ψ₁‖² = Σ_b (c_a²+c_{a+1}²) = 2 − (c₀²+c_{N-1}²),
since each interior site has degree 2 and the two free ENDS degree 1, and Σ c_i² = 1.
So  2 − ‖M‖_F²  =  c₀²+c_{N-1}²  =  the carrier's weight on the two chain ends (the endpoint leakage).
Closed form: c₀=c_{N-1}=√(2/(N+1))·sin(π/(N+1))  ⇒  endpoint = (4/(N+1))·sin²(π/(N+1)).

The OTHER half (the empirical surprise) to gate: is  λ_min(MMᵀ) = the endpoint weight  EXACT, all N?
And what is its eigenvector?

Added 2026-06-20 (F124 review): the conserved envelope Q (Part 2's bulk-cancellation heart) and the
frame-theory reading (λ_min = σ_min²(M) = 1/‖S⁻¹‖, the condition number, the K-partner kernel) are
gated here too, so the proof's prose claims are reproducible, not merely asserted.
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

    # --- Part 2's conserved envelope Q (the bulk-cancellation heart), gate-checked ---
    print("\n--- conserved envelope Q_a = c_a²+c_{a+1}² − E₁·c_a·c_{a+1},  E₁ = 2cos(π/(N+1)) ---", flush=True)
    qok = True
    for N in range(3, 21):
        c = carrier(N)
        E1 = 2.0 * np.cos(np.pi / (N + 1))
        Q = c[:-1] ** 2 + c[1:] ** 2 - E1 * c[:-1] * c[1:]
        qok &= bool(np.allclose(Q, Q[0], atol=1e-12) and abs(Q[0] - c[0] ** 2) < 1e-12)
    print(f"  Q constant along the chain AND = c₀² for all N=3..20?  {'YES' if qok else 'NO'}", flush=True)
    print("    => λ_min = 2·Q₀ = 2c₀² = c₀²+c_{N-1}² = E: the bulk telescopes, only the boundary survives.", flush=True)

    # --- the frame reading, gate-checked: λ_min = σ_min²(M) = 1/‖S⁻¹‖ (S=MᵀM), plus κ and the K-partner kernel ---
    print("\n--- frame identities: λ_min(MMᵀ) = σ_min²(M) = lower frame bound 1/‖S⁻¹‖, plus κ and the kernel ---", flush=True)
    print(f"{'N':>3} {'λ_min':>11} {'σ_min²':>11} {'1/‖S⁻¹‖':>11} {'κ=λ_max/λ_min':>14} {'k=N col‖·‖':>11}", flush=True)
    fok = True
    for N in range(3, 21):
        M, c = build_M(N)
        G = M @ M.T
        evG = np.linalg.eigvalsh(G)
        lam_min, lam_max = evG[0], evG[-1]
        sigma_min_sq = np.linalg.svd(M, compute_uv=False)[-1] ** 2     # Eckart-Young: dist² to rank-deficiency
        evS = np.linalg.eigvalsh(M.T @ M)                              # frame operator S=MᵀM, one ~0 eigenvalue
        inv_frame_bound = evS[1]                                       # smallest NONZERO (evS[0] = K-partner kernel)
        kpartner_col = np.linalg.norm(M[:, N - 1])                     # k=N column ⟨ψ_N|V_b|ψ_1⟩, must be ~0
        psiN = np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * N * (np.arange(N) + 1) / (N + 1))
        kpartner_id = np.linalg.norm(psiN - (-1.0) ** np.arange(N) * c)   # ψ_N = (−1)^i ψ_1
        fok &= bool(abs(sigma_min_sq - lam_min) < 1e-10 and abs(inv_frame_bound - lam_min) < 1e-9
                    and abs(evS[0]) < 1e-10 and kpartner_col < 1e-12 and kpartner_id < 1e-12)
        if N in (4, 5, 6):
            print(f"{N:>3} {lam_min:>11.7f} {sigma_min_sq:>11.7f} {inv_frame_bound:>11.7f} "
                  f"{lam_max / lam_min:>14.4f} {kpartner_col:>11.1e}", flush=True)
    print(f"  λ_min = σ_min² = 1/‖S⁻¹‖, K-partner kernel ψ_N=(−1)^i ψ_1 null, all N=3..20?  {'YES' if fok else 'NO'}", flush=True)
    print("    κ grows (≈3.43, 5.30, 7.60 at N=4,5,6): the end-leakage E is the conditioner of the bond→mode map.", flush=True)


if __name__ == "__main__":
    main()
