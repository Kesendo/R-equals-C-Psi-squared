"""GATE-FIRST: is F124's conditioning literally the bond-recovery inverse problem?

The borrowing-a-discipline lens (signal/control) claims: M[b,k]=<psi_k|V_b|psi_1> is the
forward map bond -> mode; recovering a bond defect from the mode response is M^+, whose
conditioning is set by the singular values of M. Term-for-term:
  - sigma_min(M)^2 = lambda_min(M M^T) = E = the lower frame bound = 1/||S^-1||
  - the worst-conditioned BOND direction (smallest left-singular vector of M) = the
    staggered (-1)^b zone-boundary pattern
  - a unit bond defect's mode-response CONTRAST = its singular value; so the staggered
    defect is sqrt(kappa) = sigma_max/sigma_min times harder to SEE than a band-edge one,
    and that ratio grows ~N (since kappa ~ N^2)
  - the K-partner psi_N is the output channel no bond excites (the null column) = the
    unobservable measurement direction

Each STAGE-0 gate can FIRE if the borrowed tool does not actually map. A firing gate IS the
find (do not loosen it). Stage 0 is pure first-principles (M only, no decoder); Stage 1
(the live DefectDecoder connection) is a separate script once Stage 0 holds.
"""
import sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
np.set_printoptions(precision=5, suppress=True, linewidth=130)


def carrier(N, k=1):
    j = np.arange(N)
    return np.sqrt(2.0 / (N + 1)) * np.sin(np.pi * k * (j + 1) / (N + 1))


def build_M(N):
    """Full transition matrix M[b,k] = <psi_k|V_b|psi_1>, bonds b=0..N-2, modes k=1..N."""
    c = carrier(N, 1)
    P = np.array([carrier(N, k) for k in range(1, N + 1)])  # P[k-1] = psi_k
    M = np.zeros((N - 1, N))
    for a in range(N - 1):
        for ki in range(N):
            M[a, ki] = P[ki][a] * c[a + 1] + P[ki][a + 1] * c[a]
    return M


def main():
    GATE = {"fired": []}

    def gate(name, cond, detail=""):
        flag = "ok " if cond else "GATE-FIRE"
        if not cond:
            GATE["fired"].append(name)
        print(f"   [{flag}] {name}" + (f"   {detail}" if detail else ""))

    Ns = list(range(4, 13))
    print("=" * 96)
    print("STAGE 0 - is F124's conditioning the bond-recovery inverse problem?")
    print("=" * 96)

    # --- G1: the worst-conditioned bond direction (smallest left-singular vector) is staggered ---
    print("\nG1: smallest left-singular vector of M (worst-conditioned bond direction) = staggered (-1)^b")
    for N in Ns:
        M = build_M(N)
        U, S, Vt = np.linalg.svd(M, full_matrices=False)   # U: (N-1)x(N-1) bond space
        u_min = U[:, -1]                                     # left vector of smallest sigma
        stag = (-1.0) ** np.arange(N - 1)
        cos = abs(u_min @ stag) / (np.linalg.norm(u_min) * np.linalg.norm(stag))
        gate(f"N={N}: |cos(u_min, staggered)| ~ 1", cos > 0.9999, f"cos={cos:.6f}")

    # --- G2: sigma_min^2 = E; the ROBUST resolution law is kappa ~ N^2; the absolute floor
    #         E ~ (N+1)^-3 only ASYMPTOTICALLY (E*(N+1)^3 -> 4*pi^2). Fitting E vs N over small N
    #         gives ~-2.5 (pre-asymptotic + wrong variable); that is a finite-size artifact, not the law. ---
    print("\nG2: sigma_min^2 = E; kappa ~ N^2 (robust); E ~ (N+1)^-3 asymptotically (E*(N+1)^3 -> 4*pi^2)")
    Es, smins, smaxs, kappas = [], [], [], []
    for N in Ns:
        M = build_M(N)
        S = np.linalg.svd(M, compute_uv=False)
        E_cf = (4.0 / (N + 1)) * np.sin(np.pi / (N + 1)) ** 2
        smin, smax = S[-1], S[0]
        Es.append(E_cf); smins.append(smin); smaxs.append(smax); kappas.append((smax / smin) ** 2)
        gate(f"N={N}: sigma_min^2 = E", abs(smin ** 2 - E_cf) < 1e-10,
             f"sigma_min^2={smin**2:.6f}, E={E_cf:.6f}, kappa={kappas[-1]:.3f}, E*(N+1)^3={E_cf*(N+1)**3:.3f}")
    # robust law: kappa ~ N^2 (the ratio; (N+1)-corrections cancel) -- fit exponent ~ 2
    p_kappa = np.polyfit(np.log(Ns), np.log(kappas), 1)[0]
    gate("kappa ~ N^2 (the robust resolution-degradation law)", 1.8 < p_kappa < 2.2,
         f"fitted exponent={p_kappa:.3f}")
    # absolute floor: E ~ (N+1)^-3 in the RIGHT variable, over an EXTENDED range (N=4..60)
    Nbig = np.arange(4, 61)
    Ebig = (4.0 / (Nbig + 1)) * np.sin(np.pi / (Nbig + 1)) ** 2
    p_E = np.polyfit(np.log(Nbig + 1), np.log(Ebig), 1)[0]
    gate("E ~ (N+1)^-3 (right variable, N=4..60; fitted exponent in [-3.1, -2.9])", -3.1 < p_E < -2.9,
         f"fitted exponent={p_E:.3f} (vs the naive N-fit over N=4..12 which gives ~-2.5, a finite-size artifact)")
    # the asymptotic constant E*(N+1)^3 -> 4*pi^2 = 39.478
    gate("E*(N+1)^3 -> 4*pi^2 = 39.478 (the asymptotic floor constant)",
         abs(Ebig[-1] * (Nbig[-1] + 1) ** 3 - 4 * np.pi ** 2) < 0.2,
         f"E*(N+1)^3 at N=60 = {Ebig[-1]*(Nbig[-1]+1)**3:.4f}")

    # --- G3: the response CONTRAST ratio (staggered vs band-edge bond pattern) = sqrt(kappa) ~ N ---
    print("\nG3: mode-response contrast of a unit staggered bond defect vs a unit band-edge bond defect = sqrt(kappa) ~ N")
    sqrt_kappas = []
    for i, N in enumerate(Ns):
        M = build_M(N)
        # band-edge bond pattern = largest left-singular vector; staggered = smallest
        U, S, Vt = np.linalg.svd(M, full_matrices=False)
        d_stag = (-1.0) ** np.arange(N - 1); d_stag /= np.linalg.norm(d_stag)
        d_edge = U[:, 0]                                    # most-visible bond direction
        contrast_stag = np.linalg.norm(M.T @ d_stag)       # = sigma_min (response magnitude)
        contrast_edge = np.linalg.norm(M.T @ d_edge)       # = sigma_max
        ratio = contrast_edge / contrast_stag
        sk = np.sqrt(kappas[i])
        sqrt_kappas.append(sk)
        gate(f"N={N}: contrast ratio = sqrt(kappa)", abs(ratio - sk) < 1e-6,
             f"ratio={ratio:.4f}, sqrt(kappa)={sk:.4f}")
    p_sk = np.polyfit(np.log(Ns), np.log(sqrt_kappas), 1)[0]
    gate("sqrt(kappa) ~ N (the staggered defect is N-times harder to see)", 0.8 < p_sk < 1.3,
         f"fitted exponent={p_sk:.3f}; sqrt(kappa) at N=4..8 = "
         f"{[round(x,2) for x in sqrt_kappas[:5]]}")

    # --- G4: the METHOD the native view did not run - a noisy Monte-Carlo detection SNR test ---
    print("\nG4: Monte-Carlo - at fixed noise, the minimum detectable defect magnitude ratio (staggered/edge) = sqrt(kappa)")
    rng = np.random.default_rng(20260620)
    sigma_eta = 1e-3
    trials = 4000
    for i, N in enumerate(Ns[:5]):                          # N=4..8
        M = build_M(N)
        U, S, Vt = np.linalg.svd(M, full_matrices=False)
        d_stag = (-1.0) ** np.arange(N - 1); d_stag /= np.linalg.norm(d_stag)
        d_edge = U[:, 0]
        # detection SNR of a unit-magnitude defect = ||M^T d|| / sigma_eta = sigma_dir / sigma_eta
        # min detectable magnitude (SNR=1) ~ sigma_eta / sigma_dir; ratio staggered/edge = sigma_max/sigma_min
        # empirical check via matched-filter detection statistic over trials
        def detect_stat(d):
            r0 = M.T @ d
            hits = 0.0
            mf = r0 / np.linalg.norm(r0)
            for _ in range(trials):
                eta = rng.normal(0, sigma_eta, size=N)
                hits += (mf @ (r0 + eta))
            return hits / trials  # ~ ||r0|| = sigma_dir for unit d
        s_stag = detect_stat(d_stag)
        s_edge = detect_stat(d_edge)
        ratio = s_edge / s_stag
        sk = sqrt_kappas[i]
        gate(f"N={N}: MC matched-filter response ratio (edge/staggered) = sqrt(kappa)",
             abs(ratio - sk) / sk < 0.05, f"MC ratio={ratio:.4f}, sqrt(kappa)={sk:.4f}")

    print("\n" + "=" * 96)
    if GATE["fired"]:
        print(f"STAGE 0: {len(GATE['fired'])} GATE(S) FIRED -> {GATE['fired']}")
        print("A firing gate is the find: diagnose (is the borrowed tool wrong, or our number?), do NOT loosen.")
    else:
        print("STAGE 0: ALL GATES PASS -> F124's conditioning IS the bond-recovery inverse problem.")
        print("The staggered q=pi defect is sqrt(kappa)~N times harder to see; the worst bond direction")
        print("is the zone-boundary mode; sigma_min = sqrt(E) is the reconstruction floor. The borrowed")
        print("signal/control tool fired: it handed us the detection-SNR method, not just a renaming.")
        print("NEXT (Stage 1): does the LIVE DefectDecoder show this sqrt(kappa)~N degradation, and is")
        print("its known 'ambiguity ~1.5 at N=5' = sqrt(kappa(5))=2.30 (gate) or a different quantity?")
    print("=" * 96)


if __name__ == "__main__":
    main()
