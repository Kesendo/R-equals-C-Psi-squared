"""
Analytical w=1 Spectrum Verification.

Extract ALL w=1 Liouvillian frequencies (XY-weight 1) for the
Heisenberg chain at N=2 through N=6. Compare with candidate
dispersion relations to find the exact analytical form.

w=1 modes: exactly one qubit carries X or Y, rest I or Z.
Under uniform Z-dephasing, all w=1 modes decay at rate 2*gamma.
Their frequencies come from the Hamiltonian structure.
"""

import numpy as np
from scipy.optimize import curve_fit
from pathlib import Path


# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(op, t, n):
    r = np.eye(1, dtype=complex)
    for k in range(n):
        r = np.kron(r, op if k == t else I2)
    return r


def build_H(n, J=1.0):
    d = 2**n
    H = np.zeros((d, d), dtype=complex)
    for i in range(n - 1):
        for p in [X, Y, Z]:
            H += J * kron_at(p, i, n) @ kron_at(p, i + 1, n)
    return H


def build_L(H, gammas):
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    n = int(np.log2(d))
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
    for k in range(n):
        if gammas[k] <= 0:
            continue
        Lk = np.sqrt(gammas[k]) * kron_at(Z, k, n)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk)
        L -= 0.5 * np.kron(Id, LdL)
        L -= 0.5 * np.kron(LdL.T, Id)
    return L


def extract_w1_frequencies(N, J=1.0, gamma=0.001):
    """Extract all w=1 Liouvillian frequencies for N-qubit chain."""
    H = build_H(N, J)
    L = build_L(H, [gamma] * N)
    evals = np.linalg.eigvals(L)

    # w=1 modes have rate = 2*gamma (within tolerance)
    target_rate = 2 * gamma
    tol = 0.1 * gamma  # 10% tolerance

    w1_freqs = []
    for ev in evals:
        rate = -ev.real
        freq = abs(ev.imag)
        if abs(rate - target_rate) < tol and freq > 1e-8:
            w1_freqs.append(round(freq, 8))

    # Remove duplicates (conjugate pairs give same |freq|)
    unique = sorted(set(round(f, 6) for f in w1_freqs))
    return unique


def main():
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "analytical_spectrum_verify.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=" * 70)
    out("ANALYTICAL w=1 SPECTRUM VERIFICATION")
    out("=" * 70)

    J = 1.0

    # ================================================================
    # Step 1: Extract all w=1 frequencies for N=2 through N=6
    # ================================================================
    out("\n--- STEP 1: Extract w=1 frequencies ---")

    all_freqs = {}
    for N in range(2, 7):
        out(f"\nComputing N={N}...", )
        freqs = extract_w1_frequencies(N, J, gamma=0.001)
        all_freqs[N] = freqs
        out(f"  N={N}: {len(freqs)} distinct w=1 frequencies")
        for i, f in enumerate(freqs):
            out(f"    f[{i}] = {f:.6f}")

    # ================================================================
    # Step 2: Test candidate dispersion relations
    # ================================================================
    out("\n" + "=" * 70)
    out("STEP 2: Test candidate dispersion relations")
    out("=" * 70)

    # Candidate A: omega_k = 4J * (1 - cos(pi*k/N)), k=1..N-1
    # Candidate B: omega_k = 4J * sin(pi*k/(2N)), k=1..N-1
    # Candidate C: omega_k = 4J * (1 - cos(pi*k/(N+1))), k=1..N
    # Candidate D: omega_k = 2J * (1 - cos(pi*k/N)), k=1..N-1
    # Candidate E: omega_k = 4J * sin(pi*k/N), k=1..N-1
    # Candidate F: omega_k = 4J * |sin(pi*k/(N))|, k=1..N-1

    candidates = {
        'A: 4J(1-cos(pi*k/N))': lambda k, N: 4 * J * (1 - np.cos(np.pi * k / N)),
        'B: 4J*sin(pi*k/(2N))': lambda k, N: 4 * J * np.sin(np.pi * k / (2 * N)),
        'C: 4J(1-cos(pi*k/(N+1)))': lambda k, N: 4 * J * (1 - np.cos(np.pi * k / (N + 1))),
        'D: 2J(1-cos(pi*k/N))': lambda k, N: 2 * J * (1 - np.cos(np.pi * k / N)),
        'E: 4J*sin(pi*k/N)': lambda k, N: 4 * J * np.abs(np.sin(np.pi * k / N)),
        'F: 4J*|sin(pi*k/(N+1))|': lambda k, N: 4 * J * np.abs(np.sin(np.pi * k / (N + 1))),
    }

    out(f"\nFor each candidate: generate predicted frequencies,")
    out(f"sort, and compare with numerically extracted w=1 frequencies.")

    for name, formula in candidates.items():
        out(f"\n--- {name} ---")
        total_err = 0
        total_count = 0
        all_match = True

        for N in range(2, 7):
            num_freqs = all_freqs[N]
            n_modes = len(num_freqs)

            # Generate predicted frequencies for different k ranges
            # Try k=1..N-1 and k=1..N
            for k_max_label, k_range in [("k=1..N-1", range(1, N)),
                                          ("k=1..N", range(1, N + 1))]:
                pred = sorted(set(round(formula(k, N), 6) for k in k_range))
                if len(pred) == n_modes:
                    # Check match
                    max_err = max(abs(p - n) for p, n in zip(pred, num_freqs))
                    if max_err < 1e-4:
                        out(f"  N={N}: MATCH ({k_max_label}),"
                            f" max_err={max_err:.2e},"
                            f" {n_modes} frequencies")
                        total_err += max_err
                        total_count += 1
                        break
            else:
                # Neither k range matched
                for k_max_label, k_range in [("k=1..N-1", range(1, N))]:
                    pred = sorted(set(round(formula(k, N), 6)
                                      for k in k_range))
                    if pred:
                        diffs = []
                        for nf in num_freqs:
                            closest = min(pred, key=lambda p: abs(p - nf))
                            diffs.append(abs(closest - nf))
                        max_err = max(diffs) if diffs else 999
                        out(f"  N={N}: {len(pred)} predicted vs"
                            f" {n_modes} numerical,"
                            f" max_err={max_err:.4f}")
                        all_match = False

        if total_count > 0:
            out(f"  >> Mean max error: {total_err/total_count:.2e}"
                f" ({total_count}/{len(range(2,7))} matched)")

    # ================================================================
    # Step 3: Deep analysis of the winning formula
    # ================================================================
    out("\n" + "=" * 70)
    out("STEP 3: Detailed verification of best candidate")
    out("=" * 70)

    # Test the most promising candidates in detail
    for name, formula, k_gen in [
        ("A: 4J(1-cos(pi*k/N))", lambda k, N: 4*J*(1-np.cos(np.pi*k/N)),
         lambda N: range(1, N)),
        ("C: 4J(1-cos(pi*k/(N+1)))", lambda k, N: 4*J*(1-np.cos(np.pi*k/(N+1))),
         lambda N: range(1, N+1)),
    ]:
        out(f"\n--- Detailed: {name} ---")
        out(f"{'N':>2} {'k':>3} {'Predicted':>12} {'Numerical':>12}"
            f" {'Error':>12} {'Match':>6}")

        for N in range(2, 7):
            num_freqs = all_freqs[N]
            pred_freqs = sorted(set(round(formula(k, N), 6)
                                    for k in k_gen(N)))

            # Pair predicted with numerical by closest match
            used = [False] * len(num_freqs)
            for ki, p in enumerate(pred_freqs):
                best_j = -1
                best_err = 999
                for j, n in enumerate(num_freqs):
                    if not used[j] and abs(p - n) < best_err:
                        best_err = abs(p - n)
                        best_j = j
                if best_j >= 0:
                    used[best_j] = True
                    match = "YES" if best_err < 1e-4 else "NO"
                    out(f"{N:2d} {ki+1:3d} {p:12.6f}"
                        f" {num_freqs[best_j]:12.6f}"
                        f" {best_err:12.2e} {match:>6}")
                else:
                    out(f"{N:2d} {ki+1:3d} {p:12.6f}"
                        f" {'---':>12} {'---':>12} {'NO':>6}")

            # Any unmatched numerical frequencies?
            unmatched = [num_freqs[j] for j in range(len(num_freqs))
                         if not used[j]]
            if unmatched:
                out(f"  Unmatched numerical: {unmatched}")

    # ================================================================
    # Step 4: Summary and implications
    # ================================================================
    out("\n" + "=" * 70)
    out("STEP 4: Summary")
    out("=" * 70)

    # Check which formula gives the max frequency = 4J(1+cos(pi/N))
    out("\nMax frequency check (must equal 4J(1+cos(pi/N))):")
    for N in range(2, 7):
        expected = 4 * J * (1 + np.cos(np.pi / N))
        actual = max(all_freqs[N]) if all_freqs[N] else 0
        match = "YES" if abs(expected - actual) < 1e-4 else "NO"
        out(f"  N={N}: expected={expected:.6f}, actual={actual:.6f},"
            f" match={match}")

    # For candidate A: max at k=N-1 gives 4J(1-cos(pi*(N-1)/N))
    # = 4J(1+cos(pi/N)) since cos(pi-x) = -cos(x). CHECK!
    out("\nCandidate A max frequency:")
    for N in range(2, 7):
        k = N - 1
        pred_max = 4 * J * (1 - np.cos(np.pi * k / N))
        expected = 4 * J * (1 + np.cos(np.pi / N))
        out(f"  N={N}: 4J(1-cos(pi*{k}/{N})) = {pred_max:.6f},"
            f" 4J(1+cos(pi/{N})) = {expected:.6f},"
            f" equal: {abs(pred_max-expected)<1e-10}")

    out("\n=== DONE ===")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
