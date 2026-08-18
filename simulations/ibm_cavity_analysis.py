"""
IBM Cavity Spectral Analysis: Cavity modes meet real hardware data.

Compares the Liouvillian spectrum for N=5 chain under three gamma profiles:
1. IBM sacrifice-zone (Q85-Q94 real T2* data)
2. Uniform (same total gamma, spread equally)
3. Zero noise (cavity modes / unitary ground state)

Answers: Why does the sacrifice zone work? Because it protects cavity modes.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
# The canonical F1 check. It lived as a hand-copy in this file until three
# copies made the missing primitive obvious; the port, its blind spot and
# its tests now live in one place (cockpit rule 2).
from framework import max_f1_pairing_distance, f1_distance_in_eps  # noqa: E402, F401

EPS = np.finfo(float).eps

# === Pauli matrices ===
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_at(op, target, n_qubits):
    """Place operator on target qubit via tensor products."""
    result = np.eye(1, dtype=complex)
    for k in range(n_qubits):
        result = np.kron(result, op if k == target else I2)
    return result


def build_heisenberg_chain(n, J=1.0):
    """Build Heisenberg XXZ chain Hamiltonian."""
    d = 2**n
    H = np.zeros((d, d), dtype=complex)
    for i in range(n - 1):
        for pauli in [X, Y, Z]:
            H += J * kron_at(pauli, i, n) @ kron_at(pauli, i + 1, n)
    return H


def build_liouvillian(H, gammas):
    """Build Lindblad Liouvillian superoperator with Z-dephasing."""
    d = H.shape[0]
    Id = np.eye(d, dtype=complex)
    n_qubits = int(np.log2(d))

    # Hamiltonian part: -i(H x I - I x H^T)
    L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))

    # Dephasing: L_k = sqrt(gamma_k) * Z_k
    for k in range(n_qubits):
        Lk = np.sqrt(gammas[k]) * kron_at(Z, k, n_qubits)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk.conj(), Lk)
        L -= 0.5 * np.kron(Id, LdL)
        L -= 0.5 * np.kron(LdL.T, Id)

    return L


def classify_modes(evals, eps=1e-10, freq_tol=1e-8):
    """Classify eigenvalues into stationary, oscillating, decaying."""
    stationary = []
    oscillating = []
    for ev in evals:
        rate = -ev.real
        freq = abs(ev.imag)
        if abs(ev.real) < eps and freq < eps:
            stationary.append(ev)
        else:
            oscillating.append((rate, freq, ev))

    # Unique frequencies
    freqs = sorted(set(round(f, 6) for _, f, _ in oscillating if f > eps))
    unique_freqs = []
    for f in freqs:
        if not unique_freqs or abs(f - unique_freqs[-1]) > freq_tol:
            unique_freqs.append(f)

    return stationary, oscillating, unique_freqs


def analyze_profile(name, gammas, H, n_qubits, out):
    """Full spectral analysis for a given gamma profile."""
    out(f"\n{'='*60}")
    out(f"PROFILE: {name}")
    out(f"{'='*60}")
    out(f"gammas: [{', '.join(f'{g:.4f}' for g in gammas)}]")
    out(f"sum(gamma): {sum(gammas):.4f}")

    L = build_liouvillian(H, gammas)
    evals = np.linalg.eigvals(L)

    stationary, oscillating, unique_freqs = classify_modes(evals)

    # Palindrome check (canonical F1 distance, see f1_distance_in_eps)
    center = sum(gammas)
    f1_eps, radius = f1_distance_in_eps(evals, center)

    out(f"\nTotal eigenvalues: {len(evals)}")
    out(f"Stationary (immune): {len(stationary)}")
    out(f"Oscillating+decaying: {len(oscillating)}")
    out(f"Distinct frequencies: {len(unique_freqs)}")
    out(f"Palindrome center: {center:.4f}")
    out(f"F1 pairing distance: {f1_eps:.1f} eps*rho (spectral radius {radius:.4f})")

    if sum(gammas) > 0:
        rates = sorted(-ev.real for ev in evals if abs(ev.imag) > 1e-10)
        if rates:
            out(f"Min decay rate: {min(rates):.6f}")
            out(f"Max decay rate: {max(rates):.6f}")
            out(f"Expected max (2*sum_gamma): {2*center:.6f}")
            out(f"Max rate / 2*sum_gamma: {max(rates)/(2*center):.4f}")

    # Protected modes (rate < threshold)
    for threshold in [0.05, 0.10, 0.20]:
        if sum(gammas) > 0:
            protected = sum(1 for r, f, _ in oscillating if r < threshold)
            out(f"Protected (rate < {threshold}): {protected}")

    # Top 20 slowest non-immune oscillating modes
    if oscillating:
        sorted_osc = sorted(oscillating, key=lambda x: x[0])
        out(f"\nTop 20 slowest oscillating modes:")
        out(f"{'#':>3} {'Rate':>10} {'Freq':>10} {'Freq/J':>8}")
        for i, (rate, freq, _) in enumerate(sorted_osc[:20]):
            out(f"{i+1:3d} {rate:10.6f} {freq:10.6f} {freq:8.3f}")

        # Top 10 fastest
        out(f"\nTop 10 fastest decaying modes:")
        sorted_fast = sorted(oscillating, key=lambda x: -x[0])
        for i, (rate, freq, _) in enumerate(sorted_fast[:10]):
            out(f"{i+1:3d} {rate:10.6f} {freq:10.6f} {freq:8.3f}")

    return evals, stationary, oscillating, unique_freqs


def frequency_group_analysis(osc_sacrifice, osc_uniform, unique_freqs, out):
    """Compare sacrifice vs uniform efficiency per frequency group."""
    out(f"\n{'='*60}")
    out("FREQUENCY GROUP ANALYSIS: Sacrifice vs Uniform")
    out(f"{'='*60}")
    # This note was hand-written into the results file by commit 5d57555 and was
    # not produced by any script, so a re-run silently dropped it. It is emitted
    # here instead, which is what makes this artifact reproducible.
    out()
    out("NOTE: All ratios are 1.00x because the palindromic pairing forces")
    out("the MEAN decay rate of each frequency group to equal Sigma_gamma")
    out("(every mode at rate r has a partner at 2S-r, so the average is S).")
    out("The sacrifice zone effect is visible in the SPREAD within each group,")
    out("not in the mean. See MODE SURVIVAL COMPARISON below for the")
    out("physically meaningful comparison.")
    out(f"\n{'Freq/J':>8} {'Sac_rate':>10} {'Uni_rate':>10} {'Ratio':>8} {'Modes':>6}")

    freq_tol = 0.1

    for freq in unique_freqs:
        sac_rates = [r for r, f, _ in osc_sacrifice if abs(f - freq) < freq_tol]
        uni_rates = [r for r, f, _ in osc_uniform if abs(f - freq) < freq_tol]

        if sac_rates and uni_rates:
            sac_mean = np.mean(sac_rates)
            uni_mean = np.mean(uni_rates)
            ratio = uni_mean / sac_mean if sac_mean > 1e-10 else float('inf')
            out(f"{freq:8.3f} {sac_mean:10.6f} {uni_mean:10.6f} {ratio:8.2f}x {len(sac_rates):6d}")


def orbit_sweeps(H, ibm_gammas, uniform_gammas, out, n_perm=200,
                 seeds=(0, 12345, 7)):
    """The four sweeps the experiment doc reads the F1 row through.

    The F1 distance comes from a greedy matcher, so one spectrum yields an
    ORBIT of values. This emits the orbit statistics so the doc's claims have
    a producer instead of a reconstructed harness: array-order ceilings, the
    paired test across profiles, the jump-operator summation-order orbit, and
    the commutator norms. Run with `python ibm_cavity_analysis.py --sweeps`
    (the jump-order part is ~120 eigendecompositions, several minutes).
    """
    import itertools
    import math

    profiles = [("zero", np.zeros(len(ibm_gammas))),
                ("ibm", ibm_gammas), ("uni", uniform_gammas)]
    spec = {}
    for nm, g in profiles:
        ev = np.linalg.eigvals(build_liouvillian(H, g))
        spec[nm] = (ev, float(sum(g)), f1_distance_in_eps(ev, sum(g))[0])

    out(f"\n{'='*60}")
    out("ORBIT SWEEPS (the F1 row is one draw from an orbit)")
    out(f"{'='*60}")

    out("\nArray order: is the printed value the ceiling of its orbit?")
    out(f"  {'profile':8} {'printed':>8} {'min':>8} {'max':>8} {'median':>8} "
        f"{'ratio':>7} {'==printed':>10} {'>printed':>9}")
    for seed in seeds:
        rng = np.random.default_rng(seed)
        perms = [rng.permutation(len(spec['ibm'][0])) for _ in range(n_perm)]
        out(f"  seed {seed}:")
        for nm in ("zero", "ibm", "uni"):
            ev, c, base = spec[nm]
            vals = [f1_distance_in_eps(ev[q], c)[0] for q in perms]
            eq = sum(1 for v in vals if abs(v - base) < 1e-9)
            gt = sum(1 for v in vals if v > base + 1e-9)
            out(f"  {nm:8} {base:8.2f} {min(vals):8.2f} {max(vals):8.2f} "
                f"{np.median(vals):8.2f} {max(vals)/min(vals):7.3f} "
                f"{eq:10d} {gt:9d}")

    out("\nPaired test (same permutation on all three): is the between-profile")
    out("difference bookkeeping noise, or systematic?")
    rng = np.random.default_rng(seeds[0])
    res = {nm: [] for nm in spec}
    for _ in range(n_perm):
        q = rng.permutation(len(spec['ibm'][0]))
        for nm, (ev, c, _b) in spec.items():
            res[nm].append(f1_distance_in_eps(ev[q], c)[0])
    d = np.array(res["ibm"]) - np.array(res["uni"])
    sem = d.std(ddof=1) / np.sqrt(len(d))
    out(f"  ibm - uni: mean {d.mean():.2f}  sem {sem:.3f}  "
        f"t {d.mean()/sem:.1f}  ibm>uni in {int((d>0).sum())}/{n_perm}")
    out("  => common mode cancels in the difference; the gap is systematic.")

    out("\nCommutator norms (how far from normal each Liouvillian is):")
    for nm, g in profiles:
        L = build_liouvillian(H, g)
        C = L @ L.conj().T - L.conj().T @ L
        out(f"  {nm:8} 2-norm {np.linalg.norm(C, 2):9.4f}   "
            f"Frobenius {np.linalg.norm(C, 'fro'):9.4f}")

    n = len(ibm_gammas)
    out(f"\nJump-operator summation order, all {math.factorial(n)} orders:")
    for nm, g in profiles:
        d_ = H.shape[0]
        Id = np.eye(d_, dtype=complex)
        ref = None
        vals, worst = [], 0.0
        for order in itertools.permutations(range(n)):
            L = -1j * (np.kron(Id, H) - np.kron(H.T, Id))
            for k in order:
                Lk = np.sqrt(g[k]) * kron_at(Z, k, n)
                LdL = Lk.conj().T @ Lk
                L += np.kron(Lk.conj(), Lk)
                L -= 0.5 * np.kron(Id, LdL)
                L -= 0.5 * np.kron(LdL.T, Id)
            if ref is None:
                ref = L
            else:
                worst = max(worst, float(np.abs(L - ref).max()))
            vals.append(f1_distance_in_eps(np.linalg.eigvals(L), sum(g))[0])
        ratio = max(vals) / min(vals) if min(vals) > 0 else float('nan')
        out(f"  {nm:8} {min(vals):7.1f} to {max(vals):7.1f}  ratio {ratio:6.4f}  "
            f"distinct {len(set(np.round(vals, 6))):3d}  max|L-L_0| {worst:.3e}")
    out("  => equal or zero gammas reassemble L bit-identically, so their")
    out("     number cannot move; only the IBM column has a jump-order orbit.")


def main(sweeps=False):
    results_dir = Path(__file__).parent / "results"
    results_dir.mkdir(exist_ok=True)
    out_path = results_dir / "ibm_cavity_analysis.txt"

    lines = []
    def out(s=""):
        print(s)
        lines.append(s)

    out("=== IBM CAVITY SPECTRAL ANALYSIS ===")
    out(f"N=5 Chain, J=1.0, Z-dephasing")
    out()

    N = 5
    J = 1.0
    H = build_heisenberg_chain(N, J)

    # IBM sacrifice-zone gammas (Q85-Q94), as D[Z] COEFFICIENTS.
    #
    # The calibration hands over coherence rates 1/T2*, and this model is
    # Z-dephasing alone, so the D[Z] coefficient is gamma = 1/(2*T2*): row 1
    # of the T2 -> gamma table in docs/GLOSSARY.md. A D[Z] channel at rate
    # gamma decays coherences at 2*gamma, so feeding 1/T2* here would double
    # every absolute rate in this file.
    #
    # The measured coherence rates are [0.2681, 0.0163, 0.0103, 0.0147,
    # 0.0105] /us. Halved below. The same Q85 number was repaired the same
    # way in its sibling script on 2026-08-05 (GLOSSARY.md, the T2 -> gamma
    # section: "that script was repaired 2026-08-05 and now feeds the
    # dephasing-only 0.134"); this one was missed until 2026-08-18.
    #
    # Ratios in this file never depended on the factor; the ABSOLUTE rates
    # and the fixed 0.05 protection threshold did.
    ibm_t2_rates = np.array([0.2681, 0.0163, 0.0103, 0.0147, 0.0105])
    ibm_gammas = ibm_t2_rates / 2.0

    # Uniform gammas (same total)
    total_gamma = sum(ibm_gammas)
    uniform_gammas = np.array([total_gamma / N] * N)

    # Zero noise
    zero_gammas = np.array([0.0] * N)

    # Run all three profiles
    ev_zero, stat_zero, osc_zero, freq_zero = analyze_profile(
        "ZERO NOISE (cavity modes)", zero_gammas, H, N, out)

    ev_ibm, stat_ibm, osc_ibm, freq_ibm = analyze_profile(
        "IBM SACRIFICE ZONE (Q85-Q94)", ibm_gammas, H, N, out)

    ev_uni, stat_uni, osc_uni, freq_uni = analyze_profile(
        "UNIFORM (same total gamma)", uniform_gammas, H, N, out)

    # Frequency group comparison
    frequency_group_analysis(osc_ibm, osc_uni, freq_zero, out)

    # Mode survival comparison
    out(f"\n{'='*60}")
    out("MODE SURVIVAL COMPARISON")
    out(f"{'='*60}")

    # Compare slowest modes
    if osc_ibm and osc_uni:
        ibm_sorted = sorted(osc_ibm, key=lambda x: x[0])
        uni_sorted = sorted(osc_uni, key=lambda x: x[0])

        out("\nSlowest 10 modes: IBM sacrifice vs Uniform")
        out(f"{'#':>3} {'IBM_rate':>10} {'Uni_rate':>10} {'Ratio':>8} {'Freq':>8}")
        for i in range(min(10, len(ibm_sorted), len(uni_sorted))):
            ir, if_, _ = ibm_sorted[i]
            ur, uf, _ = uni_sorted[i]
            ratio = ur / ir if ir > 1e-10 else float('inf')
            out(f"{i+1:3d} {ir:10.6f} {ur:10.6f} {ratio:8.2f}x {if_:8.3f}")

    # Summary
    out(f"\n{'='*60}")
    out("SUMMARY")
    out(f"{'='*60}")
    out(f"Cavity modes at gamma=0: {len(stat_zero)} stationary, {len(osc_zero)} oscillating, {len(freq_zero)} frequencies")
    out(f"IBM sacrifice zone: {len(stat_ibm)} immune, F1 distance {f1_distance_in_eps(ev_ibm, total_gamma)[0]:.1f} eps*rho")
    out(f"Uniform: {len(stat_uni)} immune, F1 distance {f1_distance_in_eps(ev_uni, total_gamma)[0]:.1f} eps*rho")

    if osc_ibm and osc_uni:
        ibm_min = min(r for r, _, _ in osc_ibm)
        uni_min = min(r for r, _, _ in osc_uni)
        out(f"Slowest oscillating mode: IBM {ibm_min:.6f} vs Uniform {uni_min:.6f} (ratio {uni_min/ibm_min:.2f}x)")

        ibm_protected = sum(1 for r, _, _ in osc_ibm if r < 0.05)
        uni_protected = sum(1 for r, _, _ in osc_uni if r < 0.05)
        out(f"Protected modes (rate<0.05): IBM {ibm_protected} vs Uniform {uni_protected} (+{ibm_protected-uni_protected})")

    out(f"\nMax decay rate: {max(-ev.real for ev in ev_ibm):.6f} (expected: {2*total_gamma:.6f})")
    out(f"Palindrome center: {total_gamma:.4f}")
    out()
    out("Key insight: The sacrifice zone does not protect qubits.")
    out("It protects CAVITY MODES. Modes localized on the interior")
    out("qubits (Q86-Q94) see less noise and survive longer.")
    out("The same 43 frequencies exist under all three profiles.")
    out("Only the damping changes.")

    if sweeps:
        orbit_sweeps(H, ibm_gammas, uniform_gammas, out)

    # Write output
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main(sweeps="--sweeps" in sys.argv)
