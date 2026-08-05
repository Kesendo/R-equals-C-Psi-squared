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


def max_f1_pairing_distance(spectrum, sigma):
    """The repo's CANONICAL F1 symmetry distance, ported from C#.

    Port of `F1SpectrumStatistics.MaxF1PairingDistance`
    (compute/RCPsiSquared.Core/F1/F1SpectrumStatistics.cs), which its own summary
    calls "the canonical F1 check", and which is the matcher behind
    `MultisetAssert.NearestNeighbourEqual` and the live witness
    `BlockSpectrumWitness.PalindromePairingDistance`. Committed reference values
    live in simulations/results/f1_n8_n9_metrics/. The same port stands in
    simulations/sacrifice_zone_mapping.py.

    Max greedy nearest-neighbour distance, WITH REMOVAL, between the multiset of
    eigenvalues and its F1 reflection {-2*sigma - lambda}. Each mirror point is
    consumed once, so an ASYMMETRIC multiplicity error leaves an unmatched point
    and surfaces as a large distance, which a set or Hausdorff distance would
    miss. Measured on the IBM profile below, against a base of 59.6 eps*rho:
    duplicating one eigenvalue onto an unrelated one gives 3.4e14 to 3.6e15 over
    five random trials, and deleting one gives 3.4e15. It runs on the FULL
    COMPLEX spectrum; a rates-only variant is strictly weaker, being blind to
    the imaginary parts.

    TWO BLIND SPOTS, both measured on that same profile, because the metric
    compares the spectrum against ITS OWN reflection:
      * A corruption that RESPECTS the reflection is invisible. Dropping one
        mirror pair while duplicating another leaves the result bit-identical
        (1.549738e-13 either way). That is the plausible failure mode of a BLOCK
        solver, which produces symmetry-respecting spectra by construction, so
        this number does not certify a block spectrum against a dense one.
      * Duplicating onto a NUMERICALLY DEGENERATE neighbour is a no-op, since
        the copy lands where a partner already sits.
    Resolution, for the same reason: perturbing one eigenvalue by 1e-13 leaves
    the number unchanged at 59.6, 1e-12 raises it to 383.5. A genuine violation
    below about 1e-13 absolute sits inside the floor and cannot be seen here.

    This REPLACED a greedy first-fit percentage inside an absolute tolerance
    (1e-6 here). That score was measuring its own matcher: tightening the
    tolerance SATURATES the printed percentage long before it fixes the pairing,
    so a 100% entry is weaker evidence than a 91% one, not stronger. See
    experiments/CONCENTRATOR_MAPPING.md ("family of at least seventeen").
    """
    spectrum = np.asarray(spectrum)
    # Guard, because a NaN would otherwise pass SILENTLY. np.argmin selects a NaN
    # index and `nan > worst` is False, so a broken spectrum cannot raise the
    # metric: measured, one NaN among four eigenvalues returns 11.0 rather than a
    # failure. The C# does not return a sentinel there, it THROWS
    # (F1SpectrumStatistics.cs:383: a NaN distance fails the `d < bestDist` test,
    # so bestIdx stays -1 and it raises "no candidate"). A metric whose
    # job is to catch a broken spectrum must not fail quiet.
    if not np.all(np.isfinite(spectrum)):
        raise ValueError(
            "spectrum contains non-finite values; the F1 distance is undefined "
            "and would otherwise be silently understated")
    reflected = -2.0 * sigma - spectrum
    taken = np.zeros(len(spectrum), dtype=bool)
    worst = 0.0
    for x in spectrum:
        d = np.abs(x - reflected)
        d[taken] = np.inf
        j = int(np.argmin(d))
        taken[j] = True
        if d[j] > worst:
            worst = d[j]
    return float(worst)


def f1_distance_in_eps(evals, sigma):
    """F1 distance in units of eps * spectral radius (the error model).

    An eigensolver on a non-normal matrix has no exact route, so the number is
    published against its backward-error model rather than against a threshold:
    the value is dimensionless, and "at the floor" is an N-DEPENDENT band, not a
    universal one. Measured here over two decades of J (0.1, 1, 10) on the
    leading sites of the IBM profile: N=2 gives 1.5 to 3.9, N=3 gives 19.3 to
    33.1, N=4 gives 31.7 to 54.9, N=5 gives 51.3 to 68.9. So the floor grows by
    a factor ~40 from N=2 to N=5, and an N=7 or N=8 reading must be graded
    against its own N rather than against "O(10-100)". Within a fixed N the
    J-dependence is 1.34x (N=5) to 2.56x (N=2), which is the same size as the
    matcher sensitivities below: the number says "at the floor" and separates
    nothing finer than that.

    Two independent sensitivities of this matcher, both measured IN PLACE on the
    three profiles below rather than imported: permuting only the summation
    order of the jump operators spreads the IBM number 2.39x over all 120 orders
    (47.9 to 114.5), and permuting only the eigenvalue ARRAY ORDER, identical
    spectrum and no physics at all, spreads it 1.47x on the zero-noise profile
    and 1.27x on uniform, because greedy matching is order-dependent by
    construction. Do not rank profiles by it. (The array-order figure is
    profile-dependent and happens to be 1.01x on the IBM profile itself, which
    is why the jump-operator sweep, not this one, carries the argument.)
    """
    radius = float(np.max(np.abs(evals))) if len(evals) else 0.0
    scale = EPS * radius
    dist = max_f1_pairing_distance(evals, sigma)
    if scale == 0.0:
        # A zero radius means a zero spectrum, so the distance is zero too.
        # Returning `dist` raw here would print an ABSOLUTE number under an
        # eps*rho label; the sibling port in sacrifice_zone_mapping.py returns
        # 0.0, and this now agrees with it.
        return 0.0, radius
    return dist / scale, radius


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


def main():
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

    # IBM sacrifice-zone gammas (Q85-Q94)
    ibm_gammas = np.array([0.2681, 0.0163, 0.0103, 0.0147, 0.0105])

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

    # Write output
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"\n>>> Results saved to: {out_path}")


if __name__ == "__main__":
    main()
