#!/usr/bin/env python3
"""Finite thermal-channel census; no blackbody or wave mechanism is inferred."""

from pathlib import Path
import os
import sys
import numpy as np

if sys.platform == "win32":
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except Exception:
        pass

RESULTS_DIR = Path(__file__).parent / "results"
J, GAMMA, TOL_FREQ = 1.0, 0.05, 1e-6
I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
Sm = np.array([[0, 1], [0, 0]], dtype=complex)  # sigma- = |0><1|: emission
Sp = np.array([[0, 0], [1, 0]], dtype=complex)  # sigma+ = |1><0|: absorption


def kron_chain(ops):
    result = ops[0]
    for operator in ops[1:]:
        result = np.kron(result, operator)
    return result


def dissipator_term(jump, dimension):
    identity = np.eye(dimension, dtype=complex)
    norm = jump.conj().T @ jump
    return np.kron(jump, jump.conj()) - 0.5 * (
        np.kron(norm, identity) + np.kron(identity, norm.T)
    )


def build_thermal_liouvillian(n, gammas, n_bar, gamma_thermal=None):
    """Build Z dephasing plus emission/absorption for n_bar >= 0."""
    if not np.isfinite(n_bar) or n_bar < 0:
        raise ValueError("n_bar must be finite and non-negative")
    if gamma_thermal is None:
        gamma_thermal = gammas
    if len(gammas) != n or len(gamma_thermal) != n:
        raise ValueError("one dephasing and thermal rate is required per site")
    gammas = np.asarray(gammas, dtype=float)
    gamma_thermal = np.asarray(gamma_thermal, dtype=float)
    if not np.all(np.isfinite(gammas)) or np.any(gammas < 0):
        raise ValueError("dephasing rates must be finite and non-negative")
    if not np.all(np.isfinite(gamma_thermal)) or np.any(gamma_thermal < 0):
        raise ValueError("thermal rates must be finite and non-negative")
    dimension = 2**n
    identity = np.eye(dimension, dtype=complex)
    hamiltonian = np.zeros((dimension, dimension), dtype=complex)
    for site in range(n - 1):
        for pauli in (Xm, Ym, Zm):
            ops = [I2] * n
            ops[site], ops[site + 1] = pauli, pauli
            hamiltonian += J * kron_chain(ops)
    generator = -1j * (np.kron(hamiltonian, identity) - np.kron(identity, hamiltonian.T))
    for site in range(n):
        ops = [I2] * n
        ops[site] = Zm
        generator += dissipator_term(np.sqrt(gammas[site]) * kron_chain(ops), dimension)
        ops = [I2] * n
        ops[site] = Sm
        generator += dissipator_term(
            np.sqrt(gamma_thermal[site] * (n_bar + 1.0)) * kron_chain(ops), dimension
        )
        if n_bar > 0:
            ops = [I2] * n
            ops[site] = Sp
            generator += dissipator_term(
                np.sqrt(gamma_thermal[site] * n_bar) * kron_chain(ops), dimension
            )
    return generator


def spectral_row(values):
    oscillating = np.abs(values.imag) > TOL_FREQ
    rates = -values.real
    q = np.abs(values.imag[oscillating]) / np.maximum(rates[oscillating], 1e-15)
    return int(np.count_nonzero(oscillating)), float(np.mean(rates)), float(np.max(q))


def zero_temperature_emission_gate(generator, expected_rate):
    """Pin signed |1> -> |0> population flow for a one-site cold bath."""
    excited = np.array([[0.0, 0.0], [0.0, 1.0]], dtype=complex).reshape(-1)
    ground = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex).reshape(-1)
    excited_flow = (generator @ excited).reshape(2, 2)
    ground_flow = (generator @ ground).reshape(2, 2)
    signed_flow = {
        "ground": float(excited_flow[0, 0].real),
        "excited": float(excited_flow[1, 1].real),
        "cold_ground_motion": float(np.linalg.norm(ground_flow)),
    }
    tolerance = 64 * np.finfo(float).eps * max(1.0, abs(expected_rate))
    if (
        abs(signed_flow["ground"] - expected_rate) > tolerance
        or abs(signed_flow["excited"] + expected_rate) > tolerance
        or signed_flow["cold_ground_motion"] > tolerance
    ):
        raise RuntimeError("zero-temperature emission direction gate failed")
    return signed_flow


def planck_1d(omega, amplitude, t_eff):
    """1D cavity Planck spectral density, omega / (exp(omega/T) - 1)."""
    x = omega / max(t_eff, 1e-10)
    return amplitude * omega / (np.exp(np.clip(x, -50, 50)) - 1 + 1e-30)


def exponential(omega, amplitude, scale):
    return amplitude * np.exp(-omega / max(scale, 1e-10))


def r_squared(observed, predicted):
    ss_res = float(np.sum((observed - predicted) ** 2))
    ss_tot = float(np.sum((observed - np.mean(observed)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


def blackbody_falsification(log, n, gammas):
    """The two tests the blackbody reading would have to pass, and does not.

    Neither is a scope statement: both are fits to the committed spectrum, run
    against a criterion fixed before the numbers are read. R2_MIN is the bar a
    curve must clear to be called a fit at all; ALPHA_TOL is how far the
    absorption exponent may sit from the value the reading predicts.
    """
    from scipy.optimize import curve_fit

    R2_MIN = 0.90
    ALPHA_TOL = 0.5

    log("=" * 75)
    log("PLANCK FIT TO THE FREQUENCY DISTRIBUTION")
    log(f"A curve counts as a fit only at R^2 >= {R2_MIN:.2f}.")
    log("=" * 75)
    log()

    planck_verdicts = []
    for n_bar in (0.1, 0.5, 1.0, 5.0, 10.0):
        values = np.linalg.eigvals(build_thermal_liouvillian(n, gammas, n_bar, gammas))
        freqs = np.abs(values.imag[np.abs(values.imag) > TOL_FREQ])
        n_bins = min(20, len(freqs) // 3)
        if n_bins < 3:
            log(f"nbar={n_bar}: {len(freqs)} modes, too few to bin")
            continue
        counts, edges = np.histogram(freqs, bins=n_bins)
        centres = (edges[:-1] + edges[1:]) / 2
        mask = counts > 0
        observed = counts[mask].astype(float)

        try:
            popt, _ = curve_fit(planck_1d, centres[mask], observed,
                                p0=[float(max(counts)), 1.0], maxfev=5000)
            r2_planck = r_squared(observed, planck_1d(centres[mask], *popt))
            t_eff = float(popt[1])
        except Exception as error:
            raise RuntimeError(f"the Planck fit did not converge at nbar={n_bar}: {error}")
        try:
            popt_e, _ = curve_fit(exponential, centres[mask], observed,
                                  p0=[float(max(counts)), 1.0], maxfev=5000)
            r2_exp = r_squared(observed, exponential(centres[mask], *popt_e))
        except Exception as error:
            raise RuntimeError(f"the exponential fit did not converge at nbar={n_bar}: {error}")
        # A fit that failed to converge must not read like a measured non-fit.
        if not (np.isfinite(r2_planck) and np.isfinite(r2_exp)):
            raise RuntimeError(
                f"a fit returned a non-finite R^2 at nbar={n_bar}: "
                f"Planck {r2_planck}, exponential {r2_exp}")

        fits = r2_planck >= R2_MIN
        planck_verdicts.append((n_bar, fits, r2_planck, r2_exp))
        log(f"nbar={n_bar}: {len(freqs)} modes, {n_bins} bins")
        log(f"  Planck:      T_eff = {t_eff:.4f}, R^2 = {r2_planck:.4f}")
        log(f"  Exponential:                 R^2 = {r2_exp:.4f}")
        log(f"  Planck is a fit at R^2 >= {R2_MIN:.2f}: {'yes' if fits else 'no'}")
        log()

    worst_r2 = max(max(p[2], p[3]) for p in planck_verdicts)
    widest_gap = max(abs(p[2] - p[3]) for p in planck_verdicts)
    log(f"  The best R^2 either curve reaches is {worst_r2:.4f}, against the {R2_MIN:.2f} bar.")
    log(f"  The widest Planck-minus-exponential gap is {widest_gap:.4f}, which is smaller than")
    log(f"  the distance {R2_MIN - worst_r2:.4f} from the better curve to being a fit at all, so")
    log("  ranking one above the other reads a preference out of noise.")
    if any(fits for _, fits, _, _ in planck_verdicts):
        raise RuntimeError(
            "a Planck fit cleared the bar: "
            + ", ".join(f"nbar={nb} R^2={r2:.4f}" for nb, fits, r2, _ in planck_verdicts if fits))
    if widest_gap >= R2_MIN - worst_r2:
        raise RuntimeError(
            f"the Planck-exponential gap {widest_gap:.4f} is no longer small against the "
            f"distance {R2_MIN - worst_r2:.4f} to a fit; the ranking sentence no longer holds")
    log()

    log("=" * 75)
    log("ABSORPTION SCALING (the Stefan-Boltzmann test)")
    log(f"The reading predicts alpha = 4 (or 2 in 1D), accepted within {ALPHA_TOL}.")
    log("=" * 75)
    log()

    sweep = [0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0]
    baseline = float(np.mean(-np.linalg.eigvals(
        build_thermal_liouvillian(n, gammas, 0.0, gammas)).real))
    log(f"{'nbar':>8} {'mean rate':>11} {'excess':>10}")
    excesses = []
    for n_bar in sweep:
        rate = float(np.mean(-np.linalg.eigvals(
            build_thermal_liouvillian(n, gammas, n_bar, gammas)).real))
        excesses.append(rate - baseline)
        log(f"{n_bar:8.3f} {rate:11.4f} {rate - baseline:10.4f}")

    excesses = np.array(excesses)
    alpha, log_amp = np.polyfit(np.log(np.array(sweep)), np.log(excesses), 1)
    log()
    log(f"  Power law: excess = {np.exp(log_amp):.4f} * nbar^{alpha:.3f}")
    log(f"  Stefan-Boltzmann would give alpha = 4, a 1D Planck law alpha = 2.")
    log(f"  Within {ALPHA_TOL} of 4: {'yes' if abs(alpha - 4) <= ALPHA_TOL else 'no'}; "
        f"of 2: {'yes' if abs(alpha - 2) <= ALPHA_TOL else 'no'}")
    log()

    # The exponent is 1 because the law behind it is exact, so it is checked
    # exactly rather than fitted. The check does not go through the eigensolver:
    # the eigenvalues sum to the trace, so the mean rate is -Re(tr L)/d^2 and the
    # whole identity is arithmetic on the generator. That route can return 0.0,
    # and does; the fitted alpha above would have hidden an exact law behind a
    # regression, and a tolerance here would have hidden nothing but itself.
    dimension_squared = float((2 ** n) ** 2)
    def mean_rate_from_trace(n_bar):
        generator = build_thermal_liouvillian(n, gammas, n_bar, gammas)
        return -float(np.trace(generator).real) / dimension_squared

    def excess_deviation(rate, points):
        baseline = mean_rate_from_trace_at(rate, 0.0)
        return float(max(
            abs((mean_rate_from_trace_at(rate, nb) - baseline) - n * rate * nb)
            for nb in points))

    def mean_rate_from_trace_at(rate, n_bar):
        generator = build_thermal_liouvillian(n, [rate] * n, n_bar, [rate] * n)
        return -float(np.trace(generator).real) / dimension_squared

    # The identity is exact, so it is asserted exactly, on inputs that can break
    # it. The generator builds each jump operator as sqrt(rate*(nbar+1)) and the
    # dissipator squares it back, so what survives is that roundtrip and nothing
    # physical: at a dyadic rate the roundtrip is exact and the deviation is
    # 0.0 to the bit, while the census rate 0.05 carries a couple of ulp. The
    # gate runs on the dyadic rate, where a real error cannot hide behind
    # rounding, and the census residual is reported rather than tolerated.
    dyadic_rate, dyadic_points = 0.25, (0.25, 0.5, 1.0, 2.0, 4.0)
    dyadic_worst = excess_deviation(dyadic_rate, dyadic_points)
    census_worst = excess_deviation(gammas[0], sweep)
    log(f"  Exact law: excess = N*gamma_T*nbar, read off the trace, not the spectrum.")
    log(f"  At the dyadic rate gamma_T={dyadic_rate}: deviation {dyadic_worst:.1f} (exact).")
    log(f"  At the census rate gamma_T={gammas[0]}: deviation {census_worst:.3e},")
    log("  which is the sqrt-then-square roundtrip in the jump operators, not the law.")
    if dyadic_worst != 0.0:
        raise RuntimeError(
            "absorption excess is not exactly N*gamma_T*nbar at a dyadic rate; "
            f"worst deviation {dyadic_worst:.3e}")
    log(f"  The excess goes as nbar^{alpha:.3f}, not as the 4 a Stefan-Boltzmann reading")
    log("  predicts nor the 2 a 1D Planck law does. The blackbody reading asks for a power")
    log("  the generator does not carry.")
    if abs(alpha - 4) <= ALPHA_TOL or abs(alpha - 2) <= ALPHA_TOL:
        raise RuntimeError(
            f"the absorption exponent {alpha:.3f} now matches a blackbody prediction; "
            "the falsification above no longer holds")
    log()
    log("  Both arms have now been checked and each raises if it changes its mind. This")
    log("  is the numerical half of the scope statement: the spectrum is not a radiation")
    log("  intensity distribution, and it also does not behave like one when asked.")
    log()
    return planck_verdicts, alpha


def render_report():
    lines = []
    def log(message=""):
        lines.append(message)

    n = 4
    gammas = [GAMMA] * n
    log("THERMAL-CHANNEL NUMERICAL CENSUS")
    log("sigma- = |0><1| emission; sigma+ = |1><0| absorption")
    log("L(nbar)=L_H+L_Z+gamma_T[(nbar+1)D[sigma-]+nbar D[sigma+]]")
    log(f"N={n}, J={J}, gamma_Z=gamma_T={GAMMA}, imag tolerance={TOL_FREQ:.0e}")
    log("Counts are numerical at the stated tolerance; they are not EP certificates.")
    log()
    log(f"{'nbar':>10} {'osc':>6} {'osc%':>8} {'mean decay':>12} {'Qmax':>10}")
    for n_bar in (0.0, 1e-9, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0):
        values = np.linalg.eigvals(build_thermal_liouvillian(n, gammas, n_bar, gammas))
        osc, mean_decay, qmax = spectral_row(values)
        log(f"{n_bar:10.3g} {osc:6d} {100*osc/len(values):7.3f}% {mean_decay:12.6f} {qmax:10.3f}")
    cold = build_thermal_liouvillian(1, [0.0], 0.0, [GAMMA])
    near = build_thermal_liouvillian(1, [0.0], 1e-9, [GAMMA])
    cold_flow = zero_temperature_emission_gate(cold, GAMMA)
    log()
    log(f"negative control: ||L(0)||_F={np.linalg.norm(cold):.6e} > 0 (spontaneous emission present)")
    log(
        "signed cold-bath gate: excited |1><1| gives "
        f"dP0/dt={cold_flow['ground']:+.6e}, "
        f"dP1/dt={cold_flow['excited']:+.6e}; "
        f"ground-state motion={cold_flow['cold_ground_motion']:.2e}"
    )
    log(f"continuity: ||L(1e-9)-L(0)||_F={np.linalg.norm(near-cold):.6e}")
    log()
    log()
    blackbody_falsification(log, n, gammas)
    log("Scope: finite spectral counts only; no thermal-photon, standing-wave, or phase-transition mechanism is inferred.")
    return "\n".join(lines) + "\n"


def main():
    output = render_report()
    print(output, end="")
    path = RESULTS_DIR / "thermal_blackbody.txt"
    path.write_text(output, encoding="utf-8")
    print(f"\nResults saved to: {path}")


if __name__ == "__main__":
    main()
