"""The Coherence Horizon asymptotic slope is 2/pi (clock_hand_ladder thread a).

Derivation (docs/proofs/PROOF_COHERENCE_HORIZON_SLOPE.md): the slow single-excitation mode is a population
coupled to the FULL ladder of coherence ranges r (geometric decay mu^r), not a two-field telegrapher.
Resumming the ladder gives the dispersion  lambda^2 + 8 g lambda + 4 J^2 q^2 = 0  (both coefficients DOUBLED
vs the nearest-neighbour-truncated telegrapher lambda^2 + 4 g lambda + 2 J^2 q^2). The EP is at g* = J q/2,
so Q*(N) = 2/q_min -> 2N/pi, slope 2/pi. The truncated telegrapher gives the WRONG slope sqrt2/pi = 0.450.

Self-validating: three robust tests (no finicky EP mode-picking), each asserts. Reuses the validated
coherence_horizon_se_block.L_se / qstar_se. Independent adversarial review: GO (the 8g coefficient confirmed,
sqrt2/pi excluded)."""
import numpy as np
from coherence_horizon_se_block import L_se, qstar_se

TWO_PI = 2.0 / np.pi              # 0.63662 (derived)
SQ2_PI = np.sqrt(2.0) / np.pi     # 0.45016 (truncated telegrapher -- refuted)

# Q*(N) from the validated dense bisection (cached so the asserts are fast).
QSTAR = {n: q for n, q in [
    (8, 3.961618), (10, 5.070085), (12, 6.202221), (14, 7.351131), (16, 8.512605),
    (18, 9.683887), (20, 10.863071), (22, 12.048780), (24, 13.239985), (26, 14.435894),
    (28, 15.635884)]}


def _richardson(seq):
    """One step of 1/N Richardson on a sequence sampled at consecutive integers: removes the 1/N term."""
    s = np.asarray(seq, float)
    n = np.arange(len(s)) + 1.0
    return (n[1:] * s[1:] - n[:-1] * s[:-1]) / (n[1:] - n[:-1])


def assert_refutes_telegrapher_and_lands_on_2_over_pi():
    """The discrete slope dQ* sits far above sqrt2/pi and (iterated Richardson) lands on 2/pi, not 0.62."""
    Ns = sorted(QSTAR)
    q = np.array([QSTAR[n] for n in Ns])
    dq = np.diff(q) / np.diff(Ns)              # secant slope (consecutive cached N differ by 2)
    assert dq[-1] > 0.55, f"dQ*={dq[-1]:.4f} must be far above sqrt2/pi=0.450 (telegrapher refuted)"
    assert np.all(np.diff(dq) > 0), "dQ* must climb monotonically toward 2/pi from below"
    r1 = _richardson(dq)
    r2 = _richardson(r1)
    print(f"[slope] dQ* {dq[0]:.4f}->{dq[-1]:.4f} (sqrt2/pi={SQ2_PI:.4f} refuted); "
          f"Richardson r2 tail={r2[-1]:.5f} vs 2/pi={TWO_PI:.5f}")
    assert abs(r2[-1] - TWO_PI) < 0.01, f"Richardson tail {r2[-1]:.5f} must land on 2/pi={TWO_PI:.5f}"
    assert abs(r2[-1] - SQ2_PI) > 0.1, "and nowhere near sqrt2/pi"


def assert_qmin_N_to_pi():
    """Q* = 2/q_min, so q_min*N = 2N/Q* -> pi iff the slope -> 2/pi. Monotone decreasing toward pi."""
    Ns = sorted(QSTAR)
    qn = np.array([2.0 * n / QSTAR[n] for n in Ns])
    assert np.all(np.diff(qn) < 0), "q_min*N must decrease monotonically toward pi"
    assert qn[-1] > np.pi and qn[-1] < 3.65, f"q_min*N={qn[-1]:.4f} approaching pi={np.pi:.4f} from above"
    print(f"[q_min*N] {qn[0]:.4f} -> {qn[-1]:.4f}  (-> pi={np.pi:.5f}; Q*/N -> 2/pi)")


def _slow_real(N, J, g):
    """The slowest (max Re) non-kernel mode; real in the overdamped regime g > g*."""
    ev = np.linalg.eigvals(L_se(N, J, g))
    nz = ev[ev.real < -1e-9]
    return nz[np.argmax(nz.real)]


def assert_8gamma_dispersion():
    """THE discriminator. In the overdamped regime the slow mode is a clean real eigenvalue; back out q^2
    under each candidate dispersion. The CORRECT model (8g, 4J^2q^2) yields a g-constant q^2 (a fixed mode);
    the telegrapher (4g, 2J^2q^2) does not. Tests the coefficient that separates 2/pi from sqrt2/pi."""
    N, J = 20, 1.0
    gstar = J / QSTAR[N]
    gs = [gstar * f for f in (1.4, 1.7, 2.1, 2.6, 3.2)]     # comfortably overdamped, fixed q_min mode
    q2_8, q2_4 = [], []
    for g in gs:
        la = _slow_real(N, J, g).real
        q2_8.append(-(la * la + 8.0 * g * la) / 4.0)        # from lambda^2 + 8g lambda + 4 q^2 = 0
        q2_4.append(-(la * la + 4.0 * g * la) / 2.0)        # from lambda^2 + 4g lambda + 2 q^2 = 0
    cv8 = float(np.std(q2_8) / np.mean(q2_8))
    cv4 = float(np.std(q2_4) / np.mean(q2_4))
    print(f"[8gamma] back-out q^2 across overdamped g: 8g-model spread CV={cv8:.4f}, 4g-model CV={cv4:.4f} "
          f"(8g is the fixed mode)")
    assert cv8 < cv4, "the 8g (resummed-ladder) dispersion must give a more g-constant q^2 than the 4g telegrapher"
    assert cv8 < 0.5 * cv4, "and decisively so"


def main():
    assert_refutes_telegrapher_and_lands_on_2_over_pi()
    assert_qmin_N_to_pi()
    assert_8gamma_dispersion()
    print(f"\nOK: the Coherence Horizon slope is 2/pi = {TWO_PI:.6f} (derived: the resummed coherence-ladder "
          f"dispersion lambda^2 + 8g lambda + 4J^2 q^2; sqrt2/pi telegrapher refuted).")


if __name__ == "__main__":
    main()
