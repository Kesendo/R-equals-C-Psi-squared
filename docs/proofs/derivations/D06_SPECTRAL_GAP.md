# D6: Spectral Gap and Mixing Time

**What this derivation is about:** The spectral gap (the decay rate of the slowest non-stationary mode, which sets the timescale for the system to forget its initial state) is exactly 2γ in the strong-coupling regime, and the mixing time then scales as N·ln(4)/(2γ), growing linearly with system size. The gap is NOT independent of the coupling: below an N-dependent threshold in Q = J/γ it is Zeno-suppressed and far smaller.

**Source formulas:** 3 (decay rate bounds); F1 plays no part in the derivation below
**Tier:** 1 above Q*_gap(N). Note the provenance is empirical: Q*_gap(N) is bisected
numerically, and no lower bound on min{<n_XY> > 0} is derived anywhere
**Status:** VERIFIED for Q > Q*_gap(N) (N=2-5, deviation < 1e-14); FAILS below Q*_gap

## Derivation

From F3: above the coupling threshold Q*_gap(N) the minimum nonzero decay rate is
2*gamma (w=1 modes).

    Spectral gap = min{|Re(lambda)| : lambda != 0} = 2*gamma      (Q > Q*_gap(N))

Note what does NOT establish this. F50 counts the modes sitting AT 2*gamma; that
is existence, not minimality. The minimum is a separate fact, and it is the one
that fails below threshold.

Measured Q*_gap(N) by bisection at gamma = 0.05: 0.5000 (N=2), 0.8002 (N=3),
1.3422 (N=4), 1.8194 (N=5). Q*_gap is a function of the Hamiltonian as well: the
same bisection on the XY chain gives 0.7071 (N=3), 0.9393 (N=4), 1.1861 (N=5).
Q*_gap is independent of gamma, since gap/gamma depends on Q alone.

As Q -> 0 the gap approaches

    Spectral gap -> 2*(1 - cos(pi/N)) * 2*J^2/gamma               (Q -> 0)

Zeno suppression: strong dephasing freezes the transport that would give a mode
its light. At N=3, gamma=0.3, J=0.001 the gap is 1.1e-5 of 2*gamma. This is an
asymptote, not a second regime covering (0, Q*_gap). It is good to about 1% only for
Q <~ 0.1, and near Q*_gap it fails badly and not even in one direction: at 98% of
Q*_gap the measured/predicted ratio is 1.67 (N=2), 1.43 (N=3), 0.91 (N=4), 0.77
(N=5), so the asymptote undershoots at small N and overshoots from N=4 on. In
between neither closed form applies and the gap is a number one computes.

Above threshold the gap determines:
- The slowest decoherence timescale: tau_slow = 1/(2*gamma)
- The mixing time (convergence to steady state):

    t_mix <= (1/gap) * ln(dim) = ln(4^N) / (2*gamma) = N*ln(4) / (2*gamma)

Two caveats on that line. There is no steady state to converge to: the kernel is
(N+1)-dimensional (F4's kernel-dimension extension; proof:
[kernel dimension by components](../PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md)),
so the flow converges to a manifold. And the bound is
quoted, not verified: `verify_derivations.py` computes t_mix and prints it
without comparing it to anything, so D6's VERIFIED status rests on the gap
alone.

Above Q*_gap(N), and within the number-conserving XY/Heisenberg family, the gap is
2*gamma whichever member of that family and whichever topology. The THRESHOLD is
not H-independent (see the XY values above), and off that family the value fails
outright: a generic real symmetric H has a one-dimensional kernel and no
2*gamma level at all.

## Numerical verification

| N | Gap (numerical) | Gap (formula) | Deviation |
|---|----------------|--------------|-----------|
| 2 | 0.1000000000   | 0.1000000000 | < 4e-16   |
| 3 | 0.1000000000   | 0.1000000000 | < 2e-15   |
| 4 | 0.1000000000   | 0.1000000000 | < 6e-15   |
| 5 | 0.1000000000   | 0.1000000000 | < 1e-14   |

(gamma = 0.05, J = 1.0 for all tests, i.e. Q = 20, far above every Q*_gap(N). The
same script at the repository's canonical Q = 1.5 reports "VERIFIED for N in
[2, 3, 4] (N in [5] below threshold, gap != 2*gamma there)", and at Q = 0.4
"NOT APPLICABLE (Q below every threshold)". The 1e-14 deviations above are a
strong-coupling result, not a universal one.)

Script: [`simulations/verify_derivations.py`](../../../simulations/verify_derivations.py)

## Replaces

Numerical search for the smallest nonzero eigenvalue, provided Q > Q*_gap(N). In
that regime mixing time estimates need only N and gamma; below it they need J
as well, and the bound inverts the gap and grows with it.

See [the Absorption Theorem proof](../PROOF_ABSORPTION_THEOREM.md) section 4.3
and [`simulations/absorption_ladder_regimes.py`](../../../simulations/absorption_ladder_regimes.py).
