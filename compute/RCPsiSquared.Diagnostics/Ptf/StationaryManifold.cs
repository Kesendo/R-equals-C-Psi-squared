using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>The closure guardian: count the stationary manifold of a Liouvillian and forecast
/// whether the PTF closure survives.
///
/// <para><b>The law (Edge 5, 2026-06-10)</b>: the PTF closure, the one-parameter time-rescaling
/// fit P_B(i, t) ≈ P_A(i, α_i·t) with Σ_i ln α_i ≈ 0 (the EQ-014 observable family), holds iff
/// the dynamics retains its COMPLETE set of N+1 stationary sector projectors, where the count
/// may be taken in any frame. For the uniform XY chain under Z-dephasing the baseline manifold
/// is the N+1 excitation-sector projectors (F4); a perturbation that lifts any of them off
/// λ = 0 kills the closure.</para>
///
/// <para><b>Frame independence</b>: a unitary frame change ρ → UρU† conjugates the Liouvillian
/// by the unitary U ⊗ U* on operator space, which preserves the spectrum; hence the stationary
/// count is the same in every frame. Hidden-frame conservation laws are therefore captured by
/// the bare count with no frame search: the pure XY+YX chain has a sublattice-X-flip frame that
/// restores U(1) while commuting with the Z-dephasers, and its bare count is already the full
/// N+1 (the tour-protocol hold). Conversely, a count below N+1 rules out a hidden full manifold
/// in ANY frame, not just the scanned ones.</para>
///
/// <para><b>Partial survival buys nothing</b>: the forecast is the EQUALITY test
/// (count == N+1), not a "how many projectors are left" heuristic. The XY+YX perturbation keeps
/// the two parity projectors yet breaks the closure fastest of all five tour bilinears
/// (canonical protocol, S = +11.3 vs J-defect control +0.097). The weaker-guardian hypothesis
/// (parity survival predicts the hold) was tested and REFUTED by the Edge-5 probe: parity
/// scored 5/6 on the canonical sweep, the full-manifold count 6/6.</para>
///
/// <para><b>Scope, stated honestly</b>: the law is established for the canonical PTF protocol
/// (uniform XY chain, Z-dephasing, bonding-mode initial state) and the 2026-05-16 multi-lens
/// tour cases at N = 3..5. It is a verified empirical law, not a derived theorem; no typed
/// Claim is created. This is a compute primitive of the same standing as
/// <see cref="SlowModeMixing"/>.</para>
///
/// <para><b>Sources</b>: probe <c>simulations/ptf_carrier_seam_retrodiction.py</c> (Edge 5 of
/// the 2026-06-10 PTF fresh-eyes chain); experiment
/// <c>experiments/PTF_PALINDROME_BREAKING_PERTURBATIONS.md</c> ("Edge 5 result"); reflection
/// <c>reflections/ON_HOW_THE_CARRIER_SHOWS_ITSELF.md</c> ("Seen again 2026-06-10"); the PTF
/// α-closure observable itself is EQ-014.</para></summary>
public static class StationaryManifold
{
    /// <summary>The baseline stationary tolerance. At N = 4, γ = 0.05 the gap is wide:
    /// surviving projectors sit at |λ| ≲ 2·10⁻¹⁵, lifted ones at |λ| ≳ 3·10⁻³.</summary>
    public const double DefaultTolerance = 1e-10;

    /// <summary>The forecast: the measured stationary count, the full-manifold size N+1, and
    /// the closure verdict. <see cref="PredictsClosureHolds"/> is the equality test
    /// (count == full manifold); any deficit, however small, predicts a break (partial
    /// survival buys nothing, see the class doc).</summary>
    public readonly record struct ClosureForecast(int StationaryCount, int FullManifold)
    {
        /// <summary>True iff the complete stationary manifold survives. This is the law's
        /// verdict, not a heuristic score.</summary>
        public bool PredictsClosureHolds => StationaryCount == FullManifold;
    }

    /// <summary>Count the stationary modes of <paramref name="liouvillian"/>: the number of
    /// eigenvalues with |λ| &lt; <paramref name="tolerance"/>. Frame-independent (the spectrum
    /// is invariant under unitary frame changes).</summary>
    public static int Count(ComplexMatrix liouvillian, double tolerance = DefaultTolerance)
    {
        var evd = liouvillian.Evd();
        int count = 0;
        foreach (var lambda in evd.EigenValues)
            if (lambda.Magnitude < tolerance) count++;
        return count;
    }

    /// <summary>Forecast the PTF closure for an N-site system: measure the stationary count of
    /// <paramref name="liouvillian"/> and compare it with the full manifold N+1 (the F4
    /// excitation-sector projector count of the unperturbed chain).</summary>
    public static ClosureForecast Forecast(int n, ComplexMatrix liouvillian, double tolerance = DefaultTolerance) =>
        new(Count(liouvillian, tolerance), FullManifold: n + 1);
}
