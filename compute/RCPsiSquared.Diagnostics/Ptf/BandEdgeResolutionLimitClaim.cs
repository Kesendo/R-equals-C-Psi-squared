using System;
using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Diagnostics.Foundation;
using Witness = RCPsiSquared.Diagnostics.Foundation.BandEdgeResolutionLimitWitness;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>The defect-localization RESOLUTION LIMIT of F124 (the optics/signal reading of the band-edge
/// transition invariant's conditioning). F124's transition matrix M[b,k]=⟨ψ_k|V_b|ψ_1⟩ is the forward map of
/// a bond-recovery inverse problem; read that way, F124's λ_min = E = σ_min²(M) is the worst-case reconstruction
/// floor (the lower frame bound), and the condition number sets the resolution:
///
/// <para>  κ = λ_max/λ_min ~ N²  (the noise amplification),  contrast σ_max/σ_min = √κ ~ N,</para>
///
/// so a staggered (−1)^b zone-boundary (q=π) bond defect is √κ ~ N times harder to localize than a band-edge
/// defect of the same magnitude (matched-filter detection SNR). The worst-conditioned bond direction is exactly
/// the staggered mode (F124's λ_min eigenvector), the q=π detail at the resolution cutoff (the optician's
/// diffraction limit); the reconstruction floor vanishes as σ_min = √E ~ (N+1)^(−3/2) (E·(N+1)³ → 4π²), so the
/// long chain goes singular along the K-direction.
///
/// <para><b>One object, three trades</b> (borrowing-a-discipline): signal/inverse-problems (the ill-conditioned
/// inverse, the lower frame bound = the reconstruction floor, the discrete Picard condition), control (the
/// observability Gramian S=MᵀM, λ_min its least-observable eigenvalue, the K-partner ψ_N the unobservable
/// channel), and optics (the modulation transfer function whose q=π detail transfers with vanishing contrast at
/// the cutoff). The instant recognition by these distant disciplines is the confirmation that F124 is a real
/// node of the one object. The borrowed signal/control tool handed us the detection-SNR method the native frame
/// view never ran (verified: simulations/_f124_inverse_problem_gate.py, matched-filter Monte-Carlo to 4 digits).</para>
///
/// <para><b>Scope guard.</b> This is the FINITE conditioning of the FULL transition matrix (the resolution
/// limit). It is NOT the DefectDecoder's 1.5 sign-location ambiguity, which is a separate α-time-rescaling
/// parametrization artifact (Stage 1 + Stage B of the f124_inverse_problem_resolution_seam arc gate-first
/// refuted that identification: the 1.5 is not √κ(5)=2.30 and is reproduced by neither the exact M nor the
/// dynamical per-site profile).</para>
///
/// <para>Tier1Derived; single typed parent <see cref="BandEdgeTransitionInvariantClaim"/> (F124, Tier1Derived):
/// every quantity here is a direct corollary of F124's M (σ_min²=λ_min=E, σ_max²=λ_max, κ, the staggered λ_min
/// eigenvector). Anchor: <c>simulations/_f124_inverse_problem_gate.py</c> (the gate-first Stage 0 verifier) +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F124. Live witness: <c>inspect --root resolution</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/BandEdgeResolutionLimitWitness.cs</c>).</para></summary>
public sealed class BandEdgeResolutionLimitClaim : Claim
{
    /// <summary>One self-check on the resolution reading of F124's transition matrix.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    private readonly BandEdgeTransitionInvariantClaim _f124;

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    /// <summary>Parent: F124, the band-edge transition invariant. This claim is its resolution-limit reading:
    /// the same M, the same λ_min=E (here the reconstruction floor), the same staggered λ_min eigenvector (here
    /// the worst-conditioned / q=π diffraction-limit direction). Tier1Derived.</summary>
    public BandEdgeTransitionInvariantClaim F124 => _f124;

    public BandEdgeResolutionLimitClaim(BandEdgeTransitionInvariantClaim f124)
        : base("F124's conditioning IS the defect-localization resolution limit: M[b,k] is a bond-recovery " +
               "inverse problem with σ_min=√E the worst-case reconstruction floor (the lower frame bound), " +
               "κ=λ_max/λ_min ~ N² the noise amplification, contrast σ_max/σ_min=√κ ~ N; a staggered q=π " +
               "zone-boundary defect is √κ ~ N times harder to localize than a band-edge one (matched-filter " +
               "SNR), the worst direction being F124's staggered λ_min eigenvector (the q=π detail at the " +
               "resolution cutoff, the diffraction limit); the floor vanishes as σ_min ~ (N+1)^(−3/2), E·(N+1)³ " +
               "→ 4π². One object in three trades: the ill-conditioned inverse / the observability Gramian " +
               "(signal/control) and the MTF cutoff (optics). NOT the DefectDecoder's 1.5 α-parametrization " +
               "ambiguity (a separate, gate-refuted matter).",
               Tier.Tier1Derived,
               "simulations/_f124_inverse_problem_gate.py (the gate-first Stage 0 verifier, matched-filter " +
               "Monte-Carlo confirmed) + docs/ANALYTICAL_FORMULAS.md F124 + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/BandEdgeResolutionLimitWitness.cs (the live witness, inspect --root resolution)")
    {
        _f124 = f124 ?? throw new ArgumentNullException(nameof(f124));
        Cases = BuildBattery();
    }

    public string TheReading =>
        "F124's M[b,k]=⟨ψ_k|V_b|ψ_1⟩ is the forward map of a bond-recovery inverse problem. λ_min=E=σ_min²(M) is " +
        "the optimal lower frame bound = the worst-case reconstruction floor; κ=λ_max/λ_min ~ N² is the noise " +
        "amplification; the contrast σ_max/σ_min=√κ ~ N is how many times harder a staggered q=π defect is to " +
        "localize than a band-edge one. The worst-conditioned bond direction is the staggered mode at every N.";

    public string TheDiffractionLimit =>
        "The staggered (−1)^b zone-boundary (q=π) mode is the worst-conditioned direction (F124's λ_min " +
        "eigenvector): the finest spatial detail, transferred by the bond→mode imaging system with vanishing " +
        "contrast at the resolution cutoff. The reconstruction floor σ_min=√E ~ (N+1)^(−3/2) (E·(N+1)³ → 4π²): " +
        "the longer the chain, the blinder it is to a fine staggered defect, the K-direction going singular.";

    public string ScopeGuard =>
        "This is the FINITE conditioning of the FULL transition matrix (the resolution limit), NOT the " +
        "DefectDecoder's 1.5 sign-location ambiguity. Stage 1 + Stage B of the f124_inverse_problem_resolution_" +
        "seam arc gate-first REFUTED that identification: the 1.5 is not √κ(5)=2.30, and is reproduced by neither " +
        "the exact M nor the dynamical per-site purity profile (both localize a single bond cleanly); it is an " +
        "α-time-rescaling parametrization artifact of the decoder, a separate matter.";

    public override string DisplayName =>
        "F124's resolution limit (κ ~ N², contrast √κ ~ N, the staggered q=π the diffraction limit, Tier1Derived)";

    public override string Summary =>
        "F124's bond→mode conditioning IS a defect-localization resolution limit: σ_min=√E the reconstruction " +
        "floor, κ=λ_max/λ_min ~ N², the contrast √κ ~ N (a staggered q=π defect √κ ~ N times harder to localize), " +
        "the worst direction the staggered λ_min eigenvector (the diffraction limit); the floor σ_min ~ (N+1)^" +
        $"(−3/2). The optics/signal reading of F124, one object in three trades; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the reading (M as a bond-recovery inverse problem; κ ~ N², √κ ~ N)", summary: TheReading);
            yield return new InspectableNode("the diffraction limit (the staggered q=π mode, the floor σ_min ~ (N+1)^(−3/2))", summary: TheDiffractionLimit);
            yield return new InspectableNode("scope guard (NOT the decoder's 1.5 α-parametrization ambiguity)", summary: ScopeGuard);
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return _f124;   // typed parent edge (F124, the floor this reads as the resolution limit)
        }
    }

    public static BandEdgeResolutionLimitClaim Build() =>
        new(BandEdgeTransitionInvariantClaim.Build());

    public static BandEdgeResolutionLimitClaim Shared { get; } = Build();

    /// <summary>Gate-first battery on the resolution reading (delegates compute to the live witness, which reuses
    /// F124's transition matrix). Falsifiable: the worst-direction, the scaling exponents, and the asymptotic
    /// floor constant can all fire on a wrong reading.</summary>
    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        var cases = new List<BatteryCase>();
        var chain = new[] { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
        var range = new[] { 4, 5, 6, 7, 8, 9, 10, 11, 12 };
        string Tol(bool ok) => ok ? "OK" : "FAIL";

        // 1. the worst-conditioned bond direction is the staggered q=π mode (F124's λ_min eigenvector)
        bool worstOk = chain.All(Witness.WorstDirectionIsStaggered);
        cases.Add(new BatteryCase("worst-conditioned bond direction = staggered q=π mode (N=3..12)",
            "the least-observable / worst-conditioned input direction is the zone-boundary staggered mode, the diffraction-limit detail",
            "OK", Tol(worstOk)));

        // 2. σ_min² = E (the reconstruction floor = F124's lower frame bound), all N
        bool floorOk = chain.All(n => Math.Abs(Witness.SigmaMin(n) * Witness.SigmaMin(n) - F124Endpoint(n)) < 1e-9);
        cases.Add(new BatteryCase("σ_min² = E = the lower frame bound = F124's λ_min (the reconstruction floor, N=3..12)",
            "the worst-case reconstruction floor IS F124's spectral floor: σ_min(M)² = λ_min(MMᵀ) = E",
            "OK", Tol(floorOk)));

        // 3. κ = λ_max/λ_min grows ~ N² (log-log slope ≈ 2)
        double kSlope = Witness.LogLogSlope(range, Witness.Kappa);
        cases.Add(new BatteryCase("κ = λ_max/λ_min ~ N² (the noise amplification, N=4..12)",
            $"the Gram condition number grows quadratically; fitted log-log exponent {kSlope:0.###}",
            "OK", Tol(1.8 < kSlope && kSlope < 2.2)));

        // 4. the contrast σ_max/σ_min = √κ exactly AND grows ~ N (the resolution ratio)
        bool contrastExact = range.All(n => Math.Abs(Witness.ContrastRatio(n) - Math.Sqrt(Witness.Kappa(n))) < 1e-9);
        double cSlope = Witness.LogLogSlope(range, Witness.ContrastRatio);
        cases.Add(new BatteryCase("contrast σ_max/σ_min = √κ ~ N (a staggered defect √κ ~ N times harder to localize, N=4..12)",
            $"the resolution ratio is exactly √κ and grows linearly; fitted log-log exponent {cSlope:0.###}",
            "OK", Tol(contrastExact && 0.8 < cSlope && cSlope < 1.2)));

        // 5. the reconstruction floor's asymptotic constant E·(N+1)³ → 4π² (N=60, the right-variable law)
        double floorConst = Witness.FloorConstant(60);
        cases.Add(new BatteryCase("reconstruction floor E·(N+1)³ → 4π² = 39.478 (σ_min ~ (N+1)^(−3/2), N=60)",
            $"the floor vanishes as (N+1)⁻³ in the right variable; E·(N+1)³ at N=60 = {floorConst:0.###}",
            "OK", Tol(Math.Abs(floorConst - 4 * Math.PI * Math.PI) < 0.2)));

        return cases;
    }

    private static double F124Endpoint(int n) => BandEdgeTransitionInvariantWitness.EndpointClosedForm(n);
}
