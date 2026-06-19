using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>(D) THE CLOSURE FUNCTIONAL (felt_time_dimensions arc, D follow-up), typed (Tier1Candidate):
/// the EXACT bond functional of the survivor's first-order rate shift -- amplitude-squared in the density
/// gradient. The live witness is <see cref="SurvivorDiffusionGradientWitness"/> (<c>inspect --root gradient</c>);
/// the gate-first Python verifiers are <c>simulations/_felt_time_amplitude_law.py</c> (the block-level law,
/// N=4..7) and <c>simulations/_felt_time_closure_functional.py</c> (the trajectory ground truth).
///
/// <para>The PTF painter closure Sum f(b) (step B, <see cref="StoneSurvivorClosureClaim"/>) reads the
/// survivor's first-order RATE shift dRe(b). This is its exact bond functional: the slow survivor is a
/// DENSITY / diffusion mode n(j) (NOT a current -- its hopping content Tr(M^dag H_b) is identically zero),
/// so a delta-J defect on bond b=(j,j+1) perturbs the LOCAL diffusion coefficient D_b, and the first-order
/// rate shift is the diffusion Rayleigh-quotient derivative dRe(b) = dlambda/dD_b ~ (n(j)-n(j+1))^2 -- the
/// SQUARED density GRADIENT ("amplitude^2"). It vanishes at the no-flux (reflecting) chain ends (gradient
/// -> 0) and peaks in the interior, mirror-symmetric, Q-invariant (the lowest diffusion harmonic k_min is
/// Q-fixed). Gate-first verified N=4..7: dRe/grad^2 bond-independent (CV 0.001..0.07), log-log slope dRe vs
/// |grad| = 2.0..2.2. The earlier single-particle phi*phi candidate used the WRONG standing wave
/// (single-particle, not the multi-magnon density mode): right power, wrong wave.</para>
///
/// <para>Typed parents: <see cref="AbsorptionTheoremClaim"/> (the -2gamma rate the diffusion mode decays
/// at) + <see cref="SurvivalIncompletenessMirrorClaim"/> (the (A) survivor whose density gradient this is;
/// caps the tier at Tier1Candidate, honestly). Sibling (see-cref, not a typed parent):
/// <see cref="StoneSurvivorClosureClaim"/> -- the same rate shift read at the TRAJECTORY level (B), this at
/// the EIGENVALUE level.</para></summary>
public sealed class SurvivorDiffusionGradientClaim : Claim
{
    /// <summary>Typed parent: the -2gamma rate law (the diffusion mode decays at this dephasing rate).</summary>
    public AbsorptionTheoremClaim RateLaw { get; }

    /// <summary>Typed parent: the (A) survivor (the slow density mode whose gradient this is).
    /// Tier1Candidate, so it caps this claim's tier.</summary>
    public SurvivalIncompletenessMirrorClaim Survivor { get; }

    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    // Lazy: the battery runs the witness (small block eigendecompositions); compute it only when the claim
    // is drilled into, so building the whole registry stays cheap.
    private IReadOnlyList<BatteryCase>? _cases;
    public IReadOnlyList<BatteryCase> Cases => _cases ??= BuildBattery();
    public int PassCount => Cases.Count(c => c.Passes);

    public SurvivorDiffusionGradientClaim(AbsorptionTheoremClaim rateLaw, SurvivalIncompletenessMirrorClaim survivor)
        : base("(D) THE CLOSURE FUNCTIONAL (felt_time arc D follow-up): the survivor's first-order bond rate shift " +
               "dRe(b) ~ (density-mode gradient at bond b)^2 -- the diffusion Rayleigh quotient (amplitude^2). The slow " +
               "survivor is a DENSITY/diffusion mode n(j) (its hopping content Tr(M^dag H_b)=0 identically), so a delta-J " +
               "defect perturbs the local diffusion coefficient D_b and dRe(b)=dlambda/dD_b ~ (n(j)-n(j+1))^2: ~0 at the " +
               "no-flux chain ends, peaked interior, mirror-symmetric, Q-invariant (k_min harmonic Q-fixed). Gate-first " +
               "N=4..7: dRe/grad^2 bond-independent (CV 0.001..0.07), log-log slope dRe vs |grad| = 2.0..2.2. The exact " +
               "bond functional UNDERLYING the PTF closure (B, the trajectory-level dual). The earlier single-particle " +
               "phi*phi candidate used the wrong standing wave (multi-magnon density mode, not single-particle): right " +
               "power, wrong wave. Live: inspect --root gradient.",
               Tier.Tier1Candidate,
               "simulations/_felt_time_amplitude_law.py + simulations/_felt_time_closure_functional.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/SurvivorDiffusionGradientWitness.cs")
    {
        RateLaw = rateLaw ?? throw new ArgumentNullException(nameof(rateLaw));
        Survivor = survivor ?? throw new ArgumentNullException(nameof(survivor));
    }

    public override string DisplayName =>
        "The closure functional (D): the survivor's bond rate shift is amplitude^2 in the density gradient (diffusion Rayleigh; Tier1Candidate)";

    public override string Summary =>
        "the survivor's first-order bond rate shift dRe(b) ~ (density-mode gradient)^2 -- the diffusion Rayleigh quotient. " +
        "The slow survivor is a density/diffusion mode; a delta-J bond defect perturbs the local diffusion coefficient, so " +
        "dRe ~ (n(j)-n(j+1))^2, ~0 at the no-flux ends. The exact bond functional underlying the PTF closure (B). " +
        "Tier1Candidate. Live: inspect --root gradient.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the amplitude^2 reading",
                summary: "dRe(b) is the survivor mode's first-order rate shift under a delta-J defect on bond b. Because the " +
                         "mode is a DENSITY/diffusion standing wave, dRe(b) is the diffusion Rayleigh-quotient derivative " +
                         "~ (n(j)-n(j+1))^2 (the squared gradient): LINEAR in grad^2 (bond-independent ratio), quadratic in " +
                         "the per-site amplitude, ~0 at the reflecting chain ends. The trajectory-level dual is inspect --root stone.");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return RateLaw;    // typed parent: the -2gamma rate
            yield return Survivor;   // typed parent: the (A) survivor whose density gradient this is
        }
    }

    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        var w = new SurvivorDiffusionGradientWitness(4);   // N=4: interior survivor (2,2), 3 bonds
        var b = w.Bonds;
        double endMax = Math.Max(b[0].GradSq, b[^1].GradSq);
        double interiorMax = b.Skip(1).Take(b.Count - 2).DefaultIfEmpty(b[0]).Max(x => x.GradSq);
        string slopeStr = w.PowerSlope.ToString("0.00", CultureInfo.InvariantCulture);
        string cvStr = w.RatioCv.ToString("0.000", CultureInfo.InvariantCulture);
        return new List<BatteryCase>
        {
            new("rate shift is amplitude^2 in the density gradient (log-log slope ~2)",
                $"N=4 (2,2): slope dRe vs |grad| = {slopeStr}", "p~2",
                w.PowerSlope is > 1.6 and < 2.4 ? "p~2" : "off"),
            new("dRe/grad^2 bond-independent (the diffusion Rayleigh law)",
                $"N=4 (2,2): CV(dRe/grad^2) = {cvStr}", "bond-independent",
                w.RatioCv < 0.15 ? "bond-independent" : "varies"),
            new("the density gradient is quiet at the no-flux chain ends",
                $"N=4 (2,2): end grad^2 max = {endMax.ToString("0.0000", CultureInfo.InvariantCulture)}, " +
                $"interior = {interiorMax.ToString("0.0000", CultureInfo.InvariantCulture)}", "quiet-ends",
                endMax < 0.5 * interiorMax ? "quiet-ends" : "loud-ends"),
        };
    }
}
