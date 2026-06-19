using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Ptf;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>THE STONE (felt_time_dimensions arc, step B), typed (Tier1Candidate): the TRAJECTORY-level
/// dual of the eigenvalue-level value/vector split. The live witness is
/// <see cref="StoneSurvivorClosureWitness"/> (<c>inspect --root stone</c>); the gate-first Python verifier
/// is <c>simulations/_stone_survivor_alpha_closure.py</c> (two-lens reviewed 2026-06-19).
///
/// <para>For the near-stationary MODE-ISOLATING probe rho_0 = I/d + eps*Herm(mode), the PTF painter
/// closure Sum_i ln(alpha_i) -- computed through the CANONICAL Symphony FitAlpha -- reads the chosen
/// mode's first-order RATE shift under a delta-J bond defect: OUT of the +-0.05 window AND sign-coherent
/// for the soft survivor interior (2,2) (Re lambda moves; K_decay defect-sensitive, soft darkness), and
/// IN the window for the rigid (0,1) band edge (Re = -2gamma frozen; K_decay defect-invariant). It is a
/// CONSTRUCTIVE confirmation of (A) for this probe, NOT a universal trajectory law (review-pinned):
/// probe-state-specific (a polarized survivor-dominated state holds), the rate shift certified by
/// sign-coherence, the magnitude scaling+sign.</para>
///
/// <para>Typed parents: <see cref="AbsorptionTheoremClaim"/> (the -2gamma rate, a_0) +
/// <see cref="SurvivalIncompletenessMirrorClaim"/> (the (A) value/vector survivor this confirms; caps the
/// tier at Tier1Candidate, honestly) + <see cref="ChiralMirrorTrajectoryClaim"/> (the PTF painter closure
/// law, EQ-014).</para></summary>
public sealed class StoneSurvivorClosureClaim : Claim
{
    /// <summary>Typed parent: the -2gamma rate law (the closure reads a shift in this rate).</summary>
    public AbsorptionTheoremClaim RateLaw { get; }

    /// <summary>Typed parent: the eigenvalue-level value/vector survivor (A) that this confirms at the
    /// trajectory level. Tier1Candidate, so it caps this claim's tier.</summary>
    public SurvivalIncompletenessMirrorClaim ValueVector { get; }

    /// <summary>Typed parent: the PTF painter closure law Sum_i ln(alpha_i) (EQ-014).</summary>
    public ChiralMirrorTrajectoryClaim PainterClosure { get; }

    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    // Lazy: the battery runs the full witness (a few small Liouvillian eigendecompositions); compute it
    // only when the claim is drilled into, so building the whole registry stays cheap.
    private IReadOnlyList<BatteryCase>? _cases;
    public IReadOnlyList<BatteryCase> Cases => _cases ??= BuildBattery();
    public int PassCount => Cases.Count(c => c.Passes);

    public StoneSurvivorClosureClaim(AbsorptionTheoremClaim rateLaw,
        SurvivalIncompletenessMirrorClaim valueVector, ChiralMirrorTrajectoryClaim painterClosure)
        : base("THE STONE (felt_time arc B): the TRAJECTORY-level dual of the eigenvalue value/vector split. For the " +
               "near-stationary mode-isolating probe rho_0 = I/d + eps*Herm(mode), the PTF painter closure Sum_i ln(alpha_i) " +
               "(through the CANONICAL Symphony FitAlpha) reads the chosen mode's first-order RATE shift under a delta-J bond " +
               "defect: OUT of the +-0.05 window AND sign-coherent for the soft survivor interior (2,2) (Re moves; K_decay " +
               "defect-sensitive), IN the window for the rigid (0,1) band edge (Re=-2gamma frozen; K_decay defect-invariant). " +
               "A CONSTRUCTIVE confirmation of (A) for this probe (two-lens reviewed 2026-06-19), NOT a universal trajectory " +
               "law: probe-state-specific (a polarized survivor-dominated state holds), the rate shift certified by " +
               "sign-coherence (coh>0.8), the magnitude scaling+sign. Live: inspect --root stone.",
               Tier.Tier1Candidate,
               "simulations/_stone_survivor_alpha_closure.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/StoneSurvivorClosureWitness.cs + " +
               "simulations/value_vector_felt_time.py")
    {
        RateLaw = rateLaw ?? throw new ArgumentNullException(nameof(rateLaw));
        ValueVector = valueVector ?? throw new ArgumentNullException(nameof(valueVector));
        PainterClosure = painterClosure ?? throw new ArgumentNullException(nameof(painterClosure));
    }

    public override string DisplayName =>
        "The stone: the PTF closure reads the mode's rate shift (survivor OUT+coh, band-edge IN) for the I/d+eps probe (Tier1Candidate)";

    public override string Summary =>
        "the PTF painter closure Sum_i ln(alpha_i), via the CANONICAL Symphony FitAlpha on the mode-isolating probe " +
        "I/d+eps*Herm(mode), reads the mode's first-order RATE shift: OUT + sign-coherent (rate-shift) for the soft " +
        "survivor interior (2,2), IN (frozen) for the rigid band edge (0,1) - the trajectory-level dual of (A). " +
        "Probe-state-specific (review-pinned, not a universal law), the rate shift certified by sign-coherence. " +
        $"Tier1Candidate. Live: inspect --root stone.";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the trajectory-level value/vector readout",
                summary: "Sum ln alpha = log of the net time-dilation of the painted trajectory. For the mode-isolating " +
                         "probe (I/d stationary => only the chosen mode drives the single-site purity), it equals the mode's " +
                         "first-order rate shift: nonzero (OUT) when K_decay moves (the soft survivor), ~0 (IN) when K_decay " +
                         "is frozen (the rigid band edge). Certified a real rate shift only via the SIGN-COHERENCE of the " +
                         "reliable per-site f (an asymmetric redistribution also breaks the closure). The eigenvalue dual is " +
                         "inspect --root survivor / horizon.");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return RateLaw;         // typed parent: the -2gamma rate
            yield return ValueVector;     // typed parent: the (A) value/vector survivor this confirms
            yield return PainterClosure;  // typed parent: the PTF painter closure law
        }
    }

    private static string Fmt(double v) => v.ToString("+0.000;-0.000", CultureInfo.InvariantCulture);

    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        var w = new StoneSurvivorClosureWitness(4);   // N=4: interior survivor (2,2) + band edge (0,1)
        var s = w.Survivor;
        var be = w.BandEdge;
        return new List<BatteryCase>
        {
            new("survivor closure OUT of the +-0.05 window (K_decay moves)",
                $"N=4 (2,2): Sum ln alpha = {Fmt(s.Closure)} ({s.Reliable}/{s.Sites} reliable)",
                "OUT", s.InWindow ? "IN" : "OUT"),
            new("survivor break is a SIGN-COHERENT rate shift (not a redistribution)",
                $"coh = {s.SignCoherence.ToString("0.00", CultureInfo.InvariantCulture)} (>0.8 => one-signed reliable f)",
                "rate-shift", s.IsRateShift ? "rate-shift" : "not-rate-shift"),
            new("band edge closure IN window (K_decay defect-invariant, rigid darkness)",
                $"N=4 (0,1): Sum ln alpha = {Fmt(be.Closure)}", "IN", be.InWindow ? "IN" : "OUT"),
            new("band edge is FROZEN, not a rate shift",
                $"coh = {be.SignCoherence.ToString("0.00", CultureInfo.InvariantCulture)}, in-window={be.InWindow}",
                "frozen", be.IsRateShift ? "rate-shift" : "frozen"),
        };
    }
}
