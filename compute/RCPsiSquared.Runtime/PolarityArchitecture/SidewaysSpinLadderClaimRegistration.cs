using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="SidewaysSpinLadderClaim"/> (the <c>sideways_spin_ladder</c>
/// open arc, typed 2026-08-09): the staggered second ladder S⁺(ρ) = Σ_l (−1)^l c_l†ρc_l† intertwines the
/// Liouvillian for Σ-odd real-symmetric pure hopping and any rate profile (derived, F142 Lemmas 2.1 + 2.3),
/// and the F125 fold family is the interiors of its two chains at p+q̃ = N∓1, one spin-(N−3)/2 multiplet per
/// chain with Clebsch-Gordan S⁺ transport norms (measured N=5, N=7 and N=9, gated in C# and Python, the
/// N=9 walk via LU shift-invert; the two η chains inferred from the 4N−8 accounting, untested). One typed
/// parent edge:
///
/// <list type="bullet">
///   <item><see cref="SpectatorIntertwinerClaim"/>: F125, whose W is the η sibling Φ (the same operator
///         under two names, the arc's founding identification) and whose fold family the chains carry.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 candidate. The intertwining half is derived (F142), but the multiplet
/// half (one irreducible chain per anti-diagonal, CG norms at every rung) is measured at three odd N
/// (5, 7, 9) without an all-N proof, so the claim as a whole stays a candidate. Parent Tier 1 derived.</para>
///
/// <para>Requires <see cref="SpectatorIntertwinerClaim"/> to be registered (the builder topo-resolves,
/// so the order of registration is free).</para></summary>
public static class SidewaysSpinLadderClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSidewaysSpinLadderClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<SidewaysSpinLadderClaim>(b =>
            new SidewaysSpinLadderClaim(b.Get<SpectatorIntertwinerClaim>()));
}
