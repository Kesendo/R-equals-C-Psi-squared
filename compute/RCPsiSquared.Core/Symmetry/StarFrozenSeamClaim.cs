using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The star's frozen seam (Tier 1 candidate): on the star the longest-lived coherence NEVER
/// un-freezes — its survivor is frozen (|Im| = 0, a pure decay) at every Q, for N ≥ 5 — and that threshold
/// IS the structural ceiling, read dynamically. This is the third member of a trichotomy the typed layer had
/// named only for the chain and the ring, and it ties the frozen-seam reading to an already-typed quantity
/// (the ceiling g2 = 4/(N−1)).
///
/// <para><b>The finding.</b> Watch the slowest non-kernel mode as the dephasing weakens (Q = J/γ grows): does
/// it stay frozen or acquire a frequency (un-freeze)? A single Q does not separate the topologies (below its
/// horizon every survivor is overdamped). The |Im|(Q) curve does: the CHAIN un-freezes at the SE-EP coherence
/// horizon Q*(N) (<see cref="CoherenceHorizonClaim"/>); the RING un-freezes at its handover (the oscillating
/// band edge overtakes the frozen (2,2) level-crossing seam); the STAR never un-freezes (N ≥ 5).</para>
///
/// <para><b>The mechanism and the threshold.</b> The star's survivor is the darkest [H,A]=0 commutant (1,1)
/// coherence — it lies in the ad_H kernel, so it commutes with H and cannot oscillate: frozen by construction.
/// It is the survivor exactly when it undercuts the −2γ Absorption floor, i.e. when the structural ceiling
/// g2 = 4/(N−1) ≤ 1 (N ≥ 5; the leaf-manifold ceiling of <see cref="StructuralCeilingClaim"/>, equal to the
/// SecondClockRegime bridge 4/(m+1) at m = N−2). At N=4 (g2 = 4/3 > 1) the commutant mode is brighter than the
/// floor, an oscillating band-edge mode is the slowest, and the star UN-freezes — the known (2,2)/K₄ outlier.
/// So the high-Q ceiling g2 ≤ 1 and the all-Q frozenness of the survivor are the same fact: the star's frozen
/// seam IS the structural ceiling, read dynamically (the survivor-level form of the "star has no coherence
/// horizon" corollary, PROOF_STRUCTURAL_CEILING §7).</para>
///
/// <para>Tier1Candidate: gate-verified across N=4,5,6,7,8 (simulations/star_frozen_seam.py — star frozen at
/// all Q for N ≥ 5, max |Im| &lt; 1e-15; N=4 un-freezes at 1.73) and live (the commutant ceiling is exact and
/// cheap; the |Im|(Q) sweep is the full 4^N at N≤5). Two typed parents: <see cref="StructuralCeilingClaim"/>
/// (Tier1Derived, the g2 = 4/(N−1) threshold + the commutant mechanism it reuses) and
/// <see cref="SecondClockRegimeClaim"/> (Tier1Candidate, the {0,2}/second-clock regime map whose star/GRADUAL
/// case this sharpens to a survivor-level statement). Capped at Tier1Candidate by the weaker parent
/// SecondClockRegimeClaim.</para>
///
/// <para><b>Model scope (XY ceiling vs Heisenberg survivor).</b> g2 = 4/(N−1) here is the XY-network
/// (hopping-only) value, by design: this claim's witness and <see cref="StructuralCeilingClaim"/> build
/// the star's XY hopping H with no ZZ diagonal, so for the XY star it is exactly the survivor darkness.
/// The project's canonical Heisenberg star (XX+YY+ZZ) survivor darkens further to ⟨n_XY⟩(Q→∞) = 4/N: the
/// single-excitation ZZ potential (hub −(N−1), leaves N−3) shifts which ad_H-kernel commutant is darkest,
/// from 4/(N−1) (XY) to 4/N (Heisenberg), verified N=5..8 + full-4^N at N=6 (simulations/star_survivor_heisenberg.py).
/// Model-robust (both XY and Heisenberg): frozen at every Q (N≥5), commutant only in the high-Q LIMIT
/// (‖[H,ρ]‖ ∝ 1/Q, not zero at finite Q), the (1,1) boundary sector, no horizon. Only the darkness VALUE
/// (and the XY g2≤1 threshold / N=4 outlier logic) is model-specific.</para>
///
/// <para>Live witness: <c>inspect --root starseam</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/StarFrozenSeamWitness.cs</c>).</para></summary>
public sealed class StarFrozenSeamClaim : Claim
{
    /// <summary>Parent: the structural ceiling. The star's leaf-manifold ceiling g2 = 4/(N−1) is BOTH the
    /// threshold (frozen seam iff ≤ 1) AND the commutant mechanism (the darkest [H,A]=0 coherence, frozen by
    /// living in the ad_H kernel). Tier1Derived.</summary>
    public StructuralCeilingClaim Ceiling { get; }

    /// <summary>Parent: the second clock's regime map. This claim sharpens its star/GRADUAL case ("only
    /// asymptotic protection, no sharp horizon") into a survivor-level statement: the slowest mode never
    /// acquires a frequency at all. Tier1Candidate; caps this child.</summary>
    public SecondClockRegimeClaim Regime { get; }

    public StarFrozenSeamClaim(StructuralCeilingClaim ceiling, SecondClockRegimeClaim regime)
        : base("On the star the longest-lived coherence never un-freezes: its survivor is frozen (|Im|=0) at " +
               "every Q for N≥5, and that threshold IS the structural ceiling read dynamically. As Q grows the " +
               "chain un-freezes at the SE-EP coherence horizon Q*(N) and the ring at its handover (the band edge " +
               "overtakes the frozen (2,2) seam), but the star's survivor is the darkest [H,A]=0 commutant (1,1) " +
               "coherence — it commutes with H, so it cannot oscillate (frozen by construction) — and it is the " +
               "survivor exactly when it undercuts the −2γ floor, i.e. g2=4/(N−1)≤1 (N≥5). N=4 (g2=4/3>1) un-freezes " +
               "(the band edge wins, the (2,2)/K₄ outlier). The third member of the trichotomy chain(SE-EP)/ring" +
               "(level crossing)/star(commutant); the high-Q ceiling and the all-Q frozenness are one fact. " +
               "g2 is not only the threshold but the value the survivor darkness saturates at: " +
               "⟨n_XY⟩(Q→∞) = min(g2, 1) — on g2 below the floor (0.800 at N=6, 0.667 at N=7) or the −2γ floor 1 " +
               "above it (N≤5, where it un-freezes). Gate-verified N=4..8 (simulations/star_frozen_seam.py).",
               Tier.Tier1Candidate,
               "simulations/star_frozen_seam.py + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/StarFrozenSeamWitness.cs (StarFrozenSeamWitness, inspect --root starseam)")
    {
        Ceiling = ceiling ?? throw new ArgumentNullException(nameof(ceiling));
        Regime = regime ?? throw new ArgumentNullException(nameof(regime));
    }

    public override string DisplayName =>
        "The star's frozen seam: the survivor never un-freezes (N≥5), the structural ceiling read dynamically";

    public override string Summary =>
        $"the star's slowest coherence is the [H,A]=0 commutant (1,1) mode (frozen by construction); it is the " +
        $"survivor — frozen at all Q — iff g2=4/(N−1)≤1 (N≥5), un-freezing only at N=4 (g2=4/3>1, the outlier). " +
        $"Third member of chain(SE-EP un-freeze)/ring(frozen level crossing)/star(frozen commutant); the ceiling " +
        $"read dynamically ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the |Im|(Q) signature (frozen vs un-freeze)",
                summary: "a single Q does not separate the topologies (below its horizon every survivor is " +
                         "overdamped, |Im|=0). The signature is the curve: the chain un-freezes at Q*(N), the ring " +
                         "at its handover, the star never (N≥5). The star's survivor stays frozen at every Q.");
            yield return new InspectableNode("the mechanism: the commutant commutes with H",
                summary: "the survivor is the darkest [H,A]=0 commutant (1,1) coherence — it lies in the ad_H " +
                         "kernel, so −i[H,ρ]=0 and it cannot oscillate. Frozen by construction. The ring's (2,2) " +
                         "seam is frozen by a DIFFERENT route (a level crossing: two real eigenvalues coincide); " +
                         "the chain is neither (it oscillates above Q*).");
            yield return new InspectableNode("the threshold = the structural ceiling g2 = 4/(N−1) ≤ 1",
                summary: "the commutant coherence is the survivor only when it undercuts the −2γ floor, i.e. " +
                         "g2=4/(N−1)≤1 (N≥5). At N=4 (4/3>1) it is brighter than the floor → an oscillating band-edge " +
                         "mode wins → the star un-freezes. The high-Q ceiling and the all-Q frozenness are the same " +
                         "fact; the frozen seam is the ceiling read dynamically. Sharper: g2 is the survivor darkness " +
                         "itself in the high-Q limit, ⟨n_XY⟩(Q→∞)=min(g2,1) — measured 0.800 at N=6, 0.667 at N=7 " +
                         "(= g2), and pinned to the floor 1 for N≤5. The single formula min(g2,1) is exactly the " +
                         "ladder's 'un-freezes iff g2>1'.");
            yield return Ceiling;   // typed parent edge (the ceiling threshold + commutant mechanism, Tier1Derived)
            yield return Regime;    // typed parent edge (the second-clock regime map it sharpens, Tier1Candidate, caps this child)
        }
    }

    public static StarFrozenSeamClaim Build() =>
        new StarFrozenSeamClaim(StructuralCeilingClaim.Build(), SecondClockRegimeClaim.Build());

    public static StarFrozenSeamClaim Shared { get; } = Build();
}
