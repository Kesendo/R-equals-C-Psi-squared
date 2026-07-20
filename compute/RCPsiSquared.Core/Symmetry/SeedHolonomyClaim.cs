using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The eigenvector-frame holonomy around the (1,2)-block defective seed is the mod-4 memory loop
/// i⁴ = 1: encircling the exceptional point once, the two coalescing right eigenvectors' frame (transported
/// in the biorthogonal vᵀv gauge) gives a generator with eigenvalues ±i, so M₂ = −I and M₄ = +I — the frame
/// is single-valued only after four loops. Companion to the eigenVALUE swap that <c>Numerics/Monodromy.cs</c>
/// already computes: this is the eigenVECTOR-phase half nobody had tracked.
///
/// <para>It SHARES the order-4 (Z₄) structure of the algebraic i⁴ = 1 memory loop (the Pi2 Z₄ /
/// NinetyDegreeMirrorMemory) — a correspondence DECIDED 2026-07-16 (the third-clock control,
/// experiments/SEED_HOLONOMY_THIRD_CLOCK.md): the holonomy Z₄ is NOT Π's Z₄ and not the Clifford-degree
/// mod-4 — it is generic EP2 geometry. A random complex-symmetric pencil with no palindrome, chain, or
/// Clifford structure carries the identical ±i / M₂=−I / M₄=+I signature (simulations/seed_holonomy_generic.py:
/// 5/5 random pencils plus a disguised deterministic EP2, machine-zero residuals, √ε gap slope 0.500; C# pin
/// <c>EigenvectorHolonomyTests.DisguisedPiFreeEP2_SameMod4Holonomy_TheThirdClockIsGeneric</c>). What the clocks
/// share is the scalar i alone (see experiments/ONE_FOUR_THESIS.md). Gauge note: ±i / mod-4 is the
/// biorthogonal-gauge reading; the gauge-invariant content is the swap (a Hermitian gauge gives only the
/// mod-2 ±1 swap; gated as the scout's gauge-contingency control).</para>
///
/// <para>Tier1Candidate: strong, ε-stable, independently re-verified numerics (a fresh cut confirmed
/// M₁~90°-rotation, M₂=−I, M₄=I at the N=9 seed in the R=+1 sector; the live witness reproduces it, using
/// full-block isolation at N=5 and the R=+1 sector at N=9), plus a KNOWN analytic reason (a defective EP2
/// has the ±i eigenvector holonomy — Heiss/Dembowski; the seed's defectiveness is F89's Kato simple-zero,
/// Tier1Derived). Not Tier1Derived because the from-below computation here is numerical, the analytic
/// derivation is not yet written as a repo proof, and EP2-ness ∀ odd N is itself the open β-exotic item.</para></summary>
public sealed class SeedHolonomyClaim : Claim
{
    public SeedHolonomyClaim() : base(
        "The eigenvector holonomy around the (1,2)-block defective seed is the mod-4 memory loop i⁴=1 " +
        "(in the vᵀv gauge: M₁ eigenvalues ±i, M₂ = −I, M₄ = +I; confirmed at N=5 full-block and N=9 sector).",
        Tier.Tier1Candidate,
        "compute/RCPsiSquared.Core/Numerics/EigenvectorHolonomy.cs (biorthogonal vᵀv frame transport) + " +
        "compute/RCPsiSquared.Diagnostics/Foundation/SeedHolonomyWitness.cs (inspect --root holonomy) + " +
        "experiments/F89_SEED_EXISTENCE_REDUCTION.md establishes this real-axis defective seed (its census to " +
        "N=11, backed by the exact nullity identity, plus the N=5 and N=9 numerics here) + " +
        "experiments/F86_EP_THROUGH_THE_CLOCK.md is cited for the borrowed forgetting-to-remembering reading; " +
        "note its F86a-retraction denied a real-axis defective EP on this block; the F89 census locates the " +
        "isolated seed that earlier real-axis scan reported absent, and F86 is now corrected accordingly " +
        "(PROOF_F86A_EP_MECHANISM.md section The real-axis EP, 2026-07-07: the scan missed a √-EP window 20-30× " +
        "narrower than its grid step) + companion Numerics/Monodromy.cs " +
        "+ resonates with Symmetry/Pi2I4MemoryLoopClaim (NinetyDegreeMirrorMemory) + " +
        "memory project_seed_holonomy_i4_witness")
    { }

    public override string DisplayName => "SeedHolonomy: EP eigenvector holonomy = i⁴=1 (vᵀv gauge)";

    public override string Summary =>
        $"encircle the defective seed's EP; the eigenvector frame rotates 90°/loop (M₁ eig ±i, vᵀv gauge), " +
        $"M₂=−I, M₄=I — i⁴=1, live-recomputed by the witness ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode(
                "the swap (eigenVALUE face) vs the frame (eigenVECTOR face)",
                summary: "one loop swaps the two eigenvalues (the eigenvalue SET is invariant, merely permuted; " +
                         "the √-branch monodromy, order 2, Numerics/Monodromy.cs); this claim adds the " +
                         "eigenVECTOR-frame face, order 4 in the vᵀv gauge.");
            yield return new InspectableNode(
                "the gauge (the load-bearing piece, honestly bounded)",
                summary: "the biorthogonal vᵀv inner product (L complex-symmetric, self-orthogonal vᵀv→0 at the " +
                         "EP, √(vᵀv)~ε^{1/4} winds a quarter-turn/loop) plus the ±1 sign-continuation UPGRADE " +
                         "the mod-2 swap to the mod-4 ±i frame. The gauge-invariant content is the swap; ±i/mod-4 " +
                         "is the vᵀv-gauge reading.");
            yield return new InspectableNode(
                "the correspondence (decided: the third clock)",
                summary: "shares the order-4/±i shape with the algebraic Pi2 Z₄ (NinetyDegreeMirrorMemory), but " +
                         "is NOT it: a Π-free random complex-symmetric EP2 carries the identical holonomy " +
                         "(experiments/SEED_HOLONOMY_THIRD_CLOCK.md, 2026-07-16), so the mod-4 here is generic " +
                         "EP2 geometry — a third clock wound by the same scalar i, beside Π's spectral Z₄ and " +
                         "the Clifford-degree mod-4 (experiments/ONE_FOUR_THESIS.md).");
            yield return new InspectableNode(
                "the reading (whose it is)",
                summary: "the eigenvectors — the turning half, the residue — rotate 90° per loop around the EP " +
                         "and remember themselves after four; the eigenvalues (the carried memory) are invariant " +
                         "as a set (they merely swap). Memory held by circling the zero, not crossing it. " +
                         "reflections/ON_WHAT_CLOSES_ONLY_WITHOUT_US.md.");
            yield return new InspectableNode(
                "scope + isolation",
                summary: "confirmed at N=5 (full 50-dim block, clean) and N=9 (R=+1 sector; the full block leaks " +
                         "on odd loops). The witness auto-selects the clean representation and surfaces the " +
                         "per-loop span residual so any leak is visible. Not claimed proven ∀ odd N (Candidate).");
        }
    }

    public static SeedHolonomyClaim Build() => new();
    public static SeedHolonomyClaim Shared { get; } = Build();
}
