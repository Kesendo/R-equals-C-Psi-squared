using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>One cell of the F87 dissipator-resonance table: under <paramref name="DephaseLetter"/>
/// dephasing on a chain at N=4 with k=3 Z₂³-homogeneous Pauli pairs, how many of
/// <paramref name="TotalCount"/> pairs in Klein cell <paramref name="KleinIndex"/> classify
/// as F87-hard.</summary>
public sealed record KleinCellHardnessWitness(
    PauliLetter DephaseLetter,
    (int BitA, int BitB) KleinIndex,
    int HardCount,
    int TotalCount) : IInspectable
{
    public string DisplayName =>
        $"{DephaseLetter}-deph × Klein {KleinIndex}: {HardCount}/{TotalCount} hard";

    public string Summary =>
        $"under {DephaseLetter}-dephasing, {HardCount} of {TotalCount} k=3 pairs in Klein cell {KleinIndex} are F87-hard";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("dephase letter", summary: DephaseLetter.ToString());
            yield return new InspectableNode("Klein index", summary: KleinIndex.ToString());
            yield return new InspectableNode("hard count", summary: HardCount.ToString());
            yield return new InspectableNode("total count", summary: TotalCount.ToString());
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}

/// <summary>F87 dissipator-resonance law (Tier 1 derived). The F87-hardness of Klein-homogeneous
/// pairs lives exactly in the Klein cell matching the dephase letter's Klein index:
/// <list type="bullet">
///   <item>Z-dephasing (Klein (0,1)) → hardness in Klein (0,1)</item>
///   <item>X-dephasing (Klein (1,0)) → hardness in Klein (1,0)</item>
///   <item>Y-dephasing (Klein (1,1)) → hardness in Klein (1,1)</item>
/// </list>
/// The Mother sector Klein (0,0) is universally hard-free regardless of dissipator. The zeros
/// off the matched cell hold at every body count: each of the three unmatched cells is coloured
/// by a letter string that anticommutes with every jump and commutes with the cell's strings, so
/// the spectrum pairs (PROOF_F110 §2, F158); this table is the N=4 k=3 census of that.
/// SU(2)-rotation-equivalent: the three letters give bit-identical hard counts (50/76 each)
/// in their matched cells.
///
/// <para>Verified at N=4 k=3 over 294 Z₂³-homogeneous pairs (full enumeration) per
/// dephasing letter. Source: <c>simulations/klein_dissipator_resonance.py</c>.</para>
///
/// <para>Connection to the polarity-layer reading: the hardness of Klein-homogeneous pairs
/// lives INSIDE the dissipator's Klein cell, whose index depends on the letter and is (1,1), i.e. both bits, only for Y; the Z⊗N transverse-field Brecher
/// (h_y·Y or h_x·X) breaks Z⊗N from OUTSIDE the dissipator's Klein cell (bit_a-axis of the
/// polarity layer). Brecher and Hardness are the two poles of dissipator-letter resonance.</para>
/// </summary>
public sealed class DissipatorResonanceLaw : Claim
{
    public IReadOnlyList<KleinCellHardnessWitness> Witnesses { get; }

    public DissipatorResonanceLaw()
        : base("F87 dissipator-resonance law (Tier 1; SU(2)-symmetric Klein-cell alignment)",
               Tier.Tier1Derived,
               "hypotheses/THE_POLARITY_LAYER.md (Open Directions: Dissipator-resonance law CLOSED 2026-05-01) + simulations/klein_dissipator_resonance.py")
    {
        Witnesses = StandardWitnessTable;
    }

    public override string DisplayName => "F87-hardness of Klein-homogeneous pairs aligns with dephase-letter Klein index";

    public override string Summary =>
        "the F87-hardness of Klein-homogeneous pairs lives in the Klein cell matching the dephase letter (Z→(0,1), X→(1,0), Y→(1,1)); Mother (0,0) hard-free; verified N=4 k=3 over 294 pairs × 3 letters; SU(2)-symmetric";

    /// <summary>The 4×3 witness table from full N=4 k=3 enumeration. Hard counts:
    /// diagonal (matched cells) = 50/76; off-diagonal = 0/76; Mother (0,0) = 0/66 across
    /// all letters. Source: <c>simulations/klein_dissipator_resonance.py</c>.</summary>
    public static IReadOnlyList<KleinCellHardnessWitness> StandardWitnessTable { get; } = new[]
    {
        // Mother (0,0) — universally hard-free
        new KleinCellHardnessWitness(PauliLetter.Z, (0, 0), HardCount: 0, TotalCount: 66),
        new KleinCellHardnessWitness(PauliLetter.X, (0, 0), HardCount: 0, TotalCount: 66),
        new KleinCellHardnessWitness(PauliLetter.Y, (0, 0), HardCount: 0, TotalCount: 66),

        // Klein (0,1) Z-like — diagonal under Z-deph
        new KleinCellHardnessWitness(PauliLetter.Z, (0, 1), HardCount: 50, TotalCount: 76),
        new KleinCellHardnessWitness(PauliLetter.X, (0, 1), HardCount: 0, TotalCount: 76),
        new KleinCellHardnessWitness(PauliLetter.Y, (0, 1), HardCount: 0, TotalCount: 76),

        // Klein (1,0) X-like — diagonal under X-deph
        new KleinCellHardnessWitness(PauliLetter.Z, (1, 0), HardCount: 0, TotalCount: 76),
        new KleinCellHardnessWitness(PauliLetter.X, (1, 0), HardCount: 50, TotalCount: 76),
        new KleinCellHardnessWitness(PauliLetter.Y, (1, 0), HardCount: 0, TotalCount: 76),

        // Klein (1,1) Y-like — diagonal under Y-deph
        new KleinCellHardnessWitness(PauliLetter.Z, (1, 1), HardCount: 0, TotalCount: 76),
        new KleinCellHardnessWitness(PauliLetter.X, (1, 1), HardCount: 0, TotalCount: 76),
        new KleinCellHardnessWitness(PauliLetter.Y, (1, 1), HardCount: 50, TotalCount: 76),
    };

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("source",
                summary: "simulations/klein_dissipator_resonance.py (294 Klein-homogeneous + Y-par-homogeneous k=3 pairs at N=4, full enumeration, classified per Pauli-letter dephasing)");
            yield return new InspectableNode("structural fact 1: Mother is universally hard-free",
                summary: "Klein (0,0) produces zero F87-hard cases regardless of dissipator letter, at every body count: its strings commute with every letter string, so a lit one colours the cell (PROOF_F109 Step 5)");
            yield return new InspectableNode("structural fact 2: SU(2) rotation",
                summary: "the three Pauli letters give bit-identical 50/76 hard counts in their matched diagonal cell; the dissipator picks which Klein axis hosts the hardness of Klein-homogeneous pairs");
            yield return new InspectableNode("connection to Z⊗N-Brecher",
                summary: "transverse-field Brecher breaks Z⊗N from OUTSIDE the dissipator's Klein cell (bit_a-axis); the F87-hardness of Klein-homogeneous pairs lives INSIDE the matched cell, whose index is the letter's own and is both bits only for Y; the two poles of dissipator-letter resonance act on the same polarity layer");
            yield return InspectableNode.Group("witness table (4×3)",
                Witnesses.Cast<IInspectable>().ToArray());
        }
    }
}
