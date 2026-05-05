using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum): elementary block-structure constants for the (n=1, n+1=2)
/// popcount coherence block of an N-qubit chain. The c=2 stratum is the smallest case where
/// an inter-channel SVD between HD=1 and HD=3 is needed (Q_EP = 2/g_eff machinery), and these
/// counts are the combinatorial scaffolding under the channel-uniform / SVD-top mode construction.
///
/// <list type="bullet">
///   <item><see cref="PnDimension"/> = N — number of HD=1 product states |10..0⟩, |01..0⟩, …</item>
///   <item><see cref="PnPlus1Dimension"/> = N(N−1)/2 — number of HD=2 product states (pairs).</item>
///   <item><see cref="HdEqualsOnePairs"/> = N(N−1) — ordered (P_n, P_{n+1}) pairs at HD=1.</item>
///   <item><see cref="HdEqualsThreePairs"/> = N(N−1)(N−2)/2 — ordered pairs at HD=3.</item>
/// </list>
///
/// <para>Tier 1 derived: the values follow directly from popcount combinatorics on N-bit
/// strings; the only assumption is c=2 (enforced by the constructor).</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (c=2).</para>
/// </summary>
public sealed class C2BlockShape : Claim
{
    public CoherenceBlock Block { get; }
    public int PnDimension { get; }
    public int PnPlus1Dimension { get; }
    public int HdEqualsOnePairs { get; }
    public int HdEqualsThreePairs { get; }

    public C2BlockShape(CoherenceBlock block)
        : base("c=2 block shape", Tier.Tier1Derived, Item1Anchors.Root)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2BlockShape applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        Block = block;
        int N = block.N;
        PnDimension = N;
        PnPlus1Dimension = N * (N - 1) / 2;
        HdEqualsOnePairs = N * (N - 1);
        HdEqualsThreePairs = N * (N - 1) * (N - 2) / 2;
    }

    public override string DisplayName =>
        $"c=2 block shape (N={Block.N})";

    public override string Summary =>
        $"|P_n|={PnDimension}, |P_{{n+1}}|={PnPlus1Dimension}, " +
        $"HD=1 pairs={HdEqualsOnePairs}, HD=3 pairs={HdEqualsThreePairs} " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return InspectableNode.RealScalar("|P_n|", PnDimension);
            yield return InspectableNode.RealScalar("|P_{n+1}|", PnPlus1Dimension);
            yield return InspectableNode.RealScalar("HD=1 pairs", HdEqualsOnePairs);
            yield return InspectableNode.RealScalar("HD=3 pairs", HdEqualsThreePairs);
        }
    }
}
