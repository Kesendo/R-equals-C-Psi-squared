using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.Item1Derivation;

/// <summary>F86 Item 1 (c=2 stratum): closed-form channel-uniform vectors |c_1⟩ and |c_3⟩
/// for the (n=1, n+1=2) popcount coherence block of an N-qubit chain. These are the two
/// channel-uniform basis members in the 4-mode reduction of the c=2 inter-channel SVD
/// (see PROOF_F86_QPEAK.md Item 1).
///
/// <para>The channel-uniform vector at HD=h is the equal-weight superposition over all
/// (p, q) basis pairs with <c>popcount(p XOR q) == h</c>, unit-normalised. At c=2:</para>
///
/// <list type="bullet">
///   <item><see cref="C1Vector"/>: weight <c>1/√(N(N−1))</c> on every (p, q) with HD=1, zero elsewhere.</item>
///   <item><see cref="C3Vector"/>: weight <c>1/√(N(N−1)(N−2)/2)</c> on every (p, q) with HD=3, zero elsewhere.</item>
/// </list>
///
/// <para>Tier 1 derived: the closed forms follow from the HD=1 and HD=3 pair counts at c=2
/// (see <see cref="C2BlockShape"/>) — no SVD or numerical decomposition needed. Verified
/// against <c>FourModeBasis</c> (which builds the same vectors via <c>HdChannelBasis</c>)
/// at machine precision across N=5..8.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1 (c=2).</para>
/// </summary>
public sealed class C2ChannelUniformAnalytical : Claim
{
    public CoherenceBlock Block { get; }
    public ComplexVector C1Vector { get; }
    public ComplexVector C3Vector { get; }

    public C2ChannelUniformAnalytical(CoherenceBlock block)
        : base("c=2 channel-uniform analytical", Tier.Tier1Derived,
               Item1Anchors.Root)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2ChannelUniformAnalytical applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        Block = block;
        C1Vector = BuildChannelUniform(block, hd: 1);
        C3Vector = BuildChannelUniform(block, hd: 3);
    }

    private static ComplexVector BuildChannelUniform(CoherenceBlock block, int hd)
    {
        var basis = block.Basis;
        var vec = ComplexVector.Build.Dense(basis.MTotal);
        int count = 0;
        var indices = new List<int>();
        foreach (int pState in basis.StatesP)
        {
            foreach (int qState in basis.StatesQ)
            {
                if (BitOperations.PopCount((uint)(pState ^ qState)) == hd)
                {
                    indices.Add(basis.FlatIndex(pState, qState));
                    count++;
                }
            }
        }

        if (count == 0)
            throw new InvalidOperationException(
                $"No (p, q) pairs with HD={hd} found for N={block.N}, n={block.LowerPopcount}.");

        double weight = 1.0 / Math.Sqrt(count);
        foreach (int i in indices)
            vec[i] = new Complex(weight, 0.0);
        return vec;
    }

    public override string DisplayName =>
        $"c=2 channel-uniform analytical (N={Block.N})";

    public override string Summary =>
        $"|c_1⟩ on {Block.N * (Block.N - 1)} HD=1 pairs, " +
        $"|c_3⟩ on {Block.N * (Block.N - 1) * (Block.N - 2) / 2} HD=3 pairs " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return InspectableNode.RealScalar("|c_1⟩ weight", 1.0 / Math.Sqrt(Block.N * (Block.N - 1)));
            yield return InspectableNode.RealScalar("|c_3⟩ weight", 1.0 / Math.Sqrt(Block.N * (Block.N - 1) * (Block.N - 2) / 2.0));
            yield return InspectableNode.RealScalar("|c_1⟩ L2 norm", C1Vector.L2Norm());
            yield return InspectableNode.RealScalar("|c_3⟩ L2 norm", C3Vector.L2Norm());
        }
    }
}
