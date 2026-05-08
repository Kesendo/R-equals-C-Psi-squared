using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class QubitNecessityPi2InheritanceTests
{
    private static QubitNecessityPi2Inheritance Build() =>
        new QubitNecessityPi2Inheritance(new Pi2DyadicLadderClaim(), new Pi2OperatorSpaceMirrorClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void TotalPauliOpsPerSite_IsExactlyFour_FromLadderTermMinusOne()
    {
        // 4 = a_{-1} on the Pi2 ladder = d² for N=1 qubit. The "operator-space-side"
        // anchor matched to per-site Pauli cardinality {I, X, Y, Z}.
        var q = Build();
        Assert.Equal(4.0, q.TotalPauliOpsPerSite, precision: 14);
    }

    [Fact]
    public void ImmuneOpsPerSite_IsExactlyTwo_FromLadderTermZero()
    {
        // 2 = a_0 on the Pi2 ladder = d. The qubit dimension is the count of immune
        // per-site operators (diagonal in the dephasing eigenbasis: I, Z for Z-dephasing).
        var q = Build();
        Assert.Equal(2.0, q.ImmuneOpsPerSite, precision: 14);
    }

    [Fact]
    public void DecayingOpsPerSite_EqualsImmuneOpsPerSite()
    {
        // The d² − 2d = 0 bijection at d=2: immune count = decaying count = d.
        // For Z-dephasing: immune = {I, Z}, decaying = {X, Y}. Both 2 elements.
        var q = Build();
        Assert.Equal(q.ImmuneOpsPerSite, q.DecayingOpsPerSite, precision: 14);
    }

    [Fact]
    public void BalancedFraction_IsExactlyOneHalf()
    {
        // ImmuneOps / TotalOps = 2/4 = 1/2 = a_2 (HalfAsStructuralFixedPoint).
        // The C=0.5 universal manifested at the per-site Pauli basis level.
        var q = Build();
        Assert.Equal(0.5, q.BalancedFraction, precision: 14);
    }

    [Fact]
    public void BijectionHolds_IsTrue()
    {
        // The d² − 2d = 0 bijection: TotalPauli = 2 · ImmuneOps AND ImmuneOps = DecayingOps.
        // Both conditions hold at d=2 by the polynomial identity.
        Assert.True(Build().BijectionHolds);
    }

    [Fact]
    public void MirrorPinnedTotalOps_AgreesWithLadderTermMinusOne()
    {
        // Cross-check: Pi2OperatorSpaceMirror's pinned d²=4 for N=1 qubit must agree
        // with the ladder's a_{-1} = 4 reading.
        var q = Build();
        Assert.Equal(q.TotalPauliOpsPerSite, q.MirrorPinnedTotalOps, precision: 14);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new QubitNecessityPi2Inheritance(null!, new Pi2OperatorSpaceMirrorClaim()));
        Assert.Throws<ArgumentNullException>(() =>
            new QubitNecessityPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }

    [Fact]
    public void Anchor_References_QubitNecessity_AndMirrorSymmetryProof_AndPi2Foundation()
    {
        var q = Build();
        Assert.Contains("QUBIT_NECESSITY.md", q.Anchor);
        Assert.Contains("MIRROR_SYMMETRY_PROOF.md", q.Anchor);
        Assert.Contains("Pi2KnowledgeBaseClaims.cs", q.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", q.Anchor);
        Assert.Contains("Pi2OperatorSpaceMirrorClaim.cs", q.Anchor);
    }
}
