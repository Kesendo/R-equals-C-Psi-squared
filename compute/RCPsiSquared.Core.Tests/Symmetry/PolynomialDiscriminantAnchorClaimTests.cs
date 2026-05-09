using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class PolynomialDiscriminantAnchorClaimTests
{
    private static PolynomialDiscriminantAnchorClaim BuildClaim()
    {
        var polynomial = new PolynomialFoundationClaim();
        var qubit = new QubitDimensionalAnchorClaim();
        var ladder = new Pi2DyadicLadderClaim();
        return new PolynomialDiscriminantAnchorClaim(polynomial, qubit, ladder);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void PolynomialCoefficients_Are1Minus2Zero()
    {
        var (a, b, c) = BuildClaim().PolynomialCoefficients;
        Assert.Equal(1.0, a, precision: 14);
        Assert.Equal(-2.0, b, precision: 14);
        Assert.Equal(0.0, c, precision: 14);
    }

    [Fact]
    public void Roots_AreZeroAndTwo()
    {
        var (root0, root2) = BuildClaim().Roots;
        Assert.Equal(0.0, root0, precision: 14);
        Assert.Equal(2.0, root2, precision: 14);
    }

    [Fact]
    public void DiscriminantViaCoefficients_IsExactlyFour()
    {
        // b² − 4ac = (-2)² − 4·1·0 = 4
        Assert.Equal(4.0, BuildClaim().DiscriminantViaCoefficients, precision: 14);
    }

    [Fact]
    public void DiscriminantViaRootSeparation_IsExactlyFour()
    {
        // (2 − 0)² = 4
        Assert.Equal(4.0, BuildClaim().DiscriminantViaRootSeparation, precision: 14);
    }

    [Fact]
    public void DiscriminantViaQubitSquared_IsExactlyFour()
    {
        // d² for d=2 = 4
        Assert.Equal(4.0, BuildClaim().DiscriminantViaQubitSquared, precision: 14);
    }

    [Fact]
    public void DiscriminantViaLadder_IsExactlyFour()
    {
        // a_{-1} = 2^(1-(-1)) = 2² = 4
        Assert.Equal(4.0, BuildClaim().DiscriminantViaLadder, precision: 14);
    }

    [Fact]
    public void LadderIndex_IsMinusOne()
    {
        Assert.Equal(-1, BuildClaim().LadderIndex);
    }

    [Fact]
    public void MirrorPartnerLadderIndex_IsThree()
    {
        // 2 − (−1) = 3 (Pi2OperatorSpaceMirror inversion symmetry n + (2−n) = 2)
        Assert.Equal(3, BuildClaim().MirrorPartnerLadderIndex);
    }

    [Fact]
    public void MirrorPartnerValue_IsExactlyOneQuarter()
    {
        // a_3 = 2^(1−3) = 2^(−2) = 1/4 (= QuarterAsBilinearMaxval)
        Assert.Equal(0.25, BuildClaim().MirrorPartnerValue, precision: 14);
    }

    [Fact]
    public void AllReadingsAgree_HoldsBitExact()
    {
        // All four discriminant readings must agree to machine precision.
        Assert.True(BuildClaim().AllReadingsAgree());
    }

    [Fact]
    public void MirrorClosureHolds_FourTimesOneQuarterIsOne()
    {
        // a_{-1} · a_3 = 4 · (1/4) = 1 (Pi2OperatorSpaceMirror inversion).
        Assert.True(BuildClaim().MirrorClosureHolds());
    }

    [Fact]
    public void Constructor_NullPolynomial_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new PolynomialDiscriminantAnchorClaim(null!, new QubitDimensionalAnchorClaim(), new Pi2DyadicLadderClaim()));
    }

    [Fact]
    public void Constructor_NullQubit_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new PolynomialDiscriminantAnchorClaim(new PolynomialFoundationClaim(), null!, new Pi2DyadicLadderClaim()));
    }

    [Fact]
    public void Constructor_NullLadder_Throws()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new PolynomialDiscriminantAnchorClaim(new PolynomialFoundationClaim(), new QubitDimensionalAnchorClaim(), null!));
    }
}
