using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F86TPeakPi2InheritanceTests
{
    private const double PinnedGammaZero = 0.05;

    private static F86TPeakPi2Inheritance Build() =>
        new F86TPeakPi2Inheritance(new Pi2DyadicLadderClaim(), new TPeakLaw(PinnedGammaZero));

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, Build().Tier);
    }

    [Fact]
    public void FourFactor_IsExactlyFour_FromLadderTermMinusOne()
    {
        // F86 t_peak's "4" denominator is exactly a_{-1} = 4 on the Pi2 dyadic ladder
        // (= d² for N=1, the operator-space side: cardinality of single-qubit Pauli basis).
        var f = Build();
        Assert.Equal(4.0, f.FourFactor, precision: 14);
    }

    [Fact]
    public void OneOverFourFactor_IsExactlyQuarter_FromLadderTermThree()
    {
        // Mirror reading: the inversion a_n · a_{2-n} = 1 means a_{-1} (= 4) and a_3
        // (= 1/4) are mirror partners. Both readings of the same fact.
        var f = Build();
        Assert.Equal(0.25, f.OneOverFourFactor, precision: 14);
        Assert.Equal(1.0, f.FourFactor * f.OneOverFourFactor, precision: 14);
    }

    [Fact]
    public void GammaZero_AgreesWithParentTPeakLaw()
    {
        var f = Build();
        Assert.Equal(PinnedGammaZero, f.GammaZero, precision: 14);
    }

    [Fact]
    public void LiveTPeak_AgreesWithParentTPeakLaw()
    {
        // Composition drift check: 1 / (FourFactor · GammaZero) bit-exact equals TPeakLaw.Value.
        var f = Build();
        var parent = new TPeakLaw(PinnedGammaZero);
        Assert.Equal(parent.Value, f.LiveTPeak, precision: 14);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        Assert.Throws<ArgumentNullException>(() =>
            new F86TPeakPi2Inheritance(null!, new TPeakLaw(PinnedGammaZero)));
        Assert.Throws<ArgumentNullException>(() =>
            new F86TPeakPi2Inheritance(new Pi2DyadicLadderClaim(), null!));
    }

    [Fact]
    public void Anchor_References_F86Proof_AndPi2Ladder()
    {
        var f = Build();
        Assert.Contains("PROOF_F86_QPEAK.md", f.Anchor);
        Assert.Contains("TPeakLaw.cs", f.Anchor);
        Assert.Contains("Pi2DyadicLadderClaim.cs", f.Anchor);
    }
}
