using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F97CardioidHalfFixedPointPi2InheritanceTests
{
    private static F97CardioidHalfFixedPointPi2Inheritance BuildClaim() =>
        new(new QuarterAsBilinearMaxvalClaim());

    [Theory]
    [InlineData(0.0)]
    [InlineData(1.5707963267948966)]
    [InlineData(3.141592653589793)]
    [InlineData(6.283185307178586)]
    public void SelectedPeriodOneRootIsMarginalAndBothRoutesAgree(double phi)
    {
        var f = BuildClaim();
        Assert.True(f.SelectedMultiplierIsMarginal(phi));
        Assert.True(f.AlgebraicIdentityHolds(phi));
        Assert.True(f.QuadraticFormulaRecoversSelectedRoot(phi));
        Assert.Equal(0.5, f.FixedPointMagnitude(phi), precision: 14);
    }

    [Fact]
    public void DenseSweepKeepsBothRoutesTogetherAndTheMagnitudeAtOneHalf()
    {
        var f = BuildClaim();
        for (int k = 0; k < 1000; k++)
        {
            double phi = 2.0 * Math.PI * k / 1000.0;
            Assert.True(f.AlgebraicIdentityHolds(phi), $"phi={phi}");
            Assert.True(f.SelectedMultiplierIsMarginal(phi), $"phi={phi}");
        }
        Assert.True(f.MagnitudeInvariantAroundCardioid());
        Assert.True(f.SquaredMagnitudeInvariantAroundCardioid());
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(1.0)]
    [InlineData(3.141592653589793)]
    public void TheClaimsOwnGateRejectsAWrongClosedForm(double phi)
    {
        // A closed form with the wrong sign on the e^(2iφ) term differs from z*(1 − z*) by
        // e^(2iφ)/2, magnitude ½ at every φ; put through the claim's own gate, it fails.
        var f = BuildClaim();
        Complex Wrong(double p)
        {
            Complex e = Complex.Exp(Complex.ImaginaryOne * p);
            return 0.5 * e + 0.25 * e * e;
        }
        Assert.False(f.AlgebraicIdentityHolds(phi, Wrong));
        Assert.True(f.AlgebraicIdentityHolds(phi, f.CardioidC));
    }

    [Theory]
    [InlineData(1e-3)]
    [InlineData(1e-5)]
    [InlineData(1e-6)]
    public void QuadraticRecoveryHoldsTowardTheCuspUnderItsErrorModel(double phi)
    {
        Assert.True(BuildClaim().QuadraticFormulaRecoversSelectedRoot(phi));
    }

    [Fact]
    public void ParametrizationPinsCuspQuarterTurnAndTail()
    {
        var f = BuildClaim();
        Assert.True(f.CuspTakesQuartersValue());
        Assert.True(f.TailAtMinusThreeQuarters());
        var atQuarterTurn = f.CardioidC(Math.PI / 2);
        Assert.Equal(0.25, atQuarterTurn.Real, precision: 14);
        Assert.Equal(0.5, atQuarterTurn.Imaginary, precision: 14);
        var z = f.CardioidFixedPoint(Math.PI / 2);
        Assert.Equal(0.0, z.Real, precision: 14);
        Assert.Equal(0.5, z.Imaginary, precision: 14);
    }

    [Fact]
    public void OtherRootIsNotGenerallyMarginal()
    {
        var f = BuildClaim();
        Complex other = f.OtherFixedPoint(Math.PI / 2);
        Assert.Equal(Math.Sqrt(5.0), F97CardioidHalfFixedPointPi2Inheritance.MultiplierMagnitude(other), precision: 14);
        Assert.NotEqual(1.0, F97CardioidHalfFixedPointPi2Inheritance.MultiplierMagnitude(other));
    }

    [Fact]
    public void OriginHasTwoDistinctFixedPointsFromTheQuadraticFormula()
    {
        Assert.True(F97CardioidHalfFixedPointPi2Inheritance.OriginFixedPointsAreDistinct());
        var (minus, plus) = F97CardioidHalfFixedPointPi2Inheritance.FixedPoints(new Complex(0.25, 0.0));
        Assert.Equal(minus, plus); // the cusp: the two fixed points meet at one half
    }

    [Fact]
    public void ArgumentUsesWrappedPrincipalBranch()
    {
        var f = BuildClaim();
        Assert.Equal(Math.PI / 2, f.FixedPointArgument(Math.PI / 2), precision: 14);
        Assert.Equal(-Math.PI / 2, f.FixedPointArgument(3 * Math.PI / 2), precision: 14);
        Assert.Equal(-1e-6, f.FixedPointArgument(2 * Math.PI - 1e-6), precision: 10);
    }

    [Fact]
    public void QuarterIsTheOneTypedParentAndThereIsNoPolarityAxis()
    {
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var f = new F97CardioidHalfFixedPointPi2Inheritance(quarter);
        Assert.Same(quarter, f.Quarter);
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(F97CardioidHalfFixedPointPi2Inheritance)));
        var parents = typeof(F97CardioidHalfFixedPointPi2Inheritance)
            .GetProperties()
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType))
            .Select(property => property.PropertyType)
            .ToArray();
        Assert.Equal(new[] { typeof(QuarterAsBilinearMaxvalClaim) }, parents);
        Assert.Throws<ArgumentNullException>(() => new F97CardioidHalfFixedPointPi2Inheritance(null!));
    }

    [Fact]
    public void TierIsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }
}
