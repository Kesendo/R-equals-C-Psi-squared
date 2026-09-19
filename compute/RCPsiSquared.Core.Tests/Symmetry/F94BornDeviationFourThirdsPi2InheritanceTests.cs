using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F94BornDeviationFourThirdsPi2InheritanceTests
{
    private static F94BornDeviationFourThirdsPi2Inheritance BuildClaim() => new();

    [Fact]
    public void NamedSetupAndExactCoefficientArePinned()
    {
        var f = BuildClaim();
        Assert.Equal(4, F94BornDeviationFourThirdsPi2Inheritance.N);
        Assert.Equal(0, F94BornDeviationFourThirdsPi2Inheritance.RetainedSiteA);
        Assert.Equal(2, F94BornDeviationFourThirdsPi2Inheritance.RetainedSiteB);
        Assert.Equal(8, F94BornDeviationFourThirdsPi2Inheritance.Sym3PartialTraceInteger);
        Assert.Equal(6, F94BornDeviationFourThirdsPi2Inheritance.TaylorThreeFactorial);
        Assert.Equal(8.0 / 6.0, f.Coefficient);
        Assert.Equal(Tier.Tier1Derived, f.Tier);
    }

    [Fact]
    public void QKAndPhysicalFormsAreTheSameLeadingTerm()
    {
        var f = BuildClaim();
        const double q = 2.5;
        const double k = 0.03;
        const double gamma = 0.4;
        double j = q * gamma;
        double t = k / gamma;
        Assert.Equal(f.DeltaDominant(q, k), f.DeltaP_Dominant(j, gamma, t), precision: 14);
        Assert.Equal((4.0 / 3.0) * q * q * k * k * k, f.DeltaDominant(q, k), precision: 15);
    }

    [Fact]
    public void StructuralIntegerChecksPass()
    {
        var f = BuildClaim();
        Assert.True(f.CoefficientAgreesWithSym3());
        Assert.True(f.CellCountsSumToSurvivingDiagrams());
        Assert.True(f.StructuralDecompositionRecoversSym3());
    }

    [Fact]
    public void ClaimIsParentlessAndHasNoPolarityAxis()
    {
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(F94BornDeviationFourThirdsPi2Inheritance)));
        var parents = typeof(F94BornDeviationFourThirdsPi2Inheritance)
            .GetProperties()
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType));
        Assert.Empty(parents);
        Assert.NotNull(typeof(F94BornDeviationFourThirdsPi2Inheritance).GetConstructor(Type.EmptyTypes));
    }

    [Theory]
    [InlineData(-1.0, 0.1)]
    [InlineData(0.1, -1.0)]
    [InlineData(double.NaN, 0.1)]
    [InlineData(0.1, double.PositiveInfinity)]
    public void QKFormRejectsInvalidInputs(double q, double k)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => BuildClaim().DeltaDominant(q, k));
    }
}
