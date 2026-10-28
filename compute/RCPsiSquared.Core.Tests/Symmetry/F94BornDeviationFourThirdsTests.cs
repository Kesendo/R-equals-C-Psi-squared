using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F94BornDeviationFourThirdsTests
{
    private static F94BornDeviationFourThirds BuildClaim() => new();

    [Fact]
    public void NamedSetupAndExactCoefficientArePinned()
    {
        var f = BuildClaim();
        Assert.Equal(4, F94BornDeviationFourThirds.N);
        Assert.Equal(0, F94BornDeviationFourThirds.RetainedSiteA);
        Assert.Equal(2, F94BornDeviationFourThirds.RetainedSiteB);
        Assert.Equal(8, F94BornDeviationFourThirds.Sym3PartialTraceInteger);
        Assert.Equal(6, F94BornDeviationFourThirds.TaylorThreeFactorial);
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
        Assert.True(f.CellCountsSumToSurvivingDiagrams());
        Assert.True(f.StructuralDecompositionRecoversSym3());
    }

    [Theory]
    [InlineData(20.0, 0.0143)]
    [InlineData(10.0, 0.005)]
    [InlineData(100.0, 0.005)]
    [InlineData(40.0, 0.005)]
    public void DeltaDominant_IsFourThirdsQSquaredKCubed(double q, double k)
    {
        // 8/6 and 4/3 round to the same binary64, and both sides multiply in the same order,
        // so the closed form and the typed route agree exactly.
        Assert.Equal((4.0 / 3.0) * q * q * k * k * k, BuildClaim().DeltaDominant(q, k));
    }

    [Theory]
    [InlineData(20.0, 0.0143)]
    [InlineData(10.0, 0.005)]
    public void C_DominantOutcome_IsOnePlusDelta(double q, double k)
    {
        // The leading-order ratio C_|00⟩ = 1 + Δ_|00⟩ of the R_i = C_i·Ψ_i² reading.
        var f = BuildClaim();
        Assert.Equal(1.0 + (4.0 / 3.0) * q * q * k * k * k, f.C_DominantOutcome(q, k));
        Assert.True(f.C_DominantOutcome(q, k) > 1.0);
        Assert.Equal(1.0, f.C_DominantOutcome(q, 0.0));
        Assert.Equal(1.0, f.C_DominantOutcome(0.0, k));
    }

    [Fact]
    public void ClaimIsParentlessAndHasNoPolarityAxis()
    {
        Assert.False(typeof(IZ2AxisClaim).IsAssignableFrom(typeof(F94BornDeviationFourThirds)));
        var parents = typeof(F94BornDeviationFourThirds)
            .GetProperties()
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType));
        Assert.Empty(parents);
        Assert.NotNull(typeof(F94BornDeviationFourThirds).GetConstructor(Type.EmptyTypes));
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
