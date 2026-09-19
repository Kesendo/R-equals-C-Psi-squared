using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>Tests for the scoped two-level Q=sqrt(3) algebra and its sole F95 parent.</summary>
public class LindbladAbsorptionMatchAtSixtyDegreesClaimTests
{
    private readonly ITestOutputHelper _out;

    public LindbladAbsorptionMatchAtSixtyDegreesClaimTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Claim_IsTier1Derived()
    {
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void QValue_IsSquareRootOfThree()
    {
        Assert.Equal(Math.Sqrt(3.0), LindbladAbsorptionMatchAtSixtyDegreesClaim.QValue, precision: 15);
    }

    [Fact]
    public void CanonicalAngleDegrees_IsSixty()
    {
        Assert.Equal(60.0, LindbladAbsorptionMatchAtSixtyDegreesClaim.CanonicalAngleDegrees);
    }

    [Fact]
    public void LindbladMagnitudeOverGamma0_IsTwo_InTheNamedPositiveGammaModel()
    {
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Equal(2.0, claim.LindbladMagnitudeOverGamma0Computed(1.0), precision: 14);
    }

    [Fact]
    public void ComputedLindbladMagnitudeOverGamma0_MatchesConstantWithinTolerance()
    {
        // Drift check: √(1 + (√3)²) = √4 = 2
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Equal(2.0, claim.LindbladMagnitudeOverGamma0Computed(0.05), precision: 14);
    }

    [Fact]
    public void ComputedF95Angle_IsSixtyDegreesWithinTolerance()
    {
        // Drift check: the positive-decay quadratic gives arctan(sqrt(3)) = 60 degrees.
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Equal(60.0, claim.F95AngleAtQSqrt3Degrees(0.05), precision: 12);
    }

    [Fact]
    public void LindbladMagnitudeSquared_EqualsFourGamma0Squared_AtAnyGamma0()
    {
        // |λ_±|² = γ₀²(1+Q²) = γ₀²(1+3) = 4γ₀²
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        foreach (var gamma0 in new[] { 0.025, 0.05, 0.10, 1.0 })
        {
            var expected = 4.0 * gamma0 * gamma0;
            Assert.Equal(expected, claim.LindbladMagnitudeSquaredAtQSqrt3(gamma0), precision: 14);
        }
    }

    [Fact]
    public void PositiveDecayRoots_AreConstructedAsMinusLambdaWithPositiveB()
    {
        const double gamma0 = 0.05;
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        var lambda = claim.LiouvillianRootsAtQSqrt3(gamma0);
        var decay = claim.PositiveDecayRootsAtQSqrt3(gamma0);

        Assert.Equal(-lambda.Plus.Real, decay.FromLambdaPlus.Real, precision: 14);
        Assert.Equal(-lambda.Plus.Imaginary, decay.FromLambdaPlus.Imaginary, precision: 14);
        Assert.Equal(-lambda.Minus.Real, decay.FromLambdaMinus.Real, precision: 14);
        Assert.Equal(-lambda.Minus.Imaginary, decay.FromLambdaMinus.Imaginary, precision: 14);
        Assert.Equal(2.0 * gamma0, (decay.FromLambdaPlus + decay.FromLambdaMinus).Real, precision: 14);
        Assert.Equal(gamma0, claim.PositiveDecayAnchorB(gamma0), precision: 14);
        Assert.True(claim.PositiveDecayAnchorB(gamma0) > 0.0);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-0.05)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void RatioAngleAndRoots_RejectNonPositiveGamma0(double gamma0)
    {
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.LindbladMagnitudeSquaredAtQSqrt3(gamma0));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.LindbladMagnitudeOverGamma0Computed(gamma0));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.F95AngleAtQSqrt3Degrees(gamma0));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.LiouvillianRootsAtQSqrt3(gamma0));
        Assert.Throws<ArgumentOutOfRangeException>(() => claim.PositiveDecayRootsAtQSqrt3(gamma0));
    }

    [Fact]
    public void SoleParent_IsF95()
    {
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Equal(Tier.Tier1Derived, claim.F95.Tier);
        Assert.NotNull(claim.F95);
        Assert.Single(claim.Children.OfType<Claim>());
    }

    [Fact]
    public void QValueMatchesQAnchorMapEntry()
    {
        // The Q=√3 anchor should appear in QAnchorMap.CanonicalAnchors as Tier1Derived.
        var map = new QAnchorMap();
        var sqrt3Anchor = map.AnkerAt(LindbladAbsorptionMatchAtSixtyDegreesClaim.QValue);
        Assert.NotNull(sqrt3Anchor);
        Assert.Equal(Tier.Tier1Derived, sqrt3Anchor!.Tier);
        Assert.Equal(QBand.Peak, sqrt3Anchor.Band);
    }

    [Fact]
    public void Render_EmitsDerivationSummary()
    {
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Contains("sqrt(3)", claim.Summary);
        Assert.Contains("60", claim.Summary);
        Assert.Contains("local algebra", claim.Summary);
        _out.WriteLine($"DisplayName: {claim.DisplayName}");
        _out.WriteLine($"Summary: {claim.Summary}");
    }

    [Fact]
    public void TypedSurface_DeclinesAbsorptionAndAncestryClaims()
    {
        var claim = LindbladAbsorptionMatchAtSixtyDegreesClaim.Build();
        var surface = $"{claim.Name} {claim.DisplayName} {claim.Summary}";

        Assert.Contains("named two-level model", surface, StringComparison.OrdinalIgnoreCase);
        Assert.Contains("not an absorption-rate identity", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("one-disagreement cell cost", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("single-site rate", surface, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("Absorption Theorem rate", surface, StringComparison.OrdinalIgnoreCase);
    }
}
