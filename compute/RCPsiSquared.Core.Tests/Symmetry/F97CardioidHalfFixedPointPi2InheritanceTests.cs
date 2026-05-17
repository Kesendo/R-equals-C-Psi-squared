using System.Numerics;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class F97CardioidHalfFixedPointPi2InheritanceTests
{
    private static F97CardioidHalfFixedPointPi2Inheritance BuildClaim() =>
        new F97CardioidHalfFixedPointPi2Inheritance(
            new HalfAsStructuralFixedPointClaim(),
            new QuarterAsBilinearMaxvalClaim(),
            new NinetyDegreeMirrorMemoryClaim(),
            new PolynomialFoundationClaim());

    [Fact]
    public void Tier_IsTier1Derived()
    {
        Assert.Equal(Tier.Tier1Derived, BuildClaim().Tier);
    }

    [Fact]
    public void B_IsExactlyOneHalf()
    {
        Assert.Equal(0.5, F97CardioidHalfFixedPointPi2Inheritance.B);
    }

    [Fact]
    public void Threshold_IsExactlyOneQuarter()
    {
        Assert.Equal(0.25, F97CardioidHalfFixedPointPi2Inheritance.Threshold);
    }

    [Fact]
    public void CardioidC_AtCusp_IsExactlyOneQuarter()
    {
        // At φ = 0: c = b·1 - b²·1 = 1/2 - 1/4 = 1/4 (real-axis cusp = F95)
        Complex c = BuildClaim().CardioidC(0.0);
        Assert.Equal(0.25, c.Real, precision: 15);
        Assert.Equal(0.0, c.Imaginary, precision: 15);
    }

    [Fact]
    public void CardioidC_AtTail_IsExactlyMinusThreeQuarters()
    {
        // At φ = π: c = b·(-1) - b²·1 = -1/2 - 1/4 = -3/4 (period-doubling boundary)
        Complex c = BuildClaim().CardioidC(Math.PI);
        Assert.Equal(-0.75, c.Real, precision: 14);
        Assert.Equal(0.0, c.Imaginary, precision: 14);
    }

    [Fact]
    public void FixedPointMagnitude_IsBitExactlyB_AtAllSampledPhi()
    {
        var f = BuildClaim();
        double[] phis = { 0.0, Math.PI / 6, Math.PI / 4, Math.PI / 3, Math.PI / 2,
                          2 * Math.PI / 3, 3 * Math.PI / 4, 5 * Math.PI / 6, Math.PI,
                          7 * Math.PI / 6, 5 * Math.PI / 4, 4 * Math.PI / 3,
                          3 * Math.PI / 2, 5 * Math.PI / 3, 7 * Math.PI / 4,
                          11 * Math.PI / 6 };
        foreach (double phi in phis)
        {
            Assert.Equal(0.5, f.FixedPointMagnitude(phi), precision: 14);
        }
    }

    [Fact]
    public void FixedPointArgument_EqualsPhiAtSampledPoints()
    {
        var f = BuildClaim();
        double[] phis = { 0.0, Math.PI / 6, Math.PI / 4, Math.PI / 3, Math.PI / 2,
                          2 * Math.PI / 3 };
        foreach (double phi in phis)
        {
            Assert.Equal(phi, f.FixedPointArgument(phi), precision: 14);
        }
    }

    [Fact]
    public void AlgebraicIdentity_HoldsAtCanonicalAngles()
    {
        var f = BuildClaim();
        double[] phis = { 0.0, Math.PI / 6, Math.PI / 4, Math.PI / 3, Math.PI / 2,
                          2 * Math.PI / 3, 3 * Math.PI / 4, 5 * Math.PI / 6, Math.PI };
        foreach (double phi in phis)
        {
            Assert.True(f.AlgebraicIdentityHolds(phi),
                $"c(φ) = z*(1 − z*) must hold bit-exact at φ = {phi}");
        }
    }

    [Fact]
    public void MagnitudeInvariantAroundCardioid_DriftCheck()
    {
        // |z*(φ)| = b for ALL φ around the cardioid
        Assert.True(BuildClaim().MagnitudeInvariantAroundCardioid());
    }

    [Fact]
    public void CuspAgreesWithF95Threshold_DriftCheck()
    {
        // At φ = 0, c = 1/4 exactly (matches F95's b² = QuarterAsBilinearMaxval)
        Assert.True(BuildClaim().CuspAgreesWithF95Threshold());
    }

    [Fact]
    public void TailAtMinusThreeQuarters_DriftCheck()
    {
        // At φ = π, c = -3/4 (boundary of the period-doubling region)
        Assert.True(BuildClaim().TailAtMinusThreeQuarters());
    }

    [Fact]
    public void CardioidC_AtPiOverTwo_HasImaginaryHalf()
    {
        // At φ = π/2: e^(iπ/2) = i, e^(iπ) = -1.
        // c = b·i - b²·(-1) = b² + i·b = 1/4 + i/2
        Complex c = BuildClaim().CardioidC(Math.PI / 2);
        Assert.Equal(0.25, c.Real, precision: 14);
        Assert.Equal(0.5, c.Imaginary, precision: 14);
    }

    [Fact]
    public void CardioidFixedPoint_AtPiOverTwo_IsImaginaryHalf()
    {
        // z*(π/2) = b·e^(iπ/2) = b·i = i/2
        Complex z = BuildClaim().CardioidFixedPoint(Math.PI / 2);
        Assert.Equal(0.0, z.Real, precision: 14);
        Assert.Equal(0.5, z.Imaginary, precision: 14);
    }

    [Fact]
    public void Constructor_RejectsNullParents()
    {
        var half = new HalfAsStructuralFixedPointClaim();
        var quarter = new QuarterAsBilinearMaxvalClaim();
        var ninety = new NinetyDegreeMirrorMemoryClaim();
        var polynomial = new PolynomialFoundationClaim();
        Assert.Throws<ArgumentNullException>(() =>
            new F97CardioidHalfFixedPointPi2Inheritance(null!, quarter, ninety, polynomial));
        Assert.Throws<ArgumentNullException>(() =>
            new F97CardioidHalfFixedPointPi2Inheritance(half, null!, ninety, polynomial));
        Assert.Throws<ArgumentNullException>(() =>
            new F97CardioidHalfFixedPointPi2Inheritance(half, quarter, null!, polynomial));
        Assert.Throws<ArgumentNullException>(() =>
            new F97CardioidHalfFixedPointPi2Inheritance(half, quarter, ninety, null!));
    }

    [Fact]
    public void TypedParents_AreExposed()
    {
        var f = BuildClaim();
        Assert.NotNull(f.Half);
        Assert.NotNull(f.Quarter);
        Assert.NotNull(f.NinetyDegree);
        Assert.NotNull(f.Polynomial);
    }

    [Fact]
    public void Anchor_References_ProofAndScripts()
    {
        var f = BuildClaim();
        Assert.Contains("PROOF_F97_CARDIOID_HALF_FIXED_POINT.md", f.Anchor);
        Assert.Contains("ANALYTICAL_FORMULAS.md F97", f.Anchor);
        Assert.Contains("_cardioid_parametrization_tier1.py", f.Anchor);
        Assert.Contains("CPSI_COMPLEX_PLANE.md", f.Anchor);
        Assert.Contains("F95AngleAtQuadraticZeroPi2Inheritance.cs", f.Anchor);
    }

    [Theory]
    [InlineData(0.0, 0.5)]
    [InlineData(0.7853981633974483, 0.5)]  // π/4
    [InlineData(1.5707963267948966, 0.5)]  // π/2
    [InlineData(3.141592653589793, 0.5)]   // π
    [InlineData(4.71238898038469, 0.5)]    // 3π/2
    public void MagnitudeOfFixedPoint_IsAlwaysOneHalf(double phi, double expected)
    {
        Assert.Equal(expected, BuildClaim().FixedPointMagnitude(phi), precision: 14);
    }
}
