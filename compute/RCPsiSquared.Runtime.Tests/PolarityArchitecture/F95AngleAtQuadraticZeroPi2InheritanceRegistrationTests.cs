using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F95AngleAtQuadraticZeroPi2InheritanceRegistrationTests
{
    private static ClaimRegistry Build() => new ClaimRegistryBuilder()
        .RegisterF95AngleAtQuadraticZeroPi2Inheritance()
        .Build();

    [Fact]
    public void RegistrationAddsClaim()
    {
        Assert.True(Build().Contains<F95AngleAtQuadraticZeroPi2Inheritance>());
    }

    [Fact]
    public void RegisteredClaimDoesNotImplementIZ2AxisClaim()
    {
        Assert.False(Build().Get<F95AngleAtQuadraticZeroPi2Inheritance>() is IZ2AxisClaim);
    }

    [Fact]
    public void RegisteredClaimHasNoAncestors()
    {
        Assert.Empty(Build().AncestorsOf<F95AngleAtQuadraticZeroPi2Inheritance>());
    }

    [Fact]
    public void RegisteredFormulaHasPositiveBBranches()
    {
        var claim = Build().Get<F95AngleAtQuadraticZeroPi2Inheritance>();
        Assert.Equal(0.0, claim.ThetaGeneral(0.25, 0.5));
        Assert.Equal(Math.PI / 6.0, claim.ThetaGeneral(1.0 / 3.0, 0.5), precision: 14);
        Assert.True(double.IsNaN(claim.ThetaGeneral(0.2, 0.5)));
    }

    [Fact]
    public void RegisteredClaimHasNoPublicClaimProperties()
    {
        var parents = typeof(F95AngleAtQuadraticZeroPi2Inheritance)
            .GetProperties()
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType));
        Assert.Empty(parents);
    }
}
