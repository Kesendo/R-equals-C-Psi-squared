using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F95AngleAtQuadraticZeroRegistrationTests
{
    private static ClaimRegistry Build() => new ClaimRegistryBuilder()
        .RegisterF95AngleAtQuadraticZero()
        .Build();

    [Fact]
    public void RegistrationAddsClaim()
    {
        Assert.True(Build().Contains<F95AngleAtQuadraticZero>());
    }

    [Fact]
    public void RegisteredClaimDoesNotImplementIZ2AxisClaim()
    {
        Assert.False(Build().Get<F95AngleAtQuadraticZero>() is IZ2AxisClaim);
    }

    [Fact]
    public void RegisteredClaimHasNoAncestors()
    {
        Assert.Empty(Build().AncestorsOf<F95AngleAtQuadraticZero>());
    }

    [Fact]
    public void RegisteredFormulaHasPositiveBBranches()
    {
        var claim = Build().Get<F95AngleAtQuadraticZero>();
        Assert.Equal(0.0, claim.ThetaGeneral(0.25, 0.5));
        Assert.Equal(Math.PI / 6.0, claim.ThetaGeneral(1.0 / 3.0, 0.5), precision: 14);
        Assert.True(double.IsNaN(claim.ThetaGeneral(0.2, 0.5)));
    }

    [Fact]
    public void RegisteredClaimHasNoPublicClaimProperties()
    {
        var parents = typeof(F95AngleAtQuadraticZero)
            .GetProperties()
            .Where(property => typeof(Claim).IsAssignableFrom(property.PropertyType));
        Assert.Empty(parents);
    }
}
