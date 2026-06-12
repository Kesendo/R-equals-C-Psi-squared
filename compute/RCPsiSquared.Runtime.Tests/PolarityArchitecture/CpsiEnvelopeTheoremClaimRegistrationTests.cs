using System.Linq;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class CpsiEnvelopeTheoremClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF25CPsiBellPlusPi2Inheritance();

    [Fact]
    public void Register_AddsClaim()
    {
        var registry = BuildBaseRegistry().RegisterCpsiEnvelopeTheoremClaim().Build();
        Assert.True(registry.Contains<CpsiEnvelopeTheoremClaim>());
    }

    [Fact]
    public void Register_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry().RegisterCpsiEnvelopeTheoremClaim().Build();
        Assert.Equal(Tier.Tier1Derived, registry.Get<CpsiEnvelopeTheoremClaim>().Tier);
    }

    [Fact]
    public void Register_AncestorsContainBothParents()
    {
        var registry = BuildBaseRegistry().RegisterCpsiEnvelopeTheoremClaim().Build();
        var ancestors = registry.AncestorsOf<CpsiEnvelopeTheoremClaim>().ToList();
        Assert.Contains(ancestors, a => a is F25CPsiBellPlusPi2Inheritance);
        Assert.Contains(ancestors, a => a is QuarterAsBilinearMaxvalClaim);
    }
}
