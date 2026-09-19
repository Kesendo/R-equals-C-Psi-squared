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
    public void Register_DirectGetDependencies_AreExactlyF25AndQuarter()
    {
        var registry = BuildBaseRegistry().RegisterCpsiEnvelopeTheoremClaim().Build();
        var directEdges = registry.EdgesInto<CpsiEnvelopeTheoremClaim>().ToList();
        var directParentTypes = directEdges.Select(edge => edge.Parent).ToHashSet();

        Assert.Equal(2, directEdges.Count);
        Assert.Equal(2, directParentTypes.Count);
        Assert.Equal(
            new[]
            {
                typeof(F25CPsiBellPlusPi2Inheritance),
                typeof(QuarterAsBilinearMaxvalClaim),
            }.OrderBy(type => type.FullName),
            directParentTypes.OrderBy(type => type.FullName));
    }

    [Fact]
    public void Register_ExposesRepairedClaimInsteadOfHistoricalAbsorber()
    {
        var claim = BuildBaseRegistry().RegisterCpsiEnvelopeTheoremClaim().Build()
            .Get<CpsiEnvelopeTheoremClaim>();
        Assert.Contains("historical Envelope package is false", claim.Name);
        Assert.Contains("remains unproved", claim.Name);
        Assert.Contains("conditional convergence implication", claim.Name);
        Assert.DoesNotContain("confirmed universally", claim.Name);
    }
}
