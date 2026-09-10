using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Runtime.Tests.F1Family;

[Trait("Category", "ROUTE_B_A2_N6_UNFOLDING")]
public class RouteBN6A2UnfoldingClaimRegistrationTests
{
    private static readonly Lazy<RCPsiSquared.Runtime.ObjectManager.ClaimRegistry> Registry =
        new(KnowledgeRegistryFactory.BuildDefault);

    private static Claim RegisteredClaim()
    {
        var registry = Registry.Value;
        var claim = registry.All().SingleOrDefault(c => c.GetType().Name == "RouteBN6A2UnfoldingClaim");
        Assert.NotNull(claim);
        return claim;
    }

    [Fact]
    public void F163HasTheExactN6CountsAndLocalResponseOrders()
    {
        var claim = RegisteredClaim();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        foreach (var (name, value) in new[]
        {
            ("Sites", 6), ("FullDimension", 90), ("SectorDimension", 45),
            ("LociPerParity", 133), ("TotalLoci", 266),
            ("OneEndOrder", 1), ("EqualEndOrder", 1), ("OppositeEndOrder", 2)
        })
        {
            var field = claim.GetType().GetField(name);
            Assert.NotNull(field);
            Assert.Equal(value, field.GetRawConstantValue());
        }
    }

    [Fact]
    public void ExactlyTheTwoUsedPremisesAreSharedTypedParents()
    {
        var registry = Registry.Value;
        var claim = RegisteredClaim();
        var parents = registry.AllEdges().Where(e => e.Child == claim.GetType()).Select(e => e.Parent).ToHashSet();
        Assert.True(parents.SetEquals(new[] { typeof(MirrorOrderSortingClaim), typeof(F89CrossFoldSimilarityClaim) }));
        Assert.DoesNotContain(typeof(RouteBN5RealQSemisimpleClaim), registry.AncestorsOf(claim.GetType()).Select(c => c.GetType()));
        foreach (var property in claim.GetType().GetProperties().Where(p => typeof(Claim).IsAssignableFrom(p.PropertyType)))
            Assert.Same(registry.Get(property.PropertyType), property.GetValue(registry.Get(claim.GetType())));
    }

    [Fact]
    public void RenderedClaimNamesTheLocalComplexScopeAndTheActualPremiseUses()
    {
        var claim = RegisteredClaim();
        string prose = string.Join("\n", new[] { claim.DisplayName, claim.Summary }
            .Concat(claim.Children.Where(c => c is not Claim).SelectMany(c => new[] { c.DisplayName, c.Summary })));
        Assert.Contains("N=6", prose);
        Assert.Contains("complex-q", prose);
        Assert.Contains("local", prose);
        Assert.Contains("epsilon", prose);
        Assert.Contains("F131", prose);
        Assert.Contains("partner", prose);
        Assert.Contains("not an all-N", prose);
        Assert.DoesNotContain("physical operating point", prose);
        Assert.DoesNotContain("every N", prose);
    }

    [Fact]
    public void EndBondEpsilonProfilesAndTwoLocalEpBranchesAreExplicit()
    {
        var claim = RegisteredClaim();
        string prose = claim.Name + "\n" + string.Join("\n", claim.Children.Where(c => c is not Claim).Select(c => c.Summary));
        Assert.Contains("(1+epsilon,1)", prose);
        Assert.Contains("(1+epsilon/2,1+epsilon/2)", prose);
        Assert.Contains("(1+epsilon/2,1-epsilon/2)", prose);
        Assert.Contains("two distinct local EP2 branches", prose);
    }
}
