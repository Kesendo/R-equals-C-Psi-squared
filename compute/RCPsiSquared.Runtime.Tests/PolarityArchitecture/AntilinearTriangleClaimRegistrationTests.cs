using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class AntilinearTriangleClaimRegistrationTests
{
    /// <summary>Minimal standalone registry: the claim plus its three typed parents and
    /// their own transitive parents (the mirror group's cube + Klein-V₄ + F114 chain, and
    /// F112's F108 Part 1/Part 2 + BitA-twin chain).</summary>
    private static ClaimRegistry BuildMinimalRegistry() =>
        new ClaimRegistryBuilder()
            .Register<KleinFourCellClaim>(_ => new KleinFourCellClaim())
            .RegisterKleinEightCellClaim()
            .RegisterPi2KleinV4DephaseSwapGroup()
            .RegisterCommutatorDConjugationSign()
            .RegisterMirrorGroupD4Claim()
            .RegisterF108Part2Pi2XEvenAlwaysPalindromic()
            .RegisterF108Part1Pi2EvenAlwaysPalindromic()
            .RegisterLindbladBitAPiBalance()
            .RegisterLindbladBitBPiBalance()
            .RegisterAntilinearTriangleClaim()
            .Build();

    [Fact]
    public void RegisterAntilinearTriangleClaim_AddsClaim()
    {
        var registry = BuildMinimalRegistry();
        Assert.True(registry.Contains<AntilinearTriangleClaim>());
    }

    [Fact]
    public void RegisterAntilinearTriangleClaim_TierIsTier1Derived()
    {
        var registry = BuildMinimalRegistry();
        Assert.Equal(Tier.Tier1Derived, registry.Get<AntilinearTriangleClaim>().Tier);
    }

    [Fact]
    public void Ancestors_ContainAllThreeTypedParents()
    {
        var registry = BuildMinimalRegistry();
        var ancestors = registry.AncestorsOf<AntilinearTriangleClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(MirrorGroupD4Claim), ancestors);
        Assert.Contains(typeof(CommutatorDConjugationSign), ancestors);
        Assert.Contains(typeof(LindbladBitBPiBalance), ancestors);
        // Transitive: the mirror group's Klein-V₄ owner of D, and F112's F108 Part 1 foundation.
        Assert.Contains(typeof(Pi2KleinV4DephaseSwapGroup), ancestors);
        Assert.Contains(typeof(F108Part1Pi2EvenAlwaysPalindromic), ancestors);
    }

    [Fact]
    public void Battery_AllPass_InRegisteredClaim()
    {
        var registry = BuildMinimalRegistry();
        var claim = registry.Get<AntilinearTriangleClaim>();
        Assert.NotEmpty(claim.Cases);
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void BuildDefault_ContainsAntilinearTriangleClaim()
    {
        // The claim is wired into the production registry via KnowledgeRegistryFactory
        // (registered after its last parent LindbladBitBPiBalance, 2026-06-11).
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<AntilinearTriangleClaim>());
        var claim = registry.Get<AntilinearTriangleClaim>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
