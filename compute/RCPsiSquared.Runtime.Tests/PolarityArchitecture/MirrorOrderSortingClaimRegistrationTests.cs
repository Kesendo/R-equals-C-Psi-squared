using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class MirrorOrderSortingClaimRegistrationTests
{
    /// <summary>Minimal standalone registry: the claim plus its five typed parents and
    /// their own transitive parents (the antilinear triangle's mirror-group + F114 + F112
    /// chain; the F91 family's sectors + F71-refinement + inventory chain; ChiralK is a
    /// root).</summary>
    private static ClaimRegistry BuildMinimalRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterChiralK()
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
            .RegisterSymmetryFamilyInventory()
            .RegisterJointPopcountSectors()
            .RegisterF71MirrorBlockRefinement()
            .RegisterF71AntiPalindromicGammaSpectralInvariance()
            .RegisterF92BondAntiPalindromicJSpectralInvariance()
            .RegisterF93DetuningAntiPalindromicSpectralInvariance()
            .RegisterMirrorOrderSortingClaim()
            .Build();

    [Fact]
    public void RegisterMirrorOrderSortingClaim_AddsClaim()
    {
        var registry = BuildMinimalRegistry();
        Assert.True(registry.Contains<MirrorOrderSortingClaim>());
    }

    [Fact]
    public void RegisterMirrorOrderSortingClaim_TierIsTier1Derived()
    {
        var registry = BuildMinimalRegistry();
        Assert.Equal(Tier.Tier1Derived, registry.Get<MirrorOrderSortingClaim>().Tier);
    }

    [Fact]
    public void Ancestors_ContainAllFiveTypedParents()
    {
        var registry = BuildMinimalRegistry();
        var ancestors = registry.AncestorsOf<MirrorOrderSortingClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(ChiralKClaim), ancestors);
        Assert.Contains(typeof(AntilinearTriangleClaim), ancestors);
        Assert.Contains(typeof(F71AntiPalindromicGammaSpectralInvariance), ancestors);
        Assert.Contains(typeof(F92BondAntiPalindromicJSpectralInvariance), ancestors);
        Assert.Contains(typeof(F93DetuningAntiPalindromicSpectralInvariance), ancestors);
        // Transitive: the triangle's mirror group.
        Assert.Contains(typeof(MirrorGroupD4Claim), ancestors);
    }

    [Fact]
    public void Battery_AllPass_InRegisteredClaim()
    {
        var registry = BuildMinimalRegistry();
        var claim = registry.Get<MirrorOrderSortingClaim>();
        Assert.NotEmpty(claim.Cases);
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }

    [Fact]
    public void BuildDefault_ContainsMirrorOrderSortingClaim()
    {
        // The claim is wired into the production registry via KnowledgeRegistryFactory
        // (registered after its last parent AntilinearTriangleClaim, 2026-07-16).
        var registry = KnowledgeRegistryFactory.BuildDefault();
        Assert.True(registry.Contains<MirrorOrderSortingClaim>());
        var claim = registry.Get<MirrorOrderSortingClaim>();
        Assert.Equal(claim.Cases.Count, claim.PassCount);
    }
}
