using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F108Part1Pi2EvenAlwaysPalindromicRegistrationTests
{
    [Fact]
    public void F108Part1Pi2EvenAlwaysPalindromic_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<F108Part1Pi2EvenAlwaysPalindromic>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.BitB, claim.Z2Axis);
    }

    [Fact]
    public void PolarityCubeMap_BitBClaims_IncludesF108Part1()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Contains(cubeMap.BitBClaims, c => c is F108Part1Pi2EvenAlwaysPalindromic);
    }

    [Fact]
    public void F108Part1Pi2EvenAlwaysPalindromic_HasNeedsDerivationBitATwin()
    {
        // The BitA twin (F108 Part 2 under X-dephasing) is NOT a mechanical
        // letter-swap mirror: it requires new structural derivation work (new
        // per-site Π operator + new bilinear set + new D[X] identity). Classified
        // NeedsDerivation so PolarityCubeMap surfaces it as substantive open work,
        // not as a low-cost fill target.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<F108Part1Pi2EvenAlwaysPalindromic>();
        Assert.Equal(BitATwinClassification.NeedsDerivation, claim.BitATwinStatus);
    }
}
