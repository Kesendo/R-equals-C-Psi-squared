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
    public void F108Part1Pi2EvenAlwaysPalindromic_HasFilledBitATwin()
    {
        // After F108 Part 2 closure (2026-05-25), F108 Part 1's BitATwin slot
        // is Filled with the typed F108 Part 2 Claim. PolarityCubeMap surfaces
        // this as a closed twin pair under BitB-axis coverage.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<F108Part1Pi2EvenAlwaysPalindromic>();
        Assert.Equal(BitATwinClassification.Filled, claim.BitATwinStatus);
        Assert.IsType<F108Part2Pi2XEvenAlwaysPalindromic>(claim.BitATwin);
    }
}
