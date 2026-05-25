using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F108Part2Pi2XEvenAlwaysPalindromicRegistrationTests
{
    [Fact]
    public void F108Part2Pi2XEvenAlwaysPalindromic_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<F108Part2Pi2XEvenAlwaysPalindromic>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.BitA, claim.Z2Axis);
    }

    [Fact]
    public void PolarityCubeMap_BitAClaims_IncludesF108Part2()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Contains(cubeMap.BitAClaims, c => c is F108Part2Pi2XEvenAlwaysPalindromic);
    }

    [Fact]
    public void F108Part2Pi2XEvenAlwaysPalindromic_HasNotApplicableForThisAxisBitATwin()
    {
        // BitA-axis Claims do not have BitATwin slots (twin concept is for BitB-axis
        // Claims pointing at BitA-axis siblings).
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<F108Part2Pi2XEvenAlwaysPalindromic>();
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, claim.BitATwinStatus);
    }
}
