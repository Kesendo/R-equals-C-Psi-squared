using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class F108Part3Pi2YEvenAlwaysPalindromicRegistrationTests
{
    [Fact]
    public void F108Part3Pi2YEvenAlwaysPalindromic_IsRegistered_AfterBuildDefault()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<F108Part3Pi2YEvenAlwaysPalindromic>();
        Assert.NotNull(claim);
        Assert.Equal(Z2Axis.BitB, claim.Z2Axis);
    }

    [Fact]
    public void PolarityCubeMap_BitBClaims_IncludesF108Part3()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var cubeMap = registry.Get<PolarityCubeMap>();
        Assert.Contains(cubeMap.BitBClaims, c => c is F108Part3Pi2YEvenAlwaysPalindromic);
    }

    [Fact]
    public void F108Part3Pi2YEvenAlwaysPalindromic_HasBitBSpecificBitATwin()
    {
        // Y-dephasing is intrinsically a bit_b-axis dephase letter (Π²_Y uses
        // bit_b same as Π²_Z); no meaningful bit_a-axis twin exists.
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<F108Part3Pi2YEvenAlwaysPalindromic>();
        Assert.Equal(BitATwinClassification.BitBSpecific, claim.BitATwinStatus);
    }
}
