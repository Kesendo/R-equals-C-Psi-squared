using System.Linq;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Wiring of <see cref="ZeroSectorImmunityClaim"/>: it lands in the default registry as
/// Tier1Derived with its four typed parents — F1PalindromeIdentity (the global palindrome it
/// refines on one block), F61 + F63 (the bit_a/bit_b sector machinery), and AbsorptionTheoremClaim
/// (the −2Σγ on the w=N corner).</summary>
public class ZeroSectorImmunityClaimRegistrationTests
{
    private static ClaimRegistry Default() => KnowledgeRegistryFactory.BuildDefault();

    [Fact]
    public void BuildDefault_ContainsClaim() =>
        Assert.True(Default().Contains<ZeroSectorImmunityClaim>());

    [Fact]
    public void TierIsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Default().Get<ZeroSectorImmunityClaim>().Tier);

    [Fact]
    public void Ancestors_ContainAllFourParents()
    {
        var ancestors = Default().AncestorsOf<ZeroSectorImmunityClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(F61BitAParityPi2Inheritance), ancestors);
        Assert.Contains(typeof(F63LCommutesPi2Pi2Inheritance), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
    }
}
