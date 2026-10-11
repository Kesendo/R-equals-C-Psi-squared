using System.Linq;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using Xunit;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

/// <summary>Wiring of <see cref="CorrelatedZDephasingPriceClaim"/>: it lands in the default registry as
/// Tier1Derived with its two typed parents, F1PalindromeIdentity (the identity it generalises) and
/// AbsorptionTheoremClaim (the cell rates the shared field reads); its closed forms read the page's
/// numbers.</summary>
public class CorrelatedZDephasingPriceClaimRegistrationTests
{
    private static ClaimRegistry Default() => KnowledgeRegistryFactory.BuildDefault();

    [Fact]
    public void BuildDefault_ContainsClaim() =>
        Assert.True(Default().Contains<CorrelatedZDephasingPriceClaim>());

    [Fact]
    public void TierIsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, Default().Get<CorrelatedZDephasingPriceClaim>().Tier);

    [Fact]
    public void Ancestors_ContainBothParents()
    {
        var ancestors = Default().AncestorsOf<CorrelatedZDephasingPriceClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F1PalindromeIdentity), ancestors);
        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
    }

    [Fact]
    public void BlockSigma_ReadsThePagesExamples_TheNegativeOfThePairCentre()
    {
        // N = 4, uniform rate gamma = 10 (sigma = 40), eta = 3/10 (kappa = 3): sigma(2,2) = 28 = 2.8 gamma,
        // sigma(0,4) = 76 = 7.6 gamma; local noise (kappa = 0) gives sigma itself on every block.
        Assert.Equal(28.0, CorrelatedZDephasingPriceClaim.BlockSigma(4, 40.0, 3.0, 2, 2));
        Assert.Equal(76.0, CorrelatedZDephasingPriceClaim.BlockSigma(4, 40.0, 3.0, 0, 4));
        Assert.Equal(40.0, CorrelatedZDephasingPriceClaim.BlockSigma(4, 40.0, 0.0, 1, 3));
    }

    [Fact]
    public void ResidualNorm_ReadsTheSharedFieldForm()
    {
        // uniform shared field: sum_{k<l} a_kl^2 = N(N-1)/2 (gamma eta)^2, so 8 4^N times it is
        // 4^(N+1) N (N-1) (gamma eta)^2; at N = 3, gamma eta = 3: 4^4 * 6 * 9 = 13824.
        Assert.Equal(13824.0, CorrelatedZDephasingPriceClaim.ResidualNormSquared(3, 3 * 9.0));
    }
}
