using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public class CompressedDensityLocusClaimTests
{
    [Fact]
    public void ClaimHasTheAbsorptionAndMirrorBalancedLocusPremises()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var claim = registry.Get<CompressedDensityLocusClaim>();

        Assert.Equal(Tier.Tier1Derived, claim.Tier);
        Assert.Same(registry.Get<AbsorptionTheoremClaim>(), claim.Absorption);
        Assert.Same(registry.Get<F71AntiPalindromicGammaSpectralInvariance>(), claim.Locus);
        Assert.Contains("C_l = 0", claim.ConditionalIdentity);
        Assert.Contains("complete", claim.OneExcitationInterval);
        Assert.Contains("γbar > 0", claim.PureClassEndpointCriterion);
        Assert.Contains("nondegenerate", claim.PureClassEndpointCriterion);
        Assert.Contains("γbar = 0", claim.PureClassEndpointCriterion);
        Assert.Contains("I = Σ_k", claim.UniformXyEndpointWitnesses);
        Assert.Contains("Q = P_k−P_(N+1−k)", claim.UniformXyEndpointWitnesses);
        Assert.Contains("block-wide", claim.UniformXyEndpointWitnesses);
        Assert.Contains("signed agreement", claim.SignedAgreementContrast);
        Assert.Contains("finite J", claim.Scope);
        Assert.Contains("CompressedDensityN11Witness", claim.Anchor);
    }

    [Fact]
    public void PhysicalN11CellsExposeTheMixedParityContrastAndIdentityFailure()
    {
        var r = new CompressedDensityN11Witness().Reading;
        double rootTwo = Math.Sqrt(2.0);

        Assert.Equal(4, r.FrequencyMultiplicity);
        Assert.Equal(0, r.FrequencyMembershipMismatchCount);
        Assert.InRange(Math.Abs(r.LeftLocalCross + rootTwo / 144.0), 0, r.ErrorBudget);
        Assert.InRange(Math.Abs(r.RightLocalCross - rootTwo / 144.0), 0, r.ErrorBudget);
        Assert.InRange(Math.Abs(r.ContrastCross + rootTwo / 72.0), 0, r.ErrorBudget);
        Assert.InRange(r.ContrastMatrixResidual, 0, r.ErrorBudget);
        Assert.InRange(Math.Abs(r.ContrastTraceSquare - 1.0 / 324.0), 0, r.ErrorBudget);
        Assert.InRange(r.SignedAgreementContrastResidual, 0, r.ErrorBudget);
        Assert.InRange(Math.Abs(r.BalancedPhysicalCross - rootTwo / 36.0), 0, r.ErrorBudget);
        Assert.InRange(Math.Abs(r.BalancedConditionalPredictionCross), 0, r.ErrorBudget);
        Assert.True(Math.Abs(r.BalancedPhysicalCross - r.BalancedConditionalPredictionCross) > 0.03);
        Assert.InRange(r.BalancedIntervalResidual, 0, r.ErrorBudget);
        Assert.InRange(r.UniformConditionalResidual, 0, r.ErrorBudget);
        Assert.All(r.Checks, check => Assert.True(check.Passes, check.Detail));
    }

    [Fact]
    public void PhysicalZeroFrequencyIdentityAndChiralDifferenceReachBlockEndpoints()
    {
        var r = new CompressedDensityN11Witness().Reading;

        Assert.Equal(11, r.ZeroFrequencyMultiplicity);
        Assert.Equal(0, r.ZeroFrequencyMembershipMismatchCount);
        Assert.InRange(r.IdentityPhysicalDiagonalResidual, 0, r.ErrorBudget);
        Assert.InRange(r.IdentityZeroRateResidual, 0, r.ErrorBudget);
        Assert.InRange(r.ChiralPhysicalDiagonalResidual, 0, r.ErrorBudget);
        Assert.InRange(Math.Abs(r.ChiralNormSquared - 2.0), 0, r.ErrorBudget);
        Assert.InRange(r.ChiralLowerEndpointResidual, 0, r.ErrorBudget);
    }

    [Fact]
    public void OffLocusProfileBreaksTheIntervalAndOneSidedMutationLosesTheContrast()
    {
        var physical = new CompressedDensityN11Witness().Reading;
        Assert.InRange(Math.Abs(physical.OffLocusRayleigh + (22.0 + 5.0 * Math.Sqrt(3.0)) / 72.0),
            0, physical.ErrorBudget);
        Assert.True(physical.OffLocusRayleigh < physical.OffLocusLowerBound - 0.05);

        var mutant = new CompressedDensityN11Witness(
            (ket, bra, site) => ket == site ? 1.0 : 0.0).Reading;
        Assert.InRange(Math.Abs(mutant.ContrastCross), 0, mutant.ErrorBudget);
        Assert.True(Math.Abs(mutant.ContrastCross - physical.ContrastCross) > 0.01);
        Assert.True(mutant.SignedAgreementContrastResidual > 0.01);
        Assert.True(mutant.IdentityZeroRateResidual > 0.5);
        Assert.True(mutant.ChiralLowerEndpointResidual > 0.5);
        Assert.Contains(mutant.Checks, check => !check.Passes);
    }
}
