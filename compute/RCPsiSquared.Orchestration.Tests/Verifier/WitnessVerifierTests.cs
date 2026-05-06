using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Orchestration.Verifier;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Orchestration.Tests.Verifier;

public class WitnessVerifierTests
{
    private static ChainSystem DefaultChain(int N = 5) =>
        new ChainSystem(N: N, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    [Fact]
    public void VerifyAll_F1Family_AllMatched()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain())
            .Build();

        var verifier = new WitnessVerifier(registry);
        var result = verifier.VerifyAll();

        Assert.Equal(0, result.Drifted);
        Assert.Equal(1, result.Total); // only PalindromeResidualScalingClaim implements IDriftCheckable
    }

    private sealed class DriftedSynthetic : Claim, IDriftCheckable
    {
        public DriftedSynthetic() : base("DriftedSynthetic", Tier.Tier1Derived,
            "compute/RCPsiSquared.Runtime.Tests/TestSupport/TestClaims.cs") { }
        public override string DisplayName => "DriftedSynthetic";
        public override string Summary => "always reports drift";

        public DriftReport Verify() => new DriftReport(
            ClaimType: typeof(DriftedSynthetic),
            IsDrift: true,
            Description: "synthetic forced drift for testing the verifier",
            Magnitude: 1.0);
    }

    [Fact]
    public void VerifyAll_SyntheticDrift_DetectsAndReports()
    {
        var registry = new ClaimRegistryBuilder()
            .Register<DriftedSynthetic>(_ => new DriftedSynthetic())
            .Build();

        var verifier = new WitnessVerifier(registry);
        var result = verifier.VerifyAll();

        Assert.Equal(1, result.Drifted);
        Assert.Single(result.Reports);
        Assert.True(result.Reports[0].IsDrift);
    }
}
