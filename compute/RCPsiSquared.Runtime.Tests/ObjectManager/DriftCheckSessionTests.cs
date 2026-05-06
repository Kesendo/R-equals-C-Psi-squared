using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Runtime.F1Family;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.Tests.ObjectManager;

public class DriftCheckSessionTests
{
    private static ChainSystem DefaultChain(int N = 5) =>
        new ChainSystem(N: N, J: 1.0, GammaZero: 0.05,
            HType: HamiltonianType.XY, Topology: TopologyKind.Chain);

    [Fact]
    public void Verify_F73_AtN5_Main_BitExact_NoDrift()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain(N: 5))
            .Build();

        var session = new DriftCheckSession(registry);
        var report = session.Verify<PalindromeResidualScalingClaim>();
        Assert.False(report.IsDrift, $"unexpected drift: {report.Description}");
    }

    [Fact]
    public void VerifyAll_F1Family_NoDrift()
    {
        var registry = new ClaimRegistryBuilder()
            .RegisterF1Family(DefaultChain(N: 5))
            .Build();

        var session = new DriftCheckSession(registry);
        var reports = session.VerifyAll();
        Assert.All(reports, r => Assert.False(r.IsDrift, r.Description));
    }
}
