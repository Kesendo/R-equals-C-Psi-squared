using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Probes;

namespace RCPsiSquared.Core.Tests.Decomposition;

public class InterChannelSvdTests
{
    [Theory]
    [InlineData(5, 1, 2.7651, 1e-3)]
    [InlineData(6, 1, 2.8020, 1e-3)]
    [InlineData(7, 1, 2.8284, 1e-3)]
    [InlineData(8, 1, 2.8393, 1e-3)]
    public void Sigma0_Matches_PythonStepI(int N, int n, double expected, double tol)
    {
        // step_i SVD top singular value σ_0 (heuristic 2-level g_eff) for c=2 chains.
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
        Assert.InRange(svd.Sigma0, expected - tol, expected + tol);
    }

    [Theory]
    [InlineData(5, 1)]
    [InlineData(6, 1)]
    [InlineData(7, 1)]
    [InlineData(8, 1)]
    public void TopSingularVectors_AreOrthogonalToDickeProbe(int N, int n)
    {
        // step_i finding: at c=2 chains, ⟨probe | u_0⟩ = ⟨probe | v_0⟩ = 0 to numerical precision.
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var svd = InterChannelSvd.Build(block, hd1: 1, hd2: 3);
        var probe = DickeBlockProbe.Build(block);
        var ovU = probe.ConjugateDotProduct(svd.U0InFullBlock);
        var ovV = probe.ConjugateDotProduct(svd.V0InFullBlock);
        Assert.True(ovU.Magnitude < 1e-10,
            $"|⟨probe | u_0⟩| should be 0 at c=2 N={N}; got {ovU.Magnitude}");
        Assert.True(ovV.Magnitude < 1e-10,
            $"|⟨probe | v_0⟩| should be 0 at c=2 N={N}; got {ovV.Magnitude}");
    }
}
