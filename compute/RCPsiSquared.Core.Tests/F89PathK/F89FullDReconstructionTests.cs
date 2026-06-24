using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>FULL Option D for path-4/5/6: reconstruct the AT factor from the rate-confined invariant
/// subspace (NO F_d import), divide it out of the live block charpoly to isolate F_d (the triple),
/// and confirm F_d equals the committed oracle (a cross-check only — the witness never imports it).</summary>
public class F89FullDReconstructionTests
{
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void ForPathK_ReconstructsAt_IsolatesFd_EqualToOracle(int k)
    {
        var block = F89PathKSeDeBlock.BuildTwoTimesSymBlock(q0: 2, nBlock: k + 1);
        var charpoly = GaussianMatrixCharpoly.Characteristic(block);

        var at = F89AtFactorReconstruction.ForPathK(k);                 // reconstructed, no oracle
        var fd = F89HbMixedIsolation.Isolate(charpoly, at, F89PathKFdOracle.Degree(k));   // F_d = C / AT

        Assert.Equal(F89PathKFdOracle.FdScaled(k), fd);                 // cross-check vs the oracle
    }
}
