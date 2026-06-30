using RCPsiSquared.Core.F89PathK;
using Xunit;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>The AT-locked invariant-subspace BASIS (the path-6/N=7 exact-residual fix). The same subspace
/// whose characteristic polynomial <see cref="F89AtFactorReconstruction.ForPathK"/> returns as the AT factor
/// (validated against the oracle by <see cref="F89FullDReconstructionTests"/>); here we pin its shape so the
/// diabolic scout can split the block exactly (AT ⊎ residual) instead of by the nearest-match partition that
/// floods at F_53 density. The numeric split (AT ⊎ residual = full spectrum) is gated in the Diagnostics
/// tests, where the EVD lives.</summary>
public class F89AtInvariantSubspaceTests
{
    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void AtInvariantSubspaceBasis_HasSymDimRows_AtDegreeColumns(int k)
    {
        var basis = F89AtFactorReconstruction.AtInvariantSubspaceBasis(k);
        Assert.Equal(F89PathKFdOracle.SymDim(k), basis.GetLength(0));    // rows = the S₂-sym block dimension
        Assert.Equal(F89PathKFdOracle.AtDegree(k), basis.GetLength(1));  // cols = the AT factor degree (8/13/22)
    }
}
