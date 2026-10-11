using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live witness for Z dephasing with correlated site fields
/// (<c>experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md</c>): with integer jumps c_μ = Σ_k B_μk Z_k
/// (Kossakowski matrix a = BᵀB + 1, every covariance nonzero) and the integer XXZ chain, the
/// palindrome residual in the Pauli basis is exactly −{K − σ, ·}, K = Σ c†c, its squared norm
/// exactly 8·4^N·Σ_{k&lt;l} a_kl²; the shared field realised as independent plus collective jumps
/// shifts every block by exactly −2κ(p − q)²; the four edge blocks pay the local bill whatever the
/// covariances. Every comparison is exact (integer and dyadic entries).</summary>
public class CorrelatedZDephasingPriceWitnessTests
{
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void PriceIdentity_HoldsExactly(int n)
    {
        var w = new CorrelatedZDephasingPriceWitness(n);
        Assert.True(w.PriceIdentityResidual == 0.0,
            $"max|M + {{K − σ, ·}}| = {w.PriceIdentityResidual:e3} must be exactly 0 at N={n}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void ScalarShift_DoesNotDoIt(int n)
    {
        var w = new CorrelatedZDephasingPriceWitness(n);
        Assert.True(w.ScalarShiftResidual > 0.5,
            $"max|M| = {w.ScalarShiftResidual:e3}: the number 2σ must fail as the shift under correlation at N={n}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void ResidualNorm_IsTheClosedForm(int n)
    {
        var w = new CorrelatedZDephasingPriceWitness(n);
        Assert.True(w.PredictedNormSquared > 0.0);
        Assert.True(w.ResidualNormSquared == w.PredictedNormSquared,
            $"‖M‖² = {w.ResidualNormSquared} must equal 8·4^N·Σa_kl² = {w.PredictedNormSquared} exactly at N={n}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void FieldControl_BreaksTheIdentity(int n)
    {
        var w = new CorrelatedZDephasingPriceWitness(n);
        Assert.True(w.FieldControlResidual > 0.5,
            $"a Z field on one site (H no longer X^N-even) must break the identity at N={n}; residual {w.FieldControlResidual:e3}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void SharedField_ShiftsEveryBlockByMinusTwoKappaDSquared(int n)
    {
        var w = new CorrelatedZDephasingPriceWitness(n);
        Assert.True(w.SharedFieldBlockLawResidual == 0.0,
            $"L_shared − L_local must be exactly −2κ(p − q)² on every block at N={n}; residual {w.SharedFieldBlockLawResidual:e3}");
        Assert.True(w.SharedFieldShiftNorm > 0.0, "the collective jump must shift some block (d ≠ 0)");
        Assert.True(w.SharedFieldCentreResidual == 0.0,
            $"(K(i) + K(j))/2 must equal the claim's BlockSigma on every block at N={n}; residual {w.SharedFieldCentreResidual:e3}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void EdgeBlocks_PayTheLocalBill_TheirImageDoesNot(int n)
    {
        var w = new CorrelatedZDephasingPriceWitness(n);
        Assert.True(w.EdgeBlockResidual == 0.0,
            $"every edge-block coherence must pay exactly −2 a_ll at N={n}; residual {w.EdgeBlockResidual:e3}");
        Assert.True(w.ImageBlockCarriesCovariance,
            $"the image (1, N) of (0, 1) must carry the covariances at N={n}");
    }

    // the four shapes: the popcounts on which K takes one value, read cell by cell (from the jumps for the witness's
    // own a), meet the claim's closed form on every pattern, over 2, 4 and 3 distinct shapes at N = 3, 4, 5, and K
    // there is the shared field's value at the mean covariance, in the integers
    [Theory]
    [InlineData(3, 2)]
    [InlineData(4, 4)]
    [InlineData(5, 3)]
    public void FourShapes_ReadCellByCell_MeetTheClosedForm(int n, int shapes)
    {
        var w = new CorrelatedZDephasingPriceWitness(n);
        Assert.Equal(4, w.FourShapes.Count);
        Assert.True(w.FourShapesAgree, string.Join("; ", w.FourShapes.Select(s =>
            $"{s.Pattern}: read {{{string.Join(",", s.Read)}}} predicted {{{string.Join(",", s.Predicted)}}}")));
        Assert.Equal(shapes, w.DistinctShapes);
        Assert.True(w.MeanCovarianceHolds, $"K on a single-value popcount must be the shared field's at the mean covariance at N={n}");
    }

    // the third moment read from L (2(L + sigma) has integer entries) is the claim's closed form, for the witness's own a
    // and for a covariance -1 between the end sites alone, where the summed covariance is negative, the H-weighted sum is
    // at most zero and the third moment is nonzero: the break the GHZ argument cannot see
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void ThirdMoment_ReadFromL_IsTheClosedForm_AndTheEndSiteControlBreaks(int n)
    {
        var w = new CorrelatedZDephasingPriceWitness(n);
        Assert.True(w.CubeTrace.Imaginary.IsZero, $"imaginary part {w.CubeTrace.Imaginary} at N={n}");
        Assert.Equal(8 * w.ThirdMomentPredicted, w.CubeTrace.Real);
        Assert.True(w.ThirdMomentNoiseTerm > 0);
        Assert.True(w.EndCovarianceCubeTrace.Imaginary.IsZero);
        Assert.Equal(8 * w.EndCovarianceThirdMomentPredicted, w.EndCovarianceCubeTrace.Real);
        Assert.True(w.EndCovarianceWeightedSum <= 0, $"weighted sum {w.EndCovarianceWeightedSum} at N={n}");
        Assert.False(w.EndCovarianceThirdMomentPredicted.IsZero);
        Assert.Empty(w.EndCovarianceFarEndPopcounts);
        Assert.True(w.ThirdMomentHolds);
    }

    // the far-end modes of a negatively shared field (jumps Z_k - Z_l for every pair plus one Z_k per site) are read live
    // as eigenvectors of L at -2 K_p; the popcounts with K_p > sigma = N^2 are the middle ones, (N - 2p)^2 < N
    [Theory]
    [InlineData(3, new[] { 1, 2 })]
    [InlineData(4, new[] { 2 })]
    [InlineData(5, new[] { 2, 3 })]
    public void FarEndModes_AreEigenvectors_AndTheMiddleBreaks(int n, int[] middle)
    {
        var w = new CorrelatedZDephasingPriceWitness(n);
        Assert.All(w.FarEndModes, m => Assert.True(m.Eigenvector, $"P_{m.Popcount} X^N at N={n}"));
        Assert.Equal((long)n * n, w.FarEndSigma);
        foreach (var m in w.FarEndModes)
            Assert.Equal((long)n * n + n - (n - 2 * m.Popcount) * (n - 2 * m.Popcount), m.K);
        Assert.Equal(middle, w.FarEndRead);
        Assert.Equal(middle, w.FarEndPredicted);
        Assert.True(w.FarEndHolds);
    }

    [Fact]
    public void Kossakowski_HasEveryCovarianceNonzero()
    {
        var w = new CorrelatedZDephasingPriceWitness(4);
        for (int k = 0; k < 4; k++)
            for (int l = 0; l < 4; l++)
                if (k != l) Assert.NotEqual(0, w.Kossakowski[k, l]);
    }

    [Fact]
    public void Children_StartWithTheAssert()
    {
        var w = new CorrelatedZDephasingPriceWitness(3);
        var children = ((IInspectable)w).Children.ToList();
        Assert.Equal(8, children.Count);
        Assert.Contains("price identity", children[0].DisplayName);
    }

    [Fact]
    public void Constructor_RejectsBadArgs()
    {
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CorrelatedZDephasingPriceWitness(2));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CorrelatedZDephasingPriceWitness(6));
    }
}
