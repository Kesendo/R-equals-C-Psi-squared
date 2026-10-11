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
        Assert.Equal(5, children.Count);
        Assert.Contains("price identity", children[0].DisplayName);
    }

    [Fact]
    public void Constructor_RejectsBadArgs()
    {
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CorrelatedZDephasingPriceWitness(2));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => new CorrelatedZDephasingPriceWitness(6));
    }
}
