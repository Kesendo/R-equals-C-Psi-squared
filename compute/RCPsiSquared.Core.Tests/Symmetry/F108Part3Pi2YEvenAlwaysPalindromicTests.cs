using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using static RCPsiSquared.Core.Tests.Symmetry.F108TestSupport;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>F108 Part 3 closure tests: Π_5bilinear (Y-deph variant) · L · Π⁻¹ =
/// −L − 2σI bit-exactly for every Π²_Y-even non-truly Hamiltonian + Y-dephasing.
/// Mirrors the Part 1 + Part 2 test pattern; reuses shared
/// <see cref="F108TestSupport.ComputeOperatorResidual"/> with
/// <see cref="PauliLetter.Y"/>. Reproduction:
/// <c>simulations/f108_part3_y_dephasing_scan.py</c>.</summary>
public class F108Part3Pi2YEvenAlwaysPalindromicTests
{
    // ============================================================
    // Claim metadata
    // ============================================================

    [Fact]
    public void Z2Axis_IsBitB() =>
        Assert.Equal(Z2Axis.BitB, new F108Part3Pi2YEvenAlwaysPalindromic().Z2Axis);

    [Fact]
    public void BitATwin_IsNull() =>
        Assert.Null(new F108Part3Pi2YEvenAlwaysPalindromic().BitATwin);

    [Fact]
    public void BitATwinStatus_IsBitBSpecific() =>
        // Y-dephasing is intrinsically a bit_b-axis dephase (Π²_Y uses bit_b same
        // as Π²_Z per PiOperator.SquaredEigenvalue), so no meaningful bit_a-axis
        // twin exists.
        Assert.Equal(BitATwinClassification.BitBSpecific,
            new F108Part3Pi2YEvenAlwaysPalindromic().BitATwinStatus);

    [Fact]
    public void Tier_IsTier1Derived() =>
        Assert.Equal(Tier.Tier1Derived, new F108Part3Pi2YEvenAlwaysPalindromic().Tier);

    // ============================================================
    // Pi5BilinearOperator Y-deph branch algebraic properties
    // ============================================================

    [Fact]
    public void Pi5BilinearOperator_YDeph_ActOnLetter_MatchesPaperMap()
    {
        // Per PROOF_F108_PART3: I → +X, X → −I, Y → −iZ, Z → +iY.
        // Same I↔X, Y↔Z permutation as Z-deph variant; Y/Z 2-cycle phase flipped
        // from +i to −i.
        Assert.Equal((PauliLetter.X, Complex.One),
            Pi5BilinearOperator.ActOnLetter(PauliLetter.I, PauliLetter.Y));
        Assert.Equal((PauliLetter.I, -Complex.One),
            Pi5BilinearOperator.ActOnLetter(PauliLetter.X, PauliLetter.Y));
        Assert.Equal((PauliLetter.Z, new Complex(0, -1)),
            Pi5BilinearOperator.ActOnLetter(PauliLetter.Y, PauliLetter.Y));
        Assert.Equal((PauliLetter.Y, new Complex(0, 1)),
            Pi5BilinearOperator.ActOnLetter(PauliLetter.Z, PauliLetter.Y));
    }

    [Fact]
    public void Pi5BilinearOperator_YDeph_BuildFull_IsUnitary_AtN3()
    {
        var pi = Pi5BilinearOperator.BuildFull(3, PauliLetter.Y);
        long d2 = 1L << (2 * 3);
        var product = pi * pi.ConjugateTranspose();
        var identity = Matrix<Complex>.Build.DenseIdentity((int)d2);
        double err = (product - identity).FrobeniusNorm();
        Assert.True(err < 1e-12,
            $"Π_5bilinear (Y-deph) is not unitary at N=3: ‖Π Π† − I‖ = {err:E3}");
    }

    [Fact]
    public void Pi5BilinearOperator_AllThreeDephaseVariants_CacheKeyDistinct()
    {
        // Cache key (N, dephaseLetter) must distinguish Z, X, Y matrices at same N.
        var piZ = Pi5BilinearOperator.BuildFull(2, PauliLetter.Z);
        var piX = Pi5BilinearOperator.BuildFull(2, PauliLetter.X);
        var piY = Pi5BilinearOperator.BuildFull(2, PauliLetter.Y);
        Assert.True((piZ - piY).FrobeniusNorm() > 1.0,
            "Y-deph variant must differ from Z-deph variant (phase difference on Y/Z 2-cycle)");
        Assert.True((piX - piY).FrobeniusNorm() > 1.0,
            "Y-deph variant must differ from X-deph variant (different permutation)");
        Assert.True((piZ - piX).FrobeniusNorm() > 1.0,
            "Z-deph and X-deph variants must differ");
    }

    // ============================================================
    // Numerical residual: the main F108 Part 3 verification
    // ============================================================

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void YzZyChain_YDephResidualIsZero(int N)
    {
        // Canonical Π²_Y-even non-truly test: H = Σ_b (Y_b Z_{b+1} + Z_b Y_{b+1})
        // under Y-deph. Same Hamiltonian as Part 1's canonical test (since
        // Π²_Y-even = Π²_Z-even bilinear sets coincide), different dephase.
        var terms = BuildYzZyChain(N);
        double residual = ComputeOperatorResidual(N, terms, PauliLetter.Y);
        Assert.True(residual < ResidualTol,
            $"F108 Part 3 violated at N={N} on YZ+ZY chain: residual = {residual:E3}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void AllNinePi2YEvenNonTrulyPairs_YDephResidualIsZero_AtN(int N)
    {
        // Same 9 pairs as F108 Part 1 (since Π²_Y-even = Π²_Z-even bilinear set).
        var pairs = EnumeratePurePi2EvenNonTrulyPairs();
        Assert.Equal(9, pairs.Count);
        foreach (var (label, terms) in pairs)
        {
            var chainTerms = BuildChainFromBilinears(N, terms);
            double residual = ComputeOperatorResidual(N, chainTerms, PauliLetter.Y);
            Assert.True(residual < ResidualTol,
                $"F108 Part 3 violated at N={N} on '{label}': residual = {residual:E3}");
        }
    }

    [Fact]
    public void PureYDephasingDissipator_ResidualIsZero_AtN3()
    {
        // No Hamiltonian, pure D[Y]^⊗N. Per PROOF_F108_PART3, the per-site
        // dissipator identity Q · D[Y_l] · Q⁻¹ = −D[Y_l] − 2γ_l · I closes the
        // dissipator side bit-exactly.
        int N = 3;
        double residual = ComputeOperatorResidual(N, Array.Empty<PauliTerm>(), PauliLetter.Y);
        Assert.True(residual < ResidualTol,
            $"F108 Part 3 dissipator-side violated at N={N}: residual = {residual:E3}");
    }

    [Fact]
    public void MixedTrulyAndNonTrulyPi2YEven_YDephResidualIsZero_AtN3()
    {
        // H = XX + YZ per bond, under Y-dephasing. XX is truly under Y-deph;
        // YZ is non-truly. Both Π²_Y-even.
        int N = 3;
        var terms = new List<PauliTerm>();
        for (int b = 0; b < N - 1; b++)
        {
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.X, Complex.One));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Z, Complex.One));
        }
        double residual = ComputeOperatorResidual(N, terms, PauliLetter.Y);
        Assert.True(residual < ResidualTol,
            $"F108 Part 3 violated at N={N} on XX+YZ chain under Y-deph: residual = {residual:E3}");
    }
}
