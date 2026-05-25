using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>F108 Part 2 closure tests: Π_5bilinear (X-deph variant) · L · Π⁻¹ =
/// −L − 2σI bit-exactly for every Π²_X-even non-truly Hamiltonian + X-dephasing.
/// Mirrors the F108 Part 1 test structure and the Python scan in
/// <c>simulations/_f108_part2_x_dephasing_scan.py</c>.</summary>
public class F108Part2Pi2XEvenAlwaysPalindromicTests
{
    private const double Gamma = 0.05;
    private const double ResidualTol = 1e-10;

    // ============================================================
    // Claim metadata
    // ============================================================

    [Fact]
    public void Z2Axis_IsBitA()
    {
        var claim = new F108Part2Pi2XEvenAlwaysPalindromic();
        Assert.Equal(Z2Axis.BitA, claim.Z2Axis);
    }

    [Fact]
    public void BitATwin_IsNull_OnBitAAxis()
    {
        // BitA-axis Claims do not have BitATwin slots (the twin concept lives on
        // BitB-axis Claims pointing at BitA siblings).
        var claim = new F108Part2Pi2XEvenAlwaysPalindromic();
        Assert.Null(claim.BitATwin);
    }

    [Fact]
    public void BitATwinStatus_IsNotApplicableForThisAxis()
    {
        var claim = new F108Part2Pi2XEvenAlwaysPalindromic();
        Assert.Equal(BitATwinClassification.NotApplicableForThisAxis, claim.BitATwinStatus);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new F108Part2Pi2XEvenAlwaysPalindromic();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    // ============================================================
    // Bilinear classification predicates
    // ============================================================

    [Fact]
    public void IsPi2XEvenBilinear_AcceptsAllFivePi2XEvenPairs()
    {
        // Per F108 Part 2: {ZZ, XX, XY, YX, YY} are the five Π²_X-even 2-site
        // bilinears that Π_5bilinear (X-deph variant) handles with exact palindrome.
        Assert.True(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.Z, PauliLetter.Z));
        Assert.True(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.X, PauliLetter.X));
        Assert.True(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.X, PauliLetter.Y));
        Assert.True(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.Y, PauliLetter.X));
        Assert.True(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.Y, PauliLetter.Y));
    }

    [Fact]
    public void IsPi2XEvenBilinear_RejectsAllFourPi2XOddPairs()
    {
        // {XZ, YZ, ZX, ZY} are Π²_X-odd (bit_a sums to 1).
        Assert.False(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.X, PauliLetter.Z));
        Assert.False(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.Y, PauliLetter.Z));
        Assert.False(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.Z, PauliLetter.X));
        Assert.False(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.Z, PauliLetter.Y));
    }

    [Fact]
    public void IsPi2XEvenBilinear_RejectsIdentityContainingPairs()
    {
        Assert.False(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.I, PauliLetter.X));
        Assert.False(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.Z, PauliLetter.I));
        Assert.False(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinear(PauliLetter.I, PauliLetter.I));
    }

    [Fact]
    public void IsPi2XEvenBilinearTerm_AcceptsXYTwoSiteAtN3()
    {
        var xy = PauliTerm.TwoSite(3, 0, PauliLetter.X, 1, PauliLetter.Y, Complex.One);
        Assert.True(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinearTerm(xy));
    }

    [Fact]
    public void IsPi2XEvenBilinearTerm_RejectsYZTwoSite()
    {
        var yz = PauliTerm.TwoSite(3, 0, PauliLetter.Y, 1, PauliLetter.Z, Complex.One);
        Assert.False(F108Part2Pi2XEvenAlwaysPalindromic.IsPi2XEvenBilinearTerm(yz));
    }

    // ============================================================
    // Pi5BilinearOperator X-deph branch algebraic properties
    // ============================================================

    [Fact]
    public void Pi5BilinearOperator_XDeph_ActOnLetter_MatchesPaperMap()
    {
        // Per PROOF_F108_PART2: I → +Z, Z → −I, X → −iY, Y → +iX.
        Assert.Equal((PauliLetter.Z, Complex.One),
            Pi5BilinearOperator.ActOnLetter(PauliLetter.I, PauliLetter.X));
        Assert.Equal((PauliLetter.I, -Complex.One),
            Pi5BilinearOperator.ActOnLetter(PauliLetter.Z, PauliLetter.X));
        Assert.Equal((PauliLetter.Y, new Complex(0, -1)),
            Pi5BilinearOperator.ActOnLetter(PauliLetter.X, PauliLetter.X));
        Assert.Equal((PauliLetter.X, new Complex(0, 1)),
            Pi5BilinearOperator.ActOnLetter(PauliLetter.Y, PauliLetter.X));
    }

    [Fact]
    public void Pi5BilinearOperator_YDeph_ThrowsNotImplemented()
    {
        // F108 Part 3 has no covering Claim; Y-deph variant should throw with a
        // clear pointer to the open work.
        Assert.Throws<NotImplementedException>(() =>
            Pi5BilinearOperator.ActOnLetter(PauliLetter.I, PauliLetter.Y));
        Assert.Throws<NotImplementedException>(() =>
            Pi5BilinearOperator.BuildFull(2, PauliLetter.Y));
    }

    [Fact]
    public void Pi5BilinearOperator_XDeph_BuildFull_IsUnitary_AtN3()
    {
        var pi = Pi5BilinearOperator.BuildFull(3, PauliLetter.X);
        long d2 = 1L << (2 * 3);
        var product = pi * pi.ConjugateTranspose();
        var identity = Matrix<Complex>.Build.DenseIdentity((int)d2);
        double err = (product - identity).FrobeniusNorm();
        Assert.True(err < 1e-12,
            $"Π_5bilinear (X-deph) is not unitary at N=3: ‖Π Π† − I‖ = {err:E3}");
    }

    [Fact]
    public void Pi5BilinearOperator_XDeph_CacheKeyDistinctFromZDeph()
    {
        // Option C: the cache key (N, dephaseLetter) must distinguish Z and X
        // matrices at the same N. Verify by checking the two matrices differ.
        var piZ = Pi5BilinearOperator.BuildFull(2, PauliLetter.Z);
        var piX = Pi5BilinearOperator.BuildFull(2, PauliLetter.X);
        double diff = (piZ - piX).FrobeniusNorm();
        Assert.True(diff > 1.0,
            $"Z-deph and X-deph Π_5bilinear should be distinct matrices at N=2; ‖Π_Z − Π_X‖ = {diff:E3}");
    }

    // ============================================================
    // Numerical residual: the main F108 Part 2 verification
    // ============================================================

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void XyChain_XDephResidualIsZero(int N)
    {
        // Canonical Π²_X-even non-truly test: H = Σ_b (X_b Y_{b+1}) under X-deph.
        // Single bilinear XY on every bond; Π²_X-even (bit_a = 1+1 = 0 mod 2),
        // non-truly (#X = #Y = 1, both odd).
        var terms = BuildSingleBilinearChain(N, PauliLetter.X, PauliLetter.Y);
        double residual = ComputeXDephOperatorResidual(N, terms);
        Assert.True(residual < ResidualTol,
            $"F108 Part 2 violated at N={N} on XY chain: residual = {residual:E3}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void AllNinePi2XEvenNonTrulyPairs_XDephResidualIsZero_AtN(int N)
    {
        // Mirrors the 9-pair Phase 2 enumeration in
        // simulations/_f108_part2_x_dephasing_scan.py.
        var pairs = EnumeratePurePi2XEvenNonTrulyPairs();
        Assert.Equal(9, pairs.Count);
        foreach (var (label, terms) in pairs)
        {
            var chainTerms = BuildChainFromBilinears(N, terms);
            double residual = ComputeXDephOperatorResidual(N, chainTerms);
            Assert.True(residual < ResidualTol,
                $"F108 Part 2 violated at N={N} on '{label}': residual = {residual:E3}");
        }
    }

    [Fact]
    public void PureXDephasingDissipator_ResidualIsZero_AtN3()
    {
        // No Hamiltonian, pure D[X]^⊗N. Per PROOF_F108_PART2 Step 3, the per-site
        // dissipator identity Q·D[X]·Q⁻¹ = −D[X] − 2γ·I closes the dissipator side
        // bit-exactly.
        int N = 3;
        var emptyTerms = Array.Empty<PauliTerm>();
        double residual = ComputeXDephOperatorResidual(N, emptyTerms);
        Assert.True(residual < ResidualTol,
            $"F108 Part 2 dissipator-side violated at N={N}: residual = {residual:E3}");
    }

    [Fact]
    public void MixedTrulyAndNonTrulyPi2XEven_XDephResidualIsZero_AtN3()
    {
        // H = ZZ + XY per bond. ZZ is truly under X-deph (#X=0 even, #Y=0 even);
        // XY is non-truly (#X=1 odd, #Y=1 odd). Both Π²_X-even.
        int N = 3;
        var terms = new List<PauliTerm>();
        for (int b = 0; b < N - 1; b++)
        {
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Z, b + 1, PauliLetter.Z, Complex.One));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.Y, Complex.One));
        }
        double residual = ComputeXDephOperatorResidual(N, terms);
        Assert.True(residual < ResidualTol,
            $"F108 Part 2 violated at N={N} on ZZ+XY chain: residual = {residual:E3}");
    }

    // ============================================================
    // Helpers
    // ============================================================

    private static IReadOnlyList<PauliTerm> BuildSingleBilinearChain(
        int N, PauliLetter a, PauliLetter c)
    {
        var terms = new List<PauliTerm>();
        for (int b = 0; b < N - 1; b++)
            terms.Add(PauliTerm.TwoSite(N, b, a, b + 1, c, Complex.One));
        return terms;
    }

    private static IReadOnlyList<PauliTerm> BuildChainFromBilinears(
        int N, IReadOnlyList<(PauliLetter, PauliLetter)> bilinears)
    {
        var terms = new List<PauliTerm>();
        for (int b = 0; b < N - 1; b++)
            foreach (var (a, c) in bilinears)
                terms.Add(PauliTerm.TwoSite(N, b, a, b + 1, c, Complex.One));
        return terms;
    }

    /// <summary>Enumerate the 9 pure-Π²_X-even non-truly pairs: 2 single-bilinear
    /// (XY, YX) and 7 two-term combinations. Each pair placed on every bond.</summary>
    private static IReadOnlyList<(string Label, IReadOnlyList<(PauliLetter, PauliLetter)> Bilinears)>
        EnumeratePurePi2XEvenNonTrulyPairs()
    {
        var XY = (PauliLetter.X, PauliLetter.Y);
        var YX = (PauliLetter.Y, PauliLetter.X);
        var ZZ = (PauliLetter.Z, PauliLetter.Z);
        var XX = (PauliLetter.X, PauliLetter.X);
        var YY = (PauliLetter.Y, PauliLetter.Y);
        return new List<(string, IReadOnlyList<(PauliLetter, PauliLetter)>)>
        {
            ("XY", new[] { XY }),
            ("YX", new[] { YX }),
            ("ZZ+XY", new[] { ZZ, XY }),
            ("ZZ+YX", new[] { ZZ, YX }),
            ("XX+XY", new[] { XX, XY }),
            ("XX+YX", new[] { XX, YX }),
            ("XY+YX", new[] { XY, YX }),
            ("XY+YY", new[] { XY, YY }),
            ("YX+YY", new[] { YX, YY }),
        };
    }

    /// <summary>Compute ‖Π_5bilinear (X-deph) · L · Π⁻¹ + L + 2σ·I‖_F in the
    /// Pauli-string basis, where L uses X-dephasing on every site.</summary>
    private static double ComputeXDephOperatorResidual(int N, IReadOnlyList<PauliTerm> terms)
    {
        ComplexMatrix H;
        if (terms.Count == 0)
        {
            int d = 1 << N;
            H = Matrix<Complex>.Build.Dense(d, d);
        }
        else
        {
            H = new PauliHamiltonian(N, terms).ToMatrix();
        }

        // Build L with X-dephasing (not Z).
        var gammas = Enumerable.Repeat(Gamma, N).ToArray();
        var Lvec = PauliDephasingDissipator.Build(H, gammas, PauliLetter.X);

        var transform = PauliBasis.VecToPauliBasisTransform(N);
        double invD = 1.0 / (1 << N);
        var Lpauli = (transform.ConjugateTranspose() * Lvec * transform) * invD;

        var pi = Pi5BilinearOperator.BuildFull(N, PauliLetter.X);
        var piInv = pi.ConjugateTranspose();

        long d2 = 1L << (2 * N);
        double sigma = N * Gamma;
        var residual = pi * Lpauli * piInv + Lpauli +
            (Complex)(2.0 * sigma) * Matrix<Complex>.Build.DenseIdentity((int)d2);
        return residual.FrobeniusNorm();
    }
}
