using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

/// <summary>F108 Part 1 closure tests: Π_5bilinear · L · Π_5bilinear⁻¹ = −L − 2σI
/// bit-exactly for every Π²-even non-truly Hamiltonian + Z-dephasing. Mirrors the
/// numerical scan in <c>simulations/_f108_part1_pi_family_scan.py</c>.</summary>
public class F108Part1Pi2EvenAlwaysPalindromicTests
{
    private const double Gamma = 0.05;
    private const double ResidualTol = 1e-10;

    // ============================================================
    // Claim metadata
    // ============================================================

    [Fact]
    public void Z2Axis_IsBitB()
    {
        var claim = new F108Part1Pi2EvenAlwaysPalindromic();
        Assert.Equal(Z2Axis.BitB, claim.Z2Axis);
    }

    [Fact]
    public void BitATwin_IsNull()
    {
        var claim = new F108Part1Pi2EvenAlwaysPalindromic();
        Assert.Null(claim.BitATwin);
    }

    [Fact]
    public void BitATwinStatus_IsNeedsDerivation()
    {
        // F108 Part 1 under X-dephasing (F108 Part 2) is NOT a mechanical
        // letter-swap mirror: it requires a new per-site Π operator, a different
        // Π²_X-even bilinear set, re-derived per-bond anti-commutation, and a new
        // D[X] per-site identity. Classified NeedsDerivation rather than
        // TrivialNotYetTyped per BitATwinClassification's enum docstring.
        var claim = new F108Part1Pi2EvenAlwaysPalindromic();
        Assert.Equal(BitATwinClassification.NeedsDerivation, claim.BitATwinStatus);
    }

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new F108Part1Pi2EvenAlwaysPalindromic();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    // ============================================================
    // Bilinear classification predicates
    // ============================================================

    [Fact]
    public void IsPi2EvenBilinear_AcceptsAllFivePi2EvenPairs()
    {
        // Per F108 Part 1: {XX, YY, YZ, ZY, ZZ} are the five Π²_Z-even 2-site
        // bilinears that Π_5bilinear handles with exact palindrome residual.
        Assert.True(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.X, PauliLetter.X));
        Assert.True(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.Y, PauliLetter.Y));
        Assert.True(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.Y, PauliLetter.Z));
        Assert.True(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.Z, PauliLetter.Y));
        Assert.True(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.Z, PauliLetter.Z));
    }

    [Fact]
    public void IsPi2EvenBilinear_RejectsAllFourPi2OddPairs()
    {
        // {XY, XZ, YX, ZX} are Π²_Z-odd; Π_5bilinear residual is non-zero (= 8.00
        // in the 2-qubit anti-commutation check), they don't close under this Π
        // family.
        Assert.False(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.X, PauliLetter.Y));
        Assert.False(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.X, PauliLetter.Z));
        Assert.False(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.Y, PauliLetter.X));
        Assert.False(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.Z, PauliLetter.X));
    }

    [Fact]
    public void IsPi2EvenBilinear_RejectsIdentityContainingPairs()
    {
        // I-containing pairs carry no 2-body Hamiltonian; excluded by definition.
        Assert.False(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.I, PauliLetter.X));
        Assert.False(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.Z, PauliLetter.I));
        Assert.False(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinear(PauliLetter.I, PauliLetter.I));
    }

    [Fact]
    public void IsPi2EvenBilinearTerm_AcceptsYZTwoSiteAtN3()
    {
        var yz = PauliTerm.TwoSite(3, 0, PauliLetter.Y, 1, PauliLetter.Z, Complex.One);
        Assert.True(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinearTerm(yz));
    }

    [Fact]
    public void IsPi2EvenBilinearTerm_RejectsXYTwoSite()
    {
        var xy = PauliTerm.TwoSite(3, 0, PauliLetter.X, 1, PauliLetter.Y, Complex.One);
        Assert.False(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinearTerm(xy));
    }

    [Fact]
    public void IsPi2EvenBilinearTerm_RejectsSingleSiteTerm()
    {
        // Single-site terms are k=1, not 2-body.
        var xi = PauliTerm.SingleSite(3, 0, PauliLetter.X, Complex.One);
        Assert.False(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinearTerm(xi));
    }

    [Fact]
    public void IsPi2EvenBilinearHamiltonian_AcceptsYZZYChain()
    {
        // The canonical Π²-even non-truly N=3 test case from PROOF_F108_PART1.
        var terms = BuildYzZyChain(N: 3);
        Assert.True(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinearHamiltonian(terms));
    }

    [Fact]
    public void IsPi2EvenBilinearHamiltonian_RejectsMixedXYContamination()
    {
        // Adding even a single XY term breaks the guarantee.
        var terms = BuildYzZyChain(N: 3).ToList();
        terms.Add(PauliTerm.TwoSite(3, 0, PauliLetter.X, 1, PauliLetter.Y, Complex.One));
        Assert.False(F108Part1Pi2EvenAlwaysPalindromic.IsPi2EvenBilinearHamiltonian(terms));
    }

    // ============================================================
    // Pi5BilinearOperator algebraic properties
    // ============================================================

    [Fact]
    public void Pi5BilinearOperator_ActOnLetter_MatchesPaperMap()
    {
        // Per PROOF_F108_PART1: I → +X, X → −I, Y → +iZ, Z → −iY.
        Assert.Equal((PauliLetter.X, Complex.One), Pi5BilinearOperator.ActOnLetter(PauliLetter.I));
        Assert.Equal((PauliLetter.I, -Complex.One), Pi5BilinearOperator.ActOnLetter(PauliLetter.X));
        Assert.Equal((PauliLetter.Z, new Complex(0, 1)), Pi5BilinearOperator.ActOnLetter(PauliLetter.Y));
        Assert.Equal((PauliLetter.Y, new Complex(0, -1)), Pi5BilinearOperator.ActOnLetter(PauliLetter.Z));
    }

    [Fact]
    public void Pi5BilinearOperator_BuildFull_IsUnitary_AtN3()
    {
        var pi = Pi5BilinearOperator.BuildFull(3);
        long d2 = 1L << (2 * 3); // 64
        var product = pi * pi.ConjugateTranspose();
        var identity = Matrix<Complex>.Build.DenseIdentity((int)d2);
        double err = (product - identity).FrobeniusNorm();
        Assert.True(err < 1e-12, $"Π_5bilinear is not unitary at N=3: ‖Π Π† − I‖ = {err:E3}");
    }

    [Fact]
    public void Pi5BilinearOperator_PiSquared_DiagonalSignPattern_AtN1()
    {
        // Per PROOF_F108_PART1 Step 1: M² = diag(−1, −1, +1, +1) on {I, X, Y, Z}.
        // The Pi5BilinearOperator builds Π in the Pauli-string basis where the labels
        // are indexed 0=I, 1=X, 2=Z, 3=Y (per PauliLetter enum). The expected sign
        // pattern in this index order is (−1, −1, +1, +1) on (I, X, Z, Y)
        // = (−1, −1, +1, +1) since both Y and Z fall in the +1 block.
        var pi = Pi5BilinearOperator.BuildFull(1);
        var piSquared = pi * pi;
        Assert.Equal(new Complex(-1, 0), piSquared[(int)PauliLetter.I, (int)PauliLetter.I]);
        Assert.Equal(new Complex(-1, 0), piSquared[(int)PauliLetter.X, (int)PauliLetter.X]);
        Assert.Equal(new Complex(1, 0), piSquared[(int)PauliLetter.Y, (int)PauliLetter.Y]);
        Assert.Equal(new Complex(1, 0), piSquared[(int)PauliLetter.Z, (int)PauliLetter.Z]);
        // Off-diagonal entries must be zero.
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                if (i != j)
                    Assert.Equal(Complex.Zero, piSquared[i, j]);
    }

    [Fact]
    public void Pi5BilinearOperator_OrderFour_AtN1()
    {
        // M⁴ = I per PROOF_F108_PART1.
        var pi = Pi5BilinearOperator.BuildFull(1);
        var pi4 = pi * pi * pi * pi;
        var identity = Matrix<Complex>.Build.DenseIdentity(4);
        double err = (pi4 - identity).FrobeniusNorm();
        Assert.True(err < 1e-12, $"Π_5bilinear is not order 4 at N=1: ‖Π⁴ − I‖ = {err:E3}");
    }

    // ============================================================
    // Numerical residual: the main F108 Part 1 verification
    // ============================================================

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void YzZyChain_OperatorPalindromeResidualIsZero(int N)
    {
        // The canonical Π²-even non-truly test case from PROOF_F108_PART1:
        // H = Σ_b (Y_b Z_{b+1} + Z_b Y_{b+1}). Both bilinears are Π²-even, both
        // are non-truly (Klein (1,0)). Π_5bilinear must give residual = 0.
        var terms = BuildYzZyChain(N);
        double residual = ComputeOperatorResidual(N, terms);
        Assert.True(residual < ResidualTol,
            $"F108 Part 1 violated at N={N} on YZ+ZY: residual = {residual:E3}");
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void AllNinePi2EvenNonTrulyPairs_ResidualIsZero_AtN(int N)
    {
        // Mirrors the 9-pair Phase 2 enumeration in
        // simulations/_f108_part1_pi_family_scan.py. Each pair must have
        // operator-level residual = 0 under Π_5bilinear.
        var pairs = EnumeratePureP2EvenNonTrulyPairs();
        Assert.Equal(9, pairs.Count); // Sanity: exactly 9 pairs per Phase 1 of the scan.
        foreach (var (label, terms) in pairs)
        {
            var chainTerms = BuildChainFromBilinears(N, terms);
            double residual = ComputeOperatorResidual(N, chainTerms);
            Assert.True(residual < ResidualTol,
                $"F108 Part 1 violated at N={N} on '{label}': residual = {residual:E3}");
        }
    }

    [Fact]
    public void PureZDephasingDissipator_ResidualIsZero_AtN3()
    {
        // No Hamiltonian, pure D[Z]^⊗N. Per PROOF_F108_PART1 Step 2, the per-site
        // dissipator identity Q·D[Z]·Q⁻¹ = −D[Z] − 2γ·I closes the dissipator side
        // with residual = 0.
        int N = 3;
        var emptyTerms = Array.Empty<PauliTerm>();
        double residual = ComputeOperatorResidual(N, emptyTerms);
        Assert.True(residual < ResidualTol,
            $"F108 Part 1 dissipator-side violated at N={N}: residual = {residual:E3}");
    }

    [Fact]
    public void MixedTrulyAndNonTrulyPi2Even_ResidualIsZero_AtN3()
    {
        // H = XX + YZ on every bond. XX is truly (Klein (0,0)); YZ is non-truly
        // (Klein (1,0)). Both Π²-even. Π_5bilinear handles the mix exactly.
        int N = 3;
        var terms = new List<PauliTerm>();
        for (int b = 0; b < N - 1; b++)
        {
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.X, b + 1, PauliLetter.X, Complex.One));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Z, Complex.One));
        }
        double residual = ComputeOperatorResidual(N, terms);
        Assert.True(residual < ResidualTol,
            $"F108 Part 1 violated at N={N} on XX+YZ chain: residual = {residual:E3}");
    }

    // ============================================================
    // Helpers
    // ============================================================

    private static IReadOnlyList<PauliTerm> BuildYzZyChain(int N)
    {
        var terms = new List<PauliTerm>();
        for (int b = 0; b < N - 1; b++)
        {
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Y, b + 1, PauliLetter.Z, Complex.One));
            terms.Add(PauliTerm.TwoSite(N, b, PauliLetter.Z, b + 1, PauliLetter.Y, Complex.One));
        }
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

    /// <summary>Enumerate the 9 pure-Π²-even non-truly pairs from
    /// PROOF_F108_PART1: 2 single-bilinear (YZ, ZY) and 7 two-term combinations.
    /// Each pair is given as a list of bilinears that get placed on every chain
    /// bond.</summary>
    private static IReadOnlyList<(string Label, IReadOnlyList<(PauliLetter, PauliLetter)> Bilinears)>
        EnumeratePureP2EvenNonTrulyPairs()
    {
        var YZ = (PauliLetter.Y, PauliLetter.Z);
        var ZY = (PauliLetter.Z, PauliLetter.Y);
        var XX = (PauliLetter.X, PauliLetter.X);
        var YY = (PauliLetter.Y, PauliLetter.Y);
        var ZZ = (PauliLetter.Z, PauliLetter.Z);
        return new List<(string, IReadOnlyList<(PauliLetter, PauliLetter)>)>
        {
            ("YZ", new[] { YZ }),
            ("ZY", new[] { ZY }),
            ("XX+YZ", new[] { XX, YZ }),
            ("XX+ZY", new[] { XX, ZY }),
            ("YY+YZ", new[] { YY, YZ }),
            ("YY+ZY", new[] { YY, ZY }),
            ("YZ+ZY", new[] { YZ, ZY }),
            ("YZ+ZZ", new[] { YZ, ZZ }),
            ("ZY+ZZ", new[] { ZY, ZZ }),
        };
    }

    /// <summary>Compute ‖Π_5bilinear · L · Π_5bilinear⁻¹ + L + 2σ·I‖_F in the
    /// Pauli-string basis. Returns 0 (machine precision) when F108 Part 1 holds.</summary>
    private static double ComputeOperatorResidual(int N, IReadOnlyList<PauliTerm> terms)
    {
        // Build H in the standard 2^N × 2^N basis.
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

        // Build L in vec basis (d² × d²).
        var gammas = Enumerable.Repeat(Gamma, N).ToArray();
        var Lvec = PauliDephasingDissipator.BuildZ(H, gammas);

        // Convert L to the Pauli-string basis: L_pauli = T† L T / 2^N.
        // PiOperator and Pi5BilinearOperator both live in this Pauli-string basis.
        var transform = PauliBasis.VecToPauliBasisTransform(N);
        double invD = 1.0 / (1 << N);
        var Lpauli = (transform.ConjugateTranspose() * Lvec * transform) * invD;

        // Π_5bilinear in Pauli basis.
        var pi = Pi5BilinearOperator.BuildFull(N);
        var piInv = pi.ConjugateTranspose();

        // F108 Part 1 residual: Π · L · Π⁻¹ + L + 2σ I.
        long d2 = 1L << (2 * N);
        double sigma = N * Gamma;
        var residual = pi * Lpauli * piInv + Lpauli +
            (Complex)(2.0 * sigma) * Matrix<Complex>.Build.DenseIdentity((int)d2);
        return residual.FrobeniusNorm();
    }
}
