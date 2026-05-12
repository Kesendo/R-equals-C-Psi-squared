using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F89PathK;

public class F89PathKLibTests
{
    /// <summary>F89PathKConvention.JConventionFactor must agree with the F90 bridge identity:
    /// F89 uses H = J·(XX+YY) (no 1/2) while F86 uses H = (J/2)·(XX+YY); the documented
    /// scale factor is 2.</summary>
    [Fact]
    public void Convention_JFactor_MatchesF90Bridge()
    {
        Assert.Equal(2.0, F89PathKConvention.JConventionFactor);
        Assert.Equal(F90F86C2BridgeIdentity.JConventionFactor, F89PathKConvention.JConventionFactor);
    }

    /// <summary>BuildBlockH at n_block=2 reproduces 2·H_F86 at the same J value, since F89's
    /// H = J·(XX+YY) equals 2 times F86's H = (J/2)·(XX+YY).</summary>
    [Fact]
    public void BuildBlockH_NBlock2_MatchesF86XYChainTimesTwo()
    {
        const double J = 0.7;
        var hF89 = F89BlockHamiltonian.BuildBlockH(J, nBlock: 2);
        var hF86 = PauliHamiltonian.XYChain(N: 2, J: J).ToMatrix();
        var diff = hF89 - 2.0 * hF86;
        Assert.True(diff.FrobeniusNorm() < 1e-12,
            $"F89 H should equal 2·F86 H at same J; Frobenius diff = {diff.FrobeniusNorm()}");
    }

    /// <summary>BuildBlockL at n_block=3 reproduces dρ/dt for an arbitrary Hermitian ρ via
    /// the column-major vec convention. Cross-validated against direct -i[H,ρ] + Σ γ(ZρZ - ρ).
    /// The hand-built H uses F89 convention (no 1/2 factor).</summary>
    [Fact]
    public void BuildBlockL_NBlock3_MatchesDirectDRhoDt_ColumnMajor()
    {
        const int nBlock = 3;
        const double J = 1.0;
        const double gamma = 0.5;
        int d = 1 << nBlock;

        var L = F89BlockLiouvillian.BuildBlockL(J, gamma, nBlock);
        var H = F89BlockHamiltonian.BuildBlockH(J, nBlock);
        var rho = ArbitraryHermitian(d, seed: 42);

        // Direct dρ/dt = -i[H, ρ] + Σ_l γ (Z_l ρ Z_l - ρ)
        var dRhoDirect = -Complex.ImaginaryOne * (H * rho - rho * H);
        var idMat = Matrix<Complex>.Build.DenseIdentity(d);
        for (int l = 0; l < nBlock; l++)
        {
            var Zl = PauliString.SiteOp(nBlock, l, PauliLetter.Z);
            dRhoDirect = dRhoDirect + (Complex)gamma * (Zl * rho * Zl - rho);
        }

        // L · vec_C(ρ) should equal vec_C(dρ/dt) where vec_C[b·d + a] = ρ[a, b].
        var vecRho = VecColumnMajor(rho);
        var vecResult = L * vecRho;
        var dRhoFromL = UnvecColumnMajor(vecResult, d);

        var diff = dRhoFromL - dRhoDirect;
        Assert.True(diff.FrobeniusNorm() < 1e-10,
            $"BuildBlockL action should match direct dρ/dt; Frobenius diff = {diff.FrobeniusNorm()}");
    }

    /// <summary>BuildBlockL also matches the row-major LindbladianBuilder action when both are
    /// applied to their matching vec convention. Confirms F89 (column-major) and Core's
    /// LindbladianBuilder (row-major) produce identical physical dynamics modulo vec form.</summary>
    [Fact]
    public void BuildBlockL_AgreesWithLindbladianBuilder_OnDRhoDt()
    {
        const int nBlock = 3;
        const double J = 1.0;
        const double gamma = 0.5;
        int d = 1 << nBlock;

        var lF89 = F89BlockLiouvillian.BuildBlockL(J, gamma, nBlock);
        var H = F89BlockHamiltonian.BuildBlockH(J, nBlock);
        var cOps = new ComplexMatrix[nBlock];
        for (int l = 0; l < nBlock; l++)
        {
            // Lindblad form γ(Z_l ρ Z_l - ρ) corresponds to c_l = √γ · Z_l (since Z_l† Z_l = I).
            cOps[l] = Math.Sqrt(gamma) * PauliString.SiteOp(nBlock, l, PauliLetter.Z);
        }
        var lRowMajor = LindbladianBuilder.Build(H, cOps);

        var rho = ArbitraryHermitian(d, seed: 7);

        // Apply both, un-vec via matching conventions, compare resulting dρ/dt matrices.
        var dRhoF89 = UnvecColumnMajor(lF89 * VecColumnMajor(rho), d);
        var dRhoRow = UnvecRowMajor(lRowMajor * VecRowMajor(rho), d);

        var diff = dRhoF89 - dRhoRow;
        Assert.True(diff.FrobeniusNorm() < 1e-10,
            $"F89 (column-major) and LindbladianBuilder (row-major) should give same dρ/dt; " +
            $"Frobenius diff = {diff.FrobeniusNorm()}");
    }

    /// <summary>ComputeRhoBlockZero at n_block=3, N=5: returned matrix is Hermitian and has
    /// the closed-form trace expected from the path-k partial trace structure. Trace pulled
    /// from direct evaluation of Tr(ρ_block) = Σ_a ρ_block[a, a]; for ρ_cc construction, this
    /// is zero because both terms are pure off-diagonal in the popcount basis.</summary>
    [Fact]
    public void ComputeRhoBlockZero_IsHermitian_AndOffDiagonal()
    {
        var rho = F89BlockInitialRho.ComputeRhoBlockZero(nBlock: 3, N: 5);
        Assert.Equal(8, rho.RowCount);
        Assert.Equal(8, rho.ColumnCount);

        var antiH = rho - rho.ConjugateTranspose();
        Assert.True(antiH.FrobeniusNorm() < 1e-12,
            $"ρ_block(0) must be Hermitian; ‖ρ - ρ†‖_F = {antiH.FrobeniusNorm()}");

        // ρ_cc terms are pure off-diagonal in popcount basis (popcount(c)=0 connects 1↔2;
        // popcount(c)=1 connects 0↔1). So trace = 0.
        Complex tr = Complex.Zero;
        for (int i = 0; i < rho.RowCount; i++) tr += rho[i, i];
        Assert.True(Math.Abs(tr.Real) < 1e-12 && Math.Abs(tr.Imaginary) < 1e-12,
            $"Tr(ρ_block(0)) must be 0 for popcount-off-diagonal ρ_cc; got {tr}");
    }

    /// <summary>BareSiteInitial01 at N=5 matches the closed form (N-1)/(2·√(N·C(N,2))).
    /// Computed independently from the formula: 4 / (2·√(50)) = 4/(2·5√2) = 2/(5√2).</summary>
    [Fact]
    public void BareSiteInitial01_N5_MatchesClosedForm()
    {
        double expected = 4.0 / (2.0 * Math.Sqrt(50.0));
        double actual = F89BareSiteInitial.BareSiteInitial01(5);
        Assert.Equal(expected, actual, 14);

        // Cross-check independent simplification 2 / (5·√2).
        double altForm = 2.0 / (5.0 * Math.Sqrt(2.0));
        Assert.Equal(altForm, actual, 14);
    }

    /// <summary>PerSiteReductionMatrix shape is (n_block, d²) and acts on vec_C(ρ_block) to
    /// reproduce ReduceBlockToSite01 site-by-site.</summary>
    [Fact]
    public void PerSiteReductionMatrix_ReproducesPerSiteReduction_NBlock3()
    {
        const int nBlock = 3;
        int d = 1 << nBlock;
        var w = F89BlockSiteReduction.PerSiteReductionMatrix(nBlock);
        Assert.Equal(nBlock, w.RowCount);
        Assert.Equal(d * d, w.ColumnCount);

        var rho = F89BlockInitialRho.ComputeRhoBlockZero(nBlock, N: 5);
        var vecRho = VecColumnMajor(rho);
        var perSite = w * vecRho;

        for (int l = 0; l < nBlock; l++)
        {
            var direct = F89BlockSiteReduction.ReduceBlockToSite01(rho, l, nBlock);
            var matrixForm = perSite[l];
            Assert.Equal(direct.Real, matrixForm.Real, 12);
            Assert.Equal(direct.Imaginary, matrixForm.Imaginary, 12);
        }
    }

    /// <summary>BlockBitPos returns descending powers of two: [2^(n-1), ..., 1].</summary>
    [Fact]
    public void BlockBitPos_NBlock4_IsDescendingPowersOfTwo()
    {
        var pos = F89BlockInitialRho.BlockBitPos(4);
        Assert.Equal(new[] { 8, 4, 2, 1 }, pos);
    }

    /// <summary>StateIdx round-trips bit list ↔ integer index: |0110⟩ → 6.</summary>
    [Fact]
    public void StateIdx_BigEndian_BasicCase()
    {
        var bits = new[] { 0, 1, 1, 0 };
        Assert.Equal(6, F89BlockInitialRho.StateIdx(bits));
    }

    // ---- Test helpers ----

    private static ComplexMatrix ArbitraryHermitian(int d, int seed)
    {
        var rng = new Random(seed);
        var raw = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                raw[i, j] = new Complex(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        var A = Matrix<Complex>.Build.DenseOfArray(raw);
        return 0.5 * (A + A.ConjugateTranspose());
    }

    /// <summary>Column-major vec: vec(M)[b·d + a] = M[a, b].</summary>
    private static MathNet.Numerics.LinearAlgebra.Vector<Complex> VecColumnMajor(ComplexMatrix M)
    {
        int d = M.RowCount;
        var v = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d * d);
        for (int b = 0; b < d; b++)
            for (int a = 0; a < d; a++)
                v[b * d + a] = M[a, b];
        return v;
    }

    private static ComplexMatrix UnvecColumnMajor(MathNet.Numerics.LinearAlgebra.Vector<Complex> v, int d)
    {
        var M = Matrix<Complex>.Build.Dense(d, d);
        for (int b = 0; b < d; b++)
            for (int a = 0; a < d; a++)
                M[a, b] = v[b * d + a];
        return M;
    }

    /// <summary>Row-major vec: vec(M)[a·d + b] = M[a, b], matching LindbladianBuilder.</summary>
    private static MathNet.Numerics.LinearAlgebra.Vector<Complex> VecRowMajor(ComplexMatrix M)
    {
        int d = M.RowCount;
        var v = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d * d);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                v[a * d + b] = M[a, b];
        return v;
    }

    private static ComplexMatrix UnvecRowMajor(MathNet.Numerics.LinearAlgebra.Vector<Complex> v, int d)
    {
        var M = Matrix<Complex>.Build.Dense(d, d);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                M[a, b] = v[a * d + b];
        return M;
    }
}
