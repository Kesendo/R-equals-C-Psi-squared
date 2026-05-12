using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F89PathK;

public class F89PathKLibTests
{
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

        var dRhoDirect = -Complex.ImaginaryOne * (H * rho - rho * H);
        for (int l = 0; l < nBlock; l++)
        {
            var Zl = PauliString.SiteOp(nBlock, l, PauliLetter.Z);
            dRhoDirect = dRhoDirect + (Complex)gamma * (Zl * rho * Zl - rho);
        }

        var vecRho = VecColumnMajor(rho);
        var vecResult = L * vecRho;
        var dRhoFromL = UnvecColumnMajor(vecResult, d);

        var diff = dRhoFromL - dRhoDirect;
        Assert.True(diff.FrobeniusNorm() < 1e-10,
            $"BuildBlockL action should match direct dρ/dt; Frobenius diff = {diff.FrobeniusNorm()}");
    }

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

        var dRhoF89 = UnvecColumnMajor(lF89 * VecColumnMajor(rho), d);
        var dRhoRow = UnvecRowMajor(lRowMajor * VecRowMajor(rho), d);

        var diff = dRhoF89 - dRhoRow;
        Assert.True(diff.FrobeniusNorm() < 1e-10,
            $"F89 (column-major) and LindbladianBuilder (row-major) should give same dρ/dt; " +
            $"Frobenius diff = {diff.FrobeniusNorm()}");
    }

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

    [Fact]
    public void BareSiteInitial01_N5_MatchesClosedForm()
    {
        // Closed form: (N-1)/(2·√(N·C(N,2))) = 4/(2·√50) = 2/(5·√2) at N=5.
        double expected = 4.0 / (2.0 * Math.Sqrt(50.0));
        double actual = F89BareSiteInitial.BareSiteInitial01(5);
        Assert.Equal(expected, actual, 14);

        double altForm = 2.0 / (5.0 * Math.Sqrt(2.0));
        Assert.Equal(altForm, actual, 14);
    }

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

    [Fact]
    public void BlockBitPos_NBlock4_IsDescendingPowersOfTwo()
    {
        var pos = F89BlockInitialRho.BlockBitPos(4);
        Assert.Equal(new[] { 8, 4, 2, 1 }, pos);
    }

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
