using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class Pi2KleinV4DephaseSwapGroupTests
{
    private const double Tol = 1e-10;
    private readonly ITestOutputHelper _out;

    public Pi2KleinV4DephaseSwapGroupTests(ITestOutputHelper output) => _out = output;

    // ------------------------------------------------------------------
    // Claim metadata
    // ------------------------------------------------------------------

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var claim = new Pi2KleinV4DephaseSwapGroup();
        Assert.Equal(Tier.Tier1Derived, claim.Tier);
    }

    [Fact]
    public void Anchor_ReferencesBothProofDocs()
    {
        var claim = new Pi2KleinV4DephaseSwapGroup();
        Assert.Contains("PROOF_D_PI_Z_EQUALS_PI_Y_UNIVERSAL_N.md", claim.Anchor);
        Assert.Contains("PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md", claim.Anchor);
    }

    // ------------------------------------------------------------------
    // Per-site basis-order sanity at N=1
    // ------------------------------------------------------------------

    [Fact]
    public void BuildD_AtN1_IsDiagOneOneOneMinusOne()
    {
        // Basis order (I, X, Z, Y) with PauliLetter packing i = a + 2·b means Y is at
        // letter index 3. d_l = diag(1, 1, 1, −1) places the −1 on the Y basis vector.
        var D = Pi2KleinV4DephaseSwapGroup.BuildD(N: 1);
        Assert.Equal(4, D.RowCount);
        Assert.Equal(4, D.ColumnCount);
        Assert.Equal(Complex.One, D[0, 0]);   // I
        Assert.Equal(Complex.One, D[1, 1]);   // X
        Assert.Equal(Complex.One, D[2, 2]);   // Z
        Assert.Equal(-Complex.One, D[3, 3]);  // Y
        AssertDiagonal(D);
    }

    [Fact]
    public void BuildH_AtN1_IsXZPermutation_FixingIAndY()
    {
        // h fixes I (index 0) and Y (index 3); swaps X (index 1) and Z (index 2).
        var H = Pi2KleinV4DephaseSwapGroup.BuildH(N: 1);
        Assert.Equal(4, H.RowCount);
        Assert.Equal(4, H.ColumnCount);
        // I → I
        Assert.Equal(Complex.One, H[0, 0]);
        // X (col 1) → Z (row 2); Z (col 2) → X (row 1)
        Assert.Equal(Complex.One, H[2, 1]);
        Assert.Equal(Complex.One, H[1, 2]);
        // Y → Y
        Assert.Equal(Complex.One, H[3, 3]);
        // All other entries zero
        for (int r = 0; r < 4; r++)
            for (int c = 0; c < 4; c++)
            {
                bool nonzero = (r == 0 && c == 0) || (r == 3 && c == 3) ||
                               (r == 1 && c == 2) || (r == 2 && c == 1);
                if (!nonzero) Assert.Equal(Complex.Zero, H[r, c]);
            }
    }

    // ------------------------------------------------------------------
    // Involution property: D² = H² = Q_zx² = I at N = 1, 2, 3
    // ------------------------------------------------------------------

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void D_IsInvolution_AtSmallN(int N)
    {
        var D = Pi2KleinV4DephaseSwapGroup.BuildD(N);
        AssertEqualsIdentity(D * D, $"D² at N={N}");
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void H_IsInvolution_AtSmallN(int N)
    {
        var H = Pi2KleinV4DephaseSwapGroup.BuildH(N);
        AssertEqualsIdentity(H * H, $"H² at N={N}");
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void Qzx_IsInvolution_AtSmallN(int N)
    {
        var Q = Pi2KleinV4DephaseSwapGroup.BuildQzx(N);
        AssertEqualsIdentity(Q * Q, $"Q_zx² at N={N}");
    }

    // ------------------------------------------------------------------
    // Klein-V₄ closure: D · Q_zx · H = I
    // ------------------------------------------------------------------

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void KleinV4_Closure_DTimesQzxTimesH_EqualsIdentity(int N)
    {
        // The three involutions plus identity form Klein V₄; the closure relation
        // (any two determine the third) reads D · Q_zx · H = I.
        var D = Pi2KleinV4DephaseSwapGroup.BuildD(N);
        var Q = Pi2KleinV4DephaseSwapGroup.BuildQzx(N);
        var H = Pi2KleinV4DephaseSwapGroup.BuildH(N);
        AssertEqualsIdentity(D * Q * H, $"D · Q_zx · H at N={N}");
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void KleinV4_PairwiseCommutativity(int N)
    {
        // Klein V₄ is abelian. All three pairs of the non-trivial generators commute.
        var D = Pi2KleinV4DephaseSwapGroup.BuildD(N);
        var Q = Pi2KleinV4DephaseSwapGroup.BuildQzx(N);
        var H = Pi2KleinV4DephaseSwapGroup.BuildH(N);
        AssertEqualsZero(D * Q - Q * D, $"[D, Q_zx] at N={N}");
        AssertEqualsZero(D * H - H * D, $"[D, H] at N={N}");
        AssertEqualsZero(Q * H - H * Q, $"[Q_zx, H] at N={N}");
    }

    // ------------------------------------------------------------------
    // Dephase-swap conjugation identities at N=2 via PiOperator.BuildFull
    // ------------------------------------------------------------------

    [Fact]
    public void D_ConjugatesPiZToPiY_AtN2()
    {
        // Welle 12 Task 1 closed-form: D · Π_Z · D = Π_Y bit-exact universal N.
        var D = Pi2KleinV4DephaseSwapGroup.BuildD(N: 2);
        var piZ = PiOperator.BuildFull(N: 2, PauliLetter.Z);
        var piY = PiOperator.BuildFull(N: 2, PauliLetter.Y);
        var lhs = D * piZ * D;
        AssertEqualsMatrix(lhs, piY, "D · Π_Z · D vs Π_Y at N=2");
    }

    [Fact]
    public void H_ConjugatesPiYToPiX_AtN2()
    {
        // Welle 12 Task 2: Y↔X swap via the pure basis-permutation H.
        var H = Pi2KleinV4DephaseSwapGroup.BuildH(N: 2);
        var piY = PiOperator.BuildFull(N: 2, PauliLetter.Y);
        var piX = PiOperator.BuildFull(N: 2, PauliLetter.X);
        var lhs = H * piY * H;
        AssertEqualsMatrix(lhs, piX, "H · Π_Y · H vs Π_X at N=2");
    }

    [Fact]
    public void H_ConjugatesPiXToPiY_AtN2()
    {
        // H² = I so H · Π_X · H = Π_Y is the same identity as H · Π_Y · H = Π_X
        // run in the other direction. Verifies the symmetric reading.
        var H = Pi2KleinV4DephaseSwapGroup.BuildH(N: 2);
        var piX = PiOperator.BuildFull(N: 2, PauliLetter.X);
        var piY = PiOperator.BuildFull(N: 2, PauliLetter.Y);
        var lhs = H * piX * H;
        AssertEqualsMatrix(lhs, piY, "H · Π_X · H vs Π_Y at N=2");
    }

    [Fact]
    public void Qzx_ConjugatesPiZToPiX_AtN2()
    {
        // Welle 12 Task 2: Z↔X swap via the composite Q_zx = H · D.
        var Q = Pi2KleinV4DephaseSwapGroup.BuildQzx(N: 2);
        var piZ = PiOperator.BuildFull(N: 2, PauliLetter.Z);
        var piX = PiOperator.BuildFull(N: 2, PauliLetter.X);
        var lhs = Q * piZ * Q;
        AssertEqualsMatrix(lhs, piX, "Q_zx · Π_Z · Q_zx vs Π_X at N=2");
    }

    // ------------------------------------------------------------------
    // Reconnaissance: print the 4×4 single-site matrices
    // ------------------------------------------------------------------

    [Fact]
    public void Reconnaissance_PrintsPerSiteMatrices()
    {
        var D = Pi2KleinV4DephaseSwapGroup.BuildD(N: 1);
        var H = Pi2KleinV4DephaseSwapGroup.BuildH(N: 1);
        var Q = Pi2KleinV4DephaseSwapGroup.BuildQzx(N: 1);
        _out.WriteLine("Per-site d_l (basis (I, X, Z, Y)):");
        PrintMatrix(D);
        _out.WriteLine("Per-site h:");
        PrintMatrix(H);
        _out.WriteLine("Per-site q_zx = h · d_l:");
        PrintMatrix(Q);
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    private static void AssertEqualsIdentity(ComplexMatrix M, string label)
    {
        int n = M.RowCount;
        var I = Matrix<Complex>.Build.SparseIdentity(n);
        var diff = (M - I).FrobeniusNorm();
        Assert.True(diff < Tol, $"{label}: ‖M − I‖_F = {diff} (tolerance {Tol})");
    }

    private static void AssertEqualsZero(ComplexMatrix M, string label)
    {
        var diff = M.FrobeniusNorm();
        Assert.True(diff < Tol, $"{label}: ‖M‖_F = {diff} (tolerance {Tol})");
    }

    private static void AssertEqualsMatrix(ComplexMatrix A, ComplexMatrix B, string label)
    {
        var diff = (A - B).FrobeniusNorm();
        Assert.True(diff < Tol, $"{label}: ‖A − B‖_F = {diff} (tolerance {Tol})");
    }

    private static void AssertDiagonal(ComplexMatrix M)
    {
        for (int r = 0; r < M.RowCount; r++)
            for (int c = 0; c < M.ColumnCount; c++)
                if (r != c) Assert.Equal(Complex.Zero, M[r, c]);
    }

    private void PrintMatrix(ComplexMatrix M)
    {
        for (int r = 0; r < M.RowCount; r++)
        {
            var row = new string[M.ColumnCount];
            for (int c = 0; c < M.ColumnCount; c++)
                row[c] = $"{M[r, c].Real:+0;-0; 0}{M[r, c].Imaginary:+0;-0; 0}i";
            _out.WriteLine("  [ " + string.Join("  ", row) + " ]");
        }
    }
}
