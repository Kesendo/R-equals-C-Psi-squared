using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Polarity;

namespace RCPsiSquared.Diagnostics.Tests.Polarity;

public class F112NonHermitianBasisEnumerationTests
{
    [Fact]
    public void BuildLHInPauliBasis_PauliStringXI_AtN2_HasExpectedShape()
    {
        // XI is the Pauli string "X ⊗ I" at N=2; its L_H = -i[XI, ·] in Pauli basis
        // must be a 16x16 complex matrix (4^N x 4^N at N=2).
        var letters = new[] { PauliLetter.X, PauliLetter.I };
        var H = PauliString.Build(letters);

        var L = F112NonHermitianBasisEnumeration.BuildLHInPauliBasis(H, N: 2);

        Assert.Equal(16, L.RowCount);
        Assert.Equal(16, L.ColumnCount);

        // Sanity: L = -i[XI, ·] for non-trivial XI must be non-zero. Catches sign flips
        // and accidental zero matrices that shape-only assertions would miss.
        Assert.True(L.FrobeniusNorm() > 0.1,
            $"L_XI should be non-zero at N=2; got Frobenius norm = {L.FrobeniusNorm()}");
    }

    [Fact]
    public void ProjectOntoPiEigenspace_OnL_YI_AtN2_GivesNonzeroMinusIProjection()
    {
        // YI = Y⊗I at N=2 is a Hermitian Pauli string with bit_b parity 1 (Y has bit_b=1);
        // its L_H,-i Π-projection lives in the -i eigenspace (non-trivial sanity check
        // that the projection machinery actually works). Cross-checked against the Python
        // reference simulations/_f112_open_identity_basis_enum.py: ||proj L_YI,-i|| = 4.0.
        // (Strings with bit_b parity 0, e.g. XI, project entirely to the +1 eigenspace
        // and yield zero here, which would fail this sanity test.)
        var letters = new[] { PauliLetter.Y, PauliLetter.I };
        var H = PauliString.Build(letters);
        var L = F112NonHermitianBasisEnumeration.BuildLHInPauliBasis(H, N: 2);
        var pi = PiOperator.BuildFull(N: 2, PauliLetter.Z);

        var lMinusI = F112NonHermitianBasisEnumeration.ProjectOntoPiEigenspace(L, pi, -Complex.ImaginaryOne);

        Assert.Equal(16, lMinusI.RowCount);
        Assert.Equal(16, lMinusI.ColumnCount);
        Assert.True(lMinusI.FrobeniusNorm() > 1e-10,
            $"L_YI,-i should be non-zero at N=2; got Frobenius norm = {lMinusI.FrobeniusNorm()}");
    }

    [Fact]
    public void FrobeniusInner_OfIdentityWithIdentity_EqualsDimension()
    {
        var I16 = Matrix<Complex>.Build.DenseIdentity(16);
        var inner = F112NonHermitianBasisEnumeration.FrobeniusInner(I16, I16);
        // ⟨I, I⟩ = Σ |I[i,j]|² = trace = 16
        Assert.Equal(16.0, inner.Real, precision: 10);
        Assert.Equal(0.0, inner.Imaginary, precision: 10);
    }

    [Fact]
    public void FrobeniusInner_OfPureImaginaryWithSelf_IsPositive()
    {
        // ⟨[[i, 0], [0, 0]], [[i, 0], [0, 0]]⟩ = conj(i) · i = (-i) · i = +1.
        // Without the Conjugate call this would return -1; this test catches a
        // missing conjugation that the all-real identity test cannot.
        var M = Matrix<Complex>.Build.Dense(2, 2);
        M[0, 0] = Complex.ImaginaryOne;
        var inner = F112NonHermitianBasisEnumeration.FrobeniusInner(M, M);
        Assert.Equal(1.0, inner.Real, precision: 10);
        Assert.Equal(0.0, inner.Imaginary, precision: 10);
    }

    [Fact]
    public void Enumerate_AtN2_All136UnorderedPairsAreBitExactZero()
    {
        // The Python anchor: at N=2, all 16 × 16 = 256 ordered pairs (136 distinct
        // upper-triangular) give Im⟨L_α,-i, L_β,-i⟩ = 0 bit-exact within tolerance.
        var result = F112NonHermitianBasisEnumeration.Enumerate(N: 2, tolerance: 1e-10);

        Assert.Equal(2, result.N);
        Assert.Equal(136, result.TotalPairs);  // 16 * 17 / 2 = 136 upper-triangular
        Assert.Equal(0, result.NonzeroCount);
        Assert.True(result.MaxImaginary < 1e-10,
            $"max |Im| at N=2 should be < 1e-10; got {result.MaxImaginary:E4}");
        Assert.Empty(result.NonzeroExamples);
        Assert.True(result.Elapsed > TimeSpan.Zero, "Stopwatch should have ticked");
    }

    [Theory]
    [InlineData(3, 2080)]   // 64 strings, 64 * 65 / 2 = 2080 unordered pairs
    [InlineData(4, 32896)]  // 256 strings, 256 * 257 / 2 = 32896 unordered pairs
    public void Enumerate_MatchesPythonAnchorAtNFromTwoToFour(int N, int expectedTotalPairs)
    {
        // Python anchor (simulations/_f112_open_identity_basis_enum.py output at N=2,3,4):
        // all pairs give Im = 0 bit-exact. This Theory replicates that anchor in C#.
        var result = F112NonHermitianBasisEnumeration.Enumerate(N, tolerance: 1e-10);

        Assert.Equal(N, result.N);
        Assert.Equal(expectedTotalPairs, result.TotalPairs);
        Assert.Equal(0, result.NonzeroCount);
        Assert.True(result.MaxImaginary < 1e-10,
            $"max |Im| at N={N} should be < 1e-10; got {result.MaxImaginary:E4}");
    }

    /// <summary>SLOW (~15-30 min, ~16 GB RAM): the N=5 enumeration that served as the
    /// empirical anchor (Welle 10b, 2026-05-26) for the universal-N closure subsequently
    /// established structurally (Welle 11, 2026-05-27; see
    /// <c>docs/proofs/PROOF_F112_NONHERMITIAN_UNIVERSAL_N.md</c>). After Welle 11 the
    /// structural proof gives Tier1Derived for all N, so this enumeration is the
    /// empirical anchor (historical numerical validation) rather than the
    /// Tier-promotion mechanism. Tagged SLOW_F112 so CI can filter
    /// via <c>--filter "Category!=SLOW_F112"</c> when needed.
    ///
    /// <para>Run explicitly: <c>dotnet test compute\RCPsiSquared.Diagnostics.Tests -c
    /// Release --filter "Category=SLOW_F112"</c>.</para></summary>
    [Fact]
    [Trait("Category", "SLOW_F112")]
    public void Enumerate_AtN5_All524800UpperTriangularPairsAreBitExactZero()
    {
        // 4^5 = 1024 Pauli strings → 1024 * 1025 / 2 = 524800 upper-triangular pairs.
        // ~16 GB working memory (1024 cached L_α,-i matrices × 16 MB each).
        const int N = 5;
        const int expectedPairs = 524_800;

        var result = F112NonHermitianBasisEnumeration.Enumerate(N, tolerance: 1e-10);

        Assert.Equal(N, result.N);
        Assert.Equal(expectedPairs, result.TotalPairs);
        Assert.True(result.NonzeroCount == 0,
            $"F112 non-Hermitian N=5 enumeration found {result.NonzeroCount} non-zero pairs " +
            $"(max |Im| = {result.MaxImaginary:E4}); first examples: " +
            string.Join("; ", result.NonzeroExamples.Take(5).Select(e => $"F({e.Alpha},{e.Beta})={e.Imag:E4}")));
        Assert.True(result.MaxImaginary < 1e-10,
            $"max |Im| at N=5 should be < 1e-10; got {result.MaxImaginary:E4}");
    }
}
