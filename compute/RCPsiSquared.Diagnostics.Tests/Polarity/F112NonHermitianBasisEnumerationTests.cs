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
}
