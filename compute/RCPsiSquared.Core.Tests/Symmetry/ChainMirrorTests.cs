using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.States;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class ChainMirrorTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void R_Squared_IsIdentity(int N)
    {
        var R = ChainMirror.Build(N);
        int d = 1 << N;
        var rr = R * R;
        var i = Matrix<Complex>.Build.DenseIdentity(d);
        Assert.True((rr - i).FrobeniusNorm() < 1e-12);
    }

    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    public void R_IsRealSymmetric(int N)
    {
        var R = ChainMirror.Build(N);
        Assert.True((R - R.Transpose()).FrobeniusNorm() < 1e-12);
        for (int i = 0; i < R.RowCount; i++)
            for (int j = 0; j < R.ColumnCount; j++)
                Assert.Equal(0.0, R[i, j].Imaginary, 12);
    }

    [Fact]
    public void Symmetric_PolarityState_IsF71_PlusEigenstate()
    {
        // |+⟩^N is invariant under bit reversal (uniform pattern is symmetric).
        var psi = PolarityState.Uniform(N: 4, sign: +1);
        Assert.Equal(+1, ChainMirror.EigenstateClass(psi));
    }

    [Fact]
    public void Antisymmetric_State_IsF71_MinusEigenstate()
    {
        // Construct ψ - R·ψ for a non-symmetric state — antisymmetric under R.
        int N = 3;
        var R = ChainMirror.Build(N);
        var seed = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(1 << N);
        seed[1] = Complex.One;  // |001⟩, mirrors to |100⟩ = index 4
        var psi = seed - R * seed;
        psi /= (Complex)Math.Sqrt(psi.ConjugateDotProduct(psi).Real);
        Assert.Equal(-1, ChainMirror.EigenstateClass(psi));
    }

    [Theory]
    [InlineData(3, 1, 1)]   // odd N: balanced k+k split for N-1=2 bonds → 1 sym + 1 asym
    [InlineData(5, 2, 2)]   // odd N=5: 4 bonds → 2 sym + 2 asym
    [InlineData(4, 2, 1)]   // even N=4: 3 bonds → 2 sym + 1 asym (one self-mirror)
    [InlineData(6, 3, 2)]   // even N=6: 5 bonds → 3 sym + 2 asym (one self-mirror)
    public void BondMirrorBasis_HasExpectedDimensions(int N, int nSym, int nAsym)
    {
        var (sym, asym) = ChainMirror.BondMirrorBasis(N);
        Assert.Equal(nSym, sym.Length);
        Assert.Equal(nAsym, asym.Length);
    }
}
