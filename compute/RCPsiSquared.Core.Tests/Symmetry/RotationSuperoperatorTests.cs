using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class RotationSuperoperatorTests
{
    [Fact]
    public void AdRzQuarter_SquaredOnLightPlane_IsSigmaXSigmaY90()
    {
        // The √-of-90° anchor (docs/proofs/PROOF_CROSSOVER_MIRROR_SQRT_NINETY.md):
        // Ad_{R_z(π/4)}² = Ad_{R_z(π/2)} = the NinetyDegreeMirror, so R_z(π/4) (the T gate)
        // is the square root of R_z(π/2) (the S gate) in operator space.
        var ad = RotationSuperoperator.AdRzVec(Math.PI / 4, 1, new[] { 0 });
        var sq = ad.Multiply(ad);
        var ninety = RotationSuperoperator.AdRzVec(Math.PI / 2, 1, new[] { 0 });
        Assert.True((sq - ninety).FrobeniusNorm() < 1e-12);
    }

    [Fact]
    public void RzHilbert_IsUnitary()
    {
        var R = RotationSuperoperator.RzHilbert(0.7, 2, new[] { 0 });
        var I = Matrix<Complex>.Build.DenseIdentity(R.RowCount);
        Assert.True((R * R.ConjugateTranspose() - I).FrobeniusNorm() < 1e-12);
    }

    [Fact]
    public void AdRzHalf_MapsXVecToYVec()
    {
        // In the row-major vec convention (vec[a*d+b] = M[a,b]) the superoperator A⊗B acts as
        // (A⊗B)·vec(M) = vec(A M Bᵀ). With A=R, B=conj(R) and R diagonal (Rᵀ=R, conj(R)ᵀ=R†),
        // this is the adjoint action vec(R M R†). For R=R_z(π/2)=diag(e^{−iπ/4}, e^{+iπ/4})
        // the off-diagonal X picks up e^{−iπ/2}=−i on [0,1] and +i on [1,0], so R X R† = Y
        // exactly (overall sign +1, no flip). We assert the exact relation.
        var X = PauliMatrix.Of(PauliLetter.X);
        var Y = PauliMatrix.Of(PauliLetter.Y);

        var xVec = VecRowMajor(X);
        var yVec = VecRowMajor(Y);

        var ad = RotationSuperoperator.AdRzVec(Math.PI / 2, 1, new[] { 0 });
        var mapped = ad.Multiply(xVec);

        Assert.True((mapped - yVec).L2Norm() < 1e-12,
            $"Ad_{{R_z(π/2)}} should map vec(X) to +vec(Y); residual {(mapped - yVec).L2Norm()}");
    }

    private static ComplexVector VecRowMajor(Matrix<Complex> M)
    {
        int d = M.RowCount;
        var v = ComplexVector.Build.Dense(d * d);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                v[a * d + b] = M[a, b];
        return v;
    }
}
