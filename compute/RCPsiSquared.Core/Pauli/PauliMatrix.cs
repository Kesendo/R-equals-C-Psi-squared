using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Pauli;

/// <summary>The four 2×2 Pauli matrices. The standard physical convention:
///
///   I = [[1, 0], [0, 1]]
///   X = [[0, 1], [1, 0]]
///   Z = [[1, 0], [0, −1]]
///   Y = [[0, −i], [i, 0]]
///
/// Algebraic identities: σ_α² = I; σ_x σ_y = i σ_z (and cyclic);
/// {σ_α, σ_β} = 2 δ_αβ I; [σ_α, σ_β] = 2i ε_{αβγ} σ_γ.
/// </summary>
public static class PauliMatrix
{
    private static readonly ComplexMatrix _i = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
    {
        { Complex.One, Complex.Zero },
        { Complex.Zero, Complex.One },
    });

    private static readonly ComplexMatrix _x = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
    {
        { Complex.Zero, Complex.One },
        { Complex.One, Complex.Zero },
    });

    private static readonly ComplexMatrix _z = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
    {
        { Complex.One, Complex.Zero },
        { Complex.Zero, -Complex.One },
    });

    private static readonly ComplexMatrix _y = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
    {
        { Complex.Zero, new Complex(0, -1) },
        { new Complex(0, 1), Complex.Zero },
    });

    public static ComplexMatrix Of(PauliLetter letter) => letter switch
    {
        PauliLetter.I => _i.Clone(),
        PauliLetter.X => _x.Clone(),
        PauliLetter.Z => _z.Clone(),
        PauliLetter.Y => _y.Clone(),
        _ => throw new ArgumentOutOfRangeException(nameof(letter)),
    };
}
