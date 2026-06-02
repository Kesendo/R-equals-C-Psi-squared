using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Continuous operator-space rotation Ad_{R_z(θ)} about the dephasing axis Z,
/// where R_z(θ) = exp(−iθZ/2). It is the continuous generalization of <see cref="PiOperator"/>'s
/// discrete letter action: where Π is an order-4 signed permutation of Pauli letters, this is a
/// one-parameter dial that turns the X-Y plane by angle θ while fixing I and Z.
///
/// On a single qubit the adjoint action R_z(θ) σ R_z(θ)† rotates X ↦ cos θ·X + sin θ·Y and
/// Y ↦ cos θ·Y − sin θ·X (I and Z fixed). Two angles are framework anchors:
///
///   θ = π/2: the σ_x↔σ_y 90° rotation X ↦ Y, Y ↦ −X, the operator-space NinetyDegreeMirror
///            (the adjoint action of the phase gate S).
///   θ = π/4: its square root, the adjoint action of the T gate. Ad_{R_z(π/4)}² = Ad_{R_z(π/2)},
///            the √-of-90° anchor proven in docs/proofs/PROOF_CROSSOVER_MIRROR_SQRT_NINETY.md.
///
/// The superoperator is built in the vec basis R⊗conj(R) to match the Liouvillian builders'
/// row-major vec convention (see <see cref="Lindblad.PauliDephasingDissipator"/>), so that
/// (R⊗conj(R))·vec(O) = vec(R O R†).
/// </summary>
public static class RotationSuperoperator
{
    /// <summary>R_z(θ) = exp(−iθZ/2) acting on the listed sites (identity elsewhere), as a
    /// 2^N × 2^N diagonal unitary. Per computational basis index b, the diagonal entry is the
    /// product over listed sites l of (e^{−iθ/2} if bit l is 0 else e^{+iθ/2}). Site 0 is the
    /// most-significant qubit (big-endian); bit l is extracted as (b >> (N-1-l)) &amp; 1.</summary>
    public static ComplexMatrix RzHilbert(double theta, int N, IReadOnlyList<int> sites)
    {
        if (N < 1) throw new ArgumentException($"N must be at least 1; got {N}", nameof(N));
        foreach (int l in sites)
            if (l < 0 || l >= N)
                throw new ArgumentException($"site {l} out of range [0, {N - 1}]", nameof(sites));

        int d = 1 << N;
        var minus = Complex.Exp(new Complex(0, -theta / 2.0)); // bit 0: eigenvalue +1 of Z
        var plus = Complex.Exp(new Complex(0, +theta / 2.0));   // bit 1: eigenvalue −1 of Z

        var R = Matrix<Complex>.Build.Sparse(d, d);
        for (int b = 0; b < d; b++)
        {
            Complex phase = Complex.One;
            foreach (int l in sites)
            {
                int bit = (b >> (N - 1 - l)) & 1;
                phase *= bit == 0 ? minus : plus;
            }
            R[b, b] = phase;
        }
        return R;
    }

    /// <summary>Ad_{R_z(θ)} in the vec basis: R⊗conj(R) with R = <see cref="RzHilbert"/>, a
    /// 4^N × 4^N adjoint-action superoperator acting on vec(O) by (R⊗conj(R))·vec(O) = vec(R O R†).
    /// This is the continuous generalization of <see cref="PiOperator.BuildFull"/>.</summary>
    public static ComplexMatrix AdRzVec(double theta, int N, IReadOnlyList<int> sites)
    {
        var R = RzHilbert(theta, N, sites);
        return R.KroneckerProduct(R.Conjugate());
    }
}
