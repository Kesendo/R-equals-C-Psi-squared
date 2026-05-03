using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Projectors onto the Π²-eigenspaces of the 4^N Pauli operator space, the
/// involutive-symmetry decomposition that is the most basic structural primitive in the
/// F-chain.
///
/// <para><b>What this is</b>. Π is the F1 Pauli-string conjugation operator (see
/// <see cref="PiOperator"/>). It is unitary and order 4, hence Π² is unitary involutive
/// (Π⁴ = I, Π² = Π⁻²). On any single Pauli string σ_α, Π² acts as the scalar (±1)
/// computed by <see cref="PiOperator.SquaredEigenvalue"/>: <c>+1</c> when the string's
/// bit_b parity (or bit_a, depending on dephase letter) is even, <c>−1</c> when odd.
/// Therefore Π² is diagonal in the Pauli-string basis, and the symmetric / antisymmetric
/// projectors <c>P_±² = (I ± Π²) / 2</c> are diagonal 4^N × 4^N matrices with 0/1 entries.</para>
///
/// <para><b>Why this is the first layer</b>. Decomposing operator space into the
/// Π²-fixed subspace (Π²-even, "self-conjugate under the framework's involution") and its
/// orthogonal complement (Π²-odd, "anti-self-conjugate") is the structural move that
/// underlies the F87 trichotomy, F1's palindrome residual decomposition, F79's 2-body
/// Π²-block structure, and F81's M = M_sym + M_anti split. Each downstream F-formula uses
/// this projection, often implicitly. Exposing it directly lets us read every later
/// construction as a refinement of "self-conjugate + anti-self-conjugate" content.</para>
///
/// <para>The dephase letter controls the Klein-cell rotation per
/// <see cref="PiOperator.SquaredEigenvalue"/>: Z-dephasing counts bit_b parity, Y-dephasing
/// counts bit_b too (different Π but same parity index, since Π_Y and Π_Z both flip bit_a),
/// X-dephasing counts bit_a parity (Π_X flips bit_b instead).</para>
/// </summary>
public static class Pi2Projection
{
    private static readonly System.Collections.Concurrent.ConcurrentDictionary<(int, PauliLetter, bool), ComplexMatrix> _cache = new();

    /// <summary>Diagonal 4^N × 4^N projector <c>P_+² = (I + Π²) / 2</c> onto the Π²-even
    /// (self-conjugate) Pauli-string subspace. Cached by (N, dephaseLetter).</summary>
    public static ComplexMatrix SymmetricProjector(int N, PauliLetter dephaseLetter = PauliLetter.Z) =>
        _cache.GetOrAdd((N, dephaseLetter, true), key => Build(key.Item1, key.Item2, even: true));

    /// <summary>Diagonal 4^N × 4^N projector <c>P_−² = (I − Π²) / 2</c> onto the Π²-odd
    /// (anti-self-conjugate) Pauli-string subspace. Cached by (N, dephaseLetter).</summary>
    public static ComplexMatrix AntisymmetricProjector(int N, PauliLetter dephaseLetter = PauliLetter.Z) =>
        _cache.GetOrAdd((N, dephaseLetter, false), key => Build(key.Item1, key.Item2, even: false));

    /// <summary>Counts of Π²-even and Π²-odd Pauli strings on N sites. The two counts sum
    /// to 4^N and (for any dephase letter) are equal: 4^N / 2 each.</summary>
    public static (long Even, long Odd) Counts(int N, PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 1; got {N}");
        long d2 = 1L << (2 * N);
        long evenCount = 0;
        for (long k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            if (PiOperator.SquaredEigenvalue(letters, dephaseLetter) == +1) evenCount++;
        }
        return (evenCount, d2 - evenCount);
    }

    /// <summary>Decompose a 2^N × 2^N operator A into its Π²-even and Π²-odd Pauli-basis
    /// components: <c>A = A_even + A_odd</c>, where A_even is the resummation over Π²-even
    /// Pauli strings and A_odd over Π²-odd ones. Both pieces are 2^N × 2^N matrices.
    ///
    /// <para>Round-trip identity: <c>A_even + A_odd = A</c> exactly (machine precision).
    /// Frobenius orthogonality: <c>tr(A_even† · A_odd) = 0</c>.</para>
    /// </summary>
    public static (ComplexMatrix Even, ComplexMatrix Odd) Split(
        ComplexMatrix A, int N, PauliLetter dephaseLetter = PauliLetter.Z)
    {
        int d = 1 << N;
        if (A.RowCount != d || A.ColumnCount != d)
            throw new ArgumentException($"expected {d}×{d} matrix for N={N}; got {A.RowCount}×{A.ColumnCount}");

        long d2 = 1L << (2 * N);
        var even = Matrix<Complex>.Build.Dense(d, d);
        var odd = Matrix<Complex>.Build.Dense(d, d);
        double inv = 1.0 / d;

        for (long k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            var sigma = PauliString.Build(letters);
            Complex coeff = (sigma * A).Trace() * inv;
            if (coeff.Magnitude < 1e-15) continue;
            if (PiOperator.SquaredEigenvalue(letters, dephaseLetter) == +1)
                even += coeff * sigma;
            else
                odd += coeff * sigma;
        }
        return (even, odd);
    }

    /// <summary>The 4-Klein-cell decomposition under both Π²-axes simultaneously: each Pauli
    /// string carries a pair of eigenvalues (Π²_Z, Π²_X) ∈ {+1, −1}². The 4 cells correspond
    /// to the four combinations.
    ///
    /// <para>This is the framework's own 2-axis structural cut, distinct from any standard
    /// Klein-Vierergruppen literature: the Z-axis of Π² counts bit_b parity (the F1 palindrome
    /// involution under Z-dephasing), the X-axis counts bit_a parity (the orthogonal involution
    /// under X-dephasing). Y-dephasing collapses onto the Z-axis (same bit_b parity).</para>
    ///
    /// <para>For 2-body bilinears the 4 cells are: (+, +) the 3 truly cases (XX, YY, ZZ),
    /// (−, +) the 2 non-truly Π²-even cases (YZ, ZY), (+, −) the bond-flip pair (XY, YX),
    /// (−, −) the bond-flip pair (XZ, ZX). The two-axis cut distinguishes truly from non-truly
    /// Π²-even, which the Z-axis alone cannot.</para>
    /// </summary>
    public static KleinDecomposition KleinSplit(ComplexMatrix A, int N)
    {
        int d = 1 << N;
        if (A.RowCount != d || A.ColumnCount != d)
            throw new ArgumentException($"expected {d}×{d} matrix for N={N}; got {A.RowCount}×{A.ColumnCount}");

        long d2 = 1L << (2 * N);
        var pp = Matrix<Complex>.Build.Dense(d, d);
        var pm = Matrix<Complex>.Build.Dense(d, d);
        var mp = Matrix<Complex>.Build.Dense(d, d);
        var mm = Matrix<Complex>.Build.Dense(d, d);
        double inv = 1.0 / d;

        for (long k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            var sigma = PauliString.Build(letters);
            Complex coeff = (sigma * A).Trace() * inv;
            if (coeff.Magnitude < 1e-15) continue;

            int eigZ = PiOperator.SquaredEigenvalue(letters, PauliLetter.Z);
            int eigX = PiOperator.SquaredEigenvalue(letters, PauliLetter.X);
            var contribution = coeff * sigma;
            if (eigZ == +1 && eigX == +1) pp += contribution;
            else if (eigZ == -1 && eigX == +1) mp += contribution; // Π²_Z=−, Π²_X=+
            else if (eigZ == +1 && eigX == -1) pm += contribution; // Π²_Z=+, Π²_X=−
            else mm += contribution;
        }
        return new KleinDecomposition(pp, mp, pm, mm);
    }

    private static ComplexMatrix Build(int N, PauliLetter dephaseLetter, bool even)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 1; got {N}");
        long d2 = 1L << (2 * N);
        var P = Matrix<Complex>.Build.Sparse((int)d2, (int)d2);
        for (long k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            int eig = PiOperator.SquaredEigenvalue(letters, dephaseLetter);
            bool include = even ? eig == +1 : eig == -1;
            if (include) P[(int)k, (int)k] = Complex.One;
        }
        return P;
    }
}

/// <summary>The 4-cell Klein decomposition of a 2^N × 2^N operator under (Π²_Z, Π²_X). The
/// four cells are labelled by the eigenvalue pair, with sum Pp + Pm + Mp + Mm = A exact.</summary>
public sealed record KleinDecomposition(
    ComplexMatrix Pp,  // (Π²_Z = +1, Π²_X = +1) — truly bilinears live here
    ComplexMatrix Mp,  // (Π²_Z = −1, Π²_X = +1) — XY, YX bond-flip pair
    ComplexMatrix Pm,  // (Π²_Z = +1, Π²_X = −1) — non-truly Π²-even (YZ, ZY)
    ComplexMatrix Mm); // (Π²_Z = −1, Π²_X = −1) — XZ, ZX bond-flip pair
