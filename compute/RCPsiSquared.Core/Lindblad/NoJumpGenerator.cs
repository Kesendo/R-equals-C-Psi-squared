using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Lindblad;

/// <summary>The physical (no-jump) generator G_H : ρ ↦ −i(Hρ − ρH†) for a non-Hermitian H,
/// as a 4^N × 4^N superoperator in the same vec convention the rest of the pipeline uses.
///
/// <para>This is the un-normalised conditional generator of post-selected dynamics: the operator
/// PT and gain-loss models are usually written with, obtained from a Lindbladian by dropping the
/// recycling term Σ_k c_k ρ c_k†. Writing H = A + iB with A and B Hermitian (every matrix splits
/// that way, uniquely) gives</para>
///
/// <code>
///   G_H = −i·(A ⊗ I − I ⊗ Aᵀ) + (B ⊗ I + I ⊗ Bᵀ)
/// </code>
///
/// <para>the first bracket the ordinary commutator −i[A, ·] and the second the transpose-paired
/// ANTI-commutator lift that carries the whole polarity break (F155). The Kronecker order matches
/// <see cref="PauliDephasingDissipator.Build"/> exactly, so a generator built here is in the same
/// stacking as every other superoperator in Core and may be handed straight to
/// <see cref="Symmetry.PalindromeResidual.Build"/>.</para>
///
/// <para><b>What it is not.</b> G_H is not trace-preserving, and the state-dependent
/// normalisation post-selection applies is nonlinear and is no part of this object. It is also
/// not the full Lindbladian: a physical gain medium pumps incoherently, and the recycling term
/// changes the asymmetry outside the σ⁻/σ⁺ family. Conflating G_H with the commutator −i[H, ·]
/// is the error <c>docs/CAUGHT_ERRORS.md</c> records in its 2026-06-20 entry and that entry's
/// 2026-08-19 sequel. (A separate 2026-08-19 entry concerns the F155 dephase siblings; that is
/// a different error, not this one.)</para>
///
/// <para>Python sibling: the generator built inline by
/// <c>simulations/f155_polarity_break_bilinear.py</c>, in the same Kronecker order.</para>
/// </summary>
public static class NoJumpGenerator
{
    /// <summary>Build G_H from the two Hermitian halves of H = A + iB.
    ///
    /// <para>Both halves are required to be EXACTLY Hermitian, compared entry by entry with no
    /// tolerance: an operator assembled from Pauli strings with real coefficients has entries in
    /// {0, ±1, ±i} times those coefficients and is exactly Hermitian in floats, so a failure here
    /// is a caller who has not split H rather than a rounding artefact.</para></summary>
    /// <param name="a">A = (H + H†)/2, the Hermitian part.</param>
    /// <param name="b">B = −i·(H − H†)/2, the anti-Hermitian part divided by i; also Hermitian.</param>
    public static ComplexMatrix Build(ComplexMatrix a, ComplexMatrix b)
    {
        if (a is null) throw new ArgumentNullException(nameof(a));
        if (b is null) throw new ArgumentNullException(nameof(b));
        if (a.RowCount != a.ColumnCount)
            throw new ArgumentException($"A is {a.RowCount}×{a.ColumnCount}, not square", nameof(a));
        if (b.RowCount != b.ColumnCount)
            throw new ArgumentException($"B is {b.RowCount}×{b.ColumnCount}, not square", nameof(b));
        if (a.RowCount != b.RowCount)
            throw new ArgumentException(
                $"A is {a.RowCount}×{a.RowCount} and B is {b.RowCount}×{b.RowCount}; "
              + "they are the two halves of one H and must match", nameof(b));

        int d = a.RowCount;
        int n = (int)Math.Round(Math.Log2(d));
        if ((1 << n) != d)
            throw new ArgumentException($"dimension {d} is not a power of 2", nameof(a));

        RequireExactlyHermitian(a, nameof(a), "A");
        RequireExactlyHermitian(b, nameof(b), "B");

        var identity = Matrix<Complex>.Build.DenseIdentity(d);
        var commutatorPart =
            -Complex.ImaginaryOne * (a.KroneckerProduct(identity) - identity.KroneckerProduct(a.Transpose()));
        var anticommutatorPart =
            b.KroneckerProduct(identity) + identity.KroneckerProduct(b.Transpose());
        return commutatorPart + anticommutatorPart;
    }

    private static void RequireExactlyHermitian(ComplexMatrix m, string parameterName, string label)
    {
        for (int i = 0; i < m.RowCount; i++)
        {
            for (int j = 0; j < m.ColumnCount; j++)
            {
                Complex upper = m[i, j];
                Complex lower = Complex.Conjugate(m[j, i]);
                if (upper != lower)
                    throw new ArgumentException(
                        $"{label} is not exactly Hermitian: entry [{i},{j}] = {upper} against the "
                      + $"conjugate of [{j},{i}] = {lower}. H = A + iB requires both halves "
                      + "Hermitian; split H with A = (H + H†)/2 and B = −i·(H − H†)/2.",
                        parameterName);
            }
        }
    }
}
