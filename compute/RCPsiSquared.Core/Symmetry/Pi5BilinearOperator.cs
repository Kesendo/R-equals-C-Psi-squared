using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Π_5bilinear: a phase-variant of the canonical Π operator that achieves
/// EXACT operator-level palindrome (Π·L·Π⁻¹ = −L − 2σ·I) for every Hamiltonian built
/// from the five Π²_D-even 2-site bilinears under D-dephasing, where D ∈ {Z, X}.
///
/// <para><b>Per-site label maps</b> (4×4 superop on Pauli basis {I, X, Y, Z}):</para>
/// <list type="bullet">
///   <item>Z-deph (F108 Part 1): I → +X, X → −I, Y → +iZ, Z → −iY. Same I↔X, Y↔Z
///         permutation as the canonical P1 Π with sign flips on the X→I and Z→Y back
///         arrows. M² = diag(−1, −1, +1, +1).</item>
///   <item>X-deph (F108 Part 2): I → +Z, Z → −I, X → −iY, Y → +iX. Same I↔Z, X↔Y
///         permutation as the canonical X-deph Π with sign flips on the Z→I and Y→X
///         back arrows. M² = diag(−1, +1, +1, −1).</item>
/// </list>
///
/// <para>M is unitary and order-4 (M⁴ = I). Π_5bilinear is a Liouville-space
/// automorphism on the operator (Pauli) algebra, NOT a Hilbert-space conjugation
/// (no 2×2 unitary U satisfies U·I·U† = X).</para>
///
/// <para><b>Scope</b>: Z- and X-dephasing only. The Y-dephasing analog is the open
/// F108 Part 3 (no covering Claim); calling with <c>PauliLetter.Y</c> throws
/// <see cref="NotImplementedException"/>. Full algebraic proofs and the mechanism
/// derivation live in <c>docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md</c>
/// and <c>docs/proofs/PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md</c>.</para></summary>
public static class Pi5BilinearOperator
{
    /// <summary>Π_5bilinear action on a single Pauli letter for the given dephase
    /// letter D ∈ {Z, X}. Throws <see cref="NotImplementedException"/> for D = Y
    /// (F108 Part 3, no covering Claim yet).</summary>
    public static (PauliLetter NewLetter, Complex Phase) ActOnLetter(
        PauliLetter letter, PauliLetter dephaseLetter = PauliLetter.Z)
    {
        return dephaseLetter switch
        {
            PauliLetter.Z => letter switch
            {
                PauliLetter.I => (PauliLetter.X, Complex.One),
                PauliLetter.X => (PauliLetter.I, -Complex.One),
                PauliLetter.Y => (PauliLetter.Z, new Complex(0, 1)),
                PauliLetter.Z => (PauliLetter.Y, new Complex(0, -1)),
                _ => throw new ArgumentException($"Unknown PauliLetter {letter}", nameof(letter)),
            },
            PauliLetter.X => letter switch
            {
                PauliLetter.I => (PauliLetter.Z, Complex.One),
                PauliLetter.Z => (PauliLetter.I, -Complex.One),
                PauliLetter.X => (PauliLetter.Y, new Complex(0, -1)),
                PauliLetter.Y => (PauliLetter.X, new Complex(0, 1)),
                _ => throw new ArgumentException($"Unknown PauliLetter {letter}", nameof(letter)),
            },
            PauliLetter.Y => throw new NotImplementedException(
                "Y-dephasing variant of Π_5bilinear is not yet typed (F108 Part 3, " +
                "no covering Claim). Use PauliLetter.Z (F108 Part 1) or PauliLetter.X " +
                "(F108 Part 2)."),
            _ => throw new ArgumentException(
                $"dephaseLetter must be X, Y, or Z; got {dephaseLetter}", nameof(dephaseLetter)),
        };
    }

    /// <summary>Π_5bilinear as a 4^N × 4^N signed permutation matrix in the Pauli-string
    /// basis. Cached by (N, dephaseLetter); treat as read-only.</summary>
    public static ComplexMatrix BuildFull(int N, PauliLetter dephaseLetter = PauliLetter.Z) =>
        _cache.GetOrAdd((N, dephaseLetter), key => BuildFullUncached(key.Item1, key.Item2));

    private static readonly System.Collections.Concurrent.ConcurrentDictionary<(int, PauliLetter), ComplexMatrix> _cache = new();

    private static ComplexMatrix BuildFullUncached(int N, PauliLetter dephaseLetter)
    {
        long d2 = 1L << (2 * N);
        var pi = Matrix<Complex>.Build.Sparse((int)d2, (int)d2);
        for (long k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            var newLetters = new PauliLetter[N];
            Complex sign = Complex.One;
            for (int i = 0; i < N; i++)
            {
                var (newL, phase) = ActOnLetter(letters[i], dephaseLetter);
                newLetters[i] = newL;
                sign *= phase;
            }
            long newK = PauliIndex.ToFlat(newLetters);
            pi[(int)newK, (int)k] = sign;
        }
        return pi;
    }
}
