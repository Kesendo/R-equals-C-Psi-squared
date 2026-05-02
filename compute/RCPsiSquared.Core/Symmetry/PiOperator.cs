using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Π conjugation operator on the 4^N Pauli-string basis.
///
/// Π is the palindrome operator for a single-letter dephasing dissipator. Per
/// dephase-letter ∈ {X, Y, Z} it flips the Klein axis that swaps the dissipator's
/// per-letter immune sector with its damped sector:
///
///   Z-dephasing (immune {I, Z}, damped {X, Y}): flip bit_a, phase  i^bit_b
///        I ↔ X (sign  1),  Y ↔ Z (sign  i)
///   X-dephasing (immune {I, X}, damped {Y, Z}): flip bit_b, phase −i^bit_a
///        I ↔ Z (sign  1),  X ↔ Y (sign −i)
///   Y-dephasing (immune {I, Y}, damped {X, Z}): flip bit_a, phase −i^bit_b
///        I ↔ X (sign  1),  Y ↔ Z (sign −i)
///
/// Π is unitary order-4. Per the existing class-AIII-chiral classification
/// (PT_SYMMETRY_ANALYSIS), Π is *linear* — distinct from Bender-Boettcher PT (which
/// requires anti-linear operators).
/// </summary>
public static class PiOperator
{
    /// <summary>Π action on a single Pauli letter for a given dephasing letter.
    /// Returns the new letter and the accompanying phase.</summary>
    public static (PauliLetter NewLetter, Complex Phase) ActOnLetter(PauliLetter letter, PauliLetter dephaseLetter = PauliLetter.Z)
    {
        int a = letter.BitA();
        int b = letter.BitB();
        return dephaseLetter switch
        {
            PauliLetter.Z => (LetterFromBits(1 - a, b), b == 1 ? new Complex(0, 1) : Complex.One),
            PauliLetter.X => (LetterFromBits(a, 1 - b), a == 1 ? new Complex(0, -1) : Complex.One),
            PauliLetter.Y => (LetterFromBits(1 - a, b), b == 1 ? new Complex(0, -1) : Complex.One),
            _ => throw new ArgumentException($"dephase_letter must be X, Y, or Z; got {dephaseLetter}"),
        };
    }

    /// <summary>Π² eigenvalue on a Pauli string. For Z- and Y-dephasing this is (−1)^Π²-parity
    /// where parity counts bit_b mod 2; for X-dephasing it counts bit_a mod 2.</summary>
    public static int SquaredEigenvalue(IReadOnlyList<PauliLetter> letters, PauliLetter dephaseLetter = PauliLetter.Z)
    {
        if (dephaseLetter == PauliLetter.I)
            throw new ArgumentException($"dephase_letter must be X, Y, or Z; got {dephaseLetter}");
        int sum = 0;
        foreach (var L in letters)
            sum += dephaseLetter == PauliLetter.X ? L.BitA() : L.BitB();
        return (sum & 1) == 0 ? +1 : -1;
    }

    /// <summary>Π in the 4^N Pauli-string basis as a 4^N × 4^N signed permutation matrix.
    /// Cached by (N, dephaseLetter). Callers should treat the returned matrix as read-only.</summary>
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

    private static PauliLetter LetterFromBits(int a, int b) => (PauliLetter)(a + 2 * b);
}
