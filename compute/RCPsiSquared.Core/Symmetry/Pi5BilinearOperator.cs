using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Π_5bilinear: a phase-variant of the canonical P1 Π operator that achieves
/// EXACT operator-level palindrome (Π·L·Π⁻¹ = −L − 2σ·I) for every Hamiltonian built
/// from the five Π²_Z-even 2-site bilinears {XX, YY, YZ, ZY, ZZ} under Z-dephasing.
///
/// <para><b>Per-site label action</b> (4×4 superoperator on Pauli basis {I, X, Y, Z}):</para>
/// <list type="bullet">
///   <item>I → +1 · X</item>
///   <item>X → −1 · I    (sign flip vs. canonical P1 Π's X → +I)</item>
///   <item>Y → +i · Z</item>
///   <item>Z → −i · Y    (sign flip vs. canonical P1 Π's Z → +iY)</item>
/// </list>
///
/// <para>Same I↔X, Y↔Z permutation as <see cref="PiOperator"/> (P1 family) but with
/// two phase flips. M² = diag(−1, −1, +1, +1) on {I, X, Y, Z}; M⁴ = I; M is order-4
/// and unitary on the 4-dim per-site Pauli basis.</para>
///
/// <para><b>Mechanism</b>: M anti-commutes with the commutator superop [B, ·] for
/// every Π²_Z-even 2-body bilinear B ∈ {XX, YY, YZ, ZY, ZZ} (verified bit-exact at
/// the 2-qubit level). Combined with the per-site dissipator identity M·D[Z]·M⁻¹ =
/// −D[Z] − 2γ·I, this closes F108 Part 1: every Π²-even H + Z-dephasing has a
/// guaranteed EXACT operator-level palindrome via Π_5bilinear, even when the
/// canonical Heisenberg Π gives M ≠ 0.</para>
///
/// <para><b>Subtlety</b>: Π_5bilinear is a Liouville-space automorphism, NOT a
/// Hilbert-space conjugation (ρ → U·ρ·U†). No 2×2 unitary U satisfies U·I·U† = X
/// since U·I·U† = I always. Π_5bilinear acts on the operator (Pauli) algebra
/// directly, permuting Pauli labels with phases.</para>
///
/// <para><b>Scope</b>: Z-dephasing only. The X- and Y-dephasing analogs require
/// different per-site permutations (P4-family-style I↔Y, X↔Z and similar) and are
/// tracked as the open BitATwin slot of
/// <see cref="F108Part1Pi2EvenAlwaysPalindromic"/>.</para>
///
/// <para>See <c>docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md</c> for
/// the full algebraic proof; numerical scan in
/// <c>simulations/_f108_part1_pi_family_scan.py</c>.</para></summary>
public static class Pi5BilinearOperator
{
    /// <summary>Π_5bilinear action on a single Pauli letter under Z-dephasing.
    /// Returns the new letter and the accompanying phase per the I↔X with phases
    /// (+1, −1) and Y↔Z with phases (+i, −i) rule.</summary>
    public static (PauliLetter NewLetter, Complex Phase) ActOnLetter(PauliLetter letter)
    {
        return letter switch
        {
            PauliLetter.I => (PauliLetter.X, Complex.One),
            PauliLetter.X => (PauliLetter.I, -Complex.One),
            PauliLetter.Y => (PauliLetter.Z, new Complex(0, 1)),
            PauliLetter.Z => (PauliLetter.Y, new Complex(0, -1)),
            _ => throw new ArgumentException($"Unknown PauliLetter {letter}", nameof(letter)),
        };
    }

    /// <summary>Π_5bilinear in the 4^N Pauli-string basis as a 4^N × 4^N signed
    /// permutation matrix. Cached by N. Callers should treat the returned matrix
    /// as read-only.</summary>
    public static ComplexMatrix BuildFull(int N) =>
        _cache.GetOrAdd(N, n => BuildFullUncached(n));

    private static readonly System.Collections.Concurrent.ConcurrentDictionary<int, ComplexMatrix> _cache = new();

    private static ComplexMatrix BuildFullUncached(int N)
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
                var (newL, phase) = ActOnLetter(letters[i]);
                newLetters[i] = newL;
                sign *= phase;
            }
            long newK = PauliIndex.ToFlat(newLetters);
            pi[(int)newK, (int)k] = sign;
        }
        return pi;
    }
}
