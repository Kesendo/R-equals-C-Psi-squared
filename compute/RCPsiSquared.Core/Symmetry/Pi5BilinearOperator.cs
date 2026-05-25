using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Π_5bilinear: a phase-variant of the canonical P1 Π operator that achieves
/// EXACT operator-level palindrome (Π·L·Π⁻¹ = −L − 2σ·I) for every Hamiltonian built
/// from the five Π²_D-even 2-site bilinears under D-dephasing, where D ∈ {Z, X}.
///
/// <para><b>Per-site label action under Z-dephasing</b> (4×4 superoperator on Pauli
/// basis {I, X, Y, Z}):</para>
/// <list type="bullet">
///   <item>I → +1 · X</item>
///   <item>X → −1 · I    (sign flip vs. canonical P1 Π's X → +I)</item>
///   <item>Y → +i · Z</item>
///   <item>Z → −i · Y    (sign flip vs. canonical P1 Π's Z → +iY)</item>
/// </list>
///
/// <para><b>Per-site label action under X-dephasing</b> (the BitA-axis twin, F108
/// Part 2; same structural pattern with I↔Z, X↔Y permutation):</para>
/// <list type="bullet">
///   <item>I → +1 · Z</item>
///   <item>Z → −1 · I    (sign flip vs. canonical X-deph Π's Z → +I)</item>
///   <item>X → −i · Y</item>
///   <item>Y → +i · X    (sign flip vs. canonical X-deph Π's Y → −iX)</item>
/// </list>
///
/// <para>For Z-dephasing: same I↔X, Y↔Z permutation as <see cref="PiOperator"/> with
/// two phase flips on the back arrows. M² = diag(−1, −1, +1, +1) on {I, X, Y, Z};
/// M⁴ = I; M is order-4 and unitary. For X-dephasing: same I↔Z, X↔Y permutation as
/// the canonical X-deph Π, with the analogous two phase flips on back arrows. M² =
/// diag(−1, +1, +1, −1); same order-4 structure.</para>
///
/// <para><b>Mechanism</b>: M anti-commutes with the commutator superop [B, ·] for
/// every Π²_D-even 2-body bilinear B (5 bilinears per dephase letter), verified
/// bit-exact at the 2-qubit level. Combined with the per-site dissipator identity
/// M·D[P]·M⁻¹ = −D[P] − 2γ·I (P being the dephase letter), this closes F108 Part 1
/// (Z-dephasing) and F108 Part 2 (X-dephasing): every Π²-even H + matching
/// dephasing has a guaranteed EXACT operator-level palindrome via Π_5bilinear.</para>
///
/// <para><b>Subtlety</b>: Π_5bilinear is a Liouville-space automorphism, NOT a
/// Hilbert-space conjugation (ρ → U·ρ·U†). No 2×2 unitary U satisfies U·I·U† = X
/// since U·I·U† = I always. Π_5bilinear acts on the operator (Pauli) algebra
/// directly, permuting Pauli labels with phases.</para>
///
/// <para><b>Scope</b>: Z- and X-dephasing only. The Y-dephasing analog is tracked
/// as an open structural-derivation slot (F108 Part 3, no covering Claim yet);
/// calling with <c>PauliLetter.Y</c> throws <see cref="NotImplementedException"/>.</para>
///
/// <para>See <c>docs/proofs/PROOF_F108_PART1_PI2_EVEN_ALWAYS_PALINDROMIC.md</c> and
/// <c>docs/proofs/PROOF_F108_PART2_PI2X_EVEN_ALWAYS_PALINDROMIC.md</c> for the full
/// algebraic proofs; numerical scans in
/// <c>simulations/_f108_part1_pi_family_scan.py</c> and
/// <c>simulations/_f108_part2_x_dephasing_scan.py</c>.</para></summary>
public static class Pi5BilinearOperator
{
    /// <summary>Π_5bilinear action on a single Pauli letter for a given dephase
    /// letter. Returns the new letter and the accompanying phase per the
    /// 5-bilinear phase-variant rule (see class docstring). Supports D ∈ {Z, X};
    /// throws <see cref="NotImplementedException"/> for D = Y.</summary>
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

    /// <summary>Π_5bilinear in the 4^N Pauli-string basis as a 4^N × 4^N signed
    /// permutation matrix. Cached by (N, dephaseLetter). Callers should treat the
    /// returned matrix as read-only.</summary>
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
