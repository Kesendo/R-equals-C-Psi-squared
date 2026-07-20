using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the direct-sum decomposition
/// (<c>docs/proofs/DIRECT_SUM_DECOMPOSITION.md</c>, typed home
/// <see cref="DirectSumDecompositionClaim"/>): L = L_even ⊕ L_odd by n_XY parity,
/// equal halves 2^(2N−1), odd-N Π sector exchange with mirror-image dynamics,
/// superselection charge P_XY = (−1)^{n_XY}.
///
/// <para>The witness BUILDS it at inspect time for a Heisenberg chain + uniform
/// Z-dephasing: transforms L to the 4^N Pauli-string basis, reads the
/// off-parity-block Frobenius norm (statement 1/4: exactly 0, i.e. [P_XY, L] = 0),
/// counts the sector dimensions (2^(2N−1) each), checks Π's sector map
/// column-by-column (odd N: every Pauli string lands in the OPPOSITE parity;
/// even N: the SAME), and reads the sector-restricted palindrome through the
/// global residual M = Π·L·Π⁻¹ + L + 2Σγ·I (statements 2/3; for odd N the
/// (odd,odd) block of M is exactly Π·L_even·Π⁻¹ + L_odd + 2Σγ·I).</para>
///
/// <para>The gate is two-sided via the two scope controls (the selective-breaking
/// cross of <c>simulations/direct_sum_scope_probe.py</c>): adding amplitude
/// damping keeps the off-parity block at exactly 0 while the palindrome residual
/// jumps to O(γ_T1) (T1 breaks the mirror, not the wall); adding a transverse
/// field breaks the off-parity block at O(h) while the palindrome residual stays
/// machine-zero (the field the wall, not the mirror). A construction error
/// (scrambled parity index, wrong Π sign) fails one of the four cells.</para></summary>
public sealed class DirectSumDecompositionWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The largest Pauli-basis dimension 4^N the live build will materialise.</summary>
    public const int MaxDim = 1024;

    public int N { get; }
    public double Gamma { get; }

    /// <summary>dim(V_even) counted from the Pauli basis. Theorem: 2^(2N−1).</summary>
    public long EvenSectorDim { get; }

    /// <summary>dim(V_odd) counted from the Pauli basis. Theorem: 2^(2N−1).</summary>
    public long OddSectorDim { get; }

    /// <summary>Off-parity-block Frobenius norm of L (Pauli basis). Theorem: exactly 0
    /// (equivalently [P_XY, L] = 0; ‖[P_XY, L]‖_F = 2 × this value).</summary>
    public double OffParityNorm { get; }

    /// <summary>Full ‖L‖_F, the non-triviality scale for <see cref="OffParityNorm"/>.</summary>
    public double FullNorm { get; }

    /// <summary>True iff every Pauli string's Π image lies in the opposite parity sector
    /// for odd N (sector exchange) resp. the same sector for even N.</summary>
    public bool PiSectorMapConsistent { get; }

    /// <summary>‖M‖_F of the global palindrome residual M = Π·L·Π⁻¹ + L + 2Σγ·I. Machine zero.</summary>
    public double PalindromeResidualNorm { get; }

    /// <summary>‖M‖_F restricted to the (even,even) block. For even N this is the
    /// self-palindromy of L_even; machine zero.</summary>
    public double EvenBlockResidualNorm { get; }

    /// <summary>‖M‖_F restricted to the (odd,odd) block. For odd N this is exactly
    /// ‖Π·L_even·Π⁻¹ + L_odd + 2Σγ·I‖ (statement 2); machine zero.</summary>
    public double OddBlockResidualNorm { get; }

    /// <summary>Control cell: off-parity norm with uniform amplitude damping added.
    /// Theorem (corrected scope): exactly 0 — T1 does not break the wall.</summary>
    public double AdOffParityNorm { get; }

    /// <summary>Control cell: palindrome residual with amplitude damping. O(γ_T1) &gt; 0 —
    /// T1 breaks the mirror.</summary>
    public double AdPalindromeNorm { get; }

    /// <summary>Control cell: off-parity norm with a transverse field h·Σ X_l. O(h) &gt; 0 —
    /// the field breaks the wall.</summary>
    public double FieldOffParityNorm { get; }

    /// <summary>Control cell: palindrome residual with the transverse field. Machine zero —
    /// the field does not break the mirror.</summary>
    public double FieldPalindromeNorm { get; }

    public DirectSumDecompositionWitness(int n = 3, double gamma = 0.05)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 2 (a chain needs a bond); got {n}");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        long dim = 1L << (2 * n);
        if (dim > MaxDim)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"4^N = {dim} exceeds the live-build guard {MaxDim} for N={n}; pick N ≤ 5.");

        N = n;
        Gamma = gamma;
        int d2 = (int)dim;
        double sigma = n * gamma;
        var gammas = Enumerable.Repeat(gamma, n).ToArray();

        // Heisenberg chain, J = 1 (spin convention J/4 per Pauli bond term).
        var terms = new List<PauliTerm>();
        for (int b = 0; b < n - 1; b++)
            foreach (var letter in new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z })
                terms.Add(PauliTerm.TwoSite(n, b, letter, b + 1, letter, new Complex(0.25, 0)));
        var h = new PauliHamiltonian(n, terms).ToMatrix();

        // Transverse-field control Hamiltonian: H + 0.3 · Σ X_l.
        var fieldTerms = new List<PauliTerm>(terms);
        for (int s = 0; s < n; s++)
            fieldTerms.Add(PauliTerm.SingleSite(n, s, PauliLetter.X, new Complex(0.3, 0)));
        var hField = new PauliHamiltonian(n, fieldTerms).ToMatrix();

        var lBase = PauliDephasingDissipator.BuildZ(h, gammas);
        var lAd = T1Dissipator.Build(h, gammas, Enumerable.Repeat(0.07, n).ToArray());
        var lField = PauliDephasingDissipator.BuildZ(hField, gammas);

        var transform = PauliBasis.VecToPauliBasisTransform(N);
        ComplexMatrix ToPauli(ComplexMatrix lVec) =>
            (transform.ConjugateTranspose() * lVec * transform) / Math.Pow(2, N);

        var lPauli = ToPauli(lBase);
        var lAdPauli = ToPauli(lAd);
        var lFieldPauli = ToPauli(lField);

        // bit_a parity per flat Pauli index; sector index sets.
        var parity = new int[d2];
        var evenIdx = new List<int>();
        var oddIdx = new List<int>();
        for (int k = 0; k < d2; k++)
        {
            parity[k] = PauliIndex.TotalBitA(PauliIndex.FromFlat(k, n)) % 2;
            (parity[k] == 0 ? evenIdx : oddIdx).Add(k);
        }
        EvenSectorDim = evenIdx.Count;
        OddSectorDim = oddIdx.Count;

        OffParityNorm = OffParityBlockNorm(lPauli, parity);
        FullNorm = lPauli.FrobeniusNorm();
        AdOffParityNorm = OffParityBlockNorm(lAdPauli, parity);
        FieldOffParityNorm = OffParityBlockNorm(lFieldPauli, parity);

        // Π sector map: every column's nonzero rows must sit in the opposite (odd N)
        // resp. the same (even N) parity sector.
        var pi = ComplexMatrix.Build.DenseOfMatrix(PiOperator.BuildFull(n));
        bool exchange = n % 2 == 1;
        bool consistent = true;
        for (int col = 0; col < d2 && consistent; col++)
            for (int row = 0; row < d2; row++)
            {
                if (pi[row, col] == Complex.Zero) continue;
                bool flips = parity[row] != parity[col];
                if (flips != exchange) { consistent = false; break; }
            }
        PiSectorMapConsistent = consistent;

        // Palindrome residuals (global + sector blocks + controls).
        var identity = ComplexMatrix.Build.SparseIdentity(d2);
        ComplexMatrix Residual(ComplexMatrix lp) =>
            pi * lp * pi.ConjugateTranspose() + lp + (Complex)(2.0 * sigma) * identity;

        var m = Residual(lPauli);
        PalindromeResidualNorm = m.FrobeniusNorm();
        EvenBlockResidualNorm = BlockNorm(m, evenIdx);
        OddBlockResidualNorm = BlockNorm(m, oddIdx);
        AdPalindromeNorm = Residual(lAdPauli).FrobeniusNorm();
        FieldPalindromeNorm = Residual(lFieldPauli).FrobeniusNorm();
    }

    private static double OffParityBlockNorm(ComplexMatrix m, int[] parity)
    {
        double s = 0.0;
        for (int r = 0; r < m.RowCount; r++)
            for (int c = 0; c < m.ColumnCount; c++)
                if (parity[r] != parity[c])
                {
                    var z = m[r, c];
                    s += z.Real * z.Real + z.Imaginary * z.Imaginary;
                }
        return Math.Sqrt(s);
    }

    private static double BlockNorm(ComplexMatrix m, List<int> idx)
    {
        double s = 0.0;
        foreach (int r in idx)
            foreach (int c in idx)
            {
                var z = m[r, c];
                s += z.Real * z.Real + z.Imaginary * z.Imaginary;
            }
        return Math.Sqrt(s);
    }

    public string DisplayName =>
        $"DirectSumDecompositionWitness (N={N}, γ={Gamma.ToString("0.###", Inv)}, sectors {EvenSectorDim}+{OddSectorDim})";

    public string Summary =>
        $"N={N}: L = L_even ⊕ L_odd live — off-parity ‖·‖ = {OffParityNorm.ToString("e2", Inv)} (of ‖L‖ = " +
        $"{FullNorm.ToString("e2", Inv)}); sectors {EvenSectorDim} + {OddSectorDim} (= 2^(2N−1) each); " +
        $"Π {(N % 2 == 1 ? "exchanges" : "preserves")} sectors ({(PiSectorMapConsistent ? "verified column-complete" : "INCONSISTENT")}); " +
        $"sector palindrome residual even/odd = {EvenBlockResidualNorm.ToString("e2", Inv)}/{OddBlockResidualNorm.ToString("e2", Inv)}; " +
        $"controls: T1 wall {AdOffParityNorm.ToString("e2", Inv)} mirror {AdPalindromeNorm.ToString("e2", Inv)}, " +
        $"field wall {FieldOffParityNorm.ToString("e2", Inv)} mirror {FieldPalindromeNorm.ToString("e2", Inv)}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "statement 1+4: the wall, live",
                summary: $"off-parity-block ‖·‖_F = {OffParityNorm.ToString("e3", Inv)} (exactly 0 ⟺ [P_XY, L] = 0, " +
                         $"‖[P_XY, L]‖ = 2× this) while ‖L‖_F = {FullNorm.ToString("e3", Inv)} > 0; " +
                         $"sector dims {EvenSectorDim} + {OddSectorDim} = {EvenSectorDim + OddSectorDim} (2^(2N−1) each)");

            yield return new InspectableNode(
                displayName: $"statements 2/3: Π sector map (N={N} {(N % 2 == 1 ? "odd → exchange" : "even → preserve")})",
                summary: $"every Pauli string's Π image lands in the {(N % 2 == 1 ? "OPPOSITE" : "same")} parity sector " +
                         $"(checked column-complete: {(PiSectorMapConsistent ? "consistent" : "VIOLATED")}); sector-restricted " +
                         $"palindrome residual (even,even) = {EvenBlockResidualNorm.ToString("e3", Inv)}, (odd,odd) = " +
                         $"{OddBlockResidualNorm.ToString("e3", Inv)}" +
                         (N % 2 == 1 ? " — the (odd,odd) block is exactly ‖Π·L_even·Π⁻¹ + L_odd + 2Σγ·I‖" : " — each sector self-palindromic"));

            yield return new InspectableNode(
                displayName: "control: T1 breaks the mirror, not the wall",
                summary: $"with uniform amplitude damping 0.07: off-parity ‖·‖ = {AdOffParityNorm.ToString("e3", Inv)} " +
                         $"(the direct sum survives EXACTLY; bilinear sandwich, σ∓ bit_a-homogeneous) while the palindrome " +
                         $"residual = {AdPalindromeNorm.ToString("e3", Inv)} > 0 (the sectors persist but stop being mirror images)");

            yield return new InspectableNode(
                displayName: "control: the field breaks the wall, not the mirror",
                summary: $"with a transverse field 0.3·Σ X_l: off-parity ‖·‖ = {FieldOffParityNorm.ToString("e3", Inv)} > 0 " +
                         $"(odd-n_XY Hamiltonian term, one-sided in the commutator) while the palindrome residual = " +
                         $"{FieldPalindromeNorm.ToString("e3", Inv)} stays machine-zero");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
