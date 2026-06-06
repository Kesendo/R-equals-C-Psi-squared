using System;
using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The DERIVED k-body hidden-Q routing soft-certifier (Stufe B). A periodic per-site product Q
/// palindromizes a Z-dephased k-body chain Hamiltonian H = Σ_windows Σ_terms T iff (a) each per-site map
/// swaps the dephasing classes {I, Z} ↔ {X, Y} (automatic for the candidate set) AND (b) for every term and
/// every window-parity the k-site anticommutator {Q_k, [T,·]_k} = 0. Condition (b) is checked on 4^k (the
/// term's span, NOT the 2^N Liouvillian) and additivity over windows gives the palindrome at EVERY N, so the
/// certificate is CONSTRUCTIVE (it exhibits Q), sound by derivation, and N-independent.
///
/// <para><b>The convention-safety point</b>: both Q_k and [T,·]_k are built as 4^k × 4^k matrices in the
/// SAME Pauli-string basis, the k-fold tensor of the single-site operator basis {I, X, Y, Z} (index I = 0,
/// X = 1, Y = 2, Z = 3; a k-site string is a base-4 index 0 .. 4^k − 1, most-significant digit = the
/// leftmost / first span site, matching the big-endian Core convention). [T,·]_k is built from Pauli
/// COMMUTATOR ALGEBRA ([T, σ_S] = 0 if T and σ_S commute, else 2·(phase)·σ_{T·S} via the Pauli-string
/// product), NOT the computational-basis T ⊗ I − I ⊗ Tᵀ (a different basis that silently breaks the
/// anticommutator check). This ordering is local to this file and intentionally independent of
/// <see cref="PauliLetter"/>'s a + 2·b packing.</para>
///
/// <para>Verified bit-exact against the spectral authority <see cref="PauliPairTrichotomy"/> on the 8
/// discrete-routable soft sets, the 4 non-local ceiling cases (XZX+XZY+YZX, YZY+XZY+YZX, IXI+IIY+YII,
/// IYI+IIX+XII; soft but no per-site Q), and the hard XXX+XXY+YXX (Python derivation 2026-06-06; the
/// decomposition residuals are 0.00e+00, the k-site residual equals the full-N residual, and the full
/// palindrome holds at N = 4, 5, 6). These 4 admit no per-site product Q at all and stay NotCertified; the
/// 2 cases once counted with them, XIX+XIY+YIX and YIY+XIY+YIX, are LOCAL (a continuous-uniform per-site Q
/// palindromizes them, verified; NotCertified only because that Q routes via continuous-sum, outside the
/// scalable strategies), see experiments/CEILING_FOUR_NONLOCAL_CASES.md.</para></summary>
public static class KBodyPalindromeRouting
{
    /// <summary>The single-site operator basis order for THIS file: I = 0, X = 1, Y = 2, Z = 3. The per-site
    /// maps below and the commutator superoperator share this 4-element order, so Q_k and [T,·]_k live in the
    /// same Pauli-string basis. (Distinct from <see cref="PauliLetter"/>'s I = 0, X = 1, Z = 2, Y = 3.)</summary>
    private const int I = 0, X = 1, Y = 2, Z = 3;

    /// <summary>The Frobenius-norm tolerance for the vanishing anticommutator ‖{Q_k, [T,·]_k}‖_F &lt; tol.
    /// The check is algebraically exact (signed permutations with ±1, ±i entries, or the dense M with
    /// 1/√2 entries); this only absorbs float round-off in the Kronecker products and matrix multiplies.</summary>
    private const double AnticommuteTolerance = 1e-9;

    /// <summary>The default and only span bound: a term of span k yields a 4^k × 4^k check (4^5 = 1024²).
    /// Terms of span &gt; 5 are out of scope (the strategy declines them upstream). Consumed by the Task 2
    /// <c>CertifyByRoutingKBody</c> span gate, not by this primitive's <see cref="Routes"/> (which gates only
    /// on span ≤ n).</summary>
    public const int MaxBody = 5;

    private static readonly Complex Im = Complex.ImaginaryOne;

    // ---------------------------------------------------------------------------------------------------
    // The four per-site palindrome maps (4 × 4 complex, rows = output I,X,Y,Z; cols = input I,X,Y,Z).
    // All four are class-swapping ({I,Z} ↔ {X,Y}), so the dissipator condition (a) is automatic; only (b)
    // (the per-term anticommutator) is checked. The matrix literals match the spec exactly.
    // ---------------------------------------------------------------------------------------------------

    /// <summary>P1: signed permutation I↔X, Y↔Z; signs {I:1, X:1, Y:i, Z:i}.</summary>
    private static readonly ComplexMatrix P1 = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
    {
        { 0,  1,  0,  0 },
        { 1,  0,  0,  0 },
        { 0,  0,  0, Im },
        { 0,  0, Im,  0 },
    });

    /// <summary>P4: signed permutation I↔Y, X↔Z; signs {I:1, X:i, Y:1, Z:i}.</summary>
    private static readonly ComplexMatrix P4 = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
    {
        { 0,  0,  1,  0 },
        { 0,  0,  0, Im },
        { 1,  0,  0,  0 },
        { 0, Im,  0,  0 },
    });

    /// <summary>M2: same permutation as P4 (I↔Y, X↔Z); signs {I:1, X:−i, Y:1, Z:−i}.</summary>
    private static readonly ComplexMatrix M2 = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
    {
        { 0,   0,  1,   0 },
        { 0,   0,  0, -Im },
        { 1,   0,  0,   0 },
        { 0, -Im,  0,   0 },
    });

    /// <summary>M: the continuous, dense class-swapping map with M² = −I. With s = 1/√2:
    /// I → −(X+Y)/√2, X → (I + iZ)/√2, Y → (I − iZ)/√2, Z → i(X − Y)/√2.</summary>
    private static readonly ComplexMatrix M = BuildM();

    private static ComplexMatrix BuildM()
    {
        double s = 1.0 / Math.Sqrt(2.0);
        Complex iS = Im * s;
        return Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            {  0,   s,   s,   0 },
            { -s,   0,   0,  iS },
            { -s,   0,   0, -iS },
            {  0,  iS, -iS,   0 },
        });
    }

    /// <summary>The four candidate per-site representatives, in routing-preference order P1, P4, M2, M.</summary>
    private static readonly ComplexMatrix[] Representatives = { P1, P4, M2, M };
    private static readonly string[] RepresentativeNames = { "P1", "P4", "M2", "M" };

    // ---------------------------------------------------------------------------------------------------
    // Single-site Pauli product algebra (standard physical convention: σ_x σ_y = i σ_z and cyclic), in this
    // file's {I,X,Y,Z} = {0,1,2,3} order. prod[a,b] = (result letter, phase) with σ_a σ_b = phase · σ_result.
    // ---------------------------------------------------------------------------------------------------
    private static readonly (int Letter, Complex Phase)[,] SiteProduct = BuildSiteProduct();

    private static (int Letter, Complex Phase)[,] BuildSiteProduct()
    {
        var p = new (int, Complex)[4, 4];
        // Identity row/column and the diagonal squares.
        for (int a = 0; a < 4; a++)
        {
            p[I, a] = (a, Complex.One);      // I · σ_a = σ_a
            p[a, I] = (a, Complex.One);      // σ_a · I = σ_a
            p[a, a] = (I, Complex.One);      // σ_a² = I
        }
        // The off-diagonal pairs among {X, Y, Z}: σ_x σ_y = i σ_z and cyclic.
        p[X, Y] = (Z,  Im); p[Y, X] = (Z, -Im);
        p[Y, Z] = (X,  Im); p[Z, Y] = (X, -Im);
        p[Z, X] = (Y,  Im); p[X, Z] = (Y, -Im);
        return p;
    }

    /// <summary>The product of two k-site Pauli strings (base-4 indices, most-significant digit = site 0):
    /// σ_a σ_b = phase · σ_result, computed site-by-site and accumulated. Returns (result string index,
    /// total phase).</summary>
    private static (int Result, Complex Phase) StringProduct(int a, int b, int k)
    {
        int result = 0;
        Complex phase = Complex.One;
        // Walk sites most-significant (site 0) to least-significant (site k-1).
        for (int site = 0; site < k; site++)
        {
            int shift = 2 * (k - 1 - site);
            int la = (a >> shift) & 0b11;
            int lb = (b >> shift) & 0b11;
            var (lr, ph) = SiteProduct[la, lb];
            result |= lr << shift;
            phase *= ph;
        }
        return (result, phase);
    }

    // ---------------------------------------------------------------------------------------------------
    // Q_k and [T,·]_k builders.
    // ---------------------------------------------------------------------------------------------------

    /// <summary>Q_k = the Kronecker product of the per-site maps the periodic pattern assigns to the span
    /// sites (offset + 0 .. offset + k − 1) mod period. <paramref name="pattern"/>[i] is an index into
    /// <see cref="Representatives"/>; site j of the term takes pattern[(offset + j) mod period].</summary>
    public static ComplexMatrix BuildQk(int[] pattern, int period, int offset, int k)
    {
        if (pattern is null) throw new ArgumentNullException(nameof(pattern));
        if (period <= 0) throw new ArgumentOutOfRangeException(nameof(period));
        if (k <= 0) throw new ArgumentOutOfRangeException(nameof(k));

        ComplexMatrix q = Representatives[pattern[((offset + 0) % period + period) % period]];
        for (int j = 1; j < k; j++)
        {
            int rep = pattern[((offset + j) % period + period) % period];
            q = q.KroneckerProduct(Representatives[rep]);
        }
        return q;
    }

    /// <summary>The k-site commutator superoperator [T,·]_k as a 4^k × 4^k matrix in the {I,X,Y,Z}^⊗k
    /// Pauli-string basis, built from Pauli COMMUTATOR ALGEBRA (not the computational-basis form). For each
    /// input string S (column), [T, σ_S] = 0 if T and σ_S commute, else (phase_TS − phase_ST)·σ_{T·S}; since
    /// T·S and S·T are the same string, the column has a single nonzero entry at row (T·S) with value
    /// 2·phase_TS when they anticommute, or is zero when they commute.</summary>
    public static ComplexMatrix BuildCommutatorSuperoperator(PauliTerm template)
    {
        if (template is null) throw new ArgumentNullException(nameof(template));
        int k = template.Letters.Count;
        if (k <= 0) throw new ArgumentException("template must have at least one letter", nameof(template));

        int t = TemplateToStringIndex(template, k);
        int dim = 1 << (2 * k);                                  // 4^k
        var sup = Matrix<Complex>.Build.Dense(dim, dim);

        for (int s = 0; s < dim; s++)
        {
            var (ts, phaseTs) = StringProduct(t, s, k);          // T · σ_S = phaseTs · σ_{ts}
            var (st, phaseSt) = StringProduct(s, t, k);          // σ_S · T = phaseSt · σ_{st}  (st == ts)
            Complex entry = phaseTs - phaseSt;                   // 0 if commute, 2·phaseTs if anticommute
            if (entry != Complex.Zero)
                sup[ts, s] = entry;                              // column S maps onto row T·S
        }
        return sup;
    }

    /// <summary>Encode a template's letters as a base-4 string index in THIS file's {I,X,Y,Z} = {0,1,2,3}
    /// order (most-significant digit = the first / leftmost span letter), translating from the Core
    /// <see cref="PauliLetter"/> packing.</summary>
    private static int TemplateToStringIndex(PauliTerm template, int k)
    {
        int idx = 0;
        for (int site = 0; site < k; site++)
            idx = (idx << 2) | LocalLetter(template.Letters[site]);
        return idx;
    }

    /// <summary>Map a Core <see cref="PauliLetter"/> to this file's {I:0, X:1, Y:2, Z:3} index.</summary>
    private static int LocalLetter(PauliLetter letter) => letter switch
    {
        PauliLetter.I => I,
        PauliLetter.X => X,
        PauliLetter.Y => Y,
        PauliLetter.Z => Z,
        _ => throw new ArgumentOutOfRangeException(nameof(letter)),
    };

    // ---------------------------------------------------------------------------------------------------
    // The per-term check and the candidate set.
    // ---------------------------------------------------------------------------------------------------

    /// <summary>The k-site per-term condition (b): build Q_k for the pattern at this offset and the template's
    /// commutator superoperator, and return true iff ‖Q_k [T,·]_k + [T,·]_k Q_k‖_F &lt; tol (the anticommutator
    /// vanishes). N-independent; the span k is the template's letter count.</summary>
    public static bool PerTermAnticommutes(PauliTerm template, int offset, int[] pattern, int period)
    {
        int k = template.Letters.Count;
        ComplexMatrix qk = BuildQk(pattern, period, offset, k);
        ComplexMatrix comm = BuildCommutatorSuperoperator(template);
        ComplexMatrix anti = qk * comm + comm * qk;
        return anti.FrobeniusNorm() < AnticommuteTolerance;
    }

    /// <summary>A candidate Q: a periodic pattern (each entry an index into <see cref="Representatives"/>)
    /// with its period.</summary>
    // NOTE: only enumerated, never compared/hashed by value: the compiler-generated value-equality is
    // reference equality on the Pattern array (two identical patterns compare unequal), so callers must not
    // use this as a dictionary key without supplying structural equality.
    public readonly record struct Candidate(int[] Pattern, int Period, string Description);

    /// <summary>The bounded candidate-Q set: every periodic pattern map[i mod P] over the four
    /// representatives {P1, P4, M2, M} for period P in 1 .. <paramref name="maxPeriod"/>. 4^P patterns at
    /// period P (4 representatives {P1, P4, M2, M}, P slots). (Default <c>maxPeriod = 2</c> for the strategy;
    /// Task 1 evaluates 3 in the sweep.)</summary>
    public static IReadOnlyList<Candidate> CandidateSet(int maxPeriod = 2)
    {
        var candidates = new List<Candidate>();
        int r = Representatives.Length;
        for (int period = 1; period <= maxPeriod; period++)
        {
            int total = 1;
            for (int i = 0; i < period; i++) total *= r;       // r^period patterns
            for (int code = 0; code < total; code++)
            {
                var pattern = new int[period];
                int c = code;
                for (int i = 0; i < period; i++) { pattern[i] = c % r; c /= r; }
                candidates.Add(new Candidate(pattern, period, DescribePattern(pattern)));
            }
        }
        return candidates;
    }

    private static string DescribePattern(int[] pattern)
    {
        var names = new string[pattern.Length];
        for (int i = 0; i < pattern.Length; i++) names[i] = RepresentativeNames[pattern[i]];
        return string.Join("⊗", names) + (pattern.Length > 1 ? $" (P={pattern.Length})" : "");
    }

    /// <summary>True iff some candidate (pattern, period P) routes EVERY term at EVERY offset 0 .. P − 1,
    /// i.e. <see cref="PerTermAnticommutes"/> holds for all (term, offset) under that one candidate. This is
    /// the DERIVED per-term k-site routing condition; the chain length <paramref name="n"/> is not used by
    /// the check (N-independence is the derived guarantee), only to confirm each template fits (span ≤ n).</summary>
    /// <param name="terms">The term templates (each <see cref="PauliTerm.Letters"/> is the term's span).</param>
    /// <param name="n">The chain length; used only to confirm every template's span ≤ n.</param>
    public static bool Routes(IReadOnlyList<PauliTerm> terms, int n) => Routes(terms, n, maxPeriod: 2);

    /// <summary>The candidate-set-parameterized routing decision (the public <see cref="Routes(IReadOnlyList{PauliTerm}, int)"/>
    /// fixes <paramref name="maxPeriod"/> = 2). Exposed so the soundness sweep can measure whether a larger
    /// period or a richer pattern family raises coverage.</summary>
    public static bool Routes(IReadOnlyList<PauliTerm> terms, int n, int maxPeriod)
    {
        if (terms is null) throw new ArgumentNullException(nameof(terms));
        if (terms.Count == 0) return true;                      // the empty Hamiltonian is trivially palindromic

        foreach (var term in terms)
            if (term.Letters.Count > n)
                return false;                                   // a template wider than the chain cannot be placed

        // Precompute each term's commutator superoperator once (offset-independent).
        var commByTerm = new ComplexMatrix[terms.Count];
        for (int i = 0; i < terms.Count; i++)
            commByTerm[i] = BuildCommutatorSuperoperator(terms[i]);

        foreach (var cand in CandidateSet(maxPeriod))
        {
            if (RoutesUnder(cand, terms, commByTerm))
                return true;
        }
        return false;
    }

    /// <summary>True iff the single candidate <paramref name="cand"/> makes the per-term anticommutator vanish
    /// for every term at every window-offset 0 .. period − 1.</summary>
    private static bool RoutesUnder(Candidate cand, IReadOnlyList<PauliTerm> terms, ComplexMatrix[] commByTerm)
    {
        for (int i = 0; i < terms.Count; i++)
        {
            int k = terms[i].Letters.Count;
            for (int offset = 0; offset < cand.Period; offset++)
            {
                ComplexMatrix qk = BuildQk(cand.Pattern, cand.Period, offset, k);
                ComplexMatrix anti = qk * commByTerm[i] + commByTerm[i] * qk;
                if (anti.FrobeniusNorm() >= AnticommuteTolerance)
                    return false;
            }
        }
        return true;
    }

    /// <summary>The description of the FIRST candidate that routes the term-set, or null if none does. A
    /// reporting helper for the soundness sweep (which Q exhibited the certificate).</summary>
    public static string? RoutingCandidate(IReadOnlyList<PauliTerm> terms, int n, int maxPeriod = 2)
    {
        if (terms is null || terms.Count == 0) return null;
        foreach (var term in terms)
            if (term.Letters.Count > n) return null;

        var commByTerm = new ComplexMatrix[terms.Count];
        for (int i = 0; i < terms.Count; i++)
            commByTerm[i] = BuildCommutatorSuperoperator(terms[i]);

        foreach (var cand in CandidateSet(maxPeriod))
            if (RoutesUnder(cand, terms, commByTerm))
                return cand.Description;
        return null;
    }
}
