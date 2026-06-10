using System;
using System.Collections.Generic;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The DERIVED k-body hidden-Q routing soft-certifier (Stufe B). A periodic per-site product Q
/// palindromizes a Z-dephased k-body chain Hamiltonian H = Σ_windows Σ_terms T iff (a) each per-site map
/// swaps the dephasing classes {I, Z} ↔ {X, Y} (automatic for the candidate set) AND (b′) at every
/// window-parity the TEMPLATE-SUMMED anticommutator {Q_k, Σ_T [T,·]_k} = 0: the window-level condition is
/// the sharp one. The PER-TERM variant (b) ({Q_k, [T,·]_k} = 0 for every template separately) is SUFFICIENT
/// but NOT necessary: a router may cancel the templates against each other inside one window
/// (cross-template), as the committed continuous-sum XIX+XIY+YIX case showed and the golden period-4 router
/// now proves discretely (docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md, F116). Either condition is checked on
/// 4^k (the term's span, NOT the 2^N Liouvillian) and additivity over windows gives the palindrome at EVERY
/// N, so the certificate is CONSTRUCTIVE (it exhibits Q), sound by derivation, and N-independent.
/// <see cref="Routes"/> carries the per-term lens; <see cref="RoutesWindowSummed"/> the window-summed one.
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
/// discrete-routable soft sets, the 2 Z-middle ceiling cases (XZX+XZY+YZX, YZY+XZY+YZX; soft, and routed by
/// the period-4 GOLDEN per-site router under the window-summed condition,
/// docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md + F116; the per-term <see cref="Routes"/> correctly returns
/// false for them, a documented coverage gap of the per-term lens, NOT non-locality), the 2 I-heavy cases
/// (IXI+IIY+YII, IYI+IIX+XII; soft, Routes returns false here, but LOCAL via the SingleSiteField strategy, a
/// site-varying single-site-field product), and the hard XXX+XXY+YXX (Python derivation 2026-06-06; the
/// decomposition residuals are 0.00e+00, the k-site residual equals the full-N residual, and the full
/// palindrome holds at N = 4, 5, 6). EVERY soft member of this family is per-site routable: the I-heavy via
/// single-site fields, the once-counted XIX+XIY+YIX, YIY+XIY+YIX via a continuous-uniform per-site Q, and
/// the Z-middle pair via the golden router (<see cref="RoutesWindowSummed"/>), with the per-term Routes
/// returning false on those six only because their routers are outside its per-term scalable strategies,
/// see experiments/CEILING_FOUR_NONLOCAL_CASES.md.</para></summary>
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

    // ---------------------------------------------------------------------------------------------------
    // The window-summed routing primitive and the golden period-4 candidates (Stufe B′, F116).
    // The per-term condition above is sufficient but not necessary; the sharp condition is the
    // TEMPLATE-SUMMED anticommutator at every window parity. The golden router passes it while failing
    // per term (cross-template cancellation inside one window, PROOF_CEILING_GOLDEN_ROUTER.md §2).
    // ---------------------------------------------------------------------------------------------------

    /// <summary>The golden ratio φ = (1+√5)/2, the frame constant of the F116 golden router (φ² = φ + 1).
    /// The frame directions a = φX + Y, b = X − φY are the two roots of the golden locus
    /// α² − αβ − β² = 0 (slopes 1/φ and −φ, frame angle tan 2θ = 2).</summary>
    private static readonly double Phi = (1.0 + Math.Sqrt(5.0)) / 2.0;

    /// <summary>One golden per-site map from its frame data g = q(I) and h = q(Z) (the {X, Y} images of the
    /// dephasing-dead letters). The {I, Z} images of the lit letters are forced by the block structure
    /// B = diag(−1, 1)·Cᵀ: q(X) = −g_X·I + h_X·Z, q(Y) = −g_Y·I + h_Y·Z. Class-swapping by construction
    /// (the {I,Z}→{I,Z} and {X,Y}→{X,Y} blocks are zero, so the dissipator condition (a) is automatic)
    /// with q² = −(2+φ)·I, a scalar times a unitary. Same basis order as <see cref="P1"/>
    /// (rows = output I,X,Y,Z; cols = input I,X,Y,Z).</summary>
    private static ComplexMatrix BuildGoldenSiteMap(Complex gx, Complex gy, Complex hx, Complex hy) =>
        Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            {   0, -gx, -gy,   0 },
            {  gx,   0,   0,  hx },
            {  gy,   0,   0,  hy },
            {   0,  hx,  hy,   0 },
        });

    /// <summary>The golden period-4 router (F116, PROOF_CEILING_GOLDEN_ROUTER.md §1): the per-site maps q_l
    /// for l = 0..3, applied as l mod 4 down the chain. g_l = q_l(I) follows the [a, a, b, b] rhythm with
    /// a = φX + Y and b = X − φY (the two golden-locus roots), and h_l = q_l(Z) = (−1)^(l+1)·i·R(g_l) with R
    /// the 90° rotation in the (X, Y) plane. The product W = ⊗_l q_{l mod 4} palindromizes the Z-middle
    /// XZX+XZY+YZX chain (W L W⁻¹ = −L − 2σ) at every N ≥ 3 and arbitrary site rates, via the window lemma
    /// (<see cref="RoutesWindowSummed"/>).</summary>
    public static IReadOnlyList<ComplexMatrix> GoldenSiteMaps { get; } = new[]
    {
        BuildGoldenSiteMap(Phi, 1.0,        Im, -Im * Phi),   // l = 0: g = a = (φ, 1),  h = −i·R(a) = (i, −iφ)
        BuildGoldenSiteMap(Phi, 1.0,       -Im,  Im * Phi),   // l = 1: g = a = (φ, 1),  h = +i·R(a) = (−i, iφ)
        BuildGoldenSiteMap(1.0, -Phi, -Im * Phi, -Im),        // l = 2: g = b = (1, −φ), h = −i·R(b) = (−iφ, −i)
        BuildGoldenSiteMap(1.0, -Phi,  Im * Phi,  Im),        // l = 3: g = b = (1, −φ), h = +i·R(b) = (iφ, i)
    };

    /// <summary>The per-site X↔Y conjugation s (I → I, X ↔ Y, Z → −Z) on the operator basis: the involution
    /// that maps each golden map to its sibling-routing twin, q′ = s·q·s.</summary>
    private static readonly ComplexMatrix XySwap = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
    {
        { 1, 0, 0,  0 },
        { 0, 0, 1,  0 },
        { 0, 1, 0,  0 },
        { 0, 0, 0, -1 },
    });

    /// <summary>The sibling candidate: the per-site X↔Y conjugation of the golden maps, q′_l = s·q_l·s.
    /// Routes the X↔Y sibling YZY+XZY+YZX; the X↔Y mirror is not a self-equivalence of the golden router,
    /// it maps one case's routers to the other's (PROOF_CEILING_GOLDEN_ROUTER.md §5). Class-swap and
    /// q′² = −(2+φ)·I are preserved by the conjugation.</summary>
    public static IReadOnlyList<ComplexMatrix> GoldenMirrorSiteMaps { get; } = BuildGoldenMirrorSiteMaps();

    private static ComplexMatrix[] BuildGoldenMirrorSiteMaps()
    {
        var mirror = new ComplexMatrix[GoldenSiteMaps.Count];
        for (int l = 0; l < mirror.Length; l++)
            mirror[l] = XySwap * GoldenSiteMaps[l] * XySwap;
        return mirror;
    }

    /// <summary>The golden pattern period: the maps repeat as q_{l mod 4} down the chain.</summary>
    public const int GoldenPeriod = 4;

    /// <summary>The certificate descriptions the window-summed router reports (mirroring
    /// <see cref="DescribePattern"/>'s naming for the per-term candidates).</summary>
    public const string GoldenDescription = "Golden[a,a,b,b] (P=4)";
    public const string GoldenMirrorDescription = "Golden-mirror[a,a,b,b] (P=4)";

    /// <summary>Q_k from an explicit pattern of per-site MATRICES (the golden maps are not members of
    /// <see cref="Representatives"/>): the Kronecker product of pattern[(offset + j) mod period] over the
    /// span sites j = 0 .. k − 1, in the same {I,X,Y,Z}^⊗k Pauli-string basis as <see cref="BuildQk"/>
    /// (most-significant digit = the leftmost span site).</summary>
    public static ComplexMatrix BuildQkFromMaps(IReadOnlyList<ComplexMatrix> pattern, int period, int offset, int k)
    {
        if (pattern is null) throw new ArgumentNullException(nameof(pattern));
        if (period <= 0 || period > pattern.Count) throw new ArgumentOutOfRangeException(nameof(period));
        if (k <= 0) throw new ArgumentOutOfRangeException(nameof(k));

        ComplexMatrix q = pattern[((offset + 0) % period + period) % period];
        for (int j = 1; j < k; j++)
            q = q.KroneckerProduct(pattern[((offset + j) % period + period) % period]);
        return q;
    }

    /// <summary>The TEMPLATE-SUMMED commutator superoperator Σ_T c_T·[T,·]_k on the shared 4^k window
    /// space: the per-window Hamiltonian action of the whole term-set. Coefficient-WEIGHTED, because the
    /// window-summed cancellation is cross-template (unlike the per-term lens, which one Q serves for any
    /// real linear combination): {Q_k, Σ_T c_T·[T,·]_k} = 0 certifies exactly the placed sum Σ_T c_T·T, not
    /// each template alone. All templates must share one span k (mixed spans have no common window space).</summary>
    public static ComplexMatrix BuildSummedCommutatorSuperoperator(IReadOnlyList<PauliTerm> templates)
    {
        if (templates is null) throw new ArgumentNullException(nameof(templates));
        if (templates.Count == 0) throw new ArgumentException("at least one template is required", nameof(templates));

        int k = templates[0].Letters.Count;
        ComplexMatrix? sum = null;
        foreach (var t in templates)
        {
            if (t.Letters.Count != k)
                throw new ArgumentException("all templates must share the same span k (mixed spans have no common window space)", nameof(templates));
            ComplexMatrix comm = BuildCommutatorSuperoperator(t) * t.Coefficient;
            sum = sum is null ? comm : sum + comm;
        }
        return sum!;
    }

    /// <summary>The k-site WINDOW-SUMMED condition (b′): build Q_k for the matrix pattern at this offset
    /// and the template-summed superoperator S = Σ_T c_T·[T,·]_k, and return true iff
    /// ‖Q_k·S + S·Q_k‖_F &lt; tol (the summed anticommutator vanishes). The per-term condition
    /// (<see cref="PerTermAnticommutes"/>) is the special case where every summand vanishes alone; this is
    /// the SHARP one (the golden router passes it at every offset while failing per term, the
    /// cross-template window lemma of PROOF_CEILING_GOLDEN_ROUTER.md §2). A single-template list
    /// degenerates to that template's per-term check under the same matrix pattern.</summary>
    public static bool PerWindowSummedAnticommutes(
        IReadOnlyList<PauliTerm> templates, int offset, IReadOnlyList<ComplexMatrix> pattern, int period)
    {
        ComplexMatrix summed = BuildSummedCommutatorSuperoperator(templates);
        ComplexMatrix qk = BuildQkFromMaps(pattern, period, offset, templates[0].Letters.Count);
        ComplexMatrix anti = qk * summed + summed * qk;
        return anti.FrobeniusNorm() < AnticommuteTolerance;
    }

    /// <summary>The window-summed routing decision (Stufe B′, F116): try the two golden period-4 candidates
    /// (<see cref="GoldenSiteMaps"/> and its X↔Y conjugate <see cref="GoldenMirrorSiteMaps"/>) and return
    /// the description of the FIRST whose summed anticommutator vanishes at EVERY window offset 0..3, or
    /// null if neither does. Each candidate map is class-swapping (condition (a), the dissipator leg,
    /// automatic) with q² = −(2+φ)·I (so W = ⊗_l q_{l mod 4} is invertible), and the vanishing summed
    /// anticommutator at all four offsets is the window lemma: additivity over windows then gives
    /// W L W⁻¹ = −L − 2σ at EVERY N ≥ k, so the certificate is CONSTRUCTIVE, sound by derivation, and
    /// N-independent (PROOF_CEILING_GOLDEN_ROUTER.md §3). This is the lens that certifies the two Z-middle
    /// ceiling cases (XZX+XZY+YZX golden, YZY+XZY+YZX golden-mirror) the per-term <see cref="Routes"/>
    /// correctly declines.
    ///
    /// <para>Gates: every template must share one span k (a mixed-span set has no common window space and
    /// is declined), with k ≤ <see cref="MaxBody"/> and k ≤ <paramref name="n"/>.</para></summary>
    public static string? RoutesWindowSummed(IReadOnlyList<PauliTerm> terms, int n)
    {
        if (terms is null || terms.Count == 0) return null;

        int k = terms[0].Letters.Count;
        foreach (var term in terms)
            if (term.Letters.Count != k)
                return null;                                    // mixed spans: no shared window space
        if (k > MaxBody || k > n) return null;

        ComplexMatrix summed = BuildSummedCommutatorSuperoperator(terms);
        foreach (var (maps, description) in new[]
        {
            (GoldenSiteMaps, GoldenDescription),
            (GoldenMirrorSiteMaps, GoldenMirrorDescription),
        })
        {
            bool allOffsets = true;
            for (int offset = 0; offset < GoldenPeriod && allOffsets; offset++)
            {
                ComplexMatrix qk = BuildQkFromMaps(maps, GoldenPeriod, offset, k);
                ComplexMatrix anti = qk * summed + summed * qk;
                if (anti.FrobeniusNorm() >= AnticommuteTolerance)
                    allOffsets = false;
            }
            if (allOffsets) return description;
        }
        return null;
    }
}
