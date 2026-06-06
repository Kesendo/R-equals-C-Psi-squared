using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The uniform Q-family that routes the hidden palindrome symmetry of a soft (or truly) two-term
/// bond bilinear, or <see cref="None_"/> when the palindrome is broken (hard).
///
/// <list type="bullet">
///   <item><see cref="P1"/>: truly. The canonical palindrome operator Π = P1 already pairs the spectrum.</item>
///   <item><see cref="Uniform"/>: soft. The two terms share a single uniform product Q-family.</item>
///   <item><see cref="Alternating"/>: soft. An XY/YX term closes via a site-parity (alternating) Q.</item>
///   <item><see cref="Continuous"/>: soft. A same-site X&amp;Y collision over a shared dark Z port closes
///         via a local continuous per-site rotation M (M² = −I).</item>
///   <item><see cref="None_"/>: hard. No Q closes the palindrome.</item>
/// </list>
/// </summary>
public enum QFamily { P1, Uniform, Alternating, Continuous, None_ }

/// <summary>The Klein-routed verdict of a two-term bond bilinear: its spectral fate, the hidden-symmetry
/// Q family, both terms' Klein indices, the global-spin-flip parity-break flag, and a structural reason.
///
/// <para>Produced by <see cref="TwoTermPalindromeRouting.Classify"/>; the field names mirror the dict
/// returned by the Python reference <c>classify_two_term_palindrome</c>.</para></summary>
public sealed record RoutingResult(
    string Term1,
    string Term2,
    (int A, int B) Klein1,
    (int A, int B) Klein2,
    bool ParityBreak,
    TrichotomyClass Fate,
    QFamily Family,
    string Reason);

/// <summary>A Liouvillian-free classifier for the palindrome fate of a Z-dephased two-term bond bilinear
/// H = Σ_bonds (t1 + t2), deciding the fate (truly / soft / hard) AND routing the hidden symmetry Q into a
/// family from the two bilinears' letters alone, with no eigensolve.
///
/// <para>Z-dephasing splits each site into a DC bus {I, Z} (immune) and an AC bus {X, Y} (damped). The
/// palindrome operator Q routes the AC bus, and per-site it comes in two crossover families: P1 routes the
/// X-channel, P4 routes the Y-channel. A single uniform product Q must serve BOTH bond terms at once, so
/// the chain is palindromic iff its two terms share a valid uniform Q-family; failing that it escapes via
/// an alternating (site-parity) Q, or, for the two same-site X&amp;Y collisions, via a continuous
/// (non-permutation) per-site rotation; failing THAT the palindrome is broken (hard).</para>
///
/// <para>This routing reproduces <see cref="PauliPairTrichotomy.Classify(Core.ChainSystems.ChainSystem,
/// IReadOnlyList{PauliTerm}, double, double, PauliLetter)"/> bit-exactly across all 36 two-term combos.
/// The routing is N-independent; the optional <c>n</c> drives only the parity-break cross-check
/// ‖[H, X^⊗N]‖ &gt; 0. Ported from <c>simulations/framework/diagnostics/q_family_routing.py</c>.</para></summary>
public static class TwoTermPalindromeRouting
{
    private const double ParityTolerance = 1e-9;

    // fam(t) = the set of uniform Q-families that close term t alone:
    //   fam(XX) = fam(YY) = {P1, P4}   both crossovers route a same-letter bilinear
    //   fam(ZZ) = {all}                ZZ is dephasing-aligned, every Q fixes it
    //   fam(YZ) = fam(ZY) = {P1}       only the X-router survives the Y·Z mix
    //   fam(XZ) = fam(ZX) = {P4}       only the Y-router survives the X·Z mix
    //   fam(XY) = fam(YX) = {}         no single uniform Q routes both X and Y
    // The three Q-symbols are represented as bits {P1, P4, M2} = {1, 2, 4}.
    private const int P1 = 1, P4 = 2, M2 = 4;
    private const int AllQ = P1 | P4 | M2;

    private static readonly Dictionary<string, int> Fam = new()
    {
        ["XX"] = P1 | P4,
        ["YY"] = P1 | P4,
        ["ZZ"] = AllQ,
        ["YZ"] = P1,
        ["ZY"] = P1,
        ["XZ"] = P4,
        ["ZX"] = P4,
        ["XY"] = 0,
        ["YX"] = 0,
    };

    /// <summary>Non-throwing lookup of a single bilinear's uniform Q-family mask: the bits {P1, P4, M2} =
    /// {1, 2, 4} of the uniform product Q-families that close <paramref name="term"/> alone (see the
    /// <see cref="Fam"/> table). Returns true with <paramref name="mask"/> = the family bits if
    /// <paramref name="term"/> is one of the nine recognized bilinears (fam(XY) = fam(YX) = 0 still returns
    /// true, recognized with an empty family); false with <paramref name="mask"/> = 0 for any other label
    /// (a k-body template like "XZX", a padded label like "XIZ", or an unknown symbol). Unlike
    /// <see cref="NormalizeTerm"/>, this never throws, so it is safe as a multi-term scope gate.</summary>
    /// <param name="term">A Pauli label to look up (only the nine 2-letter bilinears are recognized).</param>
    /// <param name="mask">The family bits P1|P4|M2 if recognized, else 0.</param>
    /// <returns>True iff <paramref name="term"/> is a recognized bilinear.</returns>
    public static bool TryGetUniformFamilyMask(string term, out int mask)
    {
        if (term is not null && Fam.TryGetValue(term, out mask))
            return true;
        mask = 0;
        return false;
    }

    // The two continuous-crossover escapes: same-site X&Y collision over a shared dark Z port. No discrete
    // signed-permutation crossover routes both bands, but a single continuous per-site rotation M (M² = −I)
    // does, so the mirror is LOCAL, not entangled.
    private static readonly HashSet<string> ContinuousPairs = new() { "XZ|YZ", "ZX|ZY" };

    private static readonly HashSet<string> Mother = new() { "XX", "YY", "ZZ" };

    /// <summary>Unordered key "a|b" for a pair of bilinears (sorted so order does not matter).</summary>
    private static string PairKey(string a, string b) =>
        string.CompareOrdinal(a, b) <= 0 ? $"{a}|{b}" : $"{b}|{a}";

    /// <summary>Classify the palindrome fate of H = Σ_bonds (term1 + term2) under Z-dephasing.
    ///
    /// <para>Decides the fate (truly / soft / hard) and the Klein-routing of the hidden symmetry Q from
    /// the two bilinears' letters alone, no eigensolve. The routing is N-independent; <paramref name="n"/>
    /// drives only the parity-break cross-check.</para></summary>
    /// <param name="term1">A two-letter Pauli label (e.g. "XY"); letters from I/X/Y/Z.</param>
    /// <param name="term2">A two-letter Pauli label.</param>
    /// <param name="n">Chain length used only for the parity-break test (default 3).</param>
    /// <returns>The <see cref="RoutingResult"/> with fate, Q family, Klein indices, parity-break flag, reason.</returns>
    /// <exception cref="System.ArgumentException">If a term is not exactly two recognized Pauli letters.</exception>
    public static RoutingResult Classify(string term1, string term2, int n = 3)
    {
        string a = NormalizeTerm(term1);
        string b = NormalizeTerm(term2);

        var (fate, family, ruleReason) = Route(a, b);
        string reason = fate == TrichotomyClass.Hard ? Q7Reason(a, b) : ruleReason;

        return new RoutingResult(
            Term1: a,
            Term2: b,
            Klein1: KleinIndex(a),
            Klein2: KleinIndex(b),
            ParityBreak: ParityBreaks(n, a, b),
            Fate: fate,
            Family: family,
            Reason: reason);
    }

    /// <summary>Validate a two-letter Pauli label: each letter must be one of I/X/Y/Z (case-sensitive,
    /// matching the Python reference, which rejects lower-case) and the term must have exactly two letters.</summary>
    private static string NormalizeTerm(string term)
    {
        if (term is null) throw new System.ArgumentNullException(nameof(term));
        if (term.Length != 2)
            throw new System.ArgumentException($"term must be two letters; got '{term}'", nameof(term));
        foreach (char ch in term)
            if (ch is not ('I' or 'X' or 'Y' or 'Z'))
                throw new System.ArgumentException($"unknown Pauli letter '{ch}' in term '{term}'", nameof(term));
        return term;
    }

    /// <summary>Klein-Vierergruppe (Z₂ × Z₂) index (bit_a, bit_b) of a two-letter bilinear, via the Core
    /// <see cref="PauliTerm.KleinIndex"/> (bit_a = #(X+Y) mod 2, bit_b = #(Y+Z) mod 2).</summary>
    private static (int A, int B) KleinIndex(string term)
    {
        var letters = PauliLabel.Parse(term);
        var t = new PauliTerm(letters, Complex.One);
        var (bitA, bitB) = t.KleinIndex;
        return (bitA, bitB);
    }

    // -----------------------------------------------------------------------------------------------
    // Routing rule (R2b): predict the fate + Q family from the letters alone.
    // -----------------------------------------------------------------------------------------------
    private static (TrichotomyClass Fate, QFamily Family, string Reason) Route(string a, string b)
    {
        int fa = Fam[a], fb = Fam[b];
        int shared = fa & fb;

        // 1. Shared uniform Q-family: the spectrum pairs under one product Q.
        if (shared != 0)
        {
            if (Mother.Contains(a) && Mother.Contains(b))
                return (TrichotomyClass.Truly, QFamily.P1,
                    "both Mother (XX/YY/ZZ); canonical Π = P1 pairs spectrum");
            return (TrichotomyClass.Soft, QFamily.Uniform,
                $"shared uniform Q-family {{{QSymbols(shared)}}}");
        }

        // 2. Alternating (site-parity) escape: an XY/YX term plus an XY/YX or ZZ partner closes via
        //    P1 ⊗ M2 ⊗ P1 ⊗ ... (alternating per-site crossover).
        bool aIsFather = a is "XY" or "YX";
        bool bIsFather = b is "XY" or "YX";
        bool aInSet = a is "XY" or "YX" or "ZZ";
        bool bInSet = b is "XY" or "YX" or "ZZ";
        if ((aIsFather || bIsFather) && aInSet && bInSet)
            return (TrichotomyClass.Soft, QFamily.Alternating,
                "XY/YX with site-parity (alternating) Q");

        // 3. Continuous-crossover escape: same-site X&Y over a shared dark Z port, closed by a continuous
        //    (non-permutation) uniform per-site rotation. Still LOCAL.
        if (ContinuousPairs.Contains(PairKey(a, b)))
            return (TrichotomyClass.Soft, QFamily.Continuous,
                "same-site X&Y collision over shared dark Z; local continuous per-site rotation M (M²=−I)");

        // 4. No Q closes it: the palindrome is broken.
        return (TrichotomyClass.Hard, QFamily.None_,
            "no shared uniform Q, no alternating/continuous escape");
    }

    /// <summary>Render a set of Q-symbol bits as a sorted "P1, P4, M2" string (matching the Python's
    /// sorted-set rendering of the shared family).</summary>
    private static string QSymbols(int bits)
    {
        var names = new List<string>();
        if ((bits & M2) != 0) names.Add("M2");
        if ((bits & P1) != 0) names.Add("P1");
        if ((bits & P4) != 0) names.Add("P4");
        names.Sort(System.StringComparer.Ordinal);
        return string.Join(", ", names);
    }

    // -----------------------------------------------------------------------------------------------
    // Second condition (R2c): WHY a hard combo is hard (structural Klein-cell reason).
    //   M   = (0,0) = {XX, YY, ZZ}   Mother
    //   F_a = (0,1) = {XY, YX}       Y-Father  (Π²-odd)
    //   C   = (1,0) = {YZ, ZY}       Child     (dephase-aligned)
    //   F_b = (1,1) = {XZ, ZX}       Z-Father  (Π²-odd)
    // -----------------------------------------------------------------------------------------------
    private static readonly Dictionary<(int A, int B), string> KleinLabels = new()
    {
        [(0, 0)] = "M",
        [(0, 1)] = "F_a",
        [(1, 0)] = "C",
        [(1, 1)] = "F_b",
    };

    /// <summary>The unordered Klein-cell label pair of the two terms, e.g. ("F_a", "M").</summary>
    private static (string, string) CellPair(string a, string b)
    {
        string la = KleinLabels[KleinIndex(a)];
        string lb = KleinLabels[KleinIndex(b)];
        return string.CompareOrdinal(la, lb) <= 0 ? (la, lb) : (lb, la);
    }

    /// <summary>True iff one single qubit is asked to carry BOTH an X and a Y (across the two bond terms).</summary>
    private static bool SameSiteXY(string a, string b)
    {
        var site0 = new HashSet<char> { a[0], b[0] };
        var site1 = new HashSet<char> { a[1], b[1] };
        return (site0.Contains('X') && site0.Contains('Y'))
            || (site1.Contains('X') && site1.Contains('Y'));
    }

    /// <summary>The structural reason a combo is hard, via Klein cells. Mirrors the Python <c>_q7_reason</c>;
    /// the hard-by-reason set equals the 14 hard combos bit-exactly. (Only consulted for hard fates.)</summary>
    private static string Q7Reason(string a, string b)
    {
        var k = CellPair(a, b);

        // Always-hard cells: two Fathers, or Child + Y-Father.
        if (k == ("F_a", "F_b"))
            return "F_a+F_b (Y-Father × Z-Father): both X and Y demanded, unroutable";
        if (k == ("C", "F_a"))
            return "C+F_a (Child × Y-Father): irreducible same-qubit X/Y demand";

        // Split cell A: Y-Father (XY/YX) + Mother (XX/YY/ZZ).
        if (k == ("F_a", "M"))
            return SameSiteXY(a, b)
                ? "F_a+M lit Mother: X and Y forced onto one qubit"
                : "F_a+M dark Mother (ZZ): dephasing-aligned, no X/Y conflict";

        // Split cell B: Z-Father (XZ/ZX) + Child (YZ/ZY).
        if (k == ("C", "F_b"))
            return ContinuousPairs.Contains(PairKey(a, b))
                ? "C+F_b same-site X&Y over shared dark Z: local continuous-rotation rescue"
                : "C+F_b lit X,Y on different sites: no Q closes it";

        // All other cells (M+M, M+C, M+F_b, C+C, F_a+F_a, ...): never hard.
        return "no irreducible same-qubit X/Y conflict";
    }

    // -----------------------------------------------------------------------------------------------
    // Parity-break cross-check (independent of the routing rule).
    // -----------------------------------------------------------------------------------------------

    /// <summary>True iff H = Σ_bonds (a + b) on the open chain bonds (k, k+1) fails to commute with the
    /// global spin-flip X^⊗N, i.e. ‖[H, X^⊗N]‖ &gt; 1e-9 (the X-parity Z₂ symmetry is broken).</summary>
    private static bool ParityBreaks(int n, string a, string b)
    {
        var letterTerms = new[]
        {
            (PauliLetterExtensions.FromSymbol(a[0]), PauliLetterExtensions.FromSymbol(a[1])),
            (PauliLetterExtensions.FromSymbol(b[0]), PauliLetterExtensions.FromSymbol(b[1])),
        };

        int d = 1 << n;
        var h = Matrix<Complex>.Build.Dense(d, d);
        for (int k = 0; k < n - 1; k++)
        {
            foreach (var (la, lb) in letterTerms)
                h += PauliString.SiteOp(n, k, la) * PauliString.SiteOp(n, k + 1, lb);
        }

        ComplexMatrix xn = PauliString.Build(Enumerable.Repeat(PauliLetter.X, n).ToArray());

        var comm = h * xn - xn * h;
        return comm.FrobeniusNorm() > ParityTolerance;
    }
}
