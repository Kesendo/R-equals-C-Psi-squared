using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F173, the signed balance of the flip masks (Tier 1 derived, <c>docs/proofs/PROOF_SIGNED_MASK_BALANCE.md</c>).
///
/// <para>What the page proves. Under Z dephasing on every site (any positive rates), for H a real combination of
/// Pauli strings in which no two strings share a nonzero X/Y mask and the odd-weight diagonal part vanishes, each edge
/// x – x ⊕ m of Theorem 2's graph carries the ratio (−1)^s of its one string, s = (#Y + #Z) mod 2 (the bit_b, the
/// Π²_Z parity), the same at every x; the graph is the Cayley graph of the placed masks, and the palindrome holds
/// exactly when that signed graph is balanced: no GF(2) relation among the masks uses an odd number of negative
/// strings (Theorem A; where it fails no component is good). Two templates t₁, t₂ of widths w₁, w₂ on the open chain,
/// masks m₁, m₂ of different shapes, g = gcd(m₁, m₂), qᵢ = mᵢ/g: every relation is r·(q₂, q₁), its negative count
/// r(1)·c with c = s₁·q₂(1) + s₂·q₁(1) mod 2, so c = 0 pairs at every N and c = 1 pairs exactly below
/// N* = max(w₁ + deg q₂, w₂ + deg q₁) (Theorem B; both signs negative is F115's valuation criterion, N* its window
/// floor). One template on the ring of N ≥ 2w − 1 sites breaks exactly when s = 1 and v_(1+x)(m) ≥ 2^(v₂(N)), the
/// valuation of x^N + 1 (Theorem C). Any number of templates of different shapes on the ring of N ≥ 2·max w − 1
/// sites: the sign count of a relation is read at x = 1, which sees x^N + 1 only through its factor (1 + x)^T,
/// T = 2^(v₂(N)); with v₀ the least valuation the ring breaks exactly when T ≤ v₀ and some sign is odd, or T &gt; v₀
/// and two templates of valuation v₀ differ in sign or one of greater valuation has an odd sign, so it reads N only
/// through whether T exceeds v₀ (Theorem D).</para>
///
/// <para>Typed parent: <see cref="PalindromeComplementConnectionClaim"/> (Theorem 2 of the complement connection),
/// which this claim reads where every ratio is a sign.</para>
///
/// <para>Live witness: <c>inspect --root signedbalance</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/SignedMaskBalanceWitness.cs</c>), which holds the claim's readings
/// (<see cref="Balanced"/>, <see cref="Chain"/>, <see cref="RingBreaks(IReadOnlyCollection{string}, int)"/>), pure GF(2) on
/// masks and signs, against <see cref="ComplementConnectionGraph"/>, Theorem 2 walked on the bitstrings with every entry a
/// Gaussian integer, and Theorem D also against the balance on rings past the bitstrings' reach.</para></summary>
public sealed class SignedMaskBalanceClaim : Claim
{
    public PalindromeComplementConnectionClaim ComplementConnection { get; }

    public SignedMaskBalanceClaim(PalindromeComplementConnectionClaim complementConnection)
        : base(
            "Under Z dephasing on every site, where every flip set carries one Pauli string and no diagonal string of odd " +
            "weight remains, the palindrome is the balance of the signed Cayley graph of the masks (each string signed by " +
            "(−1)^(#Y + #Z)); two templates on the chain break exactly from N* = max(w₁ + deg q₂, w₂ + deg q₁) when " +
            "s₁·q₂(1) + s₂·q₁(1) is odd, and never otherwise; on the ring the verdict reads N only through T = 2^(v₂(N)) and " +
            "each mask only through its (1 + x)-valuation: one template breaks exactly when its sign is odd and v(m) ≥ T, " +
            "and with v₀ the least valuation the ring breaks exactly when T ≤ v₀ and some sign is odd, or T > v₀ and two " +
            "templates of valuation v₀ differ in sign or one of greater valuation has an odd sign",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_SIGNED_MASK_BALANCE.md + docs/ANALYTICAL_FORMULAS.md (F173, the index entry)")
    {
        ComplementConnection = complementConnection ?? throw new ArgumentNullException(nameof(complementConnection));
    }

    /// <summary>The X/Y mask of a template, bit j for its letter j (identity letters inside the window kept).</summary>
    public static ulong Mask(string template)
    {
        ulong m = 0;
        for (int i = 0; i < template.Length; i++) if (template[i] is 'X' or 'Y') m |= 1UL << i;
        return m;
    }

    /// <summary>The sign of a string, (#Y + #Z) mod 2: X^⊗N·P·X^⊗N = (−1)^s·P.</summary>
    public static int Sign(string template) => template.Count(c => c is 'Y' or 'Z') % 2;

    /// <summary>Theorem A at one N: no relation among the placed masks with an odd number of negative strings.
    /// Elimination on (mask | sign); a row that reduces to the zero mask with sign 1 is such a relation. It decides the
    /// palindrome only under the theorem's assumption, one string per flip set: <see cref="OnChain"/> and
    /// <see cref="OnRing"/> refuse templates whose masks are translates of each other.</summary>
    public static bool Balanced(IEnumerable<(ulong Mask, int Sign)> placed)
    {
        var basis = new Dictionary<int, (ulong Mask, int Sign)>();
        foreach (var (mask0, sign0) in placed)
        {
            ulong m = mask0;
            int s = sign0 & 1;
            while (m != 0)
            {
                int top = Deg(m);
                if (!basis.TryGetValue(top, out var b)) { basis[top] = (m, s); break; }
                m ^= b.Mask;
                s ^= b.Sign;
            }
            if (m == 0 && s == 1) return false;
        }
        return true;
    }

    /// <summary>The placed masks of templates on the open chain of n sites, at every start 0 … n − w.</summary>
    public static IEnumerable<(ulong Mask, int Sign)> OnChain(IReadOnlyCollection<string> templates, int n)
    {
        DistinctShapes(templates);
        return templates.SelectMany(t => Enumerable.Range(0, Math.Max(0, n - t.Length + 1)).Select(p => (Mask(t) << p, Sign(t))));
    }

    // one string per flip set: no template without a flip, no two templates whose masks are translates
    static void DistinctShapes(IReadOnlyCollection<string> templates)
    {
        if (templates.Count == 0) throw new ArgumentException("no templates");
        var shapes = templates.Select(Mask).ToList();
        if (shapes.Any(m => m == 0)) throw new ArgumentException("every template must flip something");
        if (shapes.Select(Shape).Distinct().Count() < shapes.Count)
            throw new ArgumentException("two masks are translates: a flip set would carry two strings");
    }

    /// <summary>The placed masks of templates on the ring of n sites, every start, wrapped; n ≥ 2w − 1 for the widest
    /// template, below which placements can coincide and a flip set carry two strings, and n ≤ 64, the width of a
    /// mask.</summary>
    public static IEnumerable<(ulong Mask, int Sign)> OnRing(IReadOnlyCollection<string> templates, int n)
    {
        DistinctShapes(templates);
        int w = templates.Max(t => t.Length);
        if (n < 2 * w - 1) throw new ArgumentOutOfRangeException(nameof(n), n, $"the ring needs at least {2 * w - 1} sites");
        if (n > 64) throw new ArgumentOutOfRangeException(nameof(n), n, "a mask holds at most 64 sites");
        return templates.SelectMany(t => Enumerable.Range(0, n).Select(p => (Rotate(Mask(t), p, n), Sign(t))));
    }

    /// <summary>Theorem B: two templates on the open chain, whether the chain breaks for large N and from which N on
    /// (null when it never does).</summary>
    public static (bool Breaks, int? From) Chain(string t1, string t2)
    {
        ulong m1 = Mask(t1), m2 = Mask(t2);
        if (m1 == 0 || m2 == 0) throw new ArgumentException("both templates must flip something");
        if (Shape(m1) == Shape(m2)) throw new ArgumentException("the two masks are translates: a flip set would carry two strings");
        ulong g = Gcd(m1, m2), q1 = Div(m1, g), q2 = Div(m2, g);
        int c = (Sign(t1) * AtOne(q2) + Sign(t2) * AtOne(q1)) % 2;
        if (c == 0) return (false, null);
        return (true, Math.Max(t1.Length + Deg(q2), t2.Length + Deg(q1)));
    }

    /// <summary>Theorem C: one template on the ring of n ≥ 2w − 1 sites breaks exactly when its sign is odd and the
    /// (1 + x)-valuation of its mask reaches 2^(v₂(n)).</summary>
    public static bool RingBreaks(string template, int n) => RingBreaks(new[] { template }, n);

    /// <summary>Theorem D: templates of different shapes on the ring of n ≥ 2·max w − 1 sites, read through
    /// T = 2^(v₂(n)) and the (1 + x)-valuations vᵢ of the masks alone, v₀ the least: if T ≤ v₀ the ring breaks exactly
    /// when some sign is odd, otherwise exactly when two templates of valuation v₀ differ in sign or one of greater
    /// valuation has an odd sign. One template is Theorem C.</summary>
    public static bool RingBreaks(IReadOnlyCollection<string> templates, int n)
    {
        DistinctShapes(templates);
        int w = templates.Max(t => t.Length);
        if (n < 2 * w - 1) throw new ArgumentOutOfRangeException(nameof(n), n, $"the ring needs at least {2 * w - 1} sites");
        int big = 1 << BitOperations.TrailingZeroCount((uint)n);
        var read = templates.Select(t => (V: Valuation(Mask(t)), S: Sign(t))).ToList();
        int v0 = read.Min(r => r.V);
        if (v0 >= big) return read.Any(r => r.S == 1);
        return read.Where(r => r.V == v0).Select(r => r.S).Distinct().Count() > 1 || read.Any(r => r.V > v0 && r.S == 1);
    }

    /// <summary>v_(1+x)(m): how often 1 + x divides the mask read as a polynomial over GF(2).</summary>
    public static int Valuation(ulong m)
    {
        int v = 0;
        while (m != 0 && (BitOperations.PopCount(m) & 1) == 0) { m = Div(m, 0b11); v++; }
        return v;
    }

    static ulong Rotate(ulong m, int p, int n)
    {
        ulong r = 0;
        for (int j = 0; j < 64; j++) if (((m >> j) & 1) == 1) r |= 1UL << ((j + p) % n);
        return r;
    }

    static ulong Shape(ulong m) => m >> BitOperations.TrailingZeroCount(m);
    static int Deg(ulong p) => 63 - BitOperations.LeadingZeroCount(p);
    static int AtOne(ulong p) => BitOperations.PopCount(p) & 1;

    static ulong Gcd(ulong a, ulong b)
    {
        while (b != 0) (a, b) = (b, Mod(a, b));
        return a;
    }

    static ulong Mod(ulong a, ulong b)
    {
        while (a != 0 && Deg(a) >= Deg(b)) a ^= b << (Deg(a) - Deg(b));
        return a;
    }

    static ulong Div(ulong a, ulong b)
    {
        ulong q = 0;
        while (a != 0 && Deg(a) >= Deg(b)) { int s = Deg(a) - Deg(b); q |= 1UL << s; a ^= b << s; }
        if (a != 0) throw new InvalidOperationException("not divisible");
        return q;
    }

    public override string DisplayName =>
        "F173: the signed balance of the flip masks";

    public override string Summary =>
        "one string per flip set under Z dephasing: each edge a sign (−1)^(#Y + #Z), the palindrome the balance of the masks' " +
        "signed Cayley graph; two templates break exactly from N* when s₁·q₂(1) + s₂·q₁(1) is odd (F115 its all-negative case), " +
        $"one template on the ring exactly when s = 1 and v_(1+x)(m) ≥ 2^(v₂(N)), and any number on the ring read through " +
        $"whether 2^(v₂(N)) exceeds the least valuation ({Tier.Label()}; live witness: inspect --root signedbalance)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Theorem A, the balance",
                summary: "one string per flip set: the ratio H_x̄,x̄⊕m / H_x,x⊕m is (−1)^s at every x; a closed walk's holonomy is the parity " +
                         "of the negative strings in its relation, and every relation is walked from every vertex, so the palindrome " +
                         "holds exactly when every relation is balanced, and where one is not no component is good");
            yield return new InspectableNode("Theorem B, two templates on the chain",
                summary: "relations r·(q₂, q₁), negative count r(1)·c; c = 0 pairs at every N, c = 1 pairs exactly below " +
                         "N* = max(w₁ + deg q₂, w₂ + deg q₁); both signs negative: c = 1 iff the (1 + x)-valuations differ (F115)");
            yield return new InspectableNode("Theorem C, one template on the ring",
                summary: "relations u·m ≡ 0 mod x^N + 1; unbalanced exactly when s = 1 and v_(1+x)(m) ≥ 2^(v₂(N)): XY breaks at odd N, " +
                         "XIY unless 4 | N, XIIIY unless 8 | N");
            yield return new InspectableNode("Theorem D, any number of templates on the ring",
                summary: "x^N + 1 = (1 + x)^T·F with T = 2^(v₂(N)) and F(1) = 1; a relation reduced modulo (1 + x)^T keeps its sign count, " +
                         "and one there times F is one on the ring, so the verdict is read in GF(2)[y]/(y^T), y = 1 + x, where each mask is " +
                         "y^(vᵢ) times a unit, v₀ the least; if T ≤ v₀ every mask is 0 there and the ring breaks on any odd sign; otherwise, i₀ " +
                         "of valuation v₀, every relation is a combination of the relations that use template j once and i₀ by " +
                         "y^(v_j − v₀) times a unit, of sign count s_j + s_i₀·[v_j = v₀], and one at i₀ alone, of sign count 0");
            yield return new InspectableNode("the two pure cases",
                summary: "every sign negative: balanced means bipartite (F103 §7.4, F115 statement 3); every sign positive: X^⊗N " +
                         "commutes with H, a colouring (F108 Part 1); a balanced signing is the lit string X^⊗N·Z^φ commuting with H");
        }
    }
}
