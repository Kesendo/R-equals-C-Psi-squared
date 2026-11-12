using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
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
/// floor). Any number of templates on the chain: the relations are free of rank k − 1, and over a basis reduced for the
/// widths (weak Popov form) the chain breaks exactly from the least shifted degree of an unbalanced basis vector
/// (Theorem B′; c is the girth ladder's syzygy polynomial read at 1). One template on the ring of N ≥ 2w − 1 sites breaks exactly when s = 1 and v_(1+x)(m) ≥ 2^(v₂(N)), the
/// valuation of x^N + 1 (Theorem C). Any number of templates of different shapes on the ring of N ≥ 2·max w − 1
/// sites: the sign count of a relation is read at x = 1, which sees x^N + 1 only through its factor (1 + x)^T,
/// T = 2^(v₂(N)); with v₀ the least valuation the ring breaks exactly when T ≤ v₀ and some sign is odd, or T &gt; v₀
/// and two templates of valuation v₀ differ in sign or one of greater valuation has an odd sign, so it reads N only
/// through whether T exceeds v₀ (Theorem D). Where a flip set carries several strings whose ratio is a constant times a
/// character, ω_m·(−1)^⟨R_m, x⟩, the palindrome holds exactly when ⟨R_a, b⟩ = ⟨R_b, a⟩ for every two flip sets, every
/// relation has a vanishing read sum and hol = 1 on a basis of the relations (Theorem E); one flip shape on the ring is
/// read through gcd(f, x^N + 1) (Corollary E2), and two on the chain through the read sum D = q₂R₁ + q₁R₂ and the
/// holonomy of the relation r = 1 (Corollary E3), any number on the chain below the least shifted degree of a bad vector of
/// a reduced basis of the relations, its read sum nonzero or its holonomy −1 (Corollary E3′, <see cref="ChainShapes"/>); whether the long chain breaks is Theorem D's
/// verdict where T exceeds the least valuation, for any number of templates. With one string per flip set on the L × M torus the verdict is read on the 2-part ℤ_T₁ × ℤ_T₂ of the
/// translations, T₁ = 2^(v₂(L)), T₂ = 2^(v₂(M)), every mask folded there: the torus pairs exactly when a colouring of
/// period (T₁, T₂) exists, the signs a linear function of the folded masks modulo the differences of their translates
/// (Theorem F); the pairing tori form a staircase in (v₂(L), v₂(M)) with finitely many corners, one template's read off the
/// exponents of its mask expanded at x = y = 1, coarsened to bit lengths, and two templates' by divisibility in A = GF(2)[ℤ_T₁ × ℤ_T₂]; the ring's
/// one number per mask is the case of a cyclic 2-part. Two templates on the open L × M patch, their masks read in GF(2)[x, y]:
/// every relation is r·(q₂, q₁) and the degrees in x and in y add separately, so with c = s₁·q₂(1, 1) + s₂·q₁(1, 1) odd the
/// patch breaks exactly on the quadrant L ≥ L*, M ≥ M*, and never with c even (Theorem G, <see cref="Patch"/>). At strong
/// dephasing, H scaled to J·H, the two edges of the spectrum are −J²·spec(𝓛) and −2σ + J²·spec(𝓜) to O(J³), two Laplacians
/// whose kernels are Theorem 2's two counts, the far one magnetic, carrying Theorem 2's ratios as phases, with a potential where the
/// moduli differ (<see cref="StrongDephasingEdges"/>);
/// with one string per flip set the far edge misses −2σ by J² times twice the distance of the signs from the code of
/// balanced signings, each string weighted by c²/γ(m) (Theorem H, <see cref="StrongDephasingBreak"/>).</para>
///
/// <para>Typed parent: <see cref="PalindromeComplementConnectionClaim"/> (Theorem 2 of the complement connection),
/// which this claim reads where every ratio is a sign.</para>
///
/// <para>Live witness: <c>inspect --root signedbalance</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/SignedMaskBalanceWitness.cs</c>), which holds the claim's readings
/// (<see cref="Balanced"/>, <see cref="Chain(string, string)"/>, <see cref="Chain(IReadOnlyList{string}, bool)"/>,
/// <see cref="RingBreaks(IReadOnlyCollection{string}, int)"/>), pure GF(2) on
/// masks and signs, against <see cref="ComplementConnectionGraph"/>, Theorem 2 walked on the bitstrings with every entry a
/// Gaussian integer, and Theorem D also against the balance on rings past the bitstrings' reach; Theorem F's readings
/// (<see cref="TorusBreaks"/>, <see cref="TorusBreaksByDifferences"/>, <see cref="BreaksAtTwoParts"/>, <see cref="Corners"/>,
/// <see cref="SingleTemplateCorners"/>, <see cref="TwoTemplatesBreakByDivisibility"/>, <see cref="CornerBound"/>, <see cref="SufficientTorus"/>,
/// <see cref="OnLattice"/>, <see cref="LatticeBreaks"/>) and Theorem G's (<see cref="OnPatch"/>, <see cref="Patch"/>, <see cref="PatchBreaks"/>) have their own battery, <c>inspect --root torusbalance</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/TorusMaskBalanceWitness.cs</c>); Theorem H's reader and the lemma's two
/// generators are read in the first battery, at random rates.</para></summary>
public sealed class SignedMaskBalanceClaim : Claim
{
    public PalindromeComplementConnectionClaim ComplementConnection { get; }

    public SignedMaskBalanceClaim(PalindromeComplementConnectionClaim complementConnection)
        : base(
            "Under Z dephasing on every site, where every flip set carries one Pauli string and no diagonal string of odd " +
            "weight remains, the palindrome is the balance of the signed Cayley graph of the masks (each string signed by " +
            "(−1)^(#Y + #Z)); two templates on the chain break exactly from N* = max(w₁ + deg q₂, w₂ + deg q₁) when " +
            "s₁·q₂(1) + s₂·q₁(1) is odd, and never otherwise, and any number of templates on the chain break exactly from the least " +
            "shifted degree of an unbalanced vector in a basis of the relations reduced for the widths, and never if there is none; on the ring the verdict reads N only through T = 2^(v₂(N)) and " +
            "each mask only through its (1 + x)-valuation: one template breaks exactly when its sign is odd and v(m) ≥ T, " +
            "and with v₀ the least valuation the ring breaks exactly when T ≤ v₀ and some sign is odd, or T > v₀ and two " +
            "templates of valuation v₀ differ in sign or one of greater valuation has an odd sign; where a flip set carries several " +
            "strings and its ratio is ω·(−1)^⟨R, x⟩, the palindrome holds exactly when ⟨R_a, b⟩ = ⟨R_b, a⟩ for every two flip sets " +
            "and every relation u has ⊕ R_m = 0 and hol(u) = 1 (Theorem E), and one flip shape f on the ring of N ≥ 2w − 1 sites, " +
            "that condition holding, pairs exactly when h = gcd(f, x^N + 1) is 1, or divides R with 1 + x not dividing h or " +
            "hol((x^N + 1)/h) = 1, and two flip shapes on the chain, that condition holding, pair exactly below N* or where " +
            "D = q₂R₁ + q₁R₂ vanishes and the relation (q₂, q₁) has hol = 1, any number of flip shapes on the chain below the least shifted degree " +
            "of a basis vector of the relations with a nonzero read sum or with hol = −1; with one string per flip set on the L × M torus the verdict " +
            "is read on the 2-part ℤ_T₁ × ℤ_T₂ of the translations (T₁ = 2^(v₂(L)), T₂ = 2^(v₂(M))), where it pairs exactly when a colouring of period (T₁, T₂) " +
            "exists, the signs a linear function of the folded masks modulo the differences of their translates (Theorem F), and the pairing " +
            "tori form a staircase in (v₂(L), v₂(M)) with finitely many corners, bounded by the colength of the masks' ideal at x = y = 1 " +
            "where the masks share no factor there; two templates on the open L × M patch break exactly on the quadrant " +
            "L ≥ max(w_x,1 + deg_x q₂, w_x,2 + deg_x q₁), M ≥ max(w_y,1 + deg_y q₂, w_y,2 + deg_y q₁) when s₁·q₂(1, 1) + s₂·q₁(1, 1) is odd, " +
            "the gcd taken in GF(2)[x, y], and never otherwise (Theorem G); at strong dephasing, H scaled to J·H, the far edge's eigenvalue " +
            "nearest −2σ lies at 2J²·min over φ of Σ c_P²/γ(m_P) over the strings with ⟨φ, m_P⟩ ≠ s(P), plus O(J³), from it, twice the " +
            "rate-weighted distance of the signs from the code of balanced signings, the case of one string per flip set of a lemma: " +
            "the two edges are −J²·spec(𝓛) and −2σ + J²·spec(𝓜), two Laplacians whose kernels are Theorem 2's two counts (Theorem H)",
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
    public static int Sign(string template) => template.Count(c => c is 'Y' or 'Z') % 2;   // a torus template's '/' counts nothing

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
        if (templates.Any(t => t.Length > 64)) throw new ArgumentException("a template holds at most 64 letters");
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

    /// <summary>Theorem B′: any number of templates of different shapes on the open chain. A relation is a vector u of
    /// polynomials with Σ uᵢ·mᵢ = 0, fitting at N exactly when its shifted degree max(deg uᵢ + wᵢ) is at most N, its
    /// negative count λ(u) = Σ sᵢ·uᵢ(1). Over a width-shifted reduced basis b₁ … b_(k−1) of the relations (weak Popov
    /// form) the shifted degree of Σ rⱼ·bⱼ is max(deg rⱼ + sdeg bⱼ), so the chain breaks exactly from
    /// N* = min{sdeg bⱼ : λ(bⱼ) = 1} and never if every λ(bⱼ) = 0. With <paramref name="reduced"/> false the same reading
    /// is taken off the unreduced basis, the control that must fail (it can read N* too large). Two templates give
    /// Theorem B's N*.</summary>
    public static (bool Breaks, int? From) Chain(IReadOnlyList<string> templates, bool reduced = true)
    {
        var widths = templates.Select(t => t.Length).ToArray();
        var signs = templates.Select(Sign).ToArray();
        int? from = null;
        foreach (var b in ChainRelationBasis(templates, reduced))
        {
            int lambda = 0;
            for (int i = 0; i < b.Length; i++) lambda ^= signs[i] & AtOne(b[i]);
            if (lambda == 1) from = Math.Min(from ?? int.MaxValue, ShiftedDegree(b, widths).Degree);
        }
        return (from is not null, from);
    }

    /// <summary>A basis of the relations among the masks of the templates over GF(2)[x], k − 1 vectors (the masks are
    /// nonzero, so the relations are free of rank k − 1): the columns of the unimodular matrix that the Euclidean
    /// reduction of (m₁ … m_k) to its gcd leaves beside the gcd's column; with <paramref name="reduced"/> brought to weak
    /// Popov form for the shift by the widths (no two vectors with the same pivot, the last position attaining the
    /// shifted degree), which makes the basis reduced and the shifted degree predictable.</summary>
    public static IReadOnlyList<ulong[]> ChainRelationBasis(IReadOnlyList<string> templates, bool reduced = true)
    {
        DistinctShapes(templates);
        var widths = templates.Select(t => t.Length).ToArray();
        var basis = KernelBasis(templates.Select(Mask).ToArray());
        return reduced ? ShiftedReduce(basis, widths) : basis;
    }

    // a basis of {u : Σ uᵢ·vᵢ = 0} over GF(2)[x] for a vector v with a nonzero entry: the Euclidean reduction of v to its
    // gcd with the column operations recorded; the columns beside the gcd's span the kernel (the matrix is unimodular)
    static List<ulong[]> KernelBasis(ulong[] v0)
    {
        int k = v0.Length;
        var v = (ulong[])v0.Clone();
        if (v.All(x => x == 0)) throw new ArgumentException("the vector needs a nonzero entry");
        var columns = Enumerable.Range(0, k).Select(i => Enumerable.Range(0, k).Select(j => i == j ? 1UL : 0UL).ToArray()).ToArray();
        while (v.Count(x => x != 0) > 1)
        {
            int j = Enumerable.Range(0, k).Where(i => v[i] != 0).OrderBy(i => Deg(v[i])).First();
            for (int i = 0; i < k; i++)
            {
                if (i == j || v[i] == 0) continue;
                var (q, r) = DivMod(v[i], v[j]);
                v[i] = r;
                for (int row = 0; row < k; row++) columns[i][row] ^= MulChecked(q, columns[j][row]);
            }
        }
        int gcdAt = Array.FindIndex(v, x => x != 0);
        return Enumerable.Range(0, k).Where(i => i != gcdAt).Select(i => columns[i]).ToList();
    }

    // weak Popov form for the shift by the widths: where two vectors share the last position attaining the shifted degree,
    // subtract from the one of larger shifted degree x^d times the other (unimodular, the shifted degree or the pivot drops)
    static List<ulong[]> ShiftedReduce(List<ulong[]> basis, int[] widths)
    {
        basis = basis.Select(b => (ulong[])b.Clone()).ToList();
        while (true)
        {
            var read = basis.Select(b => ShiftedDegree(b, widths)).ToList();
            (int A, int C) clash = (from a in Enumerable.Range(0, basis.Count) from c in Enumerable.Range(0, basis.Count)
                                    where a < c && read[a].Pivot == read[c].Pivot select (a, c)).FirstOrDefault((-1, -1));
            if (clash.A < 0) return basis;
            var (hi, lo) = read[clash.A].Degree >= read[clash.C].Degree ? (clash.A, clash.C) : (clash.C, clash.A);
            int shift = read[hi].Degree - read[lo].Degree;
            basis[hi] = basis[hi].Zip(basis[lo], (x, y) => x ^ ShiftLeft(y, shift)).ToArray();
        }
    }

    /// <summary>max(deg uᵢ + wᵢ) over the nonzero entries, and the last position attaining it.</summary>
    static (int Degree, int Pivot) ShiftedDegree(ulong[] u, int[] widths)
    {
        int best = int.MinValue, pivot = -1;
        for (int i = 0; i < u.Length; i++)
            if (u[i] != 0 && Deg(u[i]) + widths[i] >= best) { best = Deg(u[i]) + widths[i]; pivot = i; }
        return (best, pivot);
    }

    static ulong MulChecked(ulong a, ulong b)
    {
        if (a == 0 || b == 0) return 0;
        if (Deg(a) + Deg(b) > 63) throw new OverflowException("a relation's entry passes degree 63");
        ulong r = 0;
        for (int j = 0; j <= Deg(b); j++) if (((b >> j) & 1) == 1) r ^= a << j;
        return r;
    }

    static ulong ShiftLeft(ulong a, int s)
    {
        if (a != 0 && Deg(a) + s > 63) throw new OverflowException("a relation's entry passes degree 63");
        return a << s;
    }

    static (ulong Q, ulong R) DivMod(ulong a, ulong b)
    {
        ulong q = 0;
        while (a != 0 && Deg(a) >= Deg(b)) { int s = Deg(a) - Deg(b); q |= 1UL << s; a ^= b << s; }
        return (q, a);
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

    // ---------------------------------------------------------------- Theorem F, the torus

    /// <summary>The letters of a template on the torus, written row by row with '/' between the rows: "XX/XY" carries X
    /// at (0, 0) and (1, 0), X at (0, 1) and Y at (1, 1); identity letters are dropped, a shorter row is read padded with I.</summary>
    public static IReadOnlyList<(int I, int J, char Letter)> TorusLetters(string template)
    {
        var rows = template.Split('/');
        var letters = new List<(int, int, char)>();
        for (int j = 0; j < rows.Length; j++)
            for (int i = 0; i < rows[j].Length; i++)
            {
                char c = rows[j][i];
                if (c is not ('I' or 'X' or 'Y' or 'Z')) throw new ArgumentException($"the letter {c} in {template}");
                if (c != 'I') letters.Add((i, j, c));
            }
        return letters;
    }

    /// <summary>The placed masks of templates on the L × M torus, site (i, j) at bit i + L·j (L·M ≤ 64), every template
    /// at every translate, wrapped, with its sign. Refuses a template that flips nothing, a template two of whose letters
    /// land on one site, and placements that coincide (two templates whose masks are translates, or a torus too small for
    /// a template), where a flip set would carry two strings.</summary>
    public static IReadOnlyList<(ulong Mask, int Sign)> OnTorus(IReadOnlyCollection<string> templates, int l, int m)
    {
        if (templates.Count == 0) throw new ArgumentException("no templates");
        if (l < 1 || m < 1 || (long)l * m > 64) throw new ArgumentOutOfRangeException(nameof(l), $"{l} × {m}: one to 64 sites");
        var placed = new List<(ulong, int)>();
        var seen = new HashSet<ulong>();
        foreach (var t in templates)
        {
            var letters = TorusLetters(t);
            if (!letters.Any(x => x.Letter is 'X' or 'Y')) throw new ArgumentException($"{t} flips nothing");
            if (letters.Select(x => (x.I % l, x.J % m)).Distinct().Count() < letters.Count)
                throw new ArgumentException($"{t} does not fit on the {l} × {m} torus");
            int sign = letters.Count(x => x.Letter is 'Y' or 'Z') % 2;
            for (int b = 0; b < m; b++)
                for (int a = 0; a < l; a++)
                {
                    ulong mask = 0;
                    foreach (var (i, j, c) in letters)
                        if (c is 'X' or 'Y') mask |= 1UL << ((i + a) % l + l * ((j + b) % m));
                    if (!seen.Add(mask)) throw new ArgumentException($"two placements on the {l} × {m} torus share a flip set: it would carry two strings");
                    placed.Add((mask, sign));
                }
        }
        return placed;
    }

    /// <summary>A template's mask folded onto the T₁ × T₂ torus and moved by (a, b): each flipped letter at (i, j) toggles
    /// the cell ((i + a) mod T₁, (j + b) mod T₂), so a cell hit an even number of times drops out.</summary>
    public static ulong FoldedMask(string template, int t1, int t2, int a = 0, int b = 0)
    {
        if (t1 < 1 || t2 < 1 || (long)t1 * t2 > 64) throw new ArgumentOutOfRangeException(nameof(t1), $"{t1} × {t2}: one to 64 cells");
        ulong mask = 0;
        foreach (var (i, j, c) in TorusLetters(template))
            if (c is 'X' or 'Y') mask ^= 1UL << ((i + a) % t1 + t1 * ((j + b) % t2));
        return mask;
    }

    /// <summary>Theorem F: templates on the L × M torus, placed as <see cref="OnTorus"/> places them, read on the 2-part of
    /// the translations, T₁ × T₂ with T₁ = 2^(v₂(L)), T₂ = 2^(v₂(M)): every mask folded there (<see cref="FoldedMask"/>) at
    /// every translate of the small torus, with its template's sign. The torus breaks exactly when a GF(2) relation among
    /// those folded masks uses an odd number of negative strings, that is when no colouring of period (T₁, T₂) exists.</summary>
    public static bool TorusBreaks(IReadOnlyCollection<string> templates, int l, int m)
    {
        OnTorus(templates, l, m);   // the setting: every placement its own flip set
        return BreaksAtTwoParts(templates, TwoPart(l), TwoPart(m));
    }

    /// <summary>Theorem F on the small torus alone: the verdict that every L × M torus with 2-parts t₁ × t₂ (powers of 2,
    /// t₁·t₂ ≤ 64) shares, wherever its placements do not coincide.</summary>
    public static bool BreaksAtTwoParts(IReadOnlyCollection<string> templates, int t1, int t2)
    {
        if (t1 < 1 || t2 < 1 || (t1 & (t1 - 1)) != 0 || (t2 & (t2 - 1)) != 0)
            throw new ArgumentOutOfRangeException(nameof(t1), $"{t1} × {t2}: the 2-parts are powers of 2");
        var rows = new List<(ulong, int)>();
        foreach (var t in templates)
        {
            int sign = Sign(t);
            for (int b = 0; b < t2; b++)
                for (int a = 0; a < t1; a++) rows.Add((FoldedMask(t, t1, t2, a, b), sign));
        }
        return !Balanced(rows);
    }

    /// <summary>The staircase of Theorem F: a colouring of period (T₁, T₂) is one of period (2T₁, T₂) and (T₁, 2T₂), so
    /// the exponents (e, f) at which the torus pairs (T₁ = 2^e, T₂ = 2^f) are closed upwards, and the torus pairs exactly
    /// when 2^eⱼ divides L and 2^fⱼ divides M for one of finitely many corners (eⱼ, fⱼ). Returned are the corners with
    /// 2^(e + f) ≤ 64, read on the small tori; every point below such a corner is read, so each is a corner of the whole
    /// staircase, and a corner past 64 cells is not seen.</summary>
    public static IReadOnlyList<(int E, int F)> Corners(IReadOnlyCollection<string> templates)
    {
        var pairing = new List<(int E, int F)>();
        for (int e = 0; e <= 6; e++)
            for (int f = 0; e + f <= 6; f++)
                if (!BreaksAtTwoParts(templates, 1 << e, 1 << f)) pairing.Add((e, f));
        return MinimalPoints(pairing);
    }

    /// <summary>The corners of one template read off its mask alone: written in u = 1 + x and v = 1 + y (Lucas: x^i holds
    /// u^a exactly when the bits of a lie among those of i), the folded mask is nonzero exactly when one of its monomials
    /// u^a·v^b has a &lt; T₁ and b &lt; T₂, and one nonzero folded mask carries either sign. So an odd sign pairs exactly
    /// above the minimal points (bit length of a, bit length of b), the least exponents of the mask coarsened to bit
    /// lengths; an even sign pairs on every torus, the corner (0, 0).</summary>
    public static IReadOnlyList<(int E, int F)> SingleTemplateCorners(string template)
    {
        if (Sign(template) == 0) return new[] { (0, 0) };
        return MinimalPoints(UvMonomials(template).Select(p => (BitLength(p.A), BitLength(p.B))).Distinct().ToList());
    }

    /// <summary>Two templates on the small torus t₁ × t₂, read by divisibility in the local ring A = GF(2)[ℤ_t₁ × ℤ_t₂]:
    /// an odd mask beside an even one breaks exactly when the even folded mask divides the odd one, two odd masks break
    /// exactly when one divides the other by an element of 𝔪, and two even masks never break. A multiple of m̄ is a sum of
    /// its translates, a multiple by an element of 𝔪 a sum of differences g·m̄ + m̄.</summary>
    public static bool TwoTemplatesBreakByDivisibility(string first, string second, int t1, int t2)
    {
        if (t1 < 1 || t2 < 1 || (t1 & (t1 - 1)) != 0 || (t2 & (t2 - 1)) != 0)
            throw new ArgumentOutOfRangeException(nameof(t1), $"{t1} × {t2}: A is local only on a 2-group");
        ulong f1 = FoldedMask(first, t1, t2), f2 = FoldedMask(second, t1, t2);
        int s1 = Sign(first), s2 = Sign(second);
        bool Divides(string d, ulong of, bool byM)
        {
            var rows = new List<(ulong, int)>();
            ulong at0 = FoldedMask(d, t1, t2);
            for (int b = 0; b < t2; b++)
                for (int a = 0; a < t1; a++) rows.Add((FoldedMask(d, t1, t2, a, b) ^ (byM ? at0 : 0), 0));
            rows.Add((of, 1));
            return !Balanced(rows);   // the last row reduces to the zero mask with sign 1 exactly when it lies in the span
        }
        return (s1, s2) switch
        {
            (1, 0) => Divides(second, f1, false),
            (0, 1) => Divides(first, f2, false),
            (1, 1) => Divides(second, f1, true) || Divides(first, f2, true),
            _ => false,
        };
    }

    /// <summary>A torus on which no two placements coincide: L ≥ 2w_x − 1 and M ≥ 2w_y − 1, w_x and w_y the largest
    /// extents of the templates' letters along x and y (sufficient, not necessary; the masks of two templates must not be
    /// translates of each other in the plane).</summary>
    public static (int L, int M) SufficientTorus(IReadOnlyCollection<string> templates)
    {
        int wx = 1, wy = 1;
        foreach (var t in templates)
        {
            var letters = TorusLetters(t);
            if (letters.Count == 0) continue;
            wx = Math.Max(wx, letters.Max(x => x.I) - letters.Min(x => x.I) + 1);
            wy = Math.Max(wy, letters.Max(x => x.J) - letters.Min(x => x.J) + 1);
        }
        return (2 * wx - 1, 2 * wy - 1);
    }

    /// <summary>The corner bound: in R = GF(2)[[u, v]] (u = 1 + x, v = 1 + y) the masks generate an ideal I; where they
    /// share no factor through x = y = 1 its colength ℓ = dim R/I is finite, 𝔪^ℓ ⊂ I, and n_u, n_v, the least n with u^n
    /// resp. v^n in 𝔪·I, are at most ℓ + 1. The verdict at (T₁, T₂) is the question whether a sum of masks with an odd
    /// number of negative ones lies in 𝔪·I + (u^T₁, v^T₂); once T₂ ≥ n_v the v-part lies in 𝔪·I and doubling T₂ changes
    /// nothing, likewise T₁ past n_u, so every corner has 2^(eⱼ − 1) &lt; n_u (where eⱼ ≥ 1) and 2^(fⱼ − 1) &lt; n_v (where
    /// fⱼ ≥ 1). Read in GF(2)[u, v]/(u^n, v^n), n = 8, 16 or 32: a read colength c below n is ℓ (Nakayama on
    /// 𝔪^c ⊂ I + 𝔪^(c + 1)), and then (u^n, v^n) ⊂ 𝔪^n ⊂ 𝔪·I and the read is exact; null where the read colength reaches
    /// n, which a common factor through (1, 1) always does.</summary>
    public static (int Colength, int Nu, int Nv)? CornerBound(IReadOnlyCollection<string> templates, int n = 16)
    {
        if (n is not (8 or 16 or 32)) throw new ArgumentOutOfRangeException(nameof(n), n, "read at n = 8, 16 or 32");
        int words = n * n / 64;
        var ideal = new Gf2Span(words);
        var maximalTimesIdeal = new Gf2Span(words);
        foreach (var t in templates)
        {
            var monomials = UvMonomials(t);
            for (int b = 0; b < n; b++)
                for (int a = 0; a < n; a++)
                {
                    var row = new ulong[words];
                    foreach (var (a0, b0) in monomials)
                        if (a0 + a < n && b0 + b < n) { int k = a0 + a + n * (b0 + b); row[k >> 6] ^= 1UL << (k & 63); }
                    ideal.Add((ulong[])row.Clone());
                    if (a != 0 || b != 0) maximalTimesIdeal.Add(row);
                }
        }
        int colength = n * n - ideal.Rank;
        if (colength >= n) return null;
        int Least(Func<int, int> cell)
        {
            for (int k = 1; k < n; k++)
            {
                var row = new ulong[words];
                int c = cell(k);
                row[c >> 6] = 1UL << (c & 63);
                if (maximalTimesIdeal.Contains(row)) return k;
            }
            return n;   // ℓ < n, so u^n and v^n lie in 𝔪^n ⊂ 𝔪^(ℓ + 1) ⊂ 𝔪·I
        }
        return (colength, Least(k => k), Least(k => n * k));
    }

    // the monomials u^a·v^b of a template's mask (Lucas: x^i holds u^a exactly when the bits of a lie among those of i)
    static HashSet<(int A, int B)> UvMonomials(string template)
    {
        var monomials = new HashSet<(int, int)>();
        foreach (var (i, j, c) in TorusLetters(template))
        {
            if (c is not ('X' or 'Y')) continue;
            for (int a = 0; a <= i; a++)
                if ((a & i) == a)
                    for (int b = 0; b <= j; b++)
                        if ((b & j) == b && !monomials.Add((a, b))) monomials.Remove((a, b));
        }
        return monomials;
    }

    // a GF(2) span of rows of several 64-bit words, pivots at the highest set bit
    sealed class Gf2Span(int words)
    {
        readonly Dictionary<int, ulong[]> pivots = new();
        public int Rank => pivots.Count;

        int Reduce(ulong[] row)
        {
            for (int w = words - 1; w >= 0; w--)
                while (row[w] != 0)
                {
                    int top = 64 * w + 63 - BitOperations.LeadingZeroCount(row[w]);
                    if (!pivots.TryGetValue(top, out var p)) return top;
                    for (int x = 0; x <= w; x++) row[x] ^= p[x];
                }
            return -1;
        }

        public void Add(ulong[] row) { int top = Reduce(row); if (top >= 0) pivots[top] = row; }
        public bool Contains(ulong[] row) => Reduce((ulong[])row.Clone()) < 0;
    }

    // ---------------------------------------------------------------- Theorem F beyond the square lattice

    /// <summary>A letter of a template on a lattice with Q sites in a unit cell: cell (I, J), site R of the cell.</summary>
    public readonly record struct LatticeLetter(int I, int J, int R, char Letter);

    /// <summary>The placed masks of templates on the L × M lattice of unit cells with Q sites each, site (i, j, r) at bit
    /// r + Q·(i + L·j) (L·M·Q ≤ 64), every template at every translate of the cells, wrapped, with its sign; refuses what
    /// <see cref="OnTorus"/> refuses.</summary>
    public static IReadOnlyList<(ulong Mask, int Sign)> OnLattice(IReadOnlyList<LatticeLetter[]> templates, int l, int m, int q)
    {
        if (templates.Count == 0) throw new ArgumentException("no templates");
        if (l < 1 || m < 1 || q < 1 || (long)l * m * q > 64) throw new ArgumentOutOfRangeException(nameof(l), $"{l} × {m} × {q}: one to 64 sites");
        var placed = new List<(ulong, int)>();
        var seen = new HashSet<ulong>();
        foreach (var t in templates)
        {
            if (!t.Any(x => x.Letter is 'X' or 'Y')) throw new ArgumentException("a template flips nothing");
            if (t.Any(x => x.R < 0 || x.R >= q || x.Letter is not ('X' or 'Y' or 'Z'))) throw new ArgumentException("a letter off the lattice");
            if (t.Select(x => (Mod(x.I, l), Mod(x.J, m), x.R)).Distinct().Count() < t.Length)
                throw new ArgumentException($"a template does not fit on the {l} × {m} lattice");
            int sign = t.Count(x => x.Letter is 'Y' or 'Z') % 2;
            for (int b = 0; b < m; b++)
                for (int a = 0; a < l; a++)
                {
                    ulong mask = 0;
                    foreach (var x in t)
                        if (x.Letter is 'X' or 'Y') mask |= 1UL << (x.R + q * (Mod(x.I + a, l) + l * Mod(x.J + b, m)));
                    if (!seen.Add(mask)) throw new ArgumentException($"two placements on the {l} × {m} lattice share a flip set: it would carry two strings");
                    placed.Add((mask, sign));
                }
        }
        return placed;
    }

    /// <summary>Theorem F on a lattice with Q sites in a unit cell: the masks folded onto the T₁ × T₂ cells of the 2-part,
    /// Q sites each, at every translate there; the lattice breaks exactly when a relation among them uses an odd number of
    /// negative strings, no colouring of the folded lattice existing.</summary>
    public static bool LatticeBreaks(IReadOnlyList<LatticeLetter[]> templates, int l, int m, int q)
    {
        OnLattice(templates, l, m, q);   // the setting
        int t1 = TwoPart(l), t2 = TwoPart(m);
        if ((long)t1 * t2 * q > 64) throw new ArgumentOutOfRangeException(nameof(l), "the folded lattice holds more than 64 sites");
        var rows = new List<(ulong, int)>();
        foreach (var t in templates)
        {
            int sign = t.Count(x => x.Letter is 'Y' or 'Z') % 2;
            for (int b = 0; b < t2; b++)
                for (int a = 0; a < t1; a++)
                {
                    ulong mask = 0;
                    foreach (var x in t)
                        if (x.Letter is 'X' or 'Y') mask ^= 1UL << (x.R + q * (Mod(x.I + a, t1) + t1 * Mod(x.J + b, t2)));
                    rows.Add((mask, sign));
                }
        }
        return !Balanced(rows);
    }

    static int Mod(int a, int n) => ((a % n) + n) % n;

    static int BitLength(int a) => a == 0 ? 0 : 32 - BitOperations.LeadingZeroCount((uint)a);

    static IReadOnlyList<(int E, int F)> MinimalPoints(IReadOnlyCollection<(int E, int F)> points) =>
        points.Where(p => !points.Any(o => o != p && o.E <= p.E && o.F <= p.F)).OrderBy(p => p.E).ThenBy(p => p.F).ToArray();

    /// <summary>Theorem F's second form: the torus breaks exactly when some choice of templates with an odd number of
    /// negative ones has its folded masks, one copy each, summing to a sum of differences g·m̄ⱼ + m̄ⱼ of folded translates
    /// (the ideal 𝔪·I of GF(2)[ℤ_T₁ × ℤ_T₂]); the pairing torus carries signs that are a linear function of the folded
    /// masks modulo those differences.</summary>
    public static bool TorusBreaksByDifferences(IReadOnlyCollection<string> templates, int l, int m)
    {
        OnTorus(templates, l, m);
        int t1 = TwoPart(l), t2 = TwoPart(m);
        var rows = new List<(ulong, int)>();
        foreach (var t in templates)
        {
            ulong at0 = FoldedMask(t, t1, t2);
            rows.Add((at0, Sign(t)));
            for (int b = 0; b < t2; b++)
                for (int a = 0; a < t1; a++) rows.Add((FoldedMask(t, t1, t2, a, b) ^ at0, 0));
        }
        return !Balanced(rows);
    }

    static int TwoPart(int n) => 1 << BitOperations.TrailingZeroCount((uint)n);

    // ---------------------------------------------------------------- Theorem G, the open patch

    /// <summary>A template's extents on the plane, written as for <see cref="TorusLetters"/>: one past its largest column
    /// and one past its largest row that carry a letter, identity letters before them counted.</summary>
    public static (int Wx, int Wy) Extent(string template)
    {
        var letters = TorusLetters(template);
        if (letters.Count == 0) throw new ArgumentException($"{template} has no letter");
        return (letters.Max(x => x.I) + 1, letters.Max(x => x.J) + 1);
    }

    /// <summary>The placed masks of templates on the open L × M patch, site (i, j) at bit i + L·j (L·M ≤ 64), every
    /// template at every (a, b) with a + w_x ≤ L and b + w_y ≤ M (<see cref="Extent"/>), with its sign. Refuses a template
    /// that flips nothing and two templates whose masks are translates in the plane, where a flip set would carry two
    /// strings.</summary>
    public static IReadOnlyList<(ulong Mask, int Sign)> OnPatch(IReadOnlyCollection<string> templates, int l, int m)
    {
        if (templates.Count == 0) throw new ArgumentException("no templates");
        if (l < 1 || m < 1 || (long)l * m > 64) throw new ArgumentOutOfRangeException(nameof(l), $"{l} × {m}: one to 64 sites");
        var shapes = new HashSet<string>();
        var placed = new List<(ulong, int)>();
        foreach (var t in templates)
        {
            var mask = PlaneMask(t);
            if (mask.Length == 0) throw new ArgumentException($"{t} flips nothing");
            if (!shapes.Add(string.Join(",", PlaneShape(mask)))) throw new ArgumentException("two masks are translates: a flip set would carry two strings");
            var (wx, wy) = Extent(t);
            int sign = Sign(t);
            for (int b = 0; b + wy <= m; b++)
                for (int a = 0; a + wx <= l; a++)
                {
                    ulong placedMask = 0;
                    for (int j = 0; j < mask.Length; j++)
                        for (int i = 0; i < 64; i++)
                            if (((mask[j] >> i) & 1) == 1) placedMask |= 1UL << (i + a + l * (j + b));
                    placed.Add((placedMask, sign));
                }
        }
        return placed;
    }

    /// <summary>Theorem G: two templates on the open L × M patch, their masks read in GF(2)[x, y] (the site (i, j) as
    /// x^i·y^j), g = gcd(m₁, m₂), qᵢ = mᵢ/g, c = s₁·q₂(1, 1) + s₂·q₁(1, 1) mod 2. Every relation is r·(q₂, q₁), and the
    /// degrees in x and in y add separately, so with c = 1 the patch breaks exactly on the quadrant L ≥ L*, M ≥ M*,
    /// L* = max(w_x1 + deg_x q₂, w_x2 + deg_x q₁) and M* the same along y, and with c = 0 on no patch (Corner null).
    /// Theorem B is its row M = 1 for templates ending in a letter (the chain counts identity letters at the end, the patch's extent
    /// does not). With <paramref name="divideByGcd"/> false the quotients are the masks themselves,
    /// the control that must fail.</summary>
    public static (bool Breaks, (int L, int M)? Corner) Patch(string t1, string t2, bool divideByGcd = true)
    {
        ulong[] m1 = PlaneMask(t1), m2 = PlaneMask(t2);
        if (m1.Length == 0 || m2.Length == 0) throw new ArgumentException("both templates must flip something");
        if (PlaneShape(m1).SequenceEqual(PlaneShape(m2))) throw new ArgumentException("the two masks are translates: a flip set would carry two strings");
        ulong[] q1 = m1, q2 = m2;
        if (divideByGcd)
        {
            var g = PlaneGcd(m1, m2);
            q1 = PlaneDivide(m1, g);
            q2 = PlaneDivide(m2, g);
        }
        int c = (Sign(t1) * PlaneAtOne(q2) + Sign(t2) * PlaneAtOne(q1)) % 2;
        if (c == 0) return (false, null);
        var (wx1, wy1) = Extent(t1);
        var (wx2, wy2) = Extent(t2);
        return (true, (Math.Max(wx1 + PlaneDegX(q2), wx2 + PlaneDegX(q1)), Math.Max(wy1 + q2.Length - 1, wy2 + q1.Length - 1)));
    }

    /// <summary>The gcd of two templates' masks in GF(2)[x, y] as Theorem G reads it, row j the coefficient of y^j (a
    /// polynomial in x, bit i the coefficient of x^i).</summary>
    public static IReadOnlyList<ulong> PatchGcd(string t1, string t2) => PlaneGcd(PlaneMask(t1), PlaneMask(t2));

    /// <summary>Theorem G at one patch: whether two templates break on the open L × M patch.</summary>
    public static bool PatchBreaks(string t1, string t2, int l, int m)
    {
        OnPatch(new[] { t1, t2 }, l, m);   // the setting
        var (breaks, corner) = Patch(t1, t2);
        return breaks && l >= corner!.Value.L && m >= corner.Value.M;
    }

    // a template's mask in GF(2)[x, y]: row j the coefficient of y^j, a polynomial in x (bit i), trailing zero rows dropped
    static ulong[] PlaneMask(string template)
    {
        var rows = new ulong[64];
        foreach (var (i, j, c) in TorusLetters(template))
        {
            if (i > 63 || j > 63) throw new ArgumentException($"{template} is wider than 64");
            if (c is 'X' or 'Y') rows[j] |= 1UL << i;
        }
        return PlaneTrim(rows);
    }

    static ulong[] PlaneTrim(ulong[] p)
    {
        int n = p.Length;
        while (n > 0 && p[n - 1] == 0) n--;
        return p.Take(n).ToArray();
    }

    // the mask moved to touch both axes, the form two translates share
    static ulong[] PlaneShape(ulong[] p)
    {
        int lead = 0;
        while (p[lead] == 0) lead++;
        int shift = p.Where(r => r != 0).Min(r => BitOperations.TrailingZeroCount(r));
        return p.Skip(lead).Select(r => r >> shift).ToArray();
    }

    static int PlaneAtOne(ulong[] p) => p.Sum(r => BitOperations.PopCount(r)) & 1;
    static int PlaneDegX(ulong[] p) => p.Where(r => r != 0).Max(Deg);

    // exact division in GF(2)[x][y]: the leading coefficient of g divides that of f wherever g divides f (every partial
    // remainder is a multiple of g), and no degree passes the dividend's
    static ulong[] PlaneDivide(ulong[] f, ulong[] g)
    {
        var rem = (ulong[])f.Clone();
        var q = new ulong[Math.Max(1, f.Length - g.Length + 1)];
        for (int d = f.Length - 1; d >= g.Length - 1; d--)
        {
            if (rem[d] == 0) continue;
            ulong k = Div(rem[d], g[^1]);
            q[d - g.Length + 1] = k;
            for (int j = 0; j < g.Length; j++) rem[d - g.Length + 1 + j] ^= MulChecked(k, g[j]);
        }
        if (rem.Any(r => r != 0)) throw new InvalidOperationException("not divisible");
        return PlaneTrim(q);
    }

    // the gcd in GF(2)[x, y] = GF(2)[x][y], a unique factorization domain: the gcd of the contents times the primitive
    // gcd, read by Euclid on pseudo-remainders, each made primitive (Gauss's lemma: a divisor of a primitive polynomial
    // is primitive, so taking the content out changes no common divisor). The remainders' degrees in x can pass 63 on the
    // way although the gcd divides both masks, so the coefficients are held as BigInteger, bit i the coefficient of x^i.
    static ulong[] PlaneGcd(ulong[] a, ulong[] b)
    {
        var p = BigPrimitive(a.Select(r => (BigInteger)r).ToArray());
        var q = BigPrimitive(b.Select(r => (BigInteger)r).ToArray());
        BigInteger c = BigGcd(BigContent(a.Select(r => (BigInteger)r).ToArray()), BigContent(b.Select(r => (BigInteger)r).ToArray()));
        if (p.Length < q.Length) (p, q) = (q, p);
        while (q.Length > 0)
        {
            var r = (BigInteger[])p.Clone();
            while (r.Length >= q.Length)
            {
                BigInteger lr = r[^1], lq = q[^1];
                int s = r.Length - q.Length;
                for (int j = 0; j < r.Length; j++) r[j] = BigMul(r[j], lq);
                for (int j = 0; j < q.Length; j++) r[s + j] ^= BigMul(lr, q[j]);
                r = BigPrimitive(BigTrim(r));
            }
            (p, q) = (q, r);
        }
        // the gcd divides both masks, so every coefficient of c·p fits in 64 bits
        return PlaneTrim(p.Select(r => (ulong)BigMul(c, r)).ToArray());
    }

    static BigInteger[] BigTrim(BigInteger[] p)
    {
        int n = p.Length;
        while (n > 0 && p[n - 1].IsZero) n--;
        return p.Take(n).ToArray();
    }

    static BigInteger BigContent(BigInteger[] p) => p.Aggregate(BigInteger.Zero, (g, r) => g.IsZero ? r : r.IsZero ? g : BigGcd(g, r));

    static BigInteger[] BigPrimitive(BigInteger[] p)
    {
        if (p.Length == 0) return p;
        BigInteger c = BigContent(p);
        return p.Select(r => r.IsZero ? r : BigDivMod(r, c).Q).ToArray();
    }

    static int BigDeg(BigInteger a) => a.IsZero ? -1 : (int)a.GetBitLength() - 1;

    static BigInteger BigMul(BigInteger a, BigInteger b)
    {
        BigInteger r = BigInteger.Zero;
        for (; !b.IsZero; b >>= 1, a <<= 1) if (!b.IsEven) r ^= a;
        return r;
    }

    static (BigInteger Q, BigInteger R) BigDivMod(BigInteger a, BigInteger b)
    {
        BigInteger q = BigInteger.Zero;
        int db = BigDeg(b);
        while (!a.IsZero && BigDeg(a) >= db) { int s = BigDeg(a) - db; q ^= BigInteger.One << s; a ^= b << s; }
        return (q, a);
    }

    static BigInteger BigGcd(BigInteger a, BigInteger b)
    {
        while (!b.IsZero) (a, b) = (b, BigDivMod(a, b).R);
        return a;
    }

    // ---------------------------------------------------------------- Theorem E, the character connection

    /// <summary>A Pauli string placed on at most 64 sites: X/Y mask, Y/Z mask, real coefficient.</summary>
    public readonly record struct PlacedString(ulong X, ulong Z, long C);

    /// <summary>What a flip set's ratio r(x) = g(x̄)/g(x) is: a constant i^Quarter times the character
    /// (−1)^⟨Read, x⟩, or an entry pair of different moduli (condition (i) fails, the palindrome breaks), or an
    /// entry that vanishes somewhere (an absent edge) or a ratio that is no character, both outside Theorem E.</summary>
    public enum RatioKind { Character, ModuliDiffer, AbsentEdge, NotCharacter }

    public readonly record struct FlipRatio(ulong Mask, RatioKind Kind, int Quarter, ulong Read);

    /// <summary>Templates placed with their identity letters at every start of the open chain (0 … n − w) or of the
    /// ring (every start, wrapped), n ≤ 64.</summary>
    public static IReadOnlyList<PlacedString> PlaceStrings(IEnumerable<(string T, long C)> templates, int n, bool ring)
    {
        if (n < 1 || n > 64) throw new ArgumentOutOfRangeException(nameof(n), n, "one to 64 sites");
        var placed = new List<PlacedString>();
        foreach (var (t, c) in templates)
        {
            if (t.Length > n) throw new ArgumentException($"{t} does not fit on {n} sites");
            int starts = ring ? n : n - t.Length + 1;
            for (int p = 0; p < starts; p++)
            {
                ulong x = 0, z = 0;
                for (int i = 0; i < t.Length; i++)
                {
                    ulong bit = 1UL << ((p + i) % n);
                    if (t[i] is 'X' or 'Y') x |= bit;
                    if (t[i] is 'Y' or 'Z') z |= bit;
                }
                placed.Add(new PlacedString(x, z, c));
            }
        }
        return placed;
    }

    /// <summary>The ratio of one flip set from its strings (all with the X/Y mask <paramref name="mask"/>): the entry
    /// g(x) = Σ c·(−i)^#Y·(−1)^⟨z, x⟩ depends only on the bits in the union U of the strings' z, so every x ⊆ U is
    /// read, against x̄'s bits U ⊕ x, in Gaussian integers, every comparison by cross-multiplication.</summary>
    public static FlipRatio ReadRatio(ulong mask, IReadOnlyList<PlacedString> strings)
    {
        ulong u = 0;
        foreach (var s in strings) u |= s.Z;
        if (BitOperations.PopCount(u) > 20) throw new ArgumentException("a flip set reading more than 20 sites");
        var bits = Enumerable.Range(0, 64).Where(j => ((u >> j) & 1) == 1).ToArray();
        GaussianInteger Entry(ulong x)
        {
            GaussianInteger sum = GaussianInteger.Zero;
            foreach (var s in strings)
            {
                var term = new GaussianInteger(s.C, BigInteger.Zero);
                for (int k = 0; k < BitOperations.PopCount(s.X & s.Z) % 4; k++) term = term * -GaussianInteger.I;
                if ((BitOperations.PopCount(x & s.Z) & 1) == 1) term = -term;
                sum = sum + term;
            }
            return sum;
        }
        static BigInteger Norm(GaussianInteger a) => a.Re * a.Re + a.Im * a.Im;
        static GaussianInteger Turn(GaussianInteger a, int q) { for (int k = 0; k < q; k++) a = a * GaussianInteger.I; return a; }
        int count = 1 << bits.Length;
        var a = new GaussianInteger[count];
        var b = new GaussianInteger[count];
        ulong Spread(int y) { ulong x = 0; for (int k = 0; k < bits.Length; k++) if (((y >> k) & 1) == 1) x |= 1UL << bits[k]; return x; }
        for (int y = 0; y < count; y++)
        {
            a[y] = Entry(Spread(y));
            b[y] = Entry(u ^ Spread(y));
        }
        bool absent = false, unequal = false;
        for (int y = 0; y < count; y++)
        {
            bool za = a[y] == GaussianInteger.Zero, zb = b[y] == GaussianInteger.Zero;
            if (za && zb) absent = true;
            else if (za || zb || Norm(a[y]) != Norm(b[y])) unequal = true;
        }
        if (unequal) return new FlipRatio(mask, RatioKind.ModuliDiffer, 0, 0);
        if (absent) return new FlipRatio(mask, RatioKind.AbsentEdge, 0, 0);
        int quarter = Enumerable.Range(0, 4).FirstOrDefault(q => b[0] == Turn(a[0], q), -1);
        if (quarter < 0) return new FlipRatio(mask, RatioKind.NotCharacter, 0, 0);
        ulong read = 0;
        for (int k = 0; k < bits.Length; k++)
        {
            int y = 1 << k;
            if (b[y] == Turn(a[y], quarter)) continue;
            if (b[y] == -Turn(a[y], quarter)) { read |= 1UL << bits[k]; continue; }
            return new FlipRatio(mask, RatioKind.NotCharacter, 0, 0);
        }
        for (int y = 0; y < count; y++)
        {
            var want = Turn(a[y], quarter);
            if ((BitOperations.PopCount(Spread(y) & read) & 1) == 1) want = -want;
            if (b[y] != want) return new FlipRatio(mask, RatioKind.NotCharacter, 0, 0);
        }
        return new FlipRatio(mask, RatioKind.Character, quarter, read);
    }

    /// <summary>Theorem E's ℤ/4 form. With every flip set's ratio a character i^(q_m)·(−1)^⟨R_m, x⟩ and condition (1)
    /// holding, β(a, b) = ⟨R(a), b⟩ is a symmetric bilinear form on the span W of the flip sets; walking W from 0 by
    /// Q₀(w ⊕ m) = Q₀(w) + q_m + 2·β(w, m) and R(w ⊕ m) = R(w) ⊕ R_m, the walk closes exactly when (2) (R linear on W)
    /// and (3) (hol = 1 on the relations) hold, and Q₀ is then a ℤ/4 quadratic refinement of β. Returns Q₀ and R on
    /// every element of W, or null where a ratio is no character or a condition fails. Reads the off-diagonal part
    /// only: the odd-weight diagonal, condition (ii), is <see cref="DecideCharacter(IReadOnlyList{PlacedString})"/>'s.</summary>
    public static IReadOnlyDictionary<ulong, (int Q, ulong R)>? QuadraticRefinement(IReadOnlyList<PlacedString> placed)
    {
        var flips = placed.Where(s => s.X != 0).GroupBy(s => s.X).Select(gr => ReadRatio(gr.Key, gr.ToList())).ToList();
        if (flips.Any(f => f.Kind != RatioKind.Character)) return null;
        static int Parity(ulong v) => BitOperations.PopCount(v) & 1;
        foreach (var a in flips)
            foreach (var b in flips)
                if (Parity(a.Read & b.Mask) != Parity(b.Read & a.Mask)) return null;
        var form = new Dictionary<ulong, (int Q, ulong R)> { [0] = (0, 0) };
        var queue = new Queue<ulong>();
        queue.Enqueue(0);
        while (queue.Count > 0)
        {
            ulong w = queue.Dequeue();
            var (q, r) = form[w];
            foreach (var m in flips)
            {
                var next = ((q + m.Quarter + 2 * Parity(r & m.Mask)) % 4, r ^ m.Read);
                if (form.TryGetValue(w ^ m.Mask, out var have)) { if (have != next) return null; }
                else
                {
                    if (form.Count >= 1 << 22) throw new ArgumentException("the span of the flip sets passes 2^22 elements");
                    form[w ^ m.Mask] = next;
                    queue.Enqueue(w ^ m.Mask);
                }
            }
        }
        return form;
    }

    /// <summary>The carrier of Theorem E's ℤ/4 form on n sites: X^⊗N·D with D = diag(i^(d_x)), d(x₀ ⊕ w) = −Q₀(w) −
    /// 2·⟨R(w), x₀⟩ (mod 4), x₀ the least element of its coset of W, so that the component through x₀ carries the
    /// refinement Q₀ + 2·⟨R(·), x₀⟩. With <paramref name="cosetTerm"/> false the term 2·⟨R(w), x₀⟩ is dropped, the
    /// control: Q₀ alone on every coset.</summary>
    public static int[] CarrierQuarters(IReadOnlyDictionary<ulong, (int Q, ulong R)> form, int n, bool cosetTerm = true)
    {
        if (n < 1 || n > 24) throw new ArgumentOutOfRangeException(nameof(n), n, "one to 24 sites");
        int dim = 1 << n;
        var d = new int[dim];
        var done = new bool[dim];
        for (int x0 = 0; x0 < dim; x0++)
        {
            if (done[x0]) continue;
            foreach (var (w, (q, r)) in form)
            {
                int x = x0 ^ (int)w;
                done[x] = true;
                int t = cosetTerm ? 2 * (BitOperations.PopCount(r & (ulong)x0) & 1) : 0;
                d[x] = ((-q - t) % 4 + 8) % 4;
            }
        }
        return d;
    }

    public enum CharacterVerdict { Pairs, Breaks, Outside }

    /// <summary>Theorem E at one N: the odd-weight diagonal, then every flip set's ratio (a different modulus breaks;
    /// an absent edge or a ratio that is no character leaves the setting outside), then (1) the symmetric reading on
    /// every two flip sets, (2) the read sets summed over each relation of a basis, (3) hol on that basis.</summary>
    public static CharacterVerdict DecideCharacter(IReadOnlyList<PlacedString> placed) => DecideCharacter(placed, true, true, true);

    /// <summary>The same with conditions (1), (2) and (3) each switchable, for the controls that must fail.</summary>
    public static CharacterVerdict DecideCharacter(IReadOnlyList<PlacedString> placed, bool squares, bool readSums, bool constants)
    {
        var diagonal = new Dictionary<ulong, long>();
        foreach (var s in placed.Where(s => s.X == 0)) diagonal[s.Z] = diagonal.GetValueOrDefault(s.Z) + s.C;
        if (diagonal.Any(kv => kv.Value != 0 && (BitOperations.PopCount(kv.Key) & 1) == 1)) return CharacterVerdict.Breaks;
        var flips = placed.Where(s => s.X != 0).GroupBy(s => s.X).Select(gr => ReadRatio(gr.Key, gr.ToList())).ToList();
        if (flips.Any(f => f.Kind == RatioKind.ModuliDiffer)) return CharacterVerdict.Breaks;
        if (flips.Any(f => f.Kind != RatioKind.Character)) return CharacterVerdict.Outside;
        static int Parity(ulong v) => BitOperations.PopCount(v) & 1;
        if (squares)
            foreach (var p in flips)
                foreach (var q in flips)
                    if (Parity(p.Read & q.Mask) != Parity(q.Read & p.Mask)) return CharacterVerdict.Breaks;
        // a basis of the relations: elimination on the masks, each row tagged with the flip sets it sums
        var basis = new Dictionary<int, (ulong Mask, BigInteger Tag)>();
        for (int i = 0; i < flips.Count; i++)
        {
            ulong m = flips[i].Mask;
            BigInteger tag = BigInteger.One << i;
            while (m != 0)
            {
                int top = 63 - BitOperations.LeadingZeroCount(m);
                if (!basis.TryGetValue(top, out var row)) { basis[top] = (m, tag); break; }
                m ^= row.Mask;
                tag ^= row.Tag;
            }
            if (m != 0) continue;
            ulong readSum = 0, at = 0;
            int quarter = 0;
            for (int k = 0; k < flips.Count; k++)
            {
                if (((tag >> k) & BigInteger.One).IsZero) continue;
                readSum ^= flips[k].Read;
                quarter += flips[k].Quarter + 2 * Parity(flips[k].Read & at);
                at ^= flips[k].Mask;
            }
            if (readSums && readSum != 0) return CharacterVerdict.Breaks;
            if (constants && quarter % 4 != 0) return CharacterVerdict.Breaks;
        }
        return CharacterVerdict.Pairs;
    }

    /// <summary>Which branch of Corollary E2 decided a ring.</summary>
    public enum RingShapeBranch { ModuliDiffer, SquaresFail, NoRelation, ReadSetNotDivisible, OnePlusXDoesNotDivide, HolonomyPlus, HolonomyMinus, Outside }

    /// <summary>Corollary E2: templates whose masks are translates of one shape f, on the ring of n ≥ 2w − 1 sites
    /// (n ≤ 63). The flip set at p = 0 gives ω and the read set R; with h = gcd(f, x^n + 1) and g = (x^n + 1)/h the
    /// ring pairs exactly when (1) holds and h = 1, or h divides R and either 1 + x does not divide h or hol(g) = 1. A flip set
    /// whose ratio is no character, or whose entry vanishes somewhere, gives the branch Outside with Breaks false: undecided,
    /// so the branch is read before Breaks.</summary>
    public static (bool Breaks, RingShapeBranch Branch) RingOneShape(IReadOnlyList<(string T, long C)> templates, int n,
        bool readSetTest = true, bool orderingTerm = true)
    {
        var masks = templates.Select(t => Mask(t.T)).ToList();
        if (masks.Any(m => m == 0)) throw new ArgumentException("every template must flip something");
        ulong f = Shape(masks[0]);
        if (masks.Any(m => Shape(m) != f)) throw new ArgumentException("the templates must share one flip shape");
        int w = templates.Max(t => t.T.Length);
        if (n < 2 * w - 1 || n > 63) throw new ArgumentOutOfRangeException(nameof(n), n, $"the ring runs from {2 * w - 1} to 63 sites");
        var placed = PlaceStrings(templates, n, ring: true);
        var flips = placed.Where(s => s.X != 0).GroupBy(s => s.X).ToDictionary(gr => gr.Key, gr => ReadRatio(gr.Key, gr.ToList()));
        if (flips.Values.Any(r => r.Kind == RatioKind.ModuliDiffer)) return (true, RingShapeBranch.ModuliDiffer);
        if (flips.Values.Any(r => r.Kind != RatioKind.Character)) return (false, RingShapeBranch.Outside);
        static int Parity(ulong v) => BitOperations.PopCount(v) & 1;
        foreach (var p in flips.Values)
            foreach (var q in flips.Values)
                if (Parity(p.Read & q.Mask) != Parity(q.Read & p.Mask)) return (true, RingShapeBranch.SquaresFail);
        ulong ring = (1UL << n) | 1UL;
        ulong h = Gcd(f, ring);
        if (h == 1) return (false, RingShapeBranch.NoRelation);
        var zero = flips[f];
        if (readSetTest && Mod(zero.Read, h) != 0) return (true, RingShapeBranch.ReadSetNotDivisible);
        if ((BitOperations.PopCount(h) & 1) == 1) return (false, RingShapeBranch.OnePlusXDoesNotDivide);
        ulong g = Div(ring, h), at = 0;
        int quarter = 0;
        for (int p = 0; p < n; p++)
        {
            if (((g >> p) & 1) == 0) continue;
            ulong m = Rotate(f, p, n), read = Rotate(zero.Read, p, n);
            quarter += zero.Quarter + (orderingTerm ? 2 * Parity(read & at) : 0);
            at ^= m;
        }
        return quarter % 4 == 0 ? (false, RingShapeBranch.HolonomyPlus) : (true, RingShapeBranch.HolonomyMinus);
    }

    /// <summary>Which branch of Corollary E3 decided a chain.</summary>
    public enum ChainShapeBranch { ModuliDiffer, SquaresFail, NoRelation, ReadSumNonzero, HolonomyPlus, HolonomyMinus, Outside }

    /// <summary>The two numbers Corollary E3 reads off two templates, each a combination of strings of one width on one
    /// nonzero mask, the masks of different shapes: N* = max(w₁ + deg q₂, w₂ + deg q₁) as in Theorem B, and the read sum
    /// D = q₂·R₁ + q₁·R₂ of the relation r = 1, R_i the read set of template i's flip set at p = 0. A template whose ratio
    /// is no character has no read set, and is refused.</summary>
    public static (int NStar, ulong ReadSum) ChainTwoShapeData(IReadOnlyList<(string T, long C)> first, IReadOnlyList<(string T, long C)> second)
    {
        var (w1, m1, r1) = OneShapeTemplate(first);
        var (w2, m2, r2) = OneShapeTemplate(second);
        if (Shape(m1) == Shape(m2)) throw new ArgumentException("the two templates must have masks of different shapes");
        if (r1.Kind != RatioKind.Character || r2.Kind != RatioKind.Character) throw new ArgumentException("a template whose ratio is no character has no read set");
        int nStar = TwoShapeFloor(w1, m1, w2, m2);
        if (nStar > 64) throw new ArgumentException("the relation r = 1 needs more than 64 sites");
        ulong g = Gcd(m1, m2);
        return (nStar, Mul(Div(m2, g), r1.Read) ^ Mul(Div(m1, g), r2.Read));
    }

    /// <summary>Corollary E3: two templates on the open chain of n sites (n ≤ 64), each a combination of strings of one
    /// width on one nonzero mask, the masks of different shapes, every flip set's ratio a constant times a character. Where
    /// (1) holds at n the chain pairs below N*, and from N* on exactly when D = 0 and hol(u₀) = 1, u₀ the relation that uses
    /// template 1 at the placements q₂ and template 2 at q₁. A flip set whose ratio is no character, or whose entry vanishes
    /// somewhere, gives the branch Outside with Breaks false: undecided, so the branch is read before Breaks.</summary>
    public static (bool Breaks, ChainShapeBranch Branch) ChainTwoShapes(IReadOnlyList<(string T, long C)> first, IReadOnlyList<(string T, long C)> second,
        int n, bool readSumTest = true, bool orderingTerm = true)
    {
        var (w1, m1, _) = OneShapeTemplate(first);
        var (w2, m2, _) = OneShapeTemplate(second);
        if (Shape(m1) == Shape(m2)) throw new ArgumentException("the two templates must have masks of different shapes");
        if (n < Math.Max(w1, w2) || n > 64) throw new ArgumentOutOfRangeException(nameof(n), n, "the chain must hold both templates, at most 64 sites");
        var placed = PlaceStrings(first.Concat(second), n, ring: false);
        var flips = placed.GroupBy(s => s.X).ToDictionary(gr => gr.Key, gr => ReadRatio(gr.Key, gr.ToList()));
        if (flips.Values.Any(r => r.Kind == RatioKind.ModuliDiffer)) return (true, ChainShapeBranch.ModuliDiffer);
        if (flips.Values.Any(r => r.Kind != RatioKind.Character)) return (false, ChainShapeBranch.Outside);
        static int Parity(ulong v) => BitOperations.PopCount(v) & 1;
        foreach (var p in flips.Values)
            foreach (var q in flips.Values)
                if (Parity(p.Read & q.Mask) != Parity(q.Read & p.Mask)) return (true, ChainShapeBranch.SquaresFail);
        if (n < TwoShapeFloor(w1, m1, w2, m2)) return (false, ChainShapeBranch.NoRelation);
        var (_, readSum) = ChainTwoShapeData(first, second);
        if (readSumTest && readSum != 0) return (true, ChainShapeBranch.ReadSumNonzero);
        ulong g = Gcd(m1, m2), at = 0;
        int quarter = 0;
        foreach (var (m, u) in new[] { (m1, Div(m2, g)), (m2, Div(m1, g)) })
            for (int p = 0; p < 64; p++)
            {
                if (((u >> p) & 1) == 0) continue;
                var r = flips[m << p];
                quarter += r.Quarter + (orderingTerm ? 2 * Parity(r.Read & at) : 0);
                at ^= m << p;
            }
        if (at != 0) throw new InvalidOperationException("the relation r = 1 does not close");
        if (readSumTest && quarter % 2 != 0) throw new InvalidOperationException("hol(u₀) = ±i where the read sums vanish, against the proof");
        return quarter % 4 == 0 ? (false, ChainShapeBranch.HolonomyPlus) : (true, ChainShapeBranch.HolonomyMinus);
    }

    /// <summary>Corollary E3′: any number of templates on the open chain of n sites (n ≤ 64), each a combination of
    /// strings of one width on one nonzero mask, no two masks translates, every flip set's ratio a constant times a
    /// character. Over one basis bⱼ of the relations S reduced for the widths, call bⱼ bad when its read sum
    /// ρ(bⱼ) = Σ bⱼᵢ·Rᵢ does not vanish (N₁* the least shifted degree of such) or when it vanishes and hol(bⱼ) = −1 (N₂*).
    /// Where (1) holds at n the chain pairs exactly below the least shifted degree of a bad basis vector. Two templates are
    /// Corollary E3, branch for branch. <paramref name="readSumTest"/> false (N₁* ignored, hol read on every basis vector)
    /// and <paramref name="orderingTerm"/> false (hol walked without it) are the controls.</summary>
    public static (bool Breaks, ChainShapeBranch Branch) ChainShapes(IReadOnlyList<IReadOnlyList<(string T, long C)>> templates,
        int n, bool readSumTest = true, bool orderingTerm = true)
    {
        if (templates.Count == 0) throw new ArgumentException("no templates");
        var read = templates.Select(OneShapeTemplate).ToArray();
        if (read.Select(r => Shape(r.Mask)).Distinct().Count() < read.Length)
            throw new ArgumentException("two masks are translates: a flip set would carry the strings of two templates");
        if (n < read.Max(r => r.Width) || n > 64) throw new ArgumentOutOfRangeException(nameof(n), n, "the chain must hold every template, at most 64 sites");
        var placed = PlaceStrings(templates.SelectMany(t => t), n, ring: false);
        var flips = placed.GroupBy(s => s.X).ToDictionary(gr => gr.Key, gr => ReadRatio(gr.Key, gr.ToList()));
        if (flips.Values.Any(r => r.Kind == RatioKind.ModuliDiffer)) return (true, ChainShapeBranch.ModuliDiffer);
        if (flips.Values.Any(r => r.Kind != RatioKind.Character) || read.Any(r => r.Ratio.Kind != RatioKind.Character))
            return (false, ChainShapeBranch.Outside);
        static int Parity(ulong v) => BitOperations.PopCount(v) & 1;
        foreach (var p in flips.Values)
            foreach (var q in flips.Values)
                if (Parity(p.Read & q.Mask) != Parity(q.Read & p.Mask)) return (true, ChainShapeBranch.SquaresFail);
        if (read.Length == 1) return (false, ChainShapeBranch.NoRelation);
        var widths = read.Select(r => r.Width).ToArray();
        var s = ShiftedReduce(KernelBasis(read.Select(r => r.Mask).ToArray()), widths);
        ulong Rho(ulong[] u)
        {
            ulong acc = 0;
            for (int i = 0; i < u.Length; i++) acc ^= MulChecked(u[i], read[i].Ratio.Read);
            return acc;
        }
        int Sdeg(ulong[] u) => ShiftedDegree(u, widths).Degree;
        int? n1 = readSumTest ? s.Where(b => Rho(b) != 0).Select(Sdeg).Cast<int?>().Min() : null;
        int Hol(ulong[] u)
        {
            ulong at = 0;
            int quarter = 0;
            for (int i = 0; i < u.Length; i++)
                for (int p = 0; p < 64; p++)
                {
                    if (((u[i] >> p) & 1) == 0) continue;
                    quarter += read[i].Ratio.Quarter + (orderingTerm ? 2 * Parity(ShiftLeft(read[i].Ratio.Read, p) & at) : 0);
                    at ^= ShiftLeft(read[i].Mask, p);
                }
            if (at != 0) throw new InvalidOperationException("a basis vector of the relations does not close");
            if (readSumTest && orderingTerm && quarter % 2 != 0) throw new InvalidOperationException("hol = ±i on a relation whose read sum vanishes, against the proof");
            return ((quarter % 4) + 4) % 4;
        }
        // below N₁* every relation lies in the span of the basis vectors whose read sums vanish (predictable degree), where hol
        // is a homomorphism: the first −1 there sits on one of them
        int? n2 = s.Where(b => (!readSumTest || Rho(b) == 0) && Hol(b) != 0).Select(Sdeg).Cast<int?>().Min();
        if (n1 is int a && n >= a && (n2 is not int b0 || a <= b0)) return (true, ChainShapeBranch.ReadSumNonzero);
        if (n2 is int b && n >= b) return (true, ChainShapeBranch.HolonomyMinus);
        return n < s.Min(Sdeg) ? (false, ChainShapeBranch.NoRelation) : (false, ChainShapeBranch.HolonomyPlus);
    }

    // a template of E3: strings of one width on one nonzero mask, with its flip set's ratio read at its own window
    static (int Width, ulong Mask, FlipRatio Ratio) OneShapeTemplate(IReadOnlyList<(string T, long C)> strings)
    {
        if (strings.Count == 0) throw new ArgumentException("a template needs a string");
        int w = strings[0].T.Length;
        ulong m = Mask(strings[0].T);
        if (m == 0) throw new ArgumentException("every template must flip something");
        if (strings.Any(s => s.T.Length != w || Mask(s.T) != m)) throw new ArgumentException("a template's strings share one width and one mask");
        return (w, m, ReadRatio(m, PlaceStrings(strings, w, ring: false)));
    }

    // Theorem B's N* = max(w₁ + deg q₂, w₂ + deg q₁)
    static int TwoShapeFloor(int w1, ulong m1, int w2, ulong m2)
    {
        ulong g = Gcd(m1, m2);
        return Math.Max(w1 + Deg(Div(m2, g)), w2 + Deg(Div(m1, g)));
    }

    static ulong Mul(ulong a, ulong b)
    {
        ulong r = 0;
        for (; b != 0; b >>= 1, a <<= 1) if ((b & 1) == 1) r ^= a;
        return r;
    }

    // ---------------------------------------------------------------- Theorem H, the size of the break at strong dephasing

    /// <summary>Theorem H under Theorem A's assumption: with H = J·Σ c_P·P the far edge's eigenvalue nearest −2σ lies at
    /// J²·λ + O(J³) from it, λ = 2·min over φ ∈ GF(2)^N of Σ c_P²/γ(m_P) over the strings with ⟨φ, m_P⟩ ≠ s(P), twice the
    /// distance of the sign vector from the code of balanced signings, each string weighted by c_P² over the rates on its
    /// flip set. Equal strings are collected first; two different strings on one nonzero flip set, or a diagonal part of
    /// odd weight left after collecting, are refused. Up to 24 sites (φ is enumerated).</summary>
    public static double StrongDephasingBreak(IReadOnlyList<PlacedString> placed, IReadOnlyList<double> gamma)
    {
        int n = gamma.Count;
        if (n < 1 || n > 24) throw new ArgumentOutOfRangeException(nameof(gamma), n, "one to 24 sites");
        if (gamma.Any(g => !(g > 0))) throw new ArgumentException("every rate strictly positive", nameof(gamma));
        var collected = new Dictionary<(ulong X, ulong Z), long>();
        foreach (var p in placed)
        {
            if ((p.X | p.Z) >> n != 0) throw new ArgumentException($"a string has a letter beyond site {n - 1}", nameof(placed));
            collected[(p.X, p.Z)] = (collected.TryGetValue((p.X, p.Z), out var c) ? c : 0) + p.C;
        }
        var strings = collected.Where(kv => kv.Value != 0).ToList();
        foreach (var kv in strings.Where(kv => kv.Key.X == 0))
            if (BitOperations.PopCount(kv.Key.Z) % 2 == 1) throw new ArgumentException("a diagonal string of odd weight remains", nameof(placed));
        var hops = strings.Where(kv => kv.Key.X != 0).ToList();
        if (hops.GroupBy(kv => kv.Key.X).Any(g => g.Count() > 1)) throw new ArgumentException("two strings on one flip set", nameof(placed));
        var terms = hops.Select(kv =>
        {
            double rate = 0;
            for (int l = 0; l < n; l++) if (((kv.Key.X >> l) & 1) == 1) rate += gamma[l];
            return (Mask: kv.Key.X, Sign: BitOperations.PopCount(kv.Key.Z) & 1, Weight: (double)kv.Value * kv.Value / rate);
        }).ToList();
        double best = double.PositiveInfinity;
        for (ulong phi = 0; phi < 1UL << n; phi++)
        {
            double cost = 0;
            foreach (var t in terms) if ((BitOperations.PopCount(phi & t.Mask) & 1) != t.Sign) cost += t.Weight;
            best = Math.Min(best, cost);
        }
        return 2 * best;
    }

    public override string DisplayName =>
        "F173: the signed balance of the flip masks";

    public override string Summary =>
        "one string per flip set under Z dephasing: each edge a sign (−1)^(#Y + #Z), the palindrome the balance of the masks' " +
        "signed Cayley graph; two templates break exactly from N* when s₁·q₂(1) + s₂·q₁(1) is odd (F115 its all-negative case), " +
        $"one template on the ring exactly when s = 1 and v_(1+x)(m) ≥ 2^(v₂(N)), and any number on the ring read through " +
        $"whether 2^(v₂(N)) exceeds the least valuation; several strings per flip set with a character ratio by a symmetric reading " +
        $"and the relations read against the read sets and a holonomy, one flip shape on the ring and any number on the chain in closed form; on the torus " +
        $"a colouring of period (2^(v₂(L)), 2^(v₂(M))), and two templates on the open patch break on one quadrant or on none; with one string per flip set, at strong dephasing the far edge misses −2σ by " +
        $"J² times twice the rate-weighted distance of the signs from the code of balanced signings ({Tier.Label()}; live witnesses: " +
        $"inspect --root signedbalance, inspect --root torusbalance)";

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
            yield return new InspectableNode("Theorem B′, any number of templates on the chain",
                summary: "relations u with Σ uᵢmᵢ = 0, free of rank k − 1, fitting at N when max(deg uᵢ + wᵢ) ≤ N; over a basis in weak " +
                         "Popov form for the widths the chain breaks exactly from the least shifted degree of an unbalanced basis vector, " +
                         "never if none is unbalanced; k = 2 is Theorem B, its c the girth ladder's syzygy polynomial read at 1");
            yield return new InspectableNode("Theorem C, one template on the ring",
                summary: "relations u·m ≡ 0 mod x^N + 1; unbalanced exactly when s = 1 and v_(1+x)(m) ≥ 2^(v₂(N)): XY breaks at odd N, " +
                         "XIY unless 4 | N, XIIIY unless 8 | N");
            yield return new InspectableNode("Theorem D, any number of templates on the ring",
                summary: "x^N + 1 = (1 + x)^T·F with T = 2^(v₂(N)) and F(1) = 1; a relation reduced modulo (1 + x)^T keeps its sign count, " +
                         "and one there times F is one on the ring, so the verdict is read in GF(2)[y]/(y^T), y = 1 + x, where each mask is " +
                         "y^(vᵢ) times a unit, v₀ the least; if T ≤ v₀ every mask is 0 there and the ring breaks on any odd sign; otherwise, i₀ " +
                         "of valuation v₀, every relation is a combination of the relations that use template j once and i₀ by " +
                         "y^(v_j − v₀) times a unit, of sign count s_j + s_i₀·[v_j = v₀], and one at i₀ alone, of sign count 0");
            yield return new InspectableNode("Theorem E, the character connection",
                summary: "every flip set's ratio ω_m·(−1)^⟨R_m, x⟩, no entry vanishing: swapping two steps costs (−1)^(⟨R_b, a⟩ + ⟨R_a, b⟩), " +
                         "so every cycle reduces to one relation walked once; moving its start by δ costs (−1)^⟨⊕ R_m, δ⟩, so a relation " +
                         "with a nonzero read sum leaves a component carrying −1; the rest is hol on a basis. Two strings of different " +
                         "signs have a character ratio exactly when their numbers of Y differ by an odd number and c₁² = c₂², R = z₁ ⊕ z₂; " +
                         "four strings can too, X ⊗ 2·CZ among them. The ℤ/4 form: with ω_m = i^(q_m) a good component is a ℤ/4 " +
                         "quadratic refinement of β(a, b) = ⟨R(a), b⟩ on the span W, the ordering term its polarization; with R " +
                         "linear on W the palindrome is the existence of one refinement Q₀ with Q₀(m) = q_m, carried by " +
                         "X^⊗N·i^(−Q₀(w) − 2⟨R(w), x₀⟩) on each coset x₀ + W");
            yield return new InspectableNode("Corollary E2, one flip shape on the ring",
                summary: "on the ring of N ≥ 2w − 1 sites with the symmetric reading holding: relations r·(x^N + 1)/h, h = gcd(f, x^N + 1); the read sum vanishes exactly when h | R; hol is rotation-invariant, " +
                         "so it lives on GF(2)[x]/(h, 1 + x): X·(X ± Y)·X breaks exactly at 3 | N, XXX + YYY pairs on every ring");
            yield return new InspectableNode("Corollary E3, two flip shapes on the chain",
                summary: "relations r·(q₂, q₁) as in Theorem B, read sum r·D with D = q₂R₁ + q₁R₂, u₀ the relation r = 1 and hol(r·u₀) = hol(u₀)^(r(1)); with the symmetric " +
                         "reading holding the chain pairs below N* and from N* exactly when D = 0 and hol(u₀) = 1; where that reading holds at " +
                         "every offset D = 0, so D ≠ 0 only brings the break forward to N*: XXX + YYY beside the field X + Y pairs on every " +
                         "chain, XXX − YYY beside it breaks from N = 3");
            yield return new InspectableNode("Theorem G, two templates on the open patch",
                summary: "GF(2)[x, y] is a unique factorization domain, so every relation is r·(q₂, q₁), negative count r(1, 1)·c; the degrees " +
                         "in x and in y add separately, so r = 1 is least in both at once and the breaking patches are one quadrant, " +
                         "L ≥ max(w_x,1 + deg_x q₂, w_x,2 + deg_x q₁), M ≥ max(w_y,1 + deg_y q₂, w_y,2 + deg_y q₁); XY across beside XX/X " +
                         "breaks exactly from 3 × 2; Theorem B is the row M = 1 (templates ending in a letter); for three or more the breaking patches, the union of the " +
                         "quadrants of the unbalanced relations, can have several corners (X beside XY across and Y over X: (2, 1) and (1, 2)), read patch by patch");
            yield return new InspectableNode("Theorem H, the size of the break at strong dephasing",
                summary: "H scaled to J·H with condition (ii) everywhere: the first order vanishes on both edges, and the second is " +
                         "−J²·𝓛 on the populations and −2σ + J²·𝓜 on the far cells, 𝓛 the Zeno Laplacian with weights |H_ba|²/γ(a ⊕ b) and " +
                         "𝓜 the same with each hop carrying H_ba·H_āb̄, both positive semidefinite, their kernels Theorem 2's two counts; " +
                         "one string per flip set: 𝓜 = Σ (c²/γ(m))(1 − (−1)^s·T_m), diagonal in the characters, λ_min twice the weighted " +
                         "distance of the signs from the code of balanced signings; the XY ring at one rate J²c²/γ on every odd ring, " +
                         "XZ + ZY on the chain two π-flux end squares, 2(3 − √5)·J²c²/γ");
            yield return new InspectableNode("the two pure cases",
                summary: "every sign negative: balanced means bipartite (F103 §7.4, F115 statement 3); every sign positive: X^⊗N " +
                         "commutes with H, a colouring (F108 Part 1); a balanced signing is the lit string X^⊗N·Z^φ commuting with H");
        }
    }
}
