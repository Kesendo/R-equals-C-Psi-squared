using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F169, the colour law for Heisenberg bonds off the letters
/// (<c>docs/proofs/PROOF_HEISENBERG_COLOURS_OFF_THE_LETTERS.md</c>). H = Σ J_ij·σ_i·σ_j + Σ h_l·σ_l with every
/// J_ij ≠ 0 on a connected graph, every site watched by one jump n_l·σ_l along an integer axis. A <b>colour</b> is
/// a direction m perpendicular to every axis; it is <b>admitted</b> when every field lies along it, and then
/// ⊗(m·σ) is an invertible element of F158's far kernel.
///
/// <para><b>The half-turn law.</b> (a·σ)⊗(b·σ) commutes with Σ W_kl σ^k⊗σ^l exactly when R_a·W·R_b = W, R_v the
/// half-turn about v; for W = J·1 that asks one colour on both ends of every bond.</para>
///
/// <para><b>The split.</b> [σ_i·σ_j, (a·σ_i)(b·σ_j)] = 2i(a×b)·(σ_j − σ_i) leaves one identity factor and a field's
/// commutator none, so the far equation splits into bond equations and field equations, and wherever the bond
/// equations admit only the colours the palindrome holds exactly when an admitted colour exists.</para>
///
/// <para><b>Every connected graph.</b> The bond equations admit only the colours. A bond whose ends share an axis
/// drops out and leaves a symmetric pair (deleted from a cycle, or a bridge gluing two smaller graphs). Every other
/// bond keeps its pair in span{m⊗m, p⊗p}; two such bonds at one site meeting at an angle neither straight nor right
/// leave nothing; otherwise every bond reads "both ends on its line or neither", which on a connected graph leaves
/// one word and its complement, and the site equations, whose terms are distinct words from N = 4 on, keep only the
/// word that carries every bond's line at both ends, the colour.</para></summary>
public sealed class HeisenbergColourLawClaim : Claim
{
    /// <summary>F158: the palindrome holds exactly when the far kernel holds an invertible element; a colour is one.</summary>
    public PalindromeTwoEndCountClaim TwoEndCount { get; }

    /// <summary>True where the law is proven for this graph: connected. The axes do not enter.</summary>
    public static bool Proven(int n, IReadOnlyList<(int A, int B)> edges, IReadOnlyList<(long X, long Y, long Z)> axes)
        => Components(n, edges) == 1;

    /// <summary>The dimension of F158's far kernel the law predicts on a connected graph: the admitted colours'
    /// symmetric tensors. N + 1 with all axes along one n and no field (Sym^N of the plane ⊥ n); 1 when one colour is
    /// admitted otherwise; 0 when none is and the bond equations' solutions already leave no room (the axes span space,
    /// or they span a plane and a field leaves its normal, or under one axis the fields are perpendicular to it and not
    /// parallel to one another). Null under one axis when a field has a part along it: there the law says only that
    /// the palindrome fails (F158 §(f5)), and the far kernel can still hold non-invertible elements (a uniform field
    /// along Z at even N leaves one).</summary>
    public static int? PredictedFarDimension(IReadOnlyList<(long X, long Y, long Z)> axes, IReadOnlyList<(long X, long Y, long Z)> fields)
    {
        var dirs = axes.Select(Primitive).Distinct().ToList();
        var f = fields.Where(v => v != (0, 0, 0)).Select(Primitive).Distinct().ToList();
        if (dirs.All(d => Cross(d, dirs[0]) == (0, 0, 0)))
        {
            if (f.Any(v => Dot(v, dirs[0]) != 0)) return null;
            if (f.Count > 1) return 0;
            return f.Count == 0 ? axes.Count + 1 : 1;
        }
        var normal = dirs.Select(d => Primitive(Cross(dirs[0], d))).First(v => v != (0, 0, 0));
        if (dirs.Any(d => Dot(normal, d) != 0)) return 0;
        return f.All(v => v == normal) ? 1 : 0;
    }

    public HeisenbergColourLawClaim(PalindromeTwoEndCountClaim twoEndCount)
        : base("F169, the colour law: Heisenberg bonds of nonzero weights on a connected graph, every site watched along one " +
               "integer axis, fields of any direction. A product of colours commutes with a bond of coupling matrix W exactly " +
               "when the half-turns about its ends keep W; the far equation splits into bond and field equations; wherever " +
               "the bond equations admit only the colours the palindrome holds exactly when a colour perpendicular to every " +
               "axis carries every field. The bond equations admit only the colours on every connected graph: a site whose " +
               "bonds meet at an angle neither straight nor right admits nothing, and otherwise each bond asks both ends " +
               "on its line or neither, which leaves at most two words, and the site equations keep only the colour. At axes one " +
               "turn makes letters F165's Theorems 1 to 3 give the same answer.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_HEISENBERG_COLOURS_OFF_THE_LETTERS.md + " +
               "docs/ANALYTICAL_FORMULAS.md (F169, and F138's paragraph \"Off the letters\") + " +
               "docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md (F165, Theorems 1 to 3) + " +
               "docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md (F158, section (h): jumps along any axis)")
    {
        TwoEndCount = twoEndCount ?? throw new ArgumentNullException(nameof(twoEndCount));
    }

    public override string DisplayName =>
        "F169: under Heisenberg bonds off the letters the palindrome needs a colour (the colour law)";

    public override string Summary =>
        "palindrome iff a colour perpendicular to every axis carries every field; proven on every connected graph " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("why a bond on one axis drops out",
                summary: "a bond's term at one end is its cross X_ij placed at the other end; for two ends on one axis it " +
                         "points along the axis, where nothing else reaches, so X_ij vanishes, the bond leaves every equation " +
                         "and only the symmetry of its pair remains");
            yield return new InspectableNode("why every connected graph admits only the colours",
                summary: "every other bond keeps its pair in span{m⊗m, p⊗p}; two bonds at one site at an angle neither " +
                         "straight nor right leave nothing; otherwise each bond asks both ends on its line or neither, " +
                         "which leaves at most one word and its complement, and a site's equation is a sum of distinct words from " +
                         "N = 4 on (smaller N read site by site), so only the word with every bond's line at both ends, the colour, survives");
            yield return new InspectableNode("the live witness",
                summary: "inspect --root colourlaw: the far kernel's dimension from the commutator in the Pauli-string basis, " +
                         "fields included, without the split, against the law's prediction");
        }
    }

    private static int Components(int n, IEnumerable<(int A, int B)> edges)
    {
        var parent = Enumerable.Range(0, n).ToArray();
        int Find(int x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; }
        int count = n;
        foreach (var (a, b) in edges)
        {
            int ra = Find(a), rb = Find(b);
            if (ra != rb) { parent[ra] = rb; count--; }
        }
        return count;
    }

    private static long Dot((long X, long Y, long Z) a, (long X, long Y, long Z) b) => checked(a.X * b.X + a.Y * b.Y + a.Z * b.Z);

    private static (long X, long Y, long Z) Cross((long X, long Y, long Z) a, (long X, long Y, long Z) b)
        => checked((a.Y * b.Z - a.Z * b.Y, a.Z * b.X - a.X * b.Z, a.X * b.Y - a.Y * b.X));

    private static (long X, long Y, long Z) Primitive((long X, long Y, long Z) v)
    {
        static long Gcd(long a, long b) { a = Math.Abs(a); b = Math.Abs(b); while (b != 0) (a, b) = (b, a % b); return a; }
        long g = Gcd(Gcd(v.X, v.Y), v.Z);
        if (g == 0) return v;
        long sign = (v.X != 0 ? v.X : v.Y != 0 ? v.Y : v.Z) < 0 ? -1 : 1;
        return (sign * v.X / g, sign * v.Y / g, sign * v.Z / g);
    }
}
