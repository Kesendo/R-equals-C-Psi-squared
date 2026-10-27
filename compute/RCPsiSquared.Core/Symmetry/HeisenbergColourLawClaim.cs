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
/// <para><b>Where it is proven.</b> A leaf's bond and a bond whose ends share an axis drop out of every equation and
/// leave a symmetric pair, so leaves peel and such bonds delete from cycles; the bond equations admit only the
/// colours wherever the reduced core is a single site or a single cycle (trees by induction, rings of every length:
/// a nonzero p⊗p coefficient needs a right hinge, and with every hinge right the slots alternate two letters, which
/// an odd ring cannot carry and an even ring's third slot separates), and wherever the axes are pairwise parallel or
/// perpendicular (the complement connection's Theorems 1 to 3). Where a site keeps three or more neighbours in the
/// reduced core it is measured, not proven.</para></summary>
public sealed class HeisenbergColourLawClaim : Claim
{
    /// <summary>F158: the palindrome holds exactly when the far kernel holds an invertible element; a colour is one.</summary>
    public PalindromeTwoEndCountClaim TwoEndCount { get; }

    /// <summary>F165: its Theorems 1 to 3 settle the axes one global turn makes letters.</summary>
    public PalindromeComplementConnectionClaim ComplementConnection { get; }

    /// <summary>The independent cycles of the reduced core: those of the graph with every bond whose ends share an
    /// axis removed (the number of such bonds that can go while the graph stays connected is their rank in the
    /// cographic matroid, so the order of the steps does not matter).</summary>
    public static int ReducedCoreCycles(int n, IReadOnlyList<(int A, int B)> edges, IReadOnlyList<(long X, long Y, long Z)> axes)
    {
        var kept = edges.Where(e => Primitive(axes[e.A]) != Primitive(axes[e.B])).ToList();
        return kept.Count - n + Components(n, kept);
    }

    /// <summary>True where the law is proven for this graph and these axes: connected, and a reduced core with at most
    /// one independent cycle or axes pairwise parallel or perpendicular.</summary>
    public static bool Proven(int n, IReadOnlyList<(int A, int B)> edges, IReadOnlyList<(long X, long Y, long Z)> axes)
    {
        if (Components(n, edges) != 1) return false;
        if (ReducedCoreCycles(n, edges, axes) <= 1) return true;
        var dirs = axes.Select(Primitive).Distinct().ToList();
        return dirs.All(a => dirs.All(b => Cross(a, b) == (0, 0, 0) || Dot(a, b) == 0));
    }

    /// <summary>The dimension of F158's far kernel the law predicts where it is proven: the admitted colours'
    /// symmetric tensors. N + 1 with all axes along one n and no field (Sym^N of the plane ⊥ n); 1 when one colour is
    /// admitted otherwise; 0 when none is (the axes span space, a field leaves the plane's normal, or under one axis
    /// the fields are not perpendicular to it and parallel to one another).</summary>
    public static int PredictedFarDimension(IReadOnlyList<(long X, long Y, long Z)> axes, IReadOnlyList<(long X, long Y, long Z)> fields)
    {
        var dirs = axes.Select(Primitive).Distinct().ToList();
        var f = fields.Where(v => v != (0, 0, 0)).Select(Primitive).Distinct().ToList();
        if (dirs.All(d => Cross(d, dirs[0]) == (0, 0, 0)))
        {
            if (f.Count > 1 || f.Any(v => Dot(v, dirs[0]) != 0)) return 0;
            return f.Count == 0 ? axes.Count + 1 : 1;
        }
        var normal = dirs.Select(d => Primitive(Cross(dirs[0], d))).First(v => v != (0, 0, 0));
        if (dirs.Any(d => Dot(normal, d) != 0)) return 0;
        return f.All(v => v == normal) ? 1 : 0;
    }

    public HeisenbergColourLawClaim(PalindromeTwoEndCountClaim twoEndCount, PalindromeComplementConnectionClaim complementConnection)
        : base("F169, the colour law: Heisenberg bonds of nonzero weights on a connected graph, every site watched along one " +
               "integer axis, fields of any direction. A product of colours commutes with a bond of coupling matrix W exactly " +
               "when the half-turns about its ends keep W; the far equation splits into bond and field equations; wherever " +
               "the bond equations admit only the colours the palindrome holds exactly when a colour perpendicular to every " +
               "axis carries every field. The bond equations admit only the colours wherever the reduced core (leaves " +
               "removed, bonds whose ends share an axis deleted from cycles) is a single site or a single cycle, and " +
               "wherever the axes are pairwise parallel or perpendicular (F165, Theorems 1 to 3). Open where a site keeps " +
               "three or more neighbours in the reduced core, measured on K4, the theta graph and a five-ring with a chord.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_HEISENBERG_COLOURS_OFF_THE_LETTERS.md + " +
               "docs/ANALYTICAL_FORMULAS.md (F169, and F138's paragraph \"Off the letters\") + " +
               "docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md (F165, Theorems 1 to 3) + " +
               "docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md (F158, section (h): jumps along any axis)")
    {
        TwoEndCount = twoEndCount ?? throw new ArgumentNullException(nameof(twoEndCount));
        ComplementConnection = complementConnection ?? throw new ArgumentNullException(nameof(complementConnection));
    }

    public override string DisplayName =>
        "F169: under Heisenberg bonds off the letters the palindrome needs a colour (the colour law)";

    public override string Summary =>
        "palindrome iff a colour perpendicular to every axis carries every field; proven on reduced cores that are a site " +
        $"or one cycle and at axes one turn makes letters; measured on K4 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("why a bond drops out",
                summary: "a bond's term at one end is its cross X_ij placed at the other end; for a leaf that term is the " +
                         "leaf's whole equation, and for two ends on one axis it points along the axis, where nothing else " +
                         "reaches, so X_ij vanishes either way, the bond leaves every equation and only the symmetry of its " +
                         "pair remains");
            yield return new InspectableNode("why a ring admits only the colours",
                summary: "each site has two terms, so one vanishing p⊗p coefficient kills them all; a nonzero one needs the " +
                         "planes of consecutive axes to meet at a right angle at every site, and then the slots carry two " +
                         "letters alternately, which an odd ring cannot do and which an even ring's third slot separates");
            yield return new InspectableNode("the live witness",
                summary: "inspect --root colourlaw: the far kernel's dimension from the commutator in the Pauli-string basis, " +
                         "fields included, without the split, against the law's prediction, on rows in and beyond the proven class");
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
