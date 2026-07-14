namespace RCPsiSquared.Core.Numerics;

/// <summary>F130 recomputed exactly at inspect time (the collision-decoupling law,
/// <c>docs/proofs/PROOF_F130_COLLISION_DECOUPLING.md</c>): for any two distinct mode triples
/// with equal levels S(τ) = S(σ) (zero NOT required), the whole cross block vanishes,
/// B(τ, σ) = 0, i.e. both half-Gram numbers U± are zero. The exact objects: with
/// ŝ(a) := ζ^a − ζ^{−a} = 2i·sin(aπ/n) ∈ ℤ[ζ_2n], the scaled Slater determinant
/// D̂_τ(x,y,z) = det[ŝ(kᵢ·col)] and the scaled half-Grams
///
/// <code>Ê± (τ,σ) = Σ_{c ≷ b} Ĝ_τ(b,c)·Ĝ_σ(b,c),   Ĝ_τ(b,c) = D̂_τ(b−1,b,c) + D̂_τ(b,b+1,c)</code>
///
/// satisfy Ê± = κ·U± with the common nonzero real scalar κ = (2i)⁶·‖D_τ‖‖D_σ‖, so
/// U± = 0 ⟺ Ê± = 0 exactly in ℤ[ζ_2n] (norms and scaling irrelevant for the zero claim),
/// and Lemma 3's U⁺ = ε·U⁻ carries over as Ê⁺ = ε·Ê⁻, an unconditional pin checked on every
/// pair including the controls. The chain walls (x = 0, x = n) are zeroed automatically by
/// ŝ (ζ^{kn} = (−1)^k). All exact via <see cref="CyclotomicRing"/>: no floats, no primes,
/// no sampling — a second method beside the committed float gate
/// <c>simulations/f130_collision_decoupling.py</c>.
///
/// <para>SCOPE, honestly: the four-cell derivation (Lemmas 3+4, the free-angle assembly (D)
/// + F128, the two-magnon lemma, the removable limit) lives in the proof doc; this witness
/// pins the exact vanishing at named pairs covering every cell of that proof, one unequal-level
/// nonzero control per level-sensitive cell, and the exact Lemma-3 sign identity.</para></summary>
public static class CollisionDecoupling
{
    /// <summary>One evaluated pair. For collision rows both zero flags must be true; for
    /// control rows both must be false; <paramref name="Lemma3SignOk"/> must hold on ALL rows.</summary>
    public sealed record PairResult(
        string Label,
        int N,
        int Eps,
        bool LevelsEqual,
        bool PlusZero,
        bool MinusZero,
        bool Lemma3SignOk);

    public sealed record Report(
        IReadOnlyList<PairResult> Collisions,
        IReadOnlyList<PairResult> Controls,
        bool AllCollisionsDecouple,
        bool AllControlsNonzero,
        bool AllLemma3SignsOk,
        bool LevelFlagsAsExpected);

    /// <summary>Named collision pairs, at least one per cell of the F130 proof (several for
    /// cell 2). The cell-1 row has UNEQUAL levels by design: that cell is level-free.</summary>
    private static readonly (string Label, int N, int[] T, int[] S, bool ExpectEqual)[] CollisionPairs =
    {
        ("cell 2 (disjoint ε=−1, level ≠ 0): n=15 (8,12,14)~(9,11,13)", 15, new[] { 8, 12, 14 }, new[] { 9, 11, 13 }, true),
        ("cell 2 (pentagon door 10|n): n=20 (1,7,9)~(3,5,10)", 20, new[] { 1, 7, 9 }, new[] { 3, 5, 10 }, true),
        ("cell 2 (τ non-clean, cleanliness not a hypothesis): n=12 (1,2,10)~(3,5,6)", 12, new[] { 1, 2, 10 }, new[] { 3, 5, 6 }, true),
        ("cell 4 (overlap-1 ε=−1): n=12 (6,10,11)~(7,9,10)", 12, new[] { 6, 10, 11 }, new[] { 7, 9, 10 }, true),
        ("cell 3 (overlap-1 ε=+1): n=18 (9,16,17)~(11,13,16)", 18, new[] { 9, 16, 17 }, new[] { 11, 13, 16 }, true),
        ("cell 2 (the resonant special case, level = 0): n=9 (2,4,8)~(1,5,7)", 9, new[] { 2, 4, 8 }, new[] { 1, 5, 7 }, true),
        ("cell 1 (level-free, levels UNEQUAL): n=12 (1,2,3)~(4,5,7)", 12, new[] { 1, 2, 3 }, new[] { 4, 5, 7 }, false),
    };

    /// <summary>Unequal-level nonzero controls, one per level-sensitive cell (the
    /// discrimination guard: without these a vacuously-zero Ê would pass silently).</summary>
    private static readonly (string Label, int N, int[] T, int[] S)[] ControlPairs =
    {
        ("cell 2 control (disjoint ε=−1, unequal): n=12 (1,2,3)~(4,5,6)", 12, new[] { 1, 2, 3 }, new[] { 4, 5, 6 }),
        ("cell 3 control (overlap-1 ε=+1, unequal): n=12 (1,2,3)~(1,4,5)", 12, new[] { 1, 2, 3 }, new[] { 1, 4, 5 }),
        ("cell 4 control (overlap-1 ε=−1, unequal): n=12 (1,2,3)~(1,4,6)", 12, new[] { 1, 2, 3 }, new[] { 1, 4, 6 }),
    };

    public static int Eps(int[] t, int[] s)
    {
        int sum = 0;
        foreach (int k in t) sum += k;
        foreach (int l in s) sum += l;
        return sum % 2 == 0 ? 1 : -1;
    }

    /// <summary>The scaled Slater determinant D̂_τ(x,y,z) = det[ŝ(kᵢ·col)] ∈ ℤ[ζ_2n].</summary>
    public static long[] DHat(int n, int[] triple, int x, int y, int z)
    {
        int m = 2 * n;
        int[] cols = { x, y, z };
        var cell = new long[3, 3][];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            {
                var v = CyclotomicRing.Zero(m);
                CyclotomicRing.AddRootPower(v, m, (long)triple[i] * cols[j], 1);
                CyclotomicRing.AddRootPower(v, m, -(long)triple[i] * cols[j], -1);
                cell[i, j] = v;
            }
        // cofactor expansion along the first row
        var det = CyclotomicRing.Zero(m);
        det = CyclotomicRing.Add(det, CyclotomicRing.Multiply(cell[0, 0], Minor(cell, m, 1, 2, 1, 2), m));
        det = CyclotomicRing.Add(det, CyclotomicRing.Negate(CyclotomicRing.Multiply(cell[0, 1], Minor(cell, m, 1, 2, 0, 2), m)));
        det = CyclotomicRing.Add(det, CyclotomicRing.Multiply(cell[0, 2], Minor(cell, m, 1, 2, 0, 1), m));
        return det;

        static long[] Minor(long[,][] c, int m, int r1, int r2, int c1, int c2) =>
            CyclotomicRing.Add(
                CyclotomicRing.Multiply(c[r1, c1], c[r2, c2], m),
                CyclotomicRing.Negate(CyclotomicRing.Multiply(c[r1, c2], c[r2, c1], m)));
    }

    /// <summary>Ĝ_τ(b,c) = D̂_τ(b−1,b,c) + D̂_τ(b,b+1,c); walls handled by ŝ automatically.</summary>
    public static long[] GHat(int n, int[] triple, int b, int c) =>
        CyclotomicRing.Add(DHat(n, triple, b - 1, b, c), DHat(n, triple, b, b + 1, c));

    /// <summary>The scaled half-Grams (Ê⁺, Ê⁻) over the two triangular regions c &gt; b, c &lt; b.</summary>
    public static (long[] Plus, long[] Minus) EHat(int n, int[] t, int[] s)
    {
        int m = 2 * n, sites = n - 1;
        var plus = CyclotomicRing.Zero(m);
        var minus = CyclotomicRing.Zero(m);
        for (int b = 1; b <= sites; b++)
        {
            for (int c = 1; c <= sites; c++)
            {
                if (c == b) continue;
                var prod = CyclotomicRing.Multiply(GHat(n, t, b, c), GHat(n, s, b, c), m);
                if (c > b) plus = CyclotomicRing.Add(plus, prod);
                else minus = CyclotomicRing.Add(minus, prod);
            }
        }
        return (plus, minus);
    }

    /// <summary>Evaluate one pair: the two exact zero flags, the exact level-equality flag,
    /// and the unconditional Lemma-3 sign identity Ê⁺ = ε·Ê⁻.</summary>
    public static PairResult Evaluate(string label, int n, int[] t, int[] s)
    {
        var (plus, minus) = EHat(n, t, s);
        int eps = Eps(t, s);
        bool signOk = CyclotomicRing.AreEqual(plus, eps == 1 ? minus : CyclotomicRing.Negate(minus));
        bool levelsEqual = CyclotomicRing.AreEqual(
            LevelCollisionCensus.LevelVector(n, t[0], t[1], t[2]),
            LevelCollisionCensus.LevelVector(n, s[0], s[1], s[2]));
        return new PairResult(label, n, eps, levelsEqual,
                              CyclotomicRing.IsZero(plus), CyclotomicRing.IsZero(minus), signOk);
    }

    /// <summary>The full fixed-pair certification: every collision row vanishes exactly, every
    /// control row is nonzero, every Lemma-3 sign holds, and every level flag matches the row's
    /// design (equal for the collision rows except the level-free cell-1 row; unequal controls).</summary>
    public static Report Analyze()
    {
        var collisions = new List<PairResult>();
        bool levelFlagsOk = true;
        foreach (var (label, n, t, s, expectEqual) in CollisionPairs)
        {
            var r = Evaluate(label, n, t, s);
            collisions.Add(r);
            if (r.LevelsEqual != expectEqual) levelFlagsOk = false;
        }
        var controls = new List<PairResult>();
        foreach (var (label, n, t, s) in ControlPairs)
        {
            var r = Evaluate(label, n, t, s);
            controls.Add(r);
            if (r.LevelsEqual) levelFlagsOk = false;
        }
        bool allDecouple = collisions.TrueForAll(r => r.PlusZero && r.MinusZero);
        bool allControlsNonzero = controls.TrueForAll(r => !r.PlusZero && !r.MinusZero);
        bool allSigns = collisions.TrueForAll(r => r.Lemma3SignOk) && controls.TrueForAll(r => r.Lemma3SignOk);
        return new Report(collisions, controls, allDecouple, allControlsNonzero, allSigns, levelFlagsOk);
    }
}
