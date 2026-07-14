namespace RCPsiSquared.Core.Numerics;

/// <summary>The F128 flip lemma, recomputed exactly over ℤ at inspect time
/// (docs/proofs/PROOF_F128_FLIP_SUM_FACTORIZATION.md §3, gate G1 of
/// <c>simulations/f128_flip_sum_factorization.py</c>): the integer trigonometric polynomial
///
/// <code>cos s · B · V_a V_b P̃,   B = Σ_u sin 2x_u − 2 sin 2s,   P̃ = Π sin((a_i−b_j)/2)</code>
///
/// (8640 monomials in half-angle monomial units t_u = e^{i x_u/2}, integer coefficients after
/// stripping the constant 2(2i)¹⁶) is annihilated by the (ℤ/2)⁶ signed character sum
/// Π_u (1 − flip_u): the totally-odd projection returns the ZERO polynomial. This is the
/// engine of the F128 factorization 𝔉 = −(e₁−f₁)²·𝒪[cos s·cot s·V_a V_b/P] and of the sharper
/// locus {e₁ = f₁}. Proof (B) is also recomputed, both halves: the 28 power-sum-shifted
/// exponent sets of the two Weyl alternants a_{M₁}, a_{M₂} all die structurally (a zero
/// exponent, a repeat, or a ±pair), AND each shifted alternant's odd projection is verified
/// to be the zero polynomial (the 720-permutation expansion, exact). Pure 64-bit integer
/// arithmetic, no floats, no randomness; independent of the Python gate (own dictionary
/// representation, own multiplication order).</summary>
public static class FlipLemma
{
    /// <summary>What <see cref="Analyze"/> recomputes. <paramref name="NumeratorMonomials"/> must
    /// be 8640, <paramref name="SurvivingMonomials"/> 0, <paramref name="ShiftedSets"/> 28.</summary>
    public sealed record Report(
        int NumeratorMonomials,
        int SurvivingMonomials,
        bool ProjectorSelfTestOk,
        int ShiftedSets,
        bool AllShiftedSetsDie,
        bool AllShiftedProjectionsVanish);

    // Exponents live in [−10, 10] per variable; pack (e + Offset) in 6 bits per variable.
    private const int Offset = 16;
    private const int Bits = 6;
    private const ulong VarMask = (1UL << Bits) - 1;

    private static ulong Pack(ReadOnlySpan<int> e)
    {
        ulong k = 0;
        for (int u = 0; u < 6; u++)
        {
            int shifted = e[u] + Offset;
            if (shifted < 0 || shifted > (int)VarMask)
                throw new InvalidOperationException($"exponent {e[u]} outside the packing window");
            k |= (ulong)shifted << (Bits * u);
        }
        return k;
    }

    private static int Exponent(ulong key, int u) => (int)((key >> (Bits * u)) & VarMask) - Offset;

    private static ulong FlipVar(ulong key, int u)
    {
        int e = Exponent(key, u);
        ulong cleared = key & ~(VarMask << (Bits * u));
        return cleared | (ulong)(-e + Offset) << (Bits * u);
    }

    private static void AddTo(Dictionary<ulong, long> target, ulong key, long coeff)
    {
        if (target.TryGetValue(key, out long v))
        {
            v += coeff;
            if (v == 0) target.Remove(key);
            else target[key] = v;
        }
        else if (coeff != 0) target[key] = coeff;
    }

    private static Dictionary<ulong, long> Multiply(Dictionary<ulong, long> a, (int[] E, long C)[] b)
    {
        var r = new Dictionary<ulong, long>(a.Count * b.Length);
        Span<int> e = stackalloc int[6];
        foreach (var (ka, ca) in a)
            foreach (var (eb, cb) in b)
            {
                for (int u = 0; u < 6; u++) e[u] = Exponent(ka, u) + eb[u];
                AddTo(r, Pack(e), ca * cb);
            }
        return r;
    }

    /// <summary>2i·sin((x_u − x_v)/2) = t_u/t_v − t_v/t_u as a two-term factor.</summary>
    private static (int[] E, long C)[] TwoTerm(int u, int v)
    {
        var e1 = new int[6]; e1[u] = 1; e1[v] = -1;
        var e2 = new int[6]; e2[u] = -1; e2[v] = 1;
        return new[] { (e1, 1L), (e2, -1L) };
    }

    private static Dictionary<ulong, long> BuildNumerator()
    {
        var p = new Dictionary<ulong, long> { [Pack(stackalloc int[6])] = 1 };
        for (int i = 0; i < 3; i++)                       // V_a (pairs within a)
            for (int j = i + 1; j < 3; j++)
                p = Multiply(p, TwoTerm(i, j));
        for (int i = 0; i < 3; i++)                       // V_b (pairs within b)
            for (int j = i + 1; j < 3; j++)
                p = Multiply(p, TwoTerm(3 + i, 3 + j));
        for (int i = 0; i < 3; i++)                       // P̃ (cross pairs a_i − b_j)
            for (int j = 0; j < 3; j++)
                p = Multiply(p, TwoTerm(i, 3 + j));

        var plus1 = new int[6]; var minus1 = new int[6];
        for (int u = 0; u < 6; u++) { plus1[u] = 1; minus1[u] = -1; }
        p = Multiply(p, new[] { (plus1, 1L), (minus1, 1L) });      // 2 cos s

        var bTerms = new List<(int[], long)>();                    // 2i·B
        for (int u = 0; u < 6; u++)
        {
            var e4 = new int[6]; e4[u] = 4;
            var em4 = new int[6]; em4[u] = -4;
            bTerms.Add((e4, 1L));
            bTerms.Add((em4, -1L));
        }
        var plus2 = new int[6]; var minus2 = new int[6];
        for (int u = 0; u < 6; u++) { plus2[u] = 2; minus2[u] = -2; }
        bTerms.Add((plus2, -2L));
        bTerms.Add((minus2, 2L));
        return Multiply(p, bTerms.ToArray());
    }

    /// <summary>Π_u (1 − flip_u), the 64·𝒪 signed character sum, applied variable by variable.</summary>
    private static Dictionary<ulong, long> OddProject(Dictionary<ulong, long> p)
    {
        for (int u = 0; u < 6; u++)
        {
            var r = new Dictionary<ulong, long>(p.Count * 2, comparer: null);
            foreach (var (k, c) in p) AddTo(r, k, c);
            foreach (var (k, c) in p) AddTo(r, FlipVar(k, u), -c);
            p = r;
        }
        return p;
    }

    private static bool ProjectorSelfTest()
    {
        // a generic monomial must project to exactly 64 distinct ±1 terms, not to zero
        var single = new Dictionary<ulong, long> { [Pack(new[] { 1, 2, 3, 4, 5, 6 })] = 1 };
        var proj = OddProject(single);
        if (proj.Count != 64) return false;
        foreach (long c in proj.Values)
            if (c != 1 && c != -1) return false;
        return true;
    }

    /// <summary>The 28 shifted exponent sets of proof (B): ±2 on one z-exponent (from
    /// p₂(z) − p₂(1/z)), ±1 on all (from the sheet term 2(Z − 1/Z)), applied to
    /// M₁ = (3,2,1,0,−1,−2) and M₂ = (2,1,0,−1,−2,−3).</summary>
    private static List<int[]> ShiftedExponentSets()
    {
        int[][] bases = { new[] { 3, 2, 1, 0, -1, -2 }, new[] { 2, 1, 0, -1, -2, -3 } };
        var sets = new List<int[]>();
        foreach (var m in bases)
        {
            for (int k = 0; k < 6; k++)
                foreach (int d in new[] { 2, -2 })
                {
                    var s = (int[])m.Clone();
                    s[k] += d;
                    sets.Add(s);
                }
            foreach (int d in new[] { 1, -1 })
            {
                var s = new int[6];
                for (int k = 0; k < 6; k++) s[k] = m[k] + d;
                sets.Add(s);
            }
        }
        return sets;
    }

    /// <summary>a_M = det[z_u^{m_k}] as an exponent dictionary in half-angle units (t-exponent
    /// = 2m), the 720-permutation expansion with the permutation sign.</summary>
    private static Dictionary<ulong, long> Alternant(int[] m)
    {
        var p = new Dictionary<ulong, long>(720);
        var perm = new int[6];
        var e = new int[6];
        void Recurse(int depth, uint usedMask, int sign)
        {
            if (depth == 6)
            {
                for (int u = 0; u < 6; u++) e[u] = 2 * m[perm[u]];
                AddTo(p, Pack(e), sign);
                return;
            }
            for (int k = 0; k < 6; k++)
            {
                if ((usedMask >> k & 1) != 0) continue;
                perm[depth] = k;
                // each choice that skips over an unused smaller index contributes inversions
                int inversions = 0;
                for (int t = 0; t < k; t++)
                    if ((usedMask >> t & 1) == 0) inversions++;
                Recurse(depth + 1, usedMask | (1u << k), inversions % 2 == 0 ? sign : -sign);
            }
        }
        Recurse(0, 0, 1);
        return p;
    }

    /// <summary>An exponent set kills its odd Weyl numerator det[z^m − z^{−m}] when it contains
    /// a zero (zero column), a repeat (the alternant itself is 0), or a ±pair (two columns equal
    /// up to sign).</summary>
    private static bool Dies(int[] m)
    {
        for (int i = 0; i < 6; i++)
        {
            if (m[i] == 0) return true;
            for (int j = i + 1; j < 6; j++)
                if (m[i] == m[j] || m[i] == -m[j]) return true;
        }
        return false;
    }

    /// <summary>Recompute the flip lemma from scratch: build the 8640-monomial numerator over ℤ,
    /// apply the signed character sum, and report the survivor count (must be 0), plus the
    /// projector self-test and the 28 shifted-set deaths.</summary>
    public static Report Analyze()
    {
        var numerator = BuildNumerator();
        var projected = OddProject(numerator);
        var sets = ShiftedExponentSets();
        bool allDie = true, allVanish = true;
        foreach (var m in sets)
        {
            if (!Dies(m)) allDie = false;
            if (OddProject(Alternant(m)).Count != 0) allVanish = false;
        }
        return new Report(numerator.Count, projected.Count, ProjectorSelfTest(), sets.Count,
                          allDie, allVanish);
    }
}
