namespace RCPsiSquared.Core.Numerics;

/// <summary>F129 recomputed exactly at inspect time (the level-collision law,
/// <c>docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md</c>): for the comb of XX-chain mode
/// energies cos(kπ/n), distinct CLEAN triples (no internal pair sums to n) with equal
/// level S(τ) = Σ cos(kᵢπ/n) exist only at 3|n (n ≥ 9) or 10|n (n ≥ 20). The exact object
/// is the level VECTOR 2S(τ) = Σ (ζ^{kᵢ} + ζ^{−kᵢ}) ∈ ℤ[ζ_2n], a free-basis integer vector
/// (<see cref="CyclotomicRing"/>), so byte-equality of vectors IS equality of levels over ℚ(ζ):
/// distinct vectors at a non-firing n are a PROOF of injectivity at that n, and an equal pair
/// at a firing n is an EXACT collision, no floats and no mod-p filter involved.
///
/// <para>SCOPE, honestly: this census runs n ≤ <see cref="LiveMaxN"/> live. The law's full
/// n ≤ 210 census, the Lam-Leung structural proof, and the one named corner family stay with
/// the committed gate <c>simulations/f129_level_collision_law.py</c> and the proof doc; this
/// witness layer is an independent exact re-derivation of the law's content on its own range,
/// not a re-proof of the law.</para></summary>
public static class LevelCollisionCensus
{
    /// <summary>The live census bound (the Python gate certifies n ≤ 210).</summary>
    public const int LiveMaxN = 60;

    /// <summary>What <see cref="Analyze"/> recomputes. The law holds on the range iff
    /// <paramref name="AllNonFiringInjective"/>, <paramref name="AllFiringCollide"/> and
    /// <paramref name="LawEquality"/> are all true.</summary>
    public sealed record Report(
        int MaxN,
        int CleanTriplesChecked,
        int NonFiringChecked,
        bool AllNonFiringInjective,
        int FiringChecked,
        bool AllFiringCollide,
        bool LawEquality,
        bool AnchorsExact);

    /// <summary>The predicted firing set of the law: 3|n (n ≥ 9) or 10|n (n ≥ 20).</summary>
    public static bool Fires(int n) => (n % 3 == 0 && n >= 9) || (n % 10 == 0 && n >= 20);

    /// <summary>All clean triples 1 ≤ k₁ &lt; k₂ &lt; k₃ ≤ n−1 (no internal pair sums to n).</summary>
    public static IEnumerable<(int K1, int K2, int K3)> CleanTriples(int n)
    {
        for (int a = 1; a < n - 1; a++)
            for (int b = a + 1; b < n; b++)
            {
                if (a + b == n) continue;
                for (int c = b + 1; c < n; c++)
                {
                    if (a + c == n || b + c == n) continue;
                    yield return (a, b, c);
                }
            }
    }

    /// <summary>The exact level vector 2S(τ) = Σᵢ (ζ^{kᵢ} + ζ^{2n−kᵢ}) ∈ ℤ[ζ_2n], reduced.</summary>
    public static long[] LevelVector(int n, int k1, int k2, int k3)
    {
        int m = 2 * n;
        var acc = CyclotomicRing.Zero(m);
        foreach (int k in stackalloc[] { k1, k2, k3 })
        {
            CyclotomicRing.AddRootPower(acc, m, k, 1);
            CyclotomicRing.AddRootPower(acc, m, -k, 1);
        }
        return acc;
    }

    /// <summary>The exact census over 5 ≤ n ≤ <paramref name="maxN"/> plus the two mechanism
    /// anchors, all in ℤ[ζ_2n]. <paramref name="fireOverride"/> exists ONLY for the
    /// discrimination test (a deliberately wrong firing predicate must break LawEquality);
    /// pass null for the real law.</summary>
    public static Report Analyze(int maxN = LiveMaxN, Func<int, bool>? fireOverride = null)
    {
        Func<int, bool> fires = fireOverride ?? Fires;
        int cleanChecked = 0, nonFiring = 0, firing = 0;
        bool allInjective = true, allCollide = true, lawEquality = true;
        for (int n = 5; n <= maxN; n++)
        {
            var groups = new Dictionary<string, int>();
            bool collision = false;
            foreach (var (a, b, c) in CleanTriples(n))
            {
                cleanChecked++;
                string key = CyclotomicRing.Key(LevelVector(n, a, b, c));
                if (groups.TryGetValue(key, out int count))
                {
                    groups[key] = count + 1;
                    collision = true;
                }
                else
                {
                    groups[key] = 1;
                }
            }
            if (fires(n))
            {
                firing++;
                if (!collision) allCollide = false;
            }
            else
            {
                nonFiring++;
                if (collision) allInjective = false;
            }
            if (collision != fires(n)) lawEquality = false;
        }
        return new Report(maxN, cleanChecked, nonFiring, allInjective, firing, allCollide,
                          lawEquality, CheckAnchors());
    }

    /// <summary>The two mechanism anchors of the proof's §4, exact (port of the Python
    /// gate's G5): n = 15, (8,12,14) ~ (9,11,13) at a NONZERO shared level, whose 12-root
    /// ± exponent set partitions into four rotated R₃ cycles (step 2n/3 = 10), each summing
    /// to zero; and n = 20, (1,7,9) ~ (3,5,10), whose set is an R₅ conjugate pair (step
    /// 2n/5 = 8) plus the zero mode {ζ^{10}, ζ^{30}} = {i, −i}, each summing to zero.</summary>
    private static bool CheckAnchors()
    {
        // n = 15: equal nonzero level, four R3 cycles partition the ± exponent multiset.
        {
            int n = 15, m = 30;
            var vt = LevelVector(n, 8, 12, 14);
            var vs = LevelVector(n, 9, 11, 13);
            if (!CyclotomicRing.AreEqual(vt, vs) || CyclotomicRing.IsZero(vt)) return false;
            // six-cosine set {8,12,14} ∪ {15−9, 15−11, 15−13} = {8,12,14,6,4,2}, ± exponents mod 30
            int[] six = { 8, 12, 14, 6, 4, 2 };
            int[][] cycles = { new[] { 8, 18, 28 }, new[] { 12, 22, 2 }, new[] { 14, 24, 4 }, new[] { 6, 16, 26 } };
            if (!PartitionMatches(six, m, cycles)) return false;
            foreach (var cyc in cycles)
                if (!RootSumIsZero(m, cyc)) return false;
        }
        // n = 20: R5 conjugate pair + zero mode.
        {
            int n = 20, m = 40;
            var vt = LevelVector(n, 1, 7, 9);
            var vs = LevelVector(n, 3, 5, 10);
            if (!CyclotomicRing.AreEqual(vt, vs) || CyclotomicRing.IsZero(vt)) return false;
            int[] six = { 1, 7, 9, 17, 15, 10 };
            int[][] pieces = { new[] { 1, 9, 17, 25, 33 }, new[] { 39, 31, 23, 15, 7 }, new[] { 10, 30 } };
            if (!PartitionMatches(six, m, pieces)) return false;
            foreach (var piece in pieces)
                if (!RootSumIsZero(m, piece)) return false;
        }
        return true;
    }

    private static bool PartitionMatches(int[] six, int m, int[][] pieces)
    {
        // the ± exponent multiset of the six-cosine set is {e, m−e : e ∈ six}
        var expected = new SortedDictionary<int, int>();
        foreach (int e in six)
        {
            Bump(expected, ((e % m) + m) % m);
            Bump(expected, ((-e % m) + m) % m);
        }
        var actual = new SortedDictionary<int, int>();
        foreach (var piece in pieces)
            foreach (int e in piece)
                Bump(actual, ((e % m) + m) % m);
        if (expected.Count != actual.Count) return false;
        foreach (var (k, v) in expected)
            if (!actual.TryGetValue(k, out int w) || w != v) return false;
        return true;

        static void Bump(SortedDictionary<int, int> d, int k) =>
            d[k] = d.TryGetValue(k, out int v) ? v + 1 : 1;
    }

    private static bool RootSumIsZero(int m, int[] exponents)
    {
        var acc = CyclotomicRing.Zero(m);
        foreach (int e in exponents)
            CyclotomicRing.AddRootPower(acc, m, e, 1);
        return CyclotomicRing.IsZero(acc);
    }
}
