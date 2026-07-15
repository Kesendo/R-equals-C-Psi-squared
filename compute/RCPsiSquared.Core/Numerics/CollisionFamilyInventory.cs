namespace RCPsiSquared.Core.Numerics;

/// <summary>The F129 family inventory recomputed exactly at inspect time (the counts, derived:
/// <c>docs/proofs/PROOF_F129_FAMILY_INVENTORY_COUNTS.md</c> +
/// <c>experiments/F129_FAMILY_INVENTORY.md</c>): every colliding pair of clean triples
/// decomposes into minimal vanishing pieces of the complete Poonen-Rubinstein/CDK list,
/// sorting the census into exactly THIRTEEN families A..M, each with a derived closed-form
/// count. This surface holds the closed arithmetic (doors, counts, the d-split, the M split)
/// and ties its SUMS live against the exact ℤ[ζ_2n] census of <see cref="LevelCollisionCensus"/>:
/// for each n the closed forms must reproduce the census's colliding-pair count, split by
/// overlap (disjoint vs shared-one-mode; two shared modes with equal levels are impossible,
/// since λ is injective on 1..n−1 the third modes would coincide, PROOF_F129 §2's d = 1 case).
///
/// <para>SCOPE, honestly: over the live range n ≤ <see cref="LiveMaxN"/> the census exercises
/// ELEVEN of the thirteen families (A..K; L's door is 70, M's is 105); the n = 70 capstone,
/// run as part of <see cref="Analyze"/>, pins L (the corner-closure's second mechanism,
/// zero + (R₇:R₅)) beside the pentagon C. Family M and its 40 + 0 + 60 split over the three
/// CDK order-210 types (<see cref="MSplit"/>; the middle type is impossible, the odd gon has
/// one axis-fixed vertex) ride the committed gate <c>simulations/f129_family_inventory.py</c>
/// (its I5 rebuilds the census M W-sets from committed substitution recipes at n = 105 and
/// 210) — here <see cref="MSplit"/> is an adopted constant with only the internal sum tie
/// to <see cref="Count"/>(M). The live tie certifies the SUMS (total + d-split); per-family
/// MEMBERSHIP (the first-fit piece decomposition and the labels) stays with the gate's
/// I2/I3/I5. The n = 105 capstone (m = 210, seven families co-firing, total 8858) is pinned
/// by the test suite via <see cref="CensusAt"/>, not run at inspect time.</para></summary>
public static class CollisionFamilyInventory
{
    /// <summary>The live tie bound (matches <see cref="LevelCollisionCensus.LiveMaxN"/>).</summary>
    public const int LiveMaxN = 60;

    /// <summary>The thirteen collision families of the inventory. A = four R₃ pieces;
    /// B = zero mode + R₃-conjugate-pair; C = zero mode + R₅-conjugate-pair (the pentagon);
    /// D = (R₅:3R₃) weight 8; E = R₃-pair + (R₅:R₃); F = (R₅:R₃) + conjugate; G = zero +
    /// (R₅:R₃); H = (R₇:R₃) weight 8; I = (R₇:5R₃); J = zero + (R₇:3R₃); K = (R₁₁:R₃);
    /// L = zero + (R₇:R₅); M = the weight-12 order-210 types.</summary>
    public enum Family { A, B, C, D, E, F, G, H, I, J, K, L, M }

    /// <summary>The family's door: the divisibility that admits it (each piece's ratio
    /// order must divide 2n).</summary>
    public static int Door(Family f) => f switch
    {
        Family.A => 3,
        Family.B => 6,
        Family.C => 10,
        Family.D => 15,
        Family.E => 15,
        Family.F => 15,
        Family.G => 30,
        Family.H => 21,
        Family.I => 21,
        Family.J => 42,
        Family.K => 33,
        Family.L => 70,
        Family.M => 105,
        _ => throw new ArgumentOutOfRangeException(nameof(f)),
    };

    /// <summary>True for the d = 2 families (the pair shares one mode; weight-8 W). Their
    /// doors all carry the factor 3: the F129 sub-law, visible piece by piece.</summary>
    public static bool SharesAMode(Family f) =>
        f is Family.B or Family.D or Family.G or Family.H;

    /// <summary>The derived closed-form count of colliding pairs in family f at comb size n
    /// (the counts proof §§2-6). Zero off the door; the F129 onsets (A and B at n = 6,
    /// C at n = 10) are zeros of the formulas themselves — A needs FLOOR division for that
    /// (C# / truncates toward zero and would give A(6) = 1).</summary>
    public static long Count(Family f, int n)
    {
        if (n < 5 || n % Door(f) != 0) return 0;
        switch (f)
        {
            case Family.A: { long k = FloorDiv(n - 9, 6); return (20 * k + 1) * (k + 1); }
            case Family.B: { long k = (n - 12) / 6; return 12 * (3 * k + 2) * (k + 1); }
            case Family.C: return 2L * (n - 10);
            case Family.D: return 12L * (n - 9);
            case Family.E: return n % 2 == 1 ? (20L * n - 264) / 3 : (20L * n - 324) / 3;
            case Family.F: return n % 2 == 1 ? 10L * n - 149 : 10L * n - 275;
            case Family.G: return 6L * (n - 8);
            case Family.H: return 6L * (n - 9);
            case Family.I: return 60;
            case Family.J: return 60;
            case Family.K: return 20;
            case Family.L: return 20;
            case Family.M: return 100;
            default: throw new ArgumentOutOfRangeException(nameof(f));
        }
    }

    /// <summary>Family M's sub-classification over the three CDK order-210 types (the counts
    /// proof §6): (R₇:(R₅:2R₃)) = 40, (R₇:R₃,(R₅:R₃)) = 0 — impossible at every n: two
    /// branches of distinct types would each need the odd gon's one axis-fixed vertex — and
    /// (R₇:2R₃,R₅) = 60. ADOPTED constants: the from-below reconstruction is the committed
    /// gate's I5; the only in-process tie is the sum to <see cref="Count"/>(M).</summary>
    public static (int FannedR5Pair, int MixedBranches, int FixedVertexR5) MSplit => (40, 0, 60);

    /// <summary>The closed-form sums at n, split by d. <paramref name="countsOverride"/>
    /// exists ONLY for the discrimination test (a perturbed closed form must break the tie);
    /// pass null for the real inventory.</summary>
    public static (long Total, long Disjoint, long Overlap1) ClosedFormsOf(
        int n, Func<Family, int, long>? countsOverride = null)
    {
        Func<Family, int, long> count = countsOverride ?? Count;
        long disjoint = 0, overlap1 = 0;
        foreach (Family f in Enum.GetValues<Family>())
        {
            long c = count(f, n);
            if (SharesAMode(f)) overlap1 += c; else disjoint += c;
        }
        return (disjoint + overlap1, disjoint, overlap1);
    }

    /// <summary>The exact census side at one n: colliding pairs (equal ℤ[ζ_2n] level vectors,
    /// byte-equality = ring equality), split by shared-mode count. Throws if a colliding pair
    /// shares two modes (impossible: λ injective forces the third modes equal — a trip means
    /// the reduction of PROOF_F129 §2 broke, never bin it silently).</summary>
    public static (long Total, long Disjoint, long Overlap1) CensusAt(int n)
    {
        // deliberately re-enumerates rather than reusing LevelCollisionCensus.Analyze's pass
        // (that one keeps only counts, this one needs the triples per bucket for the overlap
        // split); the duplication doubles the census share of inspect time, accepted.
        var buckets = new Dictionary<string, List<(int A, int B, int C)>>();
        foreach (var t in LevelCollisionCensus.CleanTriples(n))
        {
            string key = CyclotomicRing.Key(LevelCollisionCensus.LevelVector(n, t.K1, t.K2, t.K3));
            if (!buckets.TryGetValue(key, out var list)) buckets[key] = list = new List<(int, int, int)>();
            list.Add((t.K1, t.K2, t.K3));
        }
        long total = 0, disjoint = 0, overlap1 = 0;
        foreach (var bucket in buckets.Values)
        {
            if (bucket.Count < 2) continue;
            for (int i = 0; i < bucket.Count; i++)
                for (int j = i + 1; j < bucket.Count; j++)
                {
                    int shared = SharedModes(bucket[i], bucket[j]);
                    total++;
                    if (shared == 0) disjoint++;
                    else if (shared == 1) overlap1++;
                    else throw new InvalidOperationException(
                        $"colliding pair at n = {n} shares {shared} modes: the d = 1 impossibility of PROOF_F129 §2 broke");
                }
        }
        return (total, disjoint, overlap1);
    }

    /// <summary>What <see cref="Analyze"/> recomputes. The inventory holds on the range iff
    /// <paramref name="AllRowsTied"/>, <paramref name="Capstone70Tied"/> and
    /// <paramref name="MSplitSumsToCountM"/> are all true.</summary>
    public sealed record Report(
        int MaxN,
        int RowsChecked,
        bool AllRowsTied,
        bool Capstone70Tied,
        bool MSplitSumsToCountM);

    /// <summary>The live tie over 5 ≤ n ≤ <paramref name="maxN"/> (total AND d-split at every
    /// row, firing and non-firing n alike — the formulas vanish at the silent doors by
    /// themselves) plus the n = 70 capstone (L's door) and the M-split sum. The n = 105
    /// capstone stays with the tests and the committed gate.</summary>
    public static Report Analyze(int maxN = LiveMaxN, Func<Family, int, long>? countsOverride = null)
    {
        int rows = 0;
        bool allTied = true;
        for (int n = 5; n <= maxN; n++)
        {
            rows++;
            if (CensusAt(n) != ClosedFormsOf(n, countsOverride)) allTied = false;
        }
        bool capstone = CensusAt(70) == ClosedFormsOf(70, countsOverride);
        var (fanned, middle, fixedVertex) = MSplit;
        bool mSplit = fanned + middle + fixedVertex == (countsOverride ?? Count)(Family.M, 105);
        return new Report(maxN, rows, allTied, capstone, mSplit);
    }

    private static int SharedModes((int A, int B, int C) t, (int A, int B, int C) s)
    {
        int count = 0;
        Span<int> ts = stackalloc[] { t.A, t.B, t.C };
        Span<int> ss = stackalloc[] { s.A, s.B, s.C };
        foreach (int a in ts)
            foreach (int b in ss)
                if (a == b) count++;
        return count;
    }

    private static long FloorDiv(long a, long b) => a >= 0 ? a / b : -((-a + b - 1) / b);
}
