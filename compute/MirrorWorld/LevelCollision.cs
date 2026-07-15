namespace MirrorWorld;

/// <summary>The level-collision law (F129), adopted 2026-07-14 from
/// docs/proofs/PROOF_F129_LEVEL_COLLISION_LAW.md + simulations/f129_level_collision_law.py:
/// the six-cosine sibling of Seed's three-cosine world, one weight level up.
///
/// The comb at size n carries the mode energies lam_k = 2 cos(k pi / n), k = 1..n-1
/// (n = N_sites + 1). A triple {k1 &lt; k2 &lt; k3} is CLEAN when no internal pair sums to n
/// (that excludes the trivial balanced-pair family). Its LEVEL is lam_{k1}+lam_{k2}+lam_{k3}.
/// THE LAW: two distinct clean triples with EQUAL levels (zero not required) exist if and
/// only if 3|n (n &gt;= 9) or 10|n (n &gt;= 20); away from both families the level map is
/// INJECTIVE. Sub-law: a colliding pair sharing one mode forces 3|n; the 10|n door is
/// walked only by fully disjoint pairs (the pentagon, unchained from the 15|n that held it
/// in Seed's three-cosine world). Where Seed asks "which triples sum to ZERO" (resonance),
/// this asks "which triples sum to EACH OTHER" (coincidence): resonance is the special case.
///
/// Exactness, Seed's own convention: levels are computed in GF(p) over TWO independent
/// primes p = 1 (mod 2n) with zeta of exact order 2n. One-sided soundness as in the source
/// gate: exact equality implies equality mod every such p, so DISTINCTNESS mod one prime
/// PROVES exact distinctness (injectivity at a non-firing n is proof grade); an equal pair
/// mod both primes at a non-firing n would break the law assertion loudly. The mechanism
/// anchors (n = 15: four rotated R3 cycles; n = 20: the R5 conjugate pair + the zero mode)
/// are checked as root-sum zeros mod both primes.
///
/// THE FAMILY INVENTORY (adopted 2026-07-15 from experiments/F129_FAMILY_INVENTORY.md +
/// docs/proofs/PROOF_F129_FAMILY_INVENTORY_COUNTS.md): every colliding pair decomposes into
/// minimal vanishing pieces of the complete Poonen-Rubinstein/CDK list, and the pieces sort
/// the census into exactly THIRTEEN families A..M, each with a DERIVED closed-form count.
/// The degree counts the freedom, precisely: each free coset label, plus the shared-mode
/// line of a d = 2 family, contributes one factor linear in n -- A two labels, B one label
/// times the s-line (both quadratic); C..H one factor (linear); a rigid sporadic pinned by
/// conjugation symmetry contributes a constant (I..M). The
/// F129 onsets are ZEROS of the inventory (A(6) = B(6) = C(10) = 0), and family M splits
/// 40 + 0 + 60 over the three CDK order-210 types, the middle one structurally impossible.
/// What comes home is the closed arithmetic (doors, counts, the M split); the derivation
/// (label/orbit engines, the s-lemma) stays in the proof, and its section 8 carries the
/// honest code-trust flags.
///
/// The boundary, per Seed's header: the LAW and its COUNTS come home (arithmetic of the
/// comb, no eigensolver); the physical Gram decoupling at collisions (F130, B(tau,sigma) = 0,
/// the K_26/Slater construction) is an eigen-story and stays in the main repo, where it is
/// witnessed exactly (CollisionDecoupling, inspect --root crosstriple). The law's full
/// n &lt;= 210 census also stays with the committed gate (its one named corner family was
/// closed as EMPTY 2026-07-15 via the Poonen-Rubinstein classification; law unconditional).</summary>
public sealed class LevelCollision : GameObject
{
    /// <summary>Comb size n (= N_sites + 1); modes are 1..n-1.</summary>
    public int Ncomb { get; }

    public LevelCollision(World world, int ncomb) : base(world)
    {
        if (ncomb < 5) throw new ArgumentOutOfRangeException(nameof(ncomb), "comb needs n >= 5");
        Ncomb = ncomb;
    }

    public override IReadOnlyList<string> Own => new[] { "levels", "collisions", "inventory" };

    /// <summary>The law's firing predicate: 3|n (n &gt;= 9) or 10|n (n &gt;= 20).</summary>
    public static bool Fires(int n) => (n % 3 == 0 && n >= 9) || (n % 10 == 0 && n >= 20);

    /// <summary>All clean triples 1 &lt;= k1 &lt; k2 &lt; k3 &lt;= n-1 (no internal pair sums to n).</summary>
    public static List<(int K1, int K2, int K3)> CleanTriples(int n)
    {
        var list = new List<(int, int, int)>();
        for (int a = 1; a < n - 1; a++)
            for (int b = a + 1; b < n; b++)
            {
                if (a + b == n) continue;
                for (int c = b + 1; c < n; c++)
                {
                    if (a + c == n || b + c == n) continue;
                    list.Add((a, b, c));
                }
            }
        return list;
    }

    /// <summary>One comb census: clean-triple count, colliding pairs (equal level mod BOTH
    /// primes), split by overlap, plus one example pair (or null). Distinctness is proof
    /// grade; a candidate pair is confirmed by the law itself (see the header).</summary>
    public sealed record CombCensus(
        int N,
        bool Fires,
        int CleanCount,
        int CollidingPairs,
        int DisjointPairs,
        int Overlap1Pairs,
        (int, int, int)? ExampleA,
        (int, int, int)? ExampleB);

    public CombCensus Census() => CensusOf(Ncomb);

    public static CombCensus CensusOf(int n)
    {
        var triples = CleanTriples(n);
        var (p1, lam1) = CyclotomicLambdas(n, 0);
        var (p2, lam2) = CyclotomicLambdas(n, p1);
        // group by the level fingerprint mod both primes
        var groups = new Dictionary<(long, long), List<int>>();
        for (int i = 0; i < triples.Count; i++)
        {
            var (a, b, c) = triples[i];
            var key = ((lam1[a] + lam1[b] + lam1[c]) % p1, (lam2[a] + lam2[b] + lam2[c]) % p2);
            if (!groups.TryGetValue(key, out var bucket)) groups[key] = bucket = new List<int>();
            bucket.Add(i);
        }
        int pairs = 0, disjoint = 0, overlap1 = 0;
        (int, int, int)? exA = null, exB = null;
        foreach (var bucket in groups.Values)
        {
            if (bucket.Count < 2) continue;
            for (int i = 0; i < bucket.Count; i++)
                for (int j = i + 1; j < bucket.Count; j++)
                {
                    var t = triples[bucket[i]];
                    var s = triples[bucket[j]];
                    int shared = SharedModes(t, s);
                    pairs++;
                    if (shared == 0) disjoint++;
                    else if (shared == 1) overlap1++;
                    if (exA is null) { exA = t; exB = s; }
                }
        }
        return new CombCensus(n, Fires(n), triples.Count, pairs, disjoint, overlap1, exA, exB);
    }

    /// <summary>The law over a range: at every non-firing n zero colliding pairs (injectivity,
    /// proof grade via mod-p distinctness), at every firing n at least one; the observed
    /// firing set equals the predicted set exactly.</summary>
    public static bool LawHolds(int lo, int hi)
    {
        for (int n = Math.Max(lo, 5); n <= hi; n++)
            if ((CensusOf(n).CollidingPairs > 0) != Fires(n)) return false;
        return true;
    }

    /// <summary>The sub-law over a range: every overlap-1 colliding pair lives at 3|n
    /// (the 10|n door is exclusively disjoint).</summary>
    public static bool SubLawHolds(int lo, int hi)
    {
        for (int n = Math.Max(lo, 5); n <= hi; n++)
            if (CensusOf(n).Overlap1Pairs > 0 && n % 3 != 0) return false;
        return true;
    }

    // ---- the family inventory (adopted 2026-07-15): thirteen families, thirteen derived
    // ---- closed forms. Pieces per family (weight, CDK ratio order) and doors per the
    // ---- inventory table; counts per the counting proof's label/orbit engines. A family's
    // ---- count is 0 off its door AND at silent multiples below onset (the onsets are zeros
    // ---- of the formulas, never extra conditions).

    /// <summary>The thirteen collision families of the F129 inventory. A = four R3 pieces;
    /// B = zero mode + R3-conjugate-pair; C = zero mode + R5-conjugate-pair (the pentagon);
    /// D = (R5:3R3) weight 8; E = R3-pair + (R5:R3); F = (R5:R3) + conjugate; G = zero +
    /// (R5:R3); H = (R7:R3) weight 8; I = (R7:5R3); J = zero + (R7:3R3); K = (R11:R3);
    /// L = zero + (R7:R5) (the corner-closure's second mechanism); M = the weight-12
    /// order-210 types.</summary>
    public enum CollisionFamily { A, B, C, D, E, F, G, H, I, J, K, L, M }

    /// <summary>The family's door: the divisibility that admits it (each piece's ratio order
    /// must divide 2n).</summary>
    public static int FamilyDoor(CollisionFamily f) => f switch
    {
        CollisionFamily.A => 3,
        CollisionFamily.B => 6,
        CollisionFamily.C => 10,
        CollisionFamily.D => 15,
        CollisionFamily.E => 15,
        CollisionFamily.F => 15,
        CollisionFamily.G => 30,
        CollisionFamily.H => 21,
        CollisionFamily.I => 21,
        CollisionFamily.J => 42,
        CollisionFamily.K => 33,
        CollisionFamily.L => 70,
        CollisionFamily.M => 105,
        _ => throw new ArgumentOutOfRangeException(nameof(f)),
    };

    /// <summary>True for the d = 2 families (the pair shares one mode; weight-8 W). Their
    /// doors all carry the factor 3: the F129 sub-law, visible piece by piece.</summary>
    public static bool FamilySharesAMode(CollisionFamily f) =>
        f is CollisionFamily.B or CollisionFamily.D or CollisionFamily.G or CollisionFamily.H;

    /// <summary>The derived closed-form count of colliding pairs in family f at comb size n
    /// (PROOF_F129_FAMILY_INVENTORY_COUNTS.md sections 2-6). Zero off the door; the onsets
    /// (A and B at n = 6, C at n = 10) are zeros of the formulas themselves.</summary>
    public static long FamilyCount(CollisionFamily f, int n)
    {
        if (n < 5 || n % FamilyDoor(f) != 0) return 0;
        switch (f)
        {
            case CollisionFamily.A: { long k = FloorDiv(n - 9, 6); return (20 * k + 1) * (k + 1); }
            case CollisionFamily.B: { long k = (n - 12) / 6; return 12 * (3 * k + 2) * (k + 1); }
            case CollisionFamily.C: return 2L * (n - 10);
            case CollisionFamily.D: return 12L * (n - 9);
            case CollisionFamily.E: return n % 2 == 1 ? (20L * n - 264) / 3 : (20L * n - 324) / 3;
            case CollisionFamily.F: return n % 2 == 1 ? 10L * n - 149 : 10L * n - 275;
            case CollisionFamily.G: return 6L * (n - 8);
            case CollisionFamily.H: return 6L * (n - 9);
            case CollisionFamily.I: return 60;
            case CollisionFamily.J: return 60;
            case CollisionFamily.K: return 20;
            case CollisionFamily.L: return 20;
            case CollisionFamily.M: return 100;
            default: throw new ArgumentOutOfRangeException(nameof(f));
        }
    }

    /// <summary>Family M's sub-classification over the three CDK order-210 types (the counts
    /// proof section 6, pinned from below by the committed gate's I5): (R7:(R5:2R3)) = 40,
    /// (R7:R3,(R5:R3)) = 0 -- IMPOSSIBLE at every n: two branches of distinct types would
    /// each need the odd gon's one axis-fixed vertex -- and (R7:2R3,R5) = 60.</summary>
    public static (int FannedR5Pair, int MixedBranches, int FixedVertexR5) MSplit => (40, 0, 60);

    /// <summary>The inventory totals at n: the thirteen closed forms summed, split by d.
    /// The from-below tie: Total/Disjoint/Overlap1 equal the census counts exactly (the
    /// thirteen families partition the census, gate I3).</summary>
    public static (long Total, long Disjoint, long Overlap1) InventoryOf(int n)
    {
        long disjoint = 0, overlap1 = 0;
        foreach (CollisionFamily f in Enum.GetValues<CollisionFamily>())
        {
            long c = FamilyCount(f, n);
            if (FamilySharesAMode(f)) overlap1 += c; else disjoint += c;
        }
        return (disjoint + overlap1, disjoint, overlap1);
    }

    static long FloorDiv(long a, long b) => a >= 0 ? a / b : -((-a + b - 1) / b);

    /// <summary>The two mechanism anchors of the proof's paragraph 4, as root-sum zeros mod
    /// both primes: n = 15, (8,12,14) ~ (9,11,13) at a NONZERO shared level whose 12-root
    /// plus/minus set splits into four rotated R3 cycles (step 2n/3 = 10), each summing to
    /// zero; n = 20, (1,7,9) ~ (3,5,10): the R5 conjugate pair (step 2n/5 = 8) plus the zero
    /// mode {zeta^10, zeta^30} = {i, -i}.</summary>
    public static bool AnchorsExact()
    {
        // n = 15
        {
            if (!LevelsCollide(15, (8, 12, 14), (9, 11, 13), requireNonzero: true)) return false;
            int[][] cycles = { new[] { 8, 18, 28 }, new[] { 12, 22, 2 }, new[] { 14, 24, 4 }, new[] { 6, 16, 26 } };
            if (!PiecesPartitionAndVanish(15, new[] { 8, 12, 14, 6, 4, 2 }, cycles)) return false;
        }
        // n = 20
        {
            if (!LevelsCollide(20, (1, 7, 9), (3, 5, 10), requireNonzero: true)) return false;
            int[][] pieces = { new[] { 1, 9, 17, 25, 33 }, new[] { 39, 31, 23, 15, 7 }, new[] { 10, 30 } };
            if (!PiecesPartitionAndVanish(20, new[] { 1, 7, 9, 17, 15, 10 }, pieces)) return false;
        }
        return true;
    }

    static bool LevelsCollide(int n, (int, int, int) t, (int, int, int) s, bool requireNonzero)
    {
        var (p1, lam1) = CyclotomicLambdas(n, 0);
        var (p2, lam2) = CyclotomicLambdas(n, p1);
        long lt1 = (lam1[t.Item1] + lam1[t.Item2] + lam1[t.Item3]) % p1;
        long ls1 = (lam1[s.Item1] + lam1[s.Item2] + lam1[s.Item3]) % p1;
        long lt2 = (lam2[t.Item1] + lam2[t.Item2] + lam2[t.Item3]) % p2;
        long ls2 = (lam2[s.Item1] + lam2[s.Item2] + lam2[s.Item3]) % p2;
        if (lt1 != ls1 || lt2 != ls2) return false;
        return !requireNonzero || lt1 != 0 || lt2 != 0;
    }

    static bool PiecesPartitionAndVanish(int n, int[] six, int[][] pieces)
    {
        int order = 2 * n;
        // the plus/minus exponent multiset of the six-cosine set is {e, 2n - e : e in six}
        var expected = new SortedDictionary<int, int>();
        foreach (int e in six) { Bump(expected, Mod(e, order)); Bump(expected, Mod(-e, order)); }
        var actual = new SortedDictionary<int, int>();
        foreach (var piece in pieces)
            foreach (int e in piece)
                Bump(actual, Mod(e, order));
        if (expected.Count != actual.Count) return false;
        foreach (var (k, v) in expected)
            if (!actual.TryGetValue(k, out int w) || w != v) return false;
        // each piece is a root-sum zero mod both primes
        var (p1, _) = CyclotomicLambdas(n, 0);
        var (p2, _) = CyclotomicLambdas(n, p1);
        foreach (long p in new[] { p1, p2 })
        {
            long zeta = RootOfOrder(order, p);
            if (zeta == 0) return false;
            foreach (var piece in pieces)
            {
                long sum = 0;
                foreach (int e in piece) sum = (sum + ModPow(zeta, Mod(e, order), p)) % p;
                if (sum != 0) return false;
            }
        }
        return true;

        static void Bump(SortedDictionary<int, int> d, int k) =>
            d[k] = d.TryGetValue(k, out int v) ? v + 1 : 1;
        static int Mod(int e, int m) => ((e % m) + m) % m;
    }

    static int SharedModes((int, int, int) t, (int, int, int) s)
    {
        int count = 0;
        Span<int> ts = stackalloc[] { t.Item1, t.Item2, t.Item3 };
        Span<int> ss = stackalloc[] { s.Item1, s.Item2, s.Item3 };
        foreach (int a in ts)
            foreach (int b in ss)
                if (a == b) count++;
        return count;
    }

    // --- the GF(p) atoms, restated per the world's per-object self-containment (Seed has
    // --- its own copies; the local restatement is the sober base's style, drift caught by
    // --- the machine-zero pins on both sides).

    static (long P, long[] Lam) CyclotomicLambdas(int n, long above)
    {
        int order = 2 * n;
        for (long k = Math.Max(above, 1_000_000L) / order + 1; ; k++)
        {
            long p = order * k + 1;
            if (!IsPrime(p)) continue;
            long zeta = RootOfOrder(order, p);
            if (zeta == 0) continue;
            var lam = new long[n];
            for (int m = 1; m < n; m++)
                lam[m] = (ModPow(zeta, m, p) + ModPow(zeta, order - m, p)) % p;
            return (p, lam);
        }
    }

    static long RootOfOrder(int order, long p)
    {
        var qs = PrimeFactors(order);
        for (long x = 2; x < 500; x++)
        {
            long z = ModPow(x, (p - 1) / order, p);
            if (z == 1) continue;
            bool full = true;
            foreach (int q in qs)
                if (ModPow(z, order / q, p) == 1) { full = false; break; }
            if (full) return z;
        }
        return 0;
    }

    static int[] PrimeFactors(int m)
    {
        var qs = new List<int>();
        for (int q = 2; q * q <= m; q++)
            if (m % q == 0) { qs.Add(q); while (m % q == 0) m /= q; }
        if (m > 1) qs.Add(m);
        return qs.ToArray();
    }

    static bool IsPrime(long m)
    {
        if (m < 2) return false;
        long[] bases = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };
        foreach (long b in bases) { if (m % b == 0) return m == b; }
        long d = m - 1; int s = 0;
        while ((d & 1) == 0) { d >>= 1; s++; }
        foreach (long b in bases)
        {
            long x = ModPow(b, d, m);
            if (x == 1 || x == m - 1) continue;
            bool witness = true;
            for (int i = 1; i < s; i++)
            {
                x = x * x % m;
                if (x == m - 1) { witness = false; break; }
            }
            if (witness) return false;
        }
        return true;
    }

    static long ModPow(long b, long e, long p)
    {
        long r = 1;
        b %= p;
        while (e > 0)
        {
            if ((e & 1) == 1) r = r * b % p;
            b = b * b % p;
            e >>= 1;
        }
        return r;
    }
}
