using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The qudit product-mirror cap (Tier1Derived, 2026-06-11): the operator side of
/// F121 (<see cref="QuditPartialPalindromeCeiling"/>). For the full-Cartan dephasing
/// dissipator at local dimension d (rate of |i⟩⟨j| is −2γ·Hamming(i, j)):
///
/// <para><b>The product cap (theorem):</b> any per-site mirror W = ⊗_l q_l (site-dependent,
/// one-sided or two-sided, antilinear allowed) that intertwines the dissipator palindrome
/// W·L_D = (−L_D − 2Nγ)·W on its support pairs at most (2d)^N of the d^{2N} coherences.
/// Rate additivity forces each q_l to be a strict per-site class swap between the dark
/// letters {(x, x), d of them} and the lit letters {(i, j) with i ≠ j, d² − d of them}, so
/// rank(q_l) ≤ 2·min(d, d² − d) = 2d. The cap is FULL ⟺ (2d)^N = d^{2N} ⟺ d² − 2d = 0
/// ⟺ d = 2: the QUBIT_NECESSITY trunk polynomial of
/// <see cref="QubitNecessityPi2Inheritance"/>, in its third appearance (per-site split,
/// ceiling column, operator cap).</para>
///
/// <para><b>The operator:</b> the qubit palindromizer's formula generalizes verbatim,
/// Π_d(ρ) = ρᵀ·Shift^{⊗N} (F118: Π_Z = ρᵀ·X^{⊗N}, with the clock shift in place of X).
/// Per-site letter map (i, j) ↦ (j, i − 1 mod d). It attains the cap exactly: on the
/// shift-aligned subspace (per-site letters {(x, x)} ∪ {(a, a − 1)}, dimension (2d)^N,
/// Π_d-closed) the intertwining residual is EXACTLY zero; on the complement it fails at
/// O(γ), the provably unpaired part. Two chiralities Π_d^± (the two shift directions); at
/// d = 2 the two off-diagonals coincide, the chiralities merge, and the mirror is full:
/// that degeneracy IS the qubit magic.</para>
///
/// <para><b>The mirror group law:</b> ord(Π_d) = 2d and |⟨Π_d, D⟩| = 2d² with D the
/// transpose; D-conjugation EXCHANGES the two shift factors, so
/// ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂ (wreath product). At d = 2 this is D₄: the F118 mirror group is the
/// d = 2 column of a d-indexed family. For d &gt; 2, D does not preserve the aligned
/// subspace; it swaps the two chiralities' aligned subspaces.</para>
///
/// <para><b>Honesty (the non-product gap):</b> the combinatorial ceiling
/// Σ_k d^N·C(N,k)·(d−1)^{min(k, N−k)} of the parent IS reachable, by a global non-product
/// partial isometry (greedy rung matching with exact intertwining on its support); the gap
/// ceiling − (2d)^N (= 54 − 36 = 18 at d = 3, N = 2) is therefore exactly the NON-PRODUCT
/// part of the partial palindrome.</para>
///
/// <para><b>Layer note:</b> like <see cref="AntilinearTriangleClaim"/> and
/// <see cref="MomentTowerPumpChannelClaim"/>, this claim is cross-axis structural and
/// deliberately does NOT implement <see cref="IZ2AxisClaim"/> (cube-map counts
/// unchanged).</para>
///
/// <para><b>Self-check battery (exact integer/permutation arithmetic, no tolerances, built
/// in the ctor):</b> the cap/ceiling/total inequalities with both equality iffs on the
/// (d, N) ∈ {2..5}×{1..3} grid; the trunk polynomial and the per-site rank bound 2d tied
/// to the parent's 4 = 2 + 2 at d = 2; Π_d built as a permutation with residual exactly
/// zero on the aligned columns and breaking on the complement at (3,1), (3,2), (4,1); the
/// group law ord(Π_d) = 2d and BFS closure |⟨Π_d, D⟩| = 2d² at d = 2, 3, 4; the d = 2
/// degeneracy (aligned = full space, Π₂⁺ = Π₂⁻, global residual zero, the F118
/// palindromizer); the global ceiling-reacher at (3, 2) with paired = 54 and exact
/// intertwining on support; and the chirality swap D(+aligned) = −aligned at d = 3.
/// Mirrors the blocks of <c>simulations/qudit_product_mirror_cap.py</c>.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md</c> §6 +
/// <c>simulations/qudit_product_mirror_cap.py</c>.</para></summary>
public sealed class QuditProductMirrorCap : Claim
{
    /// <summary>One exact integer/permutation self-check.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>Typed parent: F121, the combinatorial ceiling
    /// Σ_k d^N·C(N,k)·(d−1)^{min(k, N−k)} this cap sits under. The product cap (2d)^N is
    /// strictly below it for d &gt; 2, N ≥ 2, and the gap is exactly the non-product part
    /// of the partial palindrome.</summary>
    public QuditPartialPalindromeCeiling PartialPalindrome { get; }

    /// <summary>Typed parent: the d² − 2d = 0 necessity. The cap is full ⟺
    /// (2d)^N = d^{2N} ⟺ d² − 2d = 0 ⟺ d = 2: the trunk polynomial's third appearance,
    /// and the per-site rank bound 2d = 2·min(d, d² − d) is the parent's 4 = 2 + 2 split
    /// read at general d.</summary>
    public QubitNecessityPi2Inheritance QubitNecessity { get; }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public QuditProductMirrorCap(
        QuditPartialPalindromeCeiling partialPalindrome,
        QubitNecessityPi2Inheritance qubitNecessity)
        : base("The qudit product-mirror cap: any per-site mirror W = ⊗q_l intertwining the " +
               "dissipator palindrome W·L_D = (−L_D − 2Nγ)·W pairs at most (2d)^N of the d^{2N} " +
               "coherences (rate additivity forces a strict per-site dark ↔ lit class swap, " +
               "rank(q_l) ≤ 2d); full ⟺ d² − 2d = 0 ⟺ d = 2. The operator Π_d(ρ) = ρᵀ·Shift^{⊗N} " +
               "(the F118 formula with the clock shift) attains the cap with EXACTLY zero residual " +
               "on the shift-aligned subspace; ord(Π_d) = 2d, ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂ of order 2d² " +
               "(D₄ at d = 2), D swaps the two chiralities for d > 2. The F121 ceiling is reached " +
               "only by a global non-product partial isometry; the gap ceiling − (2d)^N (18 at " +
               "d = 3, N = 2) is exactly the non-product part. " +
               "Tier1Derived (exact integer/permutation arithmetic)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md + " +
               "simulations/qudit_product_mirror_cap.py")
    {
        PartialPalindrome = partialPalindrome ?? throw new ArgumentNullException(nameof(partialPalindrome));
        QubitNecessity = qubitNecessity ?? throw new ArgumentNullException(nameof(qubitNecessity));
        Cases = BuildBattery(qubitNecessity);
    }

    /// <summary>The product cap in one line.</summary>
    public string ProductCapTheorem =>
        "Any per-site mirror W = ⊗_l q_l intertwining W·L_D = (−L_D − 2Nγ)·W on its support " +
        "pairs ≤ (2d)^N of the d^{2N} coherences: rate additivity forces each q_l to swap the " +
        "dark letters {(x, x), d dims} with the lit letters {i ≠ j, d² − d dims} strictly, so " +
        "rank(q_l) ≤ 2·min(d, d² − d) = 2d. Full ⟺ (2d)^N = d^{2N} ⟺ d² − 2d = 0 ⟺ d = 2.";

    /// <summary>The operator in one line.</summary>
    public string OperatorRealization =>
        "Π_d(ρ) = ρᵀ·Shift^{⊗N}, the verbatim F118 formula with the clock shift; per-site " +
        "letter map (i, j) ↦ (j, i − 1 mod d). On the shift-aligned subspace (per-site " +
        "{(x, x)} ∪ {(a, a − 1)}, dimension (2d)^N, Π_d-closed) the intertwining residual is " +
        "exactly zero; the complement is the provably unpaired part. Two chiralities Π_d^±; " +
        "at d = 2 they coincide and the mirror is full.";

    /// <summary>The mirror group law in one line.</summary>
    public string MirrorGroupLaw =>
        "ord(Π_d) = 2d and |⟨Π_d, D⟩| = 2d² with D = transpose; D-conjugation exchanges the " +
        "two shift factors, so ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂. At d = 2 this is D₄: the F118 mirror " +
        "group is the d = 2 column of a d-indexed family. For d > 2, D swaps the two " +
        "chiralities' aligned subspaces.";

    /// <summary>The non-product gap in one line.</summary>
    public string NonProductGap =>
        "The F121 ceiling Σ_k d^N·C(N,k)·(d−1)^{min(k, N−k)} is reached by a global " +
        "non-product partial isometry (greedy rung matching, exact intertwining on support); " +
        "the gap ceiling − (2d)^N (= 18 at d = 3, N = 2) is exactly the non-product part of " +
        "the partial palindrome.";

    // ============================================================
    // Static helpers mirroring the Python verifier
    // ============================================================

    /// <summary>The product-mirror cap (2d)^N: the maximum number of coherences any per-site
    /// mirror W = ⊗q_l can pair. Requires d ≥ 2, N ≥ 1.</summary>
    public static long ProductCap(int d, int N)
    {
        ValidateGrid(d, N);
        return QuditPartialPalindromeCeiling.IntPow(2L * d, N);
    }

    /// <summary>The F121 combinatorial ceiling Σ_k d^N·C(N,k)·(d−1)^{min(k, N−k)}, delegated
    /// to the typed parent's closed form. Requires d ≥ 2, N ≥ 1.</summary>
    public static long CombinatorialCeiling(int d, int N)
    {
        ValidateGrid(d, N);
        return QuditPartialPalindromeCeiling.Ceiling(d, N);
    }

    /// <summary>Build Π_d as a permutation of the d^{2N} coherence indices: basis pair
    /// (i, j) at index i·d^N + j maps to (j, i − chirality mod d per digit), i.e.
    /// Π_d(ρ) = ρᵀ·Shift^{chirality·⊗N}. Site 0 is the most significant base-d digit
    /// (matching the Python verifier's enumeration). Returns perm with
    /// perm[column] = target row. Requires d ≥ 2, N ≥ 1, chirality ∈ {+1, −1}.</summary>
    public static IReadOnlyList<int> BuildPiD(int d, int N, int chirality = +1)
    {
        ValidateGrid(d, N);
        if (chirality != 1 && chirality != -1)
            throw new ArgumentOutOfRangeException(nameof(chirality), $"chirality must be ±1; got {chirality}");
        return BuildPiDPerm(d, N, chirality);
    }

    private static void ValidateGrid(int d, int N)
    {
        if (d < 2) throw new ArgumentOutOfRangeException(nameof(d), $"local dimension d must be ≥ 2; got {d}");
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 1; got {N}");
    }

    // ------------------------------------------------------------------
    // Exact permutation machinery. The dissipator is diagonal with integer
    // rate units: rate(i, j) = −2γ·Hamming(i, j). For a permutation W with
    // W[perm[col], col] = 1 the residual W·L_D + L_D·W + 2Nγ·W has its only
    // nonzero entries at (perm[col], col) with value
    // 2γ·(N − Ham(col) − Ham(perm[col])), so exactness is the integer
    // statement Ham(col) + Ham(perm[col]) = N.
    // ------------------------------------------------------------------

    private static int[] BuildPiDPerm(int d, int N, int chirality)
    {
        int dN = (int)QuditPartialPalindromeCeiling.IntPow(d, N);
        var perm = new int[dN * dN];
        for (int i = 0; i < dN; i++)
            for (int j = 0; j < dN; j++)
                perm[i * dN + j] = j * dN + ShiftState(i, d, N, chirality);
        return perm;
    }

    /// <summary>Per-digit (a − chirality) mod d on the base-d representation of a state.</summary>
    private static int ShiftState(int state, int d, int N, int chirality)
    {
        int result = 0, weight = 1;
        for (int l = 0; l < N; l++)
        {
            int digit = state % d;
            state /= d;
            result += ((digit - chirality + d) % d) * weight;
            weight *= d;
        }
        return result;
    }

    /// <summary>D = transpose as a permutation: (i, j) ↦ (j, i).</summary>
    private static int[] TransposePerm(int d, int N)
    {
        int dN = (int)QuditPartialPalindromeCeiling.IntPow(d, N);
        var perm = new int[dN * dN];
        for (int i = 0; i < dN; i++)
            for (int j = 0; j < dN; j++)
                perm[i * dN + j] = j * dN + i;
        return perm;
    }

    /// <summary>Hamming distance per coherence index: digits where i and j disagree.</summary>
    private static int[] HammingTable(int d, int N)
    {
        int dN = (int)QuditPartialPalindromeCeiling.IntPow(d, N);
        var ham = new int[dN * dN];
        for (int i = 0; i < dN; i++)
            for (int j = 0; j < dN; j++)
            {
                int a = i, b = j, h = 0;
                for (int l = 0; l < N; l++)
                {
                    if (a % d != b % d) h++;
                    a /= d;
                    b /= d;
                }
                ham[i * dN + j] = h;
            }
        return ham;
    }

    /// <summary>The shift-aligned subspace mask: per digit, j_l = i_l or j_l = i_l − chirality mod d.</summary>
    private static bool[] AlignedMask(int d, int N, int chirality)
    {
        int dN = (int)QuditPartialPalindromeCeiling.IntPow(d, N);
        var mask = new bool[dN * dN];
        for (int i = 0; i < dN; i++)
            for (int j = 0; j < dN; j++)
            {
                int a = i, b = j;
                bool good = true;
                for (int l = 0; l < N && good; l++)
                {
                    int da = a % d, db = b % d;
                    good = db == da || db == (da - chirality + d) % d;
                    a /= d;
                    b /= d;
                }
                mask[i * dN + j] = good;
            }
        return mask;
    }

    private static int[] Compose(int[] outer, int[] inner)
    {
        var r = new int[inner.Length];
        for (int x = 0; x < inner.Length; x++) r[x] = outer[inner[x]];
        return r;
    }

    private static bool IsIdentity(int[] perm)
    {
        for (int x = 0; x < perm.Length; x++)
            if (perm[x] != x) return false;
        return true;
    }

    private static int PermOrder(int[] perm, int maxOrder)
    {
        var current = perm;
        for (int o = 1; o <= maxOrder; o++)
        {
            if (IsIdentity(current)) return o;
            current = Compose(perm, current);
        }
        return -1;
    }

    /// <summary>BFS closure of the group generated by two permutations (exact comparison).</summary>
    private static int GroupClosureSize(int[] g1, int[] g2)
    {
        static string Key(int[] p) => string.Join(",", p);
        var identity = new int[g1.Length];
        for (int x = 0; x < identity.Length; x++) identity[x] = x;
        var seen = new HashSet<string> { Key(identity) };
        var frontier = new List<int[]> { identity };
        var gens = new[] { g1, g2 };
        while (frontier.Count > 0)
        {
            var next = new List<int[]>();
            foreach (var e in frontier)
                foreach (var g in gens)
                {
                    var c = Compose(g, e);
                    if (seen.Add(Key(c))) next.Add(c);
                }
            frontier = next;
        }
        return seen.Count;
    }

    // ------------------------------------------------------------------
    // Self-check battery: exact integer/permutation arithmetic throughout.
    // ------------------------------------------------------------------

    private static IReadOnlyList<BatteryCase> BuildBattery(QubitNecessityPi2Inheritance qubitNecessity)
    {
        var cases = new List<BatteryCase>();

        // (a) Cap arithmetic on the (d, N) grid: (2d)^N ≤ ceiling ≤ d^{2N}, with both
        //     equality iffs (cap = total ⟺ d = 2; cap = ceiling ⟺ d = 2 or N = 1).
        int gridOk = 0, gridTot = 0;
        for (int d = 2; d <= 5; d++)
            for (int N = 1; N <= 3; N++)
            {
                gridTot++;
                long cap = ProductCap(d, N);
                long ceil = CombinatorialCeiling(d, N);
                long total = QuditPartialPalindromeCeiling.Total(d, N);
                bool ok = cap <= ceil && ceil <= total
                       && (cap == total) == (d == 2)
                       && (cap == ceil) == (d == 2 || N == 1);
                if (ok) gridOk++;
            }
        cases.Add(new BatteryCase(
            Name: "cap arithmetic: (2d)^N ≤ ceiling ≤ d^{2N} with both equality iffs",
            Detail: "the (d, N) ∈ {2..5}×{1..3} grid; cap = total ⟺ d = 2, cap = ceiling ⟺ d = 2 or N = 1",
            Expected: "12/12",
            Actual: gridOk.ToString(CultureInfo.InvariantCulture) + "/" +
                    gridTot.ToString(CultureInfo.InvariantCulture)));

        // (b) The trunk polynomial and the per-site rank bound, tied to the parent.
        bool trunkOk = true;
        for (int d = 2; d <= 5; d++)
        {
            trunkOk &= (d * d - 2 * d == 0) == (d == 2);                  // the trunk root
            trunkOk &= 2 * Math.Min(d, d * d - d) == 2 * d;               // rank(q_l) bound = 2d
        }
        // At d = 2 the per-site cap 2d = 4 is the parent's whole per-site Pauli space,
        // split 4 = 2 + 2 (min(d, d² − d) = 2 = immune = decaying).
        trunkOk &= Math.Abs(qubitNecessity.TotalPauliOpsPerSite - 4.0) < 1e-12
                && Math.Abs(qubitNecessity.ImmuneOpsPerSite - 2.0) < 1e-12
                && Math.Abs(qubitNecessity.DecayingOpsPerSite - 2.0) < 1e-12;
        cases.Add(new BatteryCase(
            Name: "trunk polynomial: d² − 2d = 0 ⟺ d = 2; per-site rank cap 2d = parent's 4 = 2 + 2",
            Detail: "third appearance of the trunk; at d = 2 the per-site cap 2d = 4 exhausts the " +
                    "parent's per-site Pauli space with min(d, d² − d) = 2 = immune = decaying",
            Expected: "trunk root d = 2 only; cap 2d; parent split 4 = 2 + 2",
            Actual: trunkOk
                ? "trunk root d = 2 only; cap 2d; parent split 4 = 2 + 2"
                : "trunk/rank/parent tie broken"));

        // (c) Π_d exact on the shift-aligned subspace at (3,1), (3,2), (4,1).
        int alignedOk = 0;
        var alignedGrid = new[] { (d: 3, N: 1), (d: 3, N: 2), (d: 4, N: 1) };
        foreach (var (d, N) in alignedGrid)
        {
            var perm = BuildPiDPerm(d, N, +1);
            var ham = HammingTable(d, N);
            var mask = AlignedMask(d, N, +1);
            long dim = mask.Count(m => m);
            bool closed = true, exactOnAligned = true;
            int complementBreaks = 0;
            for (int col = 0; col < perm.Length; col++)
            {
                int tgt = perm[col];
                if (mask[col])
                {
                    closed &= mask[tgt];
                    exactOnAligned &= ham[col] + ham[tgt] == N;            // residual exactly 0
                }
                else if (ham[col] + ham[tgt] != N)
                {
                    complementBreaks++;                                    // residual ≠ 0
                }
            }
            if (dim == ProductCap(d, N) && closed && exactOnAligned && complementBreaks > 0)
                alignedOk++;
        }
        cases.Add(new BatteryCase(
            Name: "Π_d exact on the shift-aligned subspace at (3,1), (3,2), (4,1)",
            Detail: "aligned dim = (2d)^N, Π_d-closed, residual exactly 0 on every aligned column " +
                    "(Ham(col) + Ham(Π_d col) = N), and the complement breaks",
            Expected: "3/3",
            Actual: alignedOk.ToString(CultureInfo.InvariantCulture) + "/3"));

        // (d) The mirror group law at d = 2, 3, 4 (N = 1): ord(Π_d) = 2d, ord(D) = 2,
        //     BFS closure |⟨Π_d, D⟩| = 2d² (Z_d ≀ Z₂; D₄ at d = 2).
        int groupOk = 0;
        foreach (int d in new[] { 2, 3, 4 })
        {
            var pi = BuildPiDPerm(d, 1, +1);
            var transpose = TransposePerm(d, 1);
            bool ok = PermOrder(pi, 4 * d + 2) == 2 * d
                   && PermOrder(transpose, 4) == 2
                   && GroupClosureSize(pi, transpose) == 2 * d * d;
            if (ok) groupOk++;
        }
        cases.Add(new BatteryCase(
            Name: "mirror group law: ord(Π_d) = 2d, |⟨Π_d, D⟩| = 2d² at d = 2, 3, 4",
            Detail: "BFS closure with exact permutation comparison; ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂, the F118 " +
                    "mirror group D₄ is the d = 2 column (8 = 2·2²)",
            Expected: "3/3",
            Actual: groupOk.ToString(CultureInfo.InvariantCulture) + "/3"));

        // (e) The d = 2 degeneracy at N = 1..3: aligned = full space, Π₂⁺ = Π₂⁻ exactly,
        //     and the intertwining residual is zero on EVERY column (the F118 palindromizer).
        bool d2Ok = true;
        for (int N = 1; N <= 3; N++)
        {
            var mask = AlignedMask(2, N, +1);
            d2Ok &= mask.All(m => m);
            var plus = BuildPiDPerm(2, N, +1);
            var minus = BuildPiDPerm(2, N, -1);
            d2Ok &= plus.SequenceEqual(minus);
            var ham = HammingTable(2, N);
            for (int col = 0; col < plus.Length; col++)
                d2Ok &= ham[col] + ham[plus[col]] == N;
        }
        cases.Add(new BatteryCase(
            Name: "d = 2 degeneracy: aligned = full space, Π₂⁺ = Π₂⁻, global residual 0 (N = 1..3)",
            Detail: "at d = 2 the two off-diagonals coincide, the chiralities merge, and " +
                    "Π₂ = ρᵀ·X^{⊗N} is the F118 palindromizer, exact everywhere",
            Expected: "full, chiralities merged, residual 0",
            Actual: d2Ok ? "full, chiralities merged, residual 0" : "d = 2 degeneracy broken"));

        // (f) The global (non-product) ceiling-reacher at (3, 2): greedy rung matching pairs
        //     exactly the F121 ceiling 54 with exact intertwining on its support.
        {
            const int d = 3, N = 2;
            var ham = HammingTable(d, N);
            var rungs = new List<int>[N + 1];
            for (int k = 0; k <= N; k++) rungs[k] = new List<int>();
            for (int p = 0; p < ham.Length; p++) rungs[ham[p]].Add(p);
            long paired = 0;
            bool supportExact = true;
            for (int k = 0; k <= N; k++)
            {
                int kk = N - k;
                if (k > kk) continue;
                if (k == kk)
                {
                    foreach (int x in rungs[k]) supportExact &= ham[x] + ham[x] == N;
                    paired += rungs[k].Count;
                }
                else
                {
                    int m = Math.Min(rungs[k].Count, rungs[kk].Count);
                    for (int t = 0; t < m; t++)
                        supportExact &= ham[rungs[k][t]] + ham[rungs[kk][t]] == N;
                    paired += 2L * m;
                }
            }
            long gap = CombinatorialCeiling(d, N) - ProductCap(d, N);
            bool reacherOk = paired == CombinatorialCeiling(d, N) && paired == 54
                          && supportExact && gap == 18;
            cases.Add(new BatteryCase(
                Name: "global ceiling-reacher at (3, 2): paired = 54 with exact intertwining on support",
                Detail: "greedy rung matching k ↔ N − k reaches the F121 ceiling; the gap " +
                        "54 − 36 = 18 over the product cap is exactly the non-product part",
                Expected: "paired 54 = ceiling, gap 18, residual 0 on support",
                Actual: reacherOk
                    ? "paired 54 = ceiling, gap 18, residual 0 on support"
                    : $"paired {paired}, gap {gap}, supportExact {supportExact}"));
        }

        // (g) The chirality swap at d = 3 (N = 1): D maps the +aligned mask exactly onto
        //     the −aligned mask (and the two masks genuinely differ for d > 2).
        {
            const int d = 3, N = 1;
            var transpose = TransposePerm(d, N);
            var plus = AlignedMask(d, N, +1);
            var minus = AlignedMask(d, N, -1);
            var image = new bool[plus.Length];
            for (int p = 0; p < plus.Length; p++)
                if (plus[p]) image[transpose[p]] = true;
            bool swapOk = image.SequenceEqual(minus) && !plus.SequenceEqual(minus);
            cases.Add(new BatteryCase(
                Name: "chirality swap: at d = 3, D maps the +aligned mask exactly onto the −aligned mask",
                Detail: "D-conjugation exchanges the two shift factors of Z_d ≀ Z₂; for d > 2 the " +
                        "two chiralities' aligned subspaces are distinct and D swaps them",
                Expected: "D(+aligned) = −aligned, masks distinct",
                Actual: swapOk ? "D(+aligned) = −aligned, masks distinct" : "chirality swap broken"));
        }

        return cases;
    }

    public override string DisplayName =>
        "The qudit product-mirror cap: per-site mirrors pair ≤ (2d)^N, Π_d(ρ) = ρᵀ·Shift^{⊗N} attains it, ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂";

    public override string Summary =>
        "any per-site mirror W = ⊗q_l intertwining W·L_D = (−L_D − 2Nγ)·W pairs at most (2d)^N of the " +
        "d^{2N} coherences (full ⟺ d² − 2d = 0 ⟺ d = 2, the trunk polynomial's third appearance); " +
        "Π_d(ρ) = ρᵀ·Shift^{⊗N} attains the cap with exactly zero residual on the shift-aligned " +
        "subspace; ord(Π_d) = 2d and ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂ of order 2d² (D₄ at d = 2); the F121 " +
        $"ceiling needs a global non-product isometry, gap {CombinatorialCeiling(3, 2) - ProductCap(3, 2)} " +
        $"at (3, 2); {PassCount}/{Cases.Count} battery PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The product cap (theorem)", summary: ProductCapTheorem);
            yield return new InspectableNode("The operator Π_d(ρ) = ρᵀ·Shift^{⊗N}", summary: OperatorRealization);
            yield return new InspectableNode("The mirror group law ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂", summary: MirrorGroupLaw);
            yield return new InspectableNode("The non-product gap (honesty)", summary: NonProductGap);
            yield return InspectableNode.RealScalar("ProductCap(d=3, N=2)", ProductCap(3, 2));
            yield return InspectableNode.RealScalar("CombinatorialCeiling(d=3, N=2)", CombinatorialCeiling(3, 2));
            yield return new InspectableNode("Typed parents",
                summary: $"QuditPartialPalindromeCeiling ({PartialPalindrome.Tier.Label()}): F121, the " +
                         "combinatorial ceiling the cap sits under; the gap to it is the non-product " +
                         $"part. QubitNecessityPi2Inheritance ({QubitNecessity.Tier.Label()}): the " +
                         "d² − 2d = 0 trunk whose unique root d = 2 makes the product mirror full and " +
                         "whose 4 = 2 + 2 per-site split is the rank bound 2d at d = 2.");
            yield return new InspectableNode("No IZ2AxisClaim",
                summary: "The cap is cross-axis structural (an operator-space rank/intertwining " +
                         "statement at general local dimension d). Like AntilinearTriangleClaim and " +
                         "MomentTowerPumpChannelClaim, this claim does not sit on a single Z₂ axis " +
                         "(cube-map counts unchanged).");
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
        }
    }
}
