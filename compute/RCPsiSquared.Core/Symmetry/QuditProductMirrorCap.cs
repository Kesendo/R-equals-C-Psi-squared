using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The qudit product-mirror cap (Tier1Derived, 2026-06-11): the operator side of
/// F121 (<see cref="QuditPartialPalindromeCeiling"/>). For the full-Cartan dephasing
/// dissipator at local dimension d (rate of |i⟩⟨j| is −2γ·Hamming(i, j)):
///
/// <para><b>The product cap (theorem):</b> let W = ⊗_l q_l be any per-site mirror
/// (site-dependent, one-sided or two-sided, antilinear allowed) that intertwines the dissipator
/// palindrome W·L_D = (−L_D − 2Nγ)·W. Write each q_l in blocks between the dark letters
/// {(x, x), d of them} and the lit letters {(i, j) with i ≠ j, d² − d of them}, and let a block
/// carry lit grade c = (lit out) + (lit in) ∈ {0, 1, 2}. The identity asks Σ_l c_l = N on every
/// nonzero product of blocks; varying one site with the others fixed forces every nonzero block
/// of q_l to share one grade c_l. The ranks per grade are r(0) = d, r(1) = 2d, r(2) = d² − d,
/// and Σ_l c_l = N forces #(c = 0) = #(c = 2) = m, so W pairs at most
/// <code>
///   P(d, N) = max_m (2d)^(N−2m) · (d³ − d²)^m
/// </code>
/// of the d^{2N} coherences. A (0, 2) pair of sites beats a (1, 1) pair iff d³ − d² &gt; 4d²,
/// i.e. iff d ≥ 6 (tie at d = 5): for d ≤ 5 the optimum is the strict per-site dark ↔ lit swap
/// and P = (2d)^N at every N; for d ≥ 6 it is not (P_dark⊗P_lit at d = 6, N = 2 has rank
/// 180 &gt; 144). P is full ⟺ d² − 2d = 0 ⟺ d = 2 at every d: the QUBIT_NECESSITY trunk
/// polynomial of <see cref="QubitNecessityPi2Inheritance"/>, in its third appearance (per-site
/// split, ceiling column, operator cap).</para>
///
/// <para><b>The operator:</b> the qubit palindromizer's formula generalizes verbatim,
/// Π_d(ρ) = ρᵀ·Shift^{⊗N} (F118: Π_Z = ρᵀ·X^{⊗N}, with the clock shift in place of X).
/// Per-site letter map (i, j) ↦ (j, i − 1 mod d). The full Π_d is a permutation of rank
/// d^{2N}; restricted to the shift-aligned subspace (per-site letters {(x, x)} ∪ {(a, a − 1)},
/// dimension (2d)^N, Π_d-closed) its intertwining residual is EXACTLY zero, so Π_d P_aligned
/// attains the cap for every d ≤ 5; on the complement this realization fails at O(γ). Two
/// chiralities Π_d^± (the two shift directions); at d = 2 the two off-diagonals coincide, the
/// chiralities merge, and the mirror is full: that degeneracy IS the qubit magic.</para>
///
/// <para><b>The mirror group law:</b> ord(Π_d) = 2d and |⟨Π_d, D⟩| = 2d² with D the
/// transpose; D-conjugation EXCHANGES the two shift factors, so
/// ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂ (wreath product). At d = 2 this is D₄: the F118 mirror group is the
/// d = 2 column of a d-indexed family. For d &gt; 2, D does not preserve the aligned
/// subspace; it swaps the two chiralities' aligned subspaces.</para>
///
/// <para><b>The non-product part:</b> the combinatorial ceiling
/// Σ_k d^N·C(N,k)·(d−1)^{min(k, N−k)} of the parent IS reachable, by a global non-product
/// partial isometry (greedy rung matching with exact intertwining on its support), so the gap
/// ceiling − P(d, N) (= 54 − 36 = 18 at d = 3, N = 2) is exactly the NON-PRODUCT part of the
/// partial palindrome, positive for every d ≥ 3, N ≥ 2 (divide by d^N: 2^N = Σ C(N,k) against
/// Σ C(N,k)(d−1)^{min(k,N−k)} for d ≤ 5; 2^{N−2m}(d−1)^m against at least C(N,m)(d−1)^m plus a
/// positive term, m = ⌊N/2⌋, for d ≥ 6; PROOF §6). At (d, N) = (3, 2), (4, 2), (3, 3) a translation-invariant mirror
/// ([W, T] = 0) reaches the whole ceiling (54, 128, 378): the rank of an integer-coefficient TI
/// intertwiner mod p equals the ceiling, and rank_p ≤ rank ≤ ceiling makes that exact. Whether
/// TI reaches the ceiling at every (d, N) is open.</para>
///
/// <para><b>Layer note:</b> like <see cref="AntilinearTriangleClaim"/> and
/// <see cref="MomentTowerPumpChannelClaim"/>, this claim is cross-axis structural and
/// deliberately does NOT implement <see cref="IZ2AxisClaim"/> (cube-map counts
/// unchanged).</para>
///
/// <para><b>Self-check battery (exact integer, permutation and GF(p) arithmetic, no
/// tolerances, built in the ctor):</b> the shift-rank/cap/ceiling/total inequalities with the
/// equality iffs on the (d, N) ∈ {2..5}×{1..3} grid; the trunk polynomial tied to the parent's
/// 4 = 2 + 2 at d = 2; the product lemma (integer c-type factors reach rank d, 2d, d² − d with
/// residual exactly zero, a factor mixing grades breaks the identity); the closed form P(d, N)
/// against the enumeration of every grade pattern, the d ≥ 6 threshold and the d = 6 projector
/// P_dark⊗P_lit (rank 180, residual 0); Π_d exact on the aligned columns and breaking on the
/// complement at (3,1), (3,2), (4,1); the group law at d = 2, 3, 4; the d = 2 degeneracy; the
/// global ceiling-reacher at (3, 2); the chirality swap at d = 3; the non-product part
/// ceiling − P; and the TI rank mod p at (2, 2), (3, 2), (4, 2) (the (3, 3) rank is a test).
/// Mirrors the blocks of <c>simulations/qudit_product_mirror_cap.py</c> and
/// <c>simulations/qudit_ti_intermediate.py</c>.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md</c> §6 (the cap) and §7
/// (the translation-invariant non-product part) +
/// <c>simulations/qudit_product_mirror_cap.py</c> +
/// <c>simulations/qudit_ti_intermediate.py</c>.</para></summary>
public sealed class QuditProductMirrorCap : Claim
{
    /// <summary>One exact integer/permutation self-check.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    /// <summary>Typed parent: F121, the combinatorial ceiling
    /// Σ_k d^N·C(N,k)·(d−1)^{min(k, N−k)} this cap sits under. The product cap P(d, N) is
    /// strictly below it for d &gt; 2, N ≥ 2, and the gap is exactly the non-product part
    /// of the partial palindrome.</summary>
    public QuditPartialPalindromeCeiling PartialPalindrome { get; }

    /// <summary>Typed parent: the d² − 2d = 0 necessity. The cap is full ⟺
    /// P(d, N) = d^{2N} ⟺ d² − 2d = 0 ⟺ d = 2: the trunk polynomial's third appearance,
    /// and the per-site swap rank 2d = 2·min(d, d² − d) is the parent's 4 = 2 + 2 split
    /// read at general d.</summary>
    public QubitNecessityPi2Inheritance QubitNecessity { get; }

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public QuditProductMirrorCap(
        QuditPartialPalindromeCeiling partialPalindrome,
        QubitNecessityPi2Inheritance qubitNecessity)
        : base("The qudit product-mirror cap: any per-site mirror W = ⊗q_l intertwining the " +
               "dissipator palindrome W·L_D = (−L_D − 2Nγ)·W carries one lit grade c_l ∈ {0, 1, 2} " +
               "per site (rank ≤ d, 2d, d² − d) with #(c = 0) = #(c = 2), so it pairs at most " +
               "P(d, N) = max_m (2d)^(N−2m)·(d³ − d²)^m of the d^{2N} coherences. P = (2d)^N for " +
               "d ≤ 5 (the strict per-site dark ↔ lit swap; tie at d = 5) and above it for d ≥ 6 " +
               "(P_dark⊗P_lit at d = 6, N = 2: rank 180 > 144); full ⟺ d² − 2d = 0 ⟺ d = 2. " +
               "The operator Π_d(ρ) = ρᵀ·Shift^{⊗N} (the F118 formula with the clock shift) attains " +
               "(2d)^N with EXACTLY zero residual on the shift-aligned subspace; ord(Π_d) = 2d, " +
               "⟨Π_d, D⟩ ≅ Z_d ≀ Z₂ of order 2d² (D₄ at d = 2), D swaps the two chiralities for " +
               "d > 2. The F121 ceiling is reached only by a non-product partial isometry; the gap " +
               "ceiling − P (18 at d = 3, N = 2) is exactly the non-product part, and a " +
               "translation-invariant mirror reaches the whole ceiling at (3, 2), (4, 2), (3, 3) " +
               "(rank mod p = ceiling). " +
               "Tier1Derived (exact integer, permutation and GF(p) arithmetic)",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md + " +
               "simulations/qudit_product_mirror_cap.py + " +
               "simulations/qudit_ti_intermediate.py")
    {
        PartialPalindrome = partialPalindrome ?? throw new ArgumentNullException(nameof(partialPalindrome));
        QubitNecessity = qubitNecessity ?? throw new ArgumentNullException(nameof(qubitNecessity));
        Cases = BuildBattery(qubitNecessity);
    }

    /// <summary>The product cap in one line.</summary>
    public string ProductCapTheorem =>
        "Any per-site mirror W = ⊗_l q_l intertwining W·L_D = (−L_D − 2Nγ)·W has every nonzero " +
        "block of q_l at one lit grade c_l ∈ {0, 1, 2} (the identity asks Σ c_l = N on every block " +
        "product, and varying one site pins its grade), with rank ≤ d, 2d, d² − d; Σ c_l = N " +
        "forces #(c = 0) = #(c = 2) = m, so rank W ≤ P(d, N) = max_m (2d)^(N−2m)·(d³ − d²)^m. " +
        "A (0, 2) pair beats a (1, 1) pair iff d³ − d² > 4d² ⟺ d ≥ 6: for d ≤ 5 the optimum is " +
        "the strict dark ↔ lit swap, P = (2d)^N; at d = 6, N = 2, P_dark⊗P_lit has rank " +
        "180 > 144. Full ⟺ P = d^{2N} ⟺ d² − 2d = 0 ⟺ d = 2.";

    /// <summary>The operator in one line.</summary>
    public string OperatorRealization =>
        "Π_d(ρ) = ρᵀ·Shift^{⊗N}, the verbatim F118 formula with the clock shift; per-site " +
        "letter map (i, j) ↦ (j, i − 1 mod d). The full Π_d is a permutation of rank d^{2N}; " +
        "on the shift-aligned subspace (per-site {(x, x)} ∪ {(a, a − 1)}, dimension (2d)^N, " +
        "Π_d-closed) the intertwining residual is exactly zero, so Π_d P_aligned attains the " +
        "cap for d ≤ 5, and on its complement this realization fails at O(γ). Two chiralities " +
        "Π_d^±; at d = 2 they coincide and the mirror is full.";

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
        "no per-site product reaches past P(d, N), so the gap ceiling − P (= 18 at d = 3, " +
        "N = 2) is exactly the non-product part of the partial palindrome.";

    /// <summary>The translation-invariant reach in one line.</summary>
    public string TranslationInvariantReach =>
        "At (d, N) = (3, 2), (4, 2), (3, 3) the non-product part is recovered ENTIRELY by a " +
        "translation-invariant (non-product) mirror: an integer-coefficient [W, T] = 0 " +
        "intertwiner has rank mod p equal to the full ceiling (54, 128, 378), strictly above " +
        "the product cap (36, 64, 216), and rank_p ≤ rank ≤ ceiling makes the equality exact. " +
        "At these three cases there is no intermediate layer: product, then translation-invariant " +
        "= ceiling. Whether translation invariance reaches the ceiling at every (d, N) is open.";

    // ============================================================
    // Static helpers mirroring the Python verifier
    // ============================================================

    /// <summary>Rank of one per-site factor of lit grade c: r(0) = d (dark → dark),
    /// r(1) = 2d (dark ↔ lit), r(2) = d² − d (lit → lit). Requires d ≥ 2, c ∈ {0, 1, 2}.</summary>
    public static long GradeRank(int d, int c)
    {
        if (d < 2) throw new ArgumentOutOfRangeException(nameof(d), $"local dimension d must be ≥ 2; got {d}");
        return c switch
        {
            0 => d,
            1 => 2L * d,
            2 => checked((long)d * d - d),
            _ => throw new ArgumentOutOfRangeException(nameof(c), $"lit grade must be 0, 1 or 2; got {c}"),
        };
    }

    /// <summary>The product-mirror cap P(d, N) = max_m (2d)^(N−2m)·(d³ − d²)^m: the maximum
    /// number of coherences any per-site mirror W = ⊗q_l can pair. Equals (2d)^N for every N
    /// when d ≤ 5 and exceeds it for d ≥ 6, N ≥ 2. Requires d ≥ 2, N ≥ 1.</summary>
    public static long ProductCap(int d, int N)
    {
        ValidateGrid(d, N);
        long best = 0;
        for (int m = 0; 2 * m <= N; m++)
        {
            long value = checked(QuditPartialPalindromeCeiling.IntPow(GradeRank(d, 1), N - 2 * m)
                * QuditPartialPalindromeCeiling.IntPow(checked(GradeRank(d, 0) * GradeRank(d, 2)), m));
            if (value > best) best = value;
        }
        return best;
    }

    /// <summary>The shift-aligned rank (2d)^N: the rank of Π_d P_aligned, the explicit
    /// product mirror that attains the cap for d ≤ 5. Requires d ≥ 2, N ≥ 1.</summary>
    public static long ShiftAlignedRank(int d, int N)
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

    /// <summary>The non-product part of the partial palindrome: ceiling − P(d, N), the number
    /// of coherence pairs no per-site product mirror can reach. Zero iff the product mirror is
    /// already the ceiling, i.e. d = 2 or N = 1. At (3, 2), (4, 2), (3, 3) a translation-invariant
    /// mirror recovers it entirely (see <see cref="TranslationInvariantRankModP"/>).
    /// Requires d ≥ 2, N ≥ 1.</summary>
    public static long NonProductPart(int d, int N) => CombinatorialCeiling(d, N) - ProductCap(d, N);

    /// <summary>The rank over GF(p), p = 1 000 000 007, of one translation-invariant palindrome
    /// intertwiner at (d, N): entries allowed where Hamming(out) + Hamming(in) = N (the
    /// intertwiner maps rung h to rung N − h), one integer coefficient per orbit of the joint
    /// cyclic site shift T, fixed by a hash of the orbit. rank_p ≤ rank over Q ≤ the maximal TI
    /// rank ≤ the ceiling, so a return value equal to the ceiling proves that translation
    /// invariance reaches it. Cost O(d^{6N}); intended for d^{2N} ≤ 729.</summary>
    public static int TranslationInvariantRankModP(int d, int N)
    {
        ValidateGrid(d, N);
        int dN = checked((int)QuditPartialPalindromeCeiling.IntPow(d, N));
        int dim = checked(dN * dN);
        var ham = HammingTable(d, N);
        var shifted = new int[dim];                                  // T on coherence indices
        for (int i = 0; i < dN; i++)
            for (int j = 0; j < dN; j++)
                shifted[i * dN + j] = CyclicShiftState(i, d, N) * dN + CyclicShiftState(j, d, N);
        var w = new long[dim][];
        for (int a = 0; a < dim; a++)
        {
            w[a] = new long[dim];
            for (int b = 0; b < dim; b++)
            {
                if (ham[a] + ham[b] != N) continue;
                long key = (long)a * dim + b;
                int x = a, y = b;
                for (int t = 1; t < N; t++)
                {
                    x = shifted[x];
                    y = shifted[y];
                    key = Math.Min(key, (long)x * dim + y);
                }
                w[a][b] = OrbitCoefficient(key);
            }
        }
        return RankModP(w);
    }

    private const long ModPrime = 1_000_000_007L;

    /// <summary>A nonzero residue mod p fixed by the orbit key (SplitMix64 finalizer).</summary>
    private static long OrbitCoefficient(long key)
    {
        ulong z = unchecked((ulong)key + 0x9E3779B97F4A7C15UL);
        z = unchecked((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL);
        z = unchecked((z ^ (z >> 27)) * 0x94D049BB133111EBUL);
        z ^= z >> 31;
        return 1 + (long)(z % (ulong)(ModPrime - 1));
    }

    /// <summary>Rank over GF(p) by Gauss-Jordan elimination on residues in [0, p).</summary>
    private static int RankModP(long[][] rows)
    {
        int n = rows.Length;
        int m = n == 0 ? 0 : rows[0].Length;
        int rank = 0;
        for (int col = 0; col < m && rank < n; col++)
        {
            int pivot = -1;
            for (int r = rank; r < n; r++)
                if (rows[r][col] % ModPrime != 0) { pivot = r; break; }
            if (pivot < 0) continue;
            (rows[rank], rows[pivot]) = (rows[pivot], rows[rank]);
            long inv = PowMod(rows[rank][col], ModPrime - 2);
            var prow = rows[rank];
            for (int c = col; c < m; c++) prow[c] = prow[c] * inv % ModPrime;
            for (int r = 0; r < n; r++)
            {
                if (r == rank) continue;
                long f = rows[r][col];
                if (f == 0) continue;
                var row = rows[r];
                for (int c = col; c < m; c++)
                {
                    if (prow[c] == 0) continue;
                    row[c] = (row[c] - f * prow[c] % ModPrime + ModPrime) % ModPrime;
                }
            }
            rank++;
        }
        return rank;
    }

    private static long PowMod(long b, long e)
    {
        long result = 1;
        b %= ModPrime;
        while (e > 0)
        {
            if ((e & 1) == 1) result = result * b % ModPrime;
            b = b * b % ModPrime;
            e >>= 1;
        }
        return result;
    }

    /// <summary>The cyclic site shift (t₀, t₁, …, t_{N−1}) ↦ (t₁, …, t_{N−1}, t₀) on a base-d
    /// state, site 0 the most significant digit.</summary>
    private static int CyclicShiftState(int state, int d, int N)
    {
        int top = checked((int)QuditPartialPalindromeCeiling.IntPow(d, N - 1));
        int lead = state / top;
        return (state % top) * d + lead;
    }

    /// <summary>Build Π_d as a permutation of the d^{2N} coherence indices: basis pair
    /// (i, j) at index i·d^N + j maps to (j, i − chirality mod d per digit), i.e.
    /// Π_d(ρ) = ρᵀ·Shift^{chirality·⊗N}. Site 0 is the most significant base-d digit
    /// (matching the Python verifier's enumeration). Returns perm with
    /// perm[column] = target row. Requires d ≥ 2, N ≥ 1, chirality ∈ {+1, −1};
    /// throws when d^N exceeds Int32 indexing.</summary>
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
    // rate units: rate(i, j) = −2γ·Hamming(i, j). For any W the residual
    // W·L_D + L_D·W + 2Nγ·W has entry W[out, in]·2γ·(N − Ham(in) − Ham(out)),
    // so exactness is the integer statement Ham(in) + Ham(out) = N on the
    // support of W.
    // ------------------------------------------------------------------

    private static int[] BuildPiDPerm(int d, int N, int chirality)
    {
        int dN = checked((int)QuditPartialPalindromeCeiling.IntPow(d, N));
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
        int dN = checked((int)QuditPartialPalindromeCeiling.IntPow(d, N));
        var perm = new int[dN * dN];
        for (int i = 0; i < dN; i++)
            for (int j = 0; j < dN; j++)
                perm[i * dN + j] = j * dN + i;
        return perm;
    }

    /// <summary>Hamming distance per coherence index: digits where i and j disagree.</summary>
    private static int[] HammingTable(int d, int N)
    {
        int dN = checked((int)QuditPartialPalindromeCeiling.IntPow(d, N));
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
        int dN = checked((int)QuditPartialPalindromeCeiling.IntPow(d, N));
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
    // The product lemma, from below: integer per-site factors on the d²
    // letters (letter index a·d + b for |a⟩⟨b|, dark iff a = b).
    // ------------------------------------------------------------------

    private static bool IsLit(int letter, int d) => letter / d != letter % d;

    /// <summary>A deterministic integer per-site factor whose nonzero blocks all carry lit
    /// grade in <paramref name="grades"/> (out-lit + in-lit). Entries are nonzero on every
    /// allowed position, so a single-grade factor has the generic rank of its block.</summary>
    private static long[][] GradeFactor(int d, IReadOnlyCollection<int> grades, int salt)
    {
        int n = d * d;
        var q = new long[n][];
        for (int o = 0; o < n; o++)
        {
            q[o] = new long[n];
            for (int i = 0; i < n; i++)
            {
                int c = (IsLit(o, d) ? 1 : 0) + (IsLit(i, d) ? 1 : 0);
                if (grades.Contains(c))
                    q[o][i] = OrbitCoefficient(((long)salt << 32) + (long)o * n + i);
            }
        }
        return q;
    }

    /// <summary>Exact residual test of W = q₀ ⊗ q₁ at N = 2: every nonzero entry of the product
    /// must sit where Ham(out) + Ham(in) = 2. Returns true iff the residual is exactly zero.</summary>
    private static bool ProductResidualIsZero(long[][] q0, long[][] q1, int d)
    {
        int n = d * d;
        for (int o0 = 0; o0 < n; o0++)
            for (int i0 = 0; i0 < n; i0++)
            {
                if (q0[o0][i0] == 0) continue;
                int h0 = (IsLit(o0, d) ? 1 : 0) + (IsLit(i0, d) ? 1 : 0);
                for (int o1 = 0; o1 < n; o1++)
                    for (int i1 = 0; i1 < n; i1++)
                    {
                        if (q1[o1][i1] == 0) continue;
                        int h1 = (IsLit(o1, d) ? 1 : 0) + (IsLit(i1, d) ? 1 : 0);
                        if (h0 + h1 != 2) return false;
                    }
            }
        return true;
    }

    private static long[][] CopyRows(long[][] rows) => rows.Select(r => (long[])r.Clone()).ToArray();

    /// <summary>Enumerate every grade pattern c ∈ {0, 1, 2}^N with Σc = N and return the
    /// largest product of grade ranks: the brute-force product optimum.</summary>
    private static long PatternMaximum(int d, int N)
    {
        long best = 0;
        var c = new int[N];
        void Walk(int site, int sum)
        {
            if (site == N)
            {
                if (sum != N) return;
                long r = 1;
                for (int l = 0; l < N; l++) r = checked(r * GradeRank(d, c[l]));
                if (r > best) best = r;
                return;
            }
            for (int g = 0; g <= 2; g++)
            {
                c[site] = g;
                Walk(site + 1, sum + g);
            }
        }
        Walk(0, 0);
        return best;
    }

    // ------------------------------------------------------------------
    // Self-check battery: exact integer, permutation and GF(p) arithmetic.
    // ------------------------------------------------------------------

    private static IReadOnlyList<BatteryCase> BuildBattery(QubitNecessityPi2Inheritance qubitNecessity)
    {
        var cases = new List<BatteryCase>();

        // (a) Rank arithmetic on the (d, N) grid: (2d)^N = P ≤ ceiling ≤ d^{2N} for d ≤ 5, with
        //     both equality iffs (P = total ⟺ d = 2; P = ceiling ⟺ d = 2 or N = 1).
        int gridOk = 0, gridTot = 0;
        for (int d = 2; d <= 5; d++)
            for (int N = 1; N <= 3; N++)
            {
                gridTot++;
                long shift = ShiftAlignedRank(d, N);
                long cap = ProductCap(d, N);
                long ceil = CombinatorialCeiling(d, N);
                long total = QuditPartialPalindromeCeiling.Total(d, N);
                bool ok = shift == cap && cap <= ceil && ceil <= total
                       && (cap == total) == (d == 2)
                       && (cap == ceil) == (d == 2 || N == 1);
                if (ok) gridOk++;
            }
        cases.Add(new BatteryCase(
            Name: "cap arithmetic: (2d)^N = P ≤ ceiling ≤ d^{2N} with both equality iffs",
            Detail: "the (d, N) ∈ {2..5}×{1..3} grid; P = total ⟺ d = 2, P = ceiling ⟺ d = 2 or N = 1",
            Expected: "12/12",
            Actual: gridOk.ToString(CultureInfo.InvariantCulture) + "/" +
                    gridTot.ToString(CultureInfo.InvariantCulture)));

        // (b) The trunk polynomial and the per-site swap rank, tied to the parent.
        bool trunkOk = true;
        for (int d = 2; d <= 8; d++)
        {
            trunkOk &= (d * d - 2 * d == 0) == (d == 2);                  // the trunk root
            trunkOk &= 2 * Math.Min(d, d * d - d) == GradeRank(d, 1);     // swap rank 2d
            trunkOk &= (ProductCap(d, 3) == QuditPartialPalindromeCeiling.Total(d, 3)) == (d == 2);
        }
        // At d = 2 the per-site swap rank 2d = 4 is the parent's whole per-site Pauli space,
        // split 4 = 2 + 2 (min(d, d² − d) = 2 = immune = decaying).
        trunkOk &= qubitNecessity.TotalPauliOpsPerSite == 4.0
                && qubitNecessity.ImmuneOpsPerSite == 2.0
                && qubitNecessity.DecayingOpsPerSite == 2.0;
        cases.Add(new BatteryCase(
            Name: "trunk polynomial: d² − 2d = 0 ⟺ d = 2 ⟺ P full; swap rank 2d = parent's 4 = 2 + 2",
            Detail: "third appearance of the trunk (d = 2..8, N = 3); at d = 2 the per-site swap rank " +
                    "2d = 4 exhausts the parent's per-site Pauli space with min(d, d² − d) = 2 = immune = decaying",
            Expected: "trunk root d = 2 only; swap rank 2d; parent split 4 = 2 + 2",
            Actual: trunkOk
                ? "trunk root d = 2 only; swap rank 2d; parent split 4 = 2 + 2"
                : "trunk/rank/parent tie broken"));

        // (c) The product lemma from below at d = 3..6: an integer factor of one lit grade c
        //     reaches rank d, 2d, d² − d (rank mod p ≤ rank ≤ block bound, so equality is exact);
        //     the N = 2 products of patterns (1, 1) and (0, 2) have residual exactly zero; a
        //     factor mixing grades 0 and 1 breaks the identity.
        {
            bool ranksOk = true, residualOk = true, mixedBreaks = true;
            for (int d = 3; d <= 6; d++)
            {
                for (int c = 0; c <= 2; c++)
                    ranksOk &= RankModP(CopyRows(GradeFactor(d, new[] { c }, salt: 10 * d + c))) == GradeRank(d, c);
                residualOk &= ProductResidualIsZero(GradeFactor(d, new[] { 1 }, 1), GradeFactor(d, new[] { 1 }, 2), d);
                residualOk &= ProductResidualIsZero(GradeFactor(d, new[] { 0 }, 3), GradeFactor(d, new[] { 2 }, 4), d);
                mixedBreaks &= !ProductResidualIsZero(GradeFactor(d, new[] { 0, 1 }, 5), GradeFactor(d, new[] { 1 }, 6), d);
            }
            string actual = ranksOk && residualOk && mixedBreaks
                ? "grade ranks d, 2d, d² − d; residual 0 on (1,1) and (0,2); mixed factor breaks"
                : $"ranks {ranksOk}, residual {residualOk}, mixed breaks {mixedBreaks}";
            cases.Add(new BatteryCase(
                Name: "product lemma: one lit grade per site, ranks d, 2d, d² − d, exact on every Σc = N pattern",
                Detail: "integer factors at d = 3..6, rank over GF(p); W = q₀⊗q₁ checked entry by entry " +
                        "(Ham(out) + Ham(in) = N on the support); a factor mixing grades 0 and 1 fails",
                Expected: "grade ranks d, 2d, d² − d; residual 0 on (1,1) and (0,2); mixed factor breaks",
                Actual: actual));
        }

        // (d) The closed form P(d, N) against the enumeration of every grade pattern, the d ≥ 6
        //     threshold (P = (2d)^N iff d ≤ 5 or N = 1), and the d = 6 projector control
        //     P_dark⊗P_lit: support on h = 1 only, residual 0, rank 180 > (2d)^N = 144.
        {
            bool closedOk = true;
            for (int d = 2; d <= 8; d++)
                for (int N = 1; N <= 5; N++)
                {
                    closedOk &= ProductCap(d, N) == PatternMaximum(d, N);
                    closedOk &= (ProductCap(d, N) == ShiftAlignedRank(d, N)) == (d <= 5 || N == 1);
                    closedOk &= ProductCap(d, N) <= CombinatorialCeiling(d, N);
                }
            for (int d = 2; d <= 30; d++)
                closedOk &= (checked((long)d * d * d - (long)d * d) > 4L * d * d) == (d >= 6);
            closedOk &= checked(5L * 5 * 5 - 5 * 5) == 4L * 5 * 5;                       // the tie at d = 5
            const int counterD = 6, counterN = 2;
            int counterDN = checked((int)QuditPartialPalindromeCeiling.IntPow(counterD, counterN));
            var counterHam = HammingTable(counterD, counterN);
            long projectorRank = 0;
            bool projectorResidualExact = true;
            for (int i = 0; i < counterDN; i++)
                for (int j = 0; j < counterDN; j++)
                {
                    if (i / counterD == j / counterD && i % counterD != j % counterD)
                    {
                        int col = i * counterDN + j;
                        projectorRank++;
                        projectorResidualExact &= 2 * counterHam[col] == counterN;
                    }
                }
            bool counterOk = projectorRank == 180 && projectorRank == ProductCap(counterD, counterN)
                          && projectorRank > ShiftAlignedRank(counterD, counterN) && projectorResidualExact;
            cases.Add(new BatteryCase(
                Name: "product optimum P(d, N) = max_m (2d)^(N−2m)(d³ − d²)^m; the swap is optimal iff d ≤ 5",
                Detail: "closed form = enumeration of every grade pattern at d = 2..8, N = 1..5; d³ − d² > 4d² ⟺ " +
                        "d ≥ 6 (d = 2..30, tie at d = 5); at d = 6, N = 2 the (0, 2) projector P_dark⊗P_lit is " +
                        "built column by column: support on h = 1, residual 0, rank 180 = P > 144 = (2d)^N",
                Expected: "closed form = pattern maximum; threshold d = 6; P_dark⊗P_lit rank 180 = P > 144, residual 0",
                Actual: closedOk && counterOk
                    ? "closed form = pattern maximum; threshold d = 6; P_dark⊗P_lit rank 180 = P > 144, residual 0"
                    : $"closed form {closedOk}, d = 6 projector rank {projectorRank} exact {projectorResidualExact}"));
        }

        // (e) Π_d exact on the shift-aligned subspace at (3,1), (3,2), (4,1).
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
            if (dim == ShiftAlignedRank(d, N) && dim == ProductCap(d, N) && closed && exactOnAligned && complementBreaks > 0)
                alignedOk++;
        }
        cases.Add(new BatteryCase(
            Name: "Π_d exact on the shift-aligned subspace at (3,1), (3,2), (4,1)",
            Detail: "aligned dim = (2d)^N = P, Π_d-closed, residual exactly 0 on every aligned column " +
                    "(Ham(col) + Ham(Π_d col) = N), and the complement breaks",
            Expected: "3/3",
            Actual: alignedOk.ToString(CultureInfo.InvariantCulture) + "/3"));

        // (f) The mirror group law at d = 2, 3, 4 (N = 1): ord(Π_d) = 2d, ord(D) = 2,
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

        // (g) The d = 2 degeneracy at N = 1..3: aligned = full space, Π₂⁺ = Π₂⁻ exactly,
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

        // (h) The global (non-product) ceiling-reacher at (3, 2): greedy rung matching pairs
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

        // (i) The chirality swap at d = 3 (N = 1): D maps the +aligned mask exactly onto
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

        // (j) The non-product part ceiling − P: positive exactly when the product mirror is not
        //     the ceiling (d > 2 and N ≥ 2), zero at d = 2 or N = 1; anchors 18 at (3, 2), 64 at
        //     (4, 2), 252 at (6, 2) where the cap is 180, not 144.
        {
            int gapGridOk = 0, gapGridTot = 0;
            for (int d = 2; d <= 7; d++)
                for (int N = 1; N <= 3; N++)
                {
                    gapGridTot++;
                    long gap = NonProductPart(d, N);
                    if (gap >= 0 && (gap > 0) == (d > 2 && N >= 2)) gapGridOk++;
                }
            bool anchorOk = NonProductPart(3, 2) == 18
                         && NonProductPart(2, 3) == 0
                         && NonProductPart(4, 2) == 64
                         && NonProductPart(6, 2) == 252;
            cases.Add(new BatteryCase(
                Name: "non-product part = ceiling − P(d, N)",
                Detail: "positive exactly when the product mirror is not the ceiling (d > 2 and N ≥ 2) " +
                        "on d = 2..7, N = 1..3; 18 at (3, 2), 0 at (2, 3), 64 at (4, 2), 252 at (6, 2)",
                Expected: "18/18 grid, anchors 18/0/64/252",
                Actual: (gapGridOk == gapGridTot && anchorOk)
                    ? "18/18 grid, anchors 18/0/64/252"
                    : $"{gapGridOk}/{gapGridTot} grid, anchors {(anchorOk ? "ok" : "off")}"));
        }

        // (k) Translation invariance reaches the ceiling, exactly: the integer-coefficient TI
        //     intertwiner's rank mod p equals the ceiling at (2, 2), (3, 2), (4, 2).
        {
            var grid = new[] { (d: 2, N: 2), (d: 3, N: 2), (d: 4, N: 2) };
            var ranks = grid.Select(g => TranslationInvariantRankModP(g.d, g.N)).ToArray();
            bool tiOk = grid.Zip(ranks, (g, r) => r == CombinatorialCeiling(g.d, g.N)).All(x => x)
                     && ranks[1] > ProductCap(3, 2) && ranks[2] > ProductCap(4, 2);
            cases.Add(new BatteryCase(
                Name: "translation-invariant rank mod p = ceiling at (2, 2), (3, 2), (4, 2)",
                Detail: "one integer coefficient per T-orbit of the allowed rung-h ↔ rung-(N − h) entries; " +
                        "rank_p ≤ rank ≤ ceiling, so equality is exact; strictly above P at d = 3, 4",
                Expected: "16, 54, 128 = ceiling; above P = 36, 64",
                Actual: tiOk
                    ? "16, 54, 128 = ceiling; above P = 36, 64"
                    : "ranks " + string.Join(", ", ranks.Select(r => r.ToString(CultureInfo.InvariantCulture)))));
        }

        return cases;
    }

    public override string DisplayName =>
        "The qudit product-mirror cap: per-site mirrors pair ≤ P(d, N) = (2d)^N for d ≤ 5, Π_d(ρ) = ρᵀ·Shift^{⊗N} attains it, ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂";

    public override string Summary =>
        "any per-site mirror W = ⊗q_l intertwining W·L_D = (−L_D − 2Nγ)·W pairs at most " +
        "P(d, N) = max_m (2d)^(N−2m)·(d³ − d²)^m of the d^{2N} coherences, = (2d)^N for d ≤ 5 and " +
        $"above it from d = 6 ({ProductCap(6, 2)} > {ShiftAlignedRank(6, 2)} at N = 2); full ⟺ " +
        "d² − 2d = 0 ⟺ d = 2 (the trunk polynomial's third appearance); Π_d(ρ) = ρᵀ·Shift^{⊗N} " +
        "attains (2d)^N with exactly zero residual on the shift-aligned subspace; ord(Π_d) = 2d and " +
        "⟨Π_d, D⟩ ≅ Z_d ≀ Z₂ of order 2d² (D₄ at d = 2); the F121 ceiling needs a non-product " +
        $"mirror, non-product part {NonProductPart(3, 2)} at (3, 2), reached by a translation-invariant " +
        $"one; {PassCount}/{Cases.Count} battery PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("The product cap (theorem)", summary: ProductCapTheorem);
            yield return new InspectableNode("The operator Π_d(ρ) = ρᵀ·Shift^{⊗N}", summary: OperatorRealization);
            yield return new InspectableNode("The mirror group law ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂", summary: MirrorGroupLaw);
            yield return new InspectableNode("The non-product part", summary: NonProductGap);
            yield return new InspectableNode("Translation-invariant recovery at (3, 2), (4, 2), (3, 3)", summary: TranslationInvariantReach);
            yield return InspectableNode.RealScalar("NonProductPart(d=3, N=2)", NonProductPart(3, 2));
            yield return InspectableNode.RealScalar("ProductCap(d=3, N=2)", ProductCap(3, 2));
            yield return InspectableNode.RealScalar("ProductCap(d=6, N=2)", ProductCap(6, 2));
            yield return InspectableNode.RealScalar("ShiftAlignedRank(d=6, N=2)", ShiftAlignedRank(6, 2));
            yield return InspectableNode.RealScalar("CombinatorialCeiling(d=3, N=2)", CombinatorialCeiling(3, 2));
            yield return new InspectableNode("Typed parents",
                summary: $"QuditPartialPalindromeCeiling ({PartialPalindrome.Tier.Label()}): F121, the " +
                         "combinatorial ceiling the cap sits under; the gap to it is the non-product " +
                         $"part. QubitNecessityPi2Inheritance ({QubitNecessity.Tier.Label()}): the " +
                         "d² − 2d = 0 trunk whose unique root d = 2 makes the product mirror full and " +
                         "whose 4 = 2 + 2 per-site split is the swap rank 2d at d = 2.");
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
