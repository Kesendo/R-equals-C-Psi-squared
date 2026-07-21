namespace RCPsiSquared.Core.Numerics;

/// <summary>F133, the symplectic closed form of the F128 cofactor W, recomputed exactly over ℤ
/// (and certified over GF(p)) at inspect time. Second implementation of the read-off core of
/// <c>simulations/f133_w_closed_form.py</c> / <c>docs/proofs/PROOF_F133_W_SYMPLECTIC_CLOSED_FORM.md</c>:
///
/// <code>W = −2⁹ · (Π_u sin x_u) · V_c(a) · V_c(b) · K / SP,   K = 2⁻³⁰ · Σ_λ n_λ · χ^{C₆}_λ</code>
///
/// over the six angles x = (a₁,a₂,a₃,b₁,b₂,b₃), t_u = e^{i x_u/2} the half-angle units. Three exact
/// legs plus a GF(p) certificate, all with disjoint code paths and their own controls:
/// <list type="bullet">
/// <item>the sin-s lemma over ℤ: 𝒪[(2i sin s)·(2i)¹⁵ Δ̂] == the zero dict (Δ̂ = Π_{u&lt;v}
/// sin((x_u−x_v)/2), 2i sin s = t^{1⁶} − t^{−1⁶}, 𝒪 = the (ℤ/2)⁶ signed character sum Π_u(1−flip_u)),
/// with a projector self-test (a generic monomial expands to 64 terms);</item>
/// <item>the alternant read-off, meet-in-the-middle: X = Δ̂·Π_{31 canonical sheets ≠ 1⁶} sin(s_L)
/// = P₁·P₂ (P₁ = Δ̂ · 15 sheets, 590016 monomials; P₂ = the other 16, 5817), and
/// n_raw(λ) = Σ_{ε∈{±1}⁶} sgn(ε)·[X]_{ε∘2(λ+ρ)} = 2·n_λ read off by convolution lookup over the
/// full 557-λ even-degree window (|λ| ≤ 18, ≤ 6 parts); packed 6×10-bit ulong keys, bias 512,
/// exact integer arithmetic;</item>
/// <item>the GF(p) certificate (p ≡ 1 mod 4, i = √−1): C₆·SP == 2⁻⁵·(2i)⁻⁴⁶·Σ_λ n_λ·A_{λ+ρ} at
/// random t_u, LHS the literal 64-flip sum of Δ̂/sin s times SP, RHS the symplectic alternants
/// A_μ = det[t_u^{2μ_j} − t_u^{−2μ_j}] weighted by the READ-OFF-derived n_λ (a code path disjoint
/// from both the table and the flip sum); a corruption control bumps one n_λ and must break the
/// identity.</item>
/// </list>
/// Pure 64-bit integer arithmetic in the lemma/read-off (ulong keys, long coefficients); GF(p) field
/// arithmetic reuses <see cref="CrossFormCertificate"/>'s helpers. No eigensolver, no floats in the
/// exact legs, deterministic (the GF(p) sampler is a fixed xorshift stream seeded from the prime).</summary>
public static class WSymplecticClosedForm
{
    // --- packing: six 10-bit fields in a ulong, bias 512 per field (exponents live in [−36, 36]) ---
    private const int Bits = 10;
    private const int Bias = 512;
    private const ulong FieldMask = (1UL << Bits) - 1;      // 0x3FF
    private static readonly ulong PackBase = ComputePackBase();

    private static ulong ComputePackBase()
    {
        ulong b = 0;
        for (int u = 0; u < 6; u++) b |= (ulong)Bias << (Bits * u);
        return b;
    }

    /// <summary>Pack an exponent vector into a biased ulong key (field u = e[u] + Bias). Throws if a
    /// field leaves [0, 2¹⁰): only ever called on in-range vectors (stored keys, factor shifts,
    /// read-off targets ≤ 36 in abs).</summary>
    private static ulong Pack(ReadOnlySpan<int> e)
    {
        ulong k = 0;
        for (int u = 0; u < 6; u++)
        {
            int shifted = e[u] + Bias;
            if (shifted < 0 || shifted > (int)FieldMask)
                throw new InvalidOperationException($"exponent {e[u]} outside the packing window");
            k |= (ulong)shifted << (Bits * u);
        }
        return k;
    }

    private static int Field(ulong key, int u) => (int)((key >> (Bits * u)) & FieldMask) - Bias;

    private static void AddTo(Dictionary<ulong, long> target, ulong key, long coeff)
    {
        if (coeff == 0) return;
        if (target.TryGetValue(key, out long v))
        {
            v += coeff;
            if (v == 0) target.Remove(key);
            else target[key] = v;
        }
        else target[key] = coeff;
    }

    /// <summary>Multiply a packed polynomial by the two-term factor t^vec − t^{−vec}. The shift is
    /// stored unbiased (Pack(vec) − PackBase); adding it to a biased key gives a biased key, exact
    /// because the field sums stay in range (the mod-2⁶⁴ ring makes intermediate negatives harmless).</summary>
    private static Dictionary<ulong, long> MulFactor(Dictionary<ulong, long> poly, ReadOnlySpan<int> vec)
    {
        Span<int> neg = stackalloc int[6];
        for (int u = 0; u < 6; u++) neg[u] = -vec[u];
        ulong up = Pack(vec) - PackBase;
        ulong dn = Pack(neg) - PackBase;
        var r = new Dictionary<ulong, long>(poly.Count * 2);
        foreach (var (k, c) in poly)
        {
            AddTo(r, k + up, c);
            AddTo(r, k + dn, -c);
        }
        return r;
    }

    // The 32 canonical sheets L ∈ {±1}⁶ with L₁ = +1, in itertools.product order (last coord fastest);
    // index 0 = 1⁶ is dropped to form the 31 non-trivial sheets. Kept identical to the Python gate so
    // the P₁ = first 15 / P₂ = last 16 split reproduces the committed 590016 / 5817 sizes.
    private static int[][] Sheets31()
    {
        var all = new List<int[]>();
        for (int mask = 0; mask < 32; mask++)
        {
            var L = new int[6];
            L[0] = 1;
            for (int j = 0; j < 5; j++) L[1 + j] = ((mask >> (4 - j)) & 1) == 1 ? -1 : 1;
            all.Add(L);
        }
        // all[0] is 1⁶; the rest are the 31 canonical sheets ≠ 1⁶.
        return all.GetRange(1, 31).ToArray();
    }

    // ---------------------------------------------------------------- the sin-s lemma (exact ℤ) ----

    /// <summary>What <see cref="AnalyzeSinSLemma"/> recomputes. <paramref name="DeltaMonomials"/> is the
    /// (2i)¹⁵ Δ̂ dict size (720, the A₅ Weyl denominator); <paramref name="SurvivingMonomials"/> must
    /// be 0; the projector self-test must expand a generic monomial to 64 terms.</summary>
    public sealed record SinSLemmaReport(int DeltaMonomials, int SurvivingMonomials, bool ProjectorSelfTestOk);

    /// <summary>(2i)¹⁵ Δ̂ = Π_{u&lt;v}(t_u/t_v − t_v/t_u), the full six-angle half-angle Vandermonde.</summary>
    private static Dictionary<ulong, long> DeltaDict()
    {
        var p = new Dictionary<ulong, long> { [PackBase] = 1 };
        Span<int> e = stackalloc int[6];
        for (int u = 0; u < 6; u++)
            for (int v = u + 1; v < 6; v++)
            {
                e.Clear(); e[u] = 1; e[v] = -1;
                p = MulFactor(p, e);
            }
        return p;
    }

    /// <summary>Π_u (1 − flip_u), the 64·𝒪 signed character sum, applied variable by variable.</summary>
    private static Dictionary<ulong, long> OddProject(Dictionary<ulong, long> p)
    {
        for (int u = 0; u < 6; u++)
        {
            var r = new Dictionary<ulong, long>(p.Count * 2);
            foreach (var (k, c) in p)
            {
                AddTo(r, k, c);
                int e = Field(k, u);
                ulong flipped = (k & ~(FieldMask << (Bits * u))) | (ulong)(-e + Bias) << (Bits * u);
                AddTo(r, flipped, -c);
            }
            p = r;
        }
        return p;
    }

    /// <summary>The sin-s lemma: 𝒪[(2i sin s)·(2i)¹⁵ Δ̂] is the zero polynomial exactly over ℤ, with a
    /// projector self-test (a generic monomial must expand to 64 distinct ±1 terms, not to zero).</summary>
    public static SinSLemmaReport AnalyzeSinSLemma()
    {
        var delta = DeltaDict();
        // 2i sin s = t^{1⁶} − t^{−1⁶}
        Span<int> ones = stackalloc int[6];
        for (int u = 0; u < 6; u++) ones[u] = 1;
        var product = MulFactor(delta, ones);   // (t^{1⁶} − t^{−1⁶}) · (2i)¹⁵ Δ̂
        var projected = OddProject(product);

        var single = new Dictionary<ulong, long> { [Pack(new[] { 1, 2, 3, 4, 5, 6 })] = 1 };
        var selfTest = OddProject(single);
        bool ok = selfTest.Count == 64;
        if (ok)
            foreach (long c in selfTest.Values)
                if (c != 1 && c != -1) { ok = false; break; }

        return new SinSLemmaReport(delta.Count, projected.Count, ok);
    }

    // ---------------------------------------------------------- the alternant read-off (exact ℤ) ----

    /// <summary>The meet-in-the-middle halves of X = Δ̂ · Π_{31 sheets} sin(s_L) = P₁ · P₂.</summary>
    public sealed class Halves
    {
        internal Dictionary<ulong, long> P1 = null!;
        internal (ulong Key, long Coeff)[] P2 = null!;
        public int P1Size => P1.Count;
        public int P2Size => P2.Length;
    }

    private static readonly int[] Rho = { 6, 5, 4, 3, 2, 1 };

    /// <summary>Build P₁ = (2i)¹⁵ Δ̂ · (first 15 non-trivial sheets) and P₂ = (last 16 sheets), each
    /// sheet factor multiplying by t^L − t^{−L}. Sizes reproduce the committed 590016 / 5817.</summary>
    public static Halves BuildHalves()
    {
        var sheets = Sheets31();
        var p1 = DeltaDict();
        for (int s = 0; s < 15; s++) p1 = MulFactor(p1, sheets[s]);
        var p2 = new Dictionary<ulong, long> { [PackBase] = 1 };
        for (int s = 15; s < 31; s++) p2 = MulFactor(p2, sheets[s]);
        return new Halves { P1 = p1, P2 = p2.Select(kv => (kv.Key, kv.Value)).ToArray() };
    }

    /// <summary>n_raw(λ) = Σ_{ε∈{±1}⁶} sgn(ε)·[X]_{ε∘2(λ+ρ)}, X = P₁·P₂, by convolution lookup.
    /// Equals 2·n_λ over the window. λ may have fewer than 6 parts (zero-padded).</summary>
    public static long NRaw(int[] lambda, Halves h)
    {
        Span<int> mu = stackalloc int[6];
        for (int j = 0; j < 6; j++)
        {
            int lj = j < lambda.Length ? lambda[j] : 0;
            mu[j] = 2 * (lj + Rho[j]);
        }
        var p1 = h.P1;
        var p2 = h.P2;
        Span<int> v = stackalloc int[6];
        long tot = 0;
        for (int mask = 0; mask < 64; mask++)
        {
            int sgn = 1;
            for (int u = 0; u < 6; u++)
            {
                if ((mask >> u & 1) == 0) v[u] = mu[u];
                else { v[u] = -mu[u]; sgn = -sgn; }
            }
            ulong target = Pack(v) + PackBase;    // biased-by-two; the second bias cancels the P₂ key bias
            foreach (var (k2, c2) in p2)
                if (p1.TryGetValue(target - k2, out long v1))
                    tot += sgn * c2 * v1;
        }
        return tot;
    }

    /// <summary>The 557-λ even-degree window (|λ| even ≤ 18, ≤ 6 parts, non-increasing), including
    /// the empty partition. Mirrors the Python gate's <c>window_lams(full=False)</c>.</summary>
    public static List<int[]> WindowLambdas()
    {
        var lams = new List<int[]> { Array.Empty<int>() };
        void Gen(List<int> prefix, int mx, int rem)
        {
            for (int q = Math.Min(mx, rem); q >= 1; q--)
            {
                if (prefix.Count < 6)
                {
                    prefix.Add(q);
                    if (prefix.Sum() % 2 == 0) lams.Add(prefix.ToArray());
                    Gen(prefix, q, rem - q);
                    prefix.RemoveAt(prefix.Count - 1);
                }
            }
        }
        Gen(new List<int>(), 18, 18);
        return lams;
    }

    /// <summary>What <see cref="AnalyzeReadoff"/> recomputes. Aggregates are n_λ (= n_raw/2):
    /// <paramref name="NonzeroCount"/> must be 143, <paramref name="MaxAbsN"/> 8, <paramref name="SumAbsN"/>
    /// 359, <paramref name="NAtEmpty"/> −1; every window λ must satisfy n_raw = 2·(committed n_λ).</summary>
    public sealed record ReadoffReport(
        int P1Size, int P2Size, int WindowCount, int NonzeroCount,
        long MaxAbsN, long SumAbsN, long NAtEmpty, int Mismatches, bool AllMatchTable);

    /// <summary>Full 557-λ window read-off vs the committed table (embedded below), with the aggregate
    /// checksums. Builds P₁/P₂ once; the 207M-lookup sweep is why this belongs in the tests, not the
    /// live witness.</summary>
    public static ReadoffReport AnalyzeReadoff()
    {
        var h = BuildHalves();
        var table = TableCoefficients();
        var tableMap = new Dictionary<string, long>();
        foreach (var (lam, n) in table) tableMap[Key(lam)] = n;

        int mismatches = 0, nonzero = 0;
        long maxAbs = 0, sumAbs = 0, nEmpty = 0;
        foreach (var lam in WindowLambdas())
        {
            long raw = NRaw(lam, h);
            long want = 2 * (tableMap.TryGetValue(Key(lam), out long t) ? t : 0);
            if (raw != want) mismatches++;
            if (raw != 0)
            {
                nonzero++;
                long n = raw / 2;
                maxAbs = Math.Max(maxAbs, Math.Abs(n));
                sumAbs += Math.Abs(n);
                if (lam.Length == 0) nEmpty = n;
            }
        }
        bool all = mismatches == 0 && nonzero == table.Count;
        return new ReadoffReport(h.P1Size, h.P2Size, WindowLambdas().Count, nonzero,
                                 maxAbs, sumAbs, nEmpty, mismatches, all);
    }

    private static string Key(int[] lam) => string.Join(",", lam);

    // ------------------------------------------------------------- the F134 two-row reflection ----

    /// <summary>What <see cref="AnalyzeTwoRowReflection"/> recomputes (F134): the two-row reflection
    /// law n_(j,k) = n_(10−j,k) on the embedded table (36 ordered pairs, 14 live), the l = 2 domain
    /// break (exactly 8 breaking pairs; the two nonzero l = 2 entries satisfying the reflection are
    /// the center-fixed j = 5 weights, μ₁ = 11), and, when halves are supplied, an NRaw spot check of
    /// every live pair on the read-off path (code-disjoint from the table).</summary>
    public sealed record TwoRowReflectionReport(
        int PairsChecked, int LivePairs, int Mismatches,
        int L2PairsChecked, int L2Breaks, bool L2CenterHoldsNonzero,
        int SpotPairsChecked, int SpotMismatches);

    public static TwoRowReflectionReport AnalyzeTwoRowReflection(Halves? spotHalves = null)
    {
        var tableMap = new Dictionary<string, long>();
        foreach (var (lam, n) in TableCoefficients()) tableMap[Key(lam)] = n;

        static int[] Lam(params int[] parts)
        {
            int len = parts.Length;
            while (len > 0 && parts[len - 1] == 0) len--;
            return parts.Take(len).ToArray();
        }
        long N(params int[] parts) => tableMap.TryGetValue(Key(Lam(parts)), out long n) ? n : 0;

        int pairs = 0, live = 0, mism = 0, spotPairs = 0, spotMism = 0;
        for (int k = 0; k <= 5; k++)
            for (int j = k; j <= 10 - k; j++)
            {
                long a = N(j, k), b = N(10 - j, k);
                pairs++;
                if (a != b) mism++;
                if (a != 0 || b != 0)
                {
                    live++;
                    if (spotHalves is not null)
                    {
                        spotPairs++;
                        long ra = NRaw(Lam(j, k), spotHalves), rb = NRaw(Lam(10 - j, k), spotHalves);
                        if (ra != rb || ra != 2 * a) spotMism++;
                    }
                }
            }

        int l2Pairs = 0, l2Breaks = 0;
        for (int k = 2; k <= 5; k++)
            for (int j = k; j <= 10 - k; j++)
            {
                l2Pairs++;
                if (N(j, k, 2) != N(10 - j, k, 2)) l2Breaks++;
            }
        bool centerHolds = N(5, 3, 2) != 0 && N(5, 5, 2) != 0;
        return new TwoRowReflectionReport(pairs, live, mism, l2Pairs, l2Breaks, centerHolds,
                                          spotPairs, spotMism);
    }

    // ------------------------------------------------------------------ the F139 seam identity ----

    /// <summary>What <see cref="AnalyzeSeamIdentity"/> recomputes (F139): the seam identity behind the
    /// F134 reflection law, on a code path independent of both the Python gate and the six-variable
    /// read-off. A priori from the letter system (five variables x = t′²): the cosine matrix
    /// ψ_{i,t} = R⁻[e_i(P)·e_t({C_M})] (P = {e₁..e₅, 𝟙} the 6-letter window; C_M the 15 Core cosine
    /// pairs), the one-variable polynomial Φ_k = (y²−1)φ₀ + yφ₁ + φ₂ with φ_i = Σ_t ψ_{i,t} y^{15−t},
    /// and the division Φ_k = S₁₀·(−1)^k P_k + R_k by the Chebyshev wall divisor
    /// S₁₀(y) = sin(11θ)/sinθ. Checks: the skew ψ_{6−i,t} = −ψ_{i,t}; deg Φ_k ≤ 15; the a-priori
    /// table n_(λ₁,k) = (−1)^{λ₁}·b_{λ₁+5} vs the embedded F133 table on every parity-live strip
    /// cell; the quotients equal (−1)^k P_k with deg R_k ≤ 4+k and R₀ = 0; the fence (l = 1
    /// remainders ≤ 4, l = 2 remainders ≥ 8: the breaks are remainder overflow); and a corruption
    /// control (one ψ entry bumped must break the table reproduction).</summary>
    public sealed record SeamIdentityReport(
        bool SkewOk, bool DegreeLemmaOk, int TableCellsChecked, int TableMismatches,
        bool QuotientsMatch, bool RemainderBoundsOk, bool K0RemainderZero,
        bool FenceL1Ok, bool FenceL2Overflow, bool CorruptionBroke);

    private const int Bits5 = 12;
    private const int Bias5 = 2048;

    private static ulong Pack5(ReadOnlySpan<int> e)
    {
        ulong k = 0;
        for (int u = 0; u < 5; u++)
        {
            int shifted = e[u] + Bias5;
            if (shifted < 0 || shifted >= (1 << Bits5))
                throw new InvalidOperationException($"5-var exponent {e[u]} outside the packing window");
            k |= (ulong)shifted << (Bits5 * u);
        }
        return k;
    }

    private static readonly ulong PackBase5 = Pack5(stackalloc int[5]);

    /// <summary>e-series of a letter multiset where each letter is a (possibly 2-term) Laurent
    /// polynomial: T[j] = coefficient dict of e_j after absorbing every letter.</summary>
    private static Dictionary<ulong, long>[] ESeries5((ulong Key, long Coeff)[][] letters, int maxDeg)
    {
        var T = new Dictionary<ulong, long>[maxDeg + 1];
        for (int j = 0; j <= maxDeg; j++) T[j] = new Dictionary<ulong, long>();
        T[0][PackBase5] = 1;
        foreach (var letter in letters)
            for (int j = maxDeg; j >= 1; j--)
                foreach (var (k, c) in T[j - 1].ToArray())
                    foreach (var (lk, lc) in letter)
                        AddTo(T[j], k + lk - PackBase5, c * lc);
        for (int j = 0; j <= maxDeg; j++)
            foreach (var (k, c) in T[j].ToArray())
                if (c == 0) T[j].Remove(k);
        return T;
    }

    private static (Dictionary<ulong, long>[] ECos, Dictionary<ulong, long>[] EWin) SeamLetterTables()
    {
        // Core cosine pairs: the 15 antipodal pairs of {±1}⁵ ∖ {±𝟙}, representative first coord +1.
        var cos = new List<(ulong, long)[]>();
        for (int mask = 0; mask < 16; mask++)
        {
            var m = new int[5];
            m[0] = 1;
            for (int j = 0; j < 4; j++) m[1 + j] = ((mask >> j) & 1) == 1 ? -1 : 1;
            bool isPole = m.All(v => v == 1);
            if (isPole) continue;                       // +𝟙 excluded; −𝟙 is the pair of +𝟙
            var neg = m.Select(v => -v).ToArray();
            cos.Add(new[] { (Pack5(m), 1L), (Pack5(neg), 1L) });
        }
        // The window P = {e₁..e₅, 𝟙}.
        var win = new List<(ulong, long)[]>();
        for (int v = 0; v < 5; v++)
        {
            var e = new int[5];
            e[v] = 1;
            win.Add(new[] { (Pack5(e), 1L) });
        }
        win.Add(new[] { (Pack5(new[] { 1, 1, 1, 1, 1 }), 1L) });
        return (ESeries5(cos.ToArray(), 15), ESeries5(win.ToArray(), 6));
    }

    private static readonly int[] Delta5 = { 2, 1, 0, -1, -2 };
    private static readonly int[] Rho5 = { 5, 4, 3, 2, 1 };

    /// <summary>ψ_{i,t} at ν: the signed W(C₅) read-off Σ_{ε,σ} sgn(ε)sgn(σ)·
    /// [e_i(P)·e_t({C_M})]_{ε∘μ + 𝟙 − σδ}, μ = ν + ρ₅.</summary>
    private static long Psi(int i, int t, ReadOnlySpan<int> nu,
                            Dictionary<ulong, long>[] eCos, Dictionary<ulong, long>[] eWin)
    {
        Span<int> mu = stackalloc int[5];
        for (int u = 0; u < 5; u++) mu[u] = nu[u] + Rho5[u];
        long tot = 0;
        var perm = new int[5];
        Span<int> arg = stackalloc int[5];
        // enumerate S₅ via Heap-free index loops (120 perms) × 32 sign vectors
        var perms = Permutations5();
        foreach (var (p, ps) in perms)
            for (int eps = 0; eps < 32; eps++)
            {
                int es = 1;
                for (int u = 0; u < 5; u++)
                {
                    int e = ((eps >> u) & 1) == 1 ? -1 : 1;
                    es *= e;
                    arg[u] = e * mu[u] + 1 - Delta5[p[u]];
                }
                ulong target = Pack5(arg);
                long inner = 0;
                foreach (var (wk, wc) in eWin[i])
                {
                    if (eCos[t].TryGetValue(target - wk + PackBase5, out long cc))
                        inner += wc * cc;
                }
                tot += (long)(ps * es) * inner;
            }
        return tot;
    }

    private static List<(int[] P, int Sign)>? _perms5;
    private static List<(int[] P, int Sign)> Permutations5()
    {
        if (_perms5 is not null) return _perms5;
        var result = new List<(int[], int)>();
        var idx = new[] { 0, 1, 2, 3, 4 };
        void Permute(int[] a, int k)
        {
            if (k == a.Length)
            {
                int s = 1;
                for (int x = 0; x < 5; x++)
                    for (int y = x + 1; y < 5; y++)
                        if (a[x] > a[y]) s = -s;
                result.Add(((int[])a.Clone(), s));
                return;
            }
            for (int x = k; x < a.Length; x++)
            {
                (a[k], a[x]) = (a[x], a[k]);
                Permute(a, k + 1);
                (a[k], a[x]) = (a[x], a[k]);
            }
        }
        Permute(idx, 0);
        _perms5 = result;
        return result;
    }

    /// <summary>Φ at ν: monomial coefficients (index = y-power) of (y²−1)φ₀ + yφ₁ + φ₂;
    /// optionally bump one ψ entry (the corruption control).</summary>
    private static long[] PhiPoly(ReadOnlySpan<int> nu, Dictionary<ulong, long>[] eCos,
                                  Dictionary<ulong, long>[] eWin, (int I, int T)? bump = null)
    {
        var outp = new long[18];
        for (int i = 0; i <= 2; i++)
            for (int t = 0; t <= 15; t++)
            {
                long v = Psi(i, t, nu, eCos, eWin);
                if (bump is not null && bump.Value.I == i && bump.Value.T == t) v += 1;
                if (v == 0) continue;
                int p = 15 - t;
                if (i == 0) { outp[p + 2] += v; outp[p] -= v; }
                else if (i == 1) outp[p + 1] += v;
                else outp[p] += v;
            }
        return outp;
    }

    /// <summary>Chebyshev S-polynomials in the 2cos normalization, S_m(2cosθ) = sin((m+1)θ)/sinθ.</summary>
    private static long[][] ChebyshevS(int maxM)
    {
        var s = new long[maxM + 1][];
        s[0] = new long[] { 1 };
        s[1] = new long[] { 0, 1 };
        for (int m = 2; m <= maxM; m++)
        {
            s[m] = new long[m + 1];
            for (int p = 0; p < s[m - 1].Length; p++) s[m][p + 1] += s[m - 1][p];
            for (int p = 0; p < s[m - 2].Length; p++) s[m][p] -= s[m - 2][p];
        }
        return s;
    }

    private static Dictionary<int, long> ToSBasis(long[] poly, long[][] s)
    {
        var a = (long[])poly.Clone();
        var b = new Dictionary<int, long>();
        for (int m = a.Length - 1; m >= 0; m--)
        {
            long c = a[m];
            if (c == 0) continue;
            b[m] = c;
            for (int p = 0; p < s[m].Length; p++) a[p] -= c * s[m][p];
        }
        if (a.Any(x => x != 0))
            throw new InvalidOperationException("S-basis reduction left a residue");
        return b;
    }

    private static (long[] Q, long[] R) DivMod(long[] a, long[] b)
    {
        var r = (long[])a.Clone();
        int db = b.Length - 1;
        while (db > 0 && b[db] == 0) db--;
        if (Math.Abs(b[db]) != 1)
            throw new InvalidOperationException("DivMod assumes a monic (or anti-monic) divisor");
        var q = new long[Math.Max(r.Length - db, 1)];
        for (int da = r.Length - 1; da >= db; da--)
        {
            if (r[da] == 0) continue;
            long c = r[da] / b[db];
            q[da - db] = c;
            for (int p = 0; p <= db; p++) r[p + da - db] -= c * b[p];
        }
        return (q, r);
    }

    private static int Deg(long[] a)
    {
        for (int p = a.Length - 1; p >= 0; p--)
            if (a[p] != 0) return p;
        return -1;
    }

    public static SeamIdentityReport AnalyzeSeamIdentity()
    {
        var (eCos, eWin) = SeamLetterTables();
        var s = ChebyshevS(17);
        var s10 = s[10];

        // skew ψ_{6−i,t} = −ψ_{i,t} on the whole grid, one-row k = 0..5
        bool skewOk = true;
        for (int k = 0; k <= 5 && skewOk; k++)
        {
            var nu = new[] { k, 0, 0, 0, 0 };
            for (int i = 0; i <= 6 && skewOk; i++)
                for (int t = 0; t <= 15 && skewOk; t++)
                    if (Psi(6 - i, t, nu, eCos, eWin) != -Psi(i, t, nu, eCos, eWin))
                        skewOk = false;
        }

        // Φ_k, degree lemma, a-priori table vs the embedded F133 table
        var tableMap = new Dictionary<string, long>();
        foreach (var (lam, n) in TableCoefficients())
            if (lam.Length <= 2)
                tableMap[$"{(lam.Length > 0 ? lam[0] : 0)},{(lam.Length > 1 ? lam[1] : 0)}"] = n;

        var phis = new long[6][];
        var sBasis = new Dictionary<int, long>[6];
        bool degOk = true;
        for (int k = 0; k <= 5; k++)
        {
            phis[k] = PhiPoly(new[] { k, 0, 0, 0, 0 }, eCos, eWin);
            if (Deg(phis[k]) > 15) degOk = false;
            sBasis[k] = ToSBasis(phis[k], s);
        }

        int cells = 0, mism = 0;
        for (int lam1 = 0; lam1 <= 10; lam1++)
            for (int k = 0; k <= Math.Min(lam1, 5); k++)
            {
                if ((lam1 + k) % 2 != 0) continue;
                cells++;
                long b = sBasis[k].TryGetValue(lam1 + 5, out long bv) ? bv : 0;
                long got = (lam1 % 2 == 0 ? 1 : -1) * b;
                long want = tableMap.TryGetValue($"{lam1},{k}", out long tv) ? tv : 0;
                if (got != want) mism++;
            }

        // the division: expected quotients (−1)^k·P_k in the S-basis
        var expectedQ = new Dictionary<int, long>[]
        {
            new() { [1] = 1, [5] = -1 },            // (−1)⁰ P₀ = S₁ − S₅
            new() { [0] = 2, [2] = 1, [4] = -1 },   // (−1)¹ P₁ = −(S₄ − S₂ − 2S₀) = 2S₀ + S₂ − S₄
            new() { [1] = 3 },                      // (−1)² P₂ = 3S₁
            new() { [0] = 2, [2] = 1 },             // (−1)³ P₃ = −P₃ = S₂ + 2S₀
            new() { [1] = 1 },                      // (−1)⁴ P₄ = S₁
            new(),                                  // (−1)⁵ P₅ = 0
        };
        bool quotOk = true, remOk = true, k0Zero = true;
        for (int k = 0; k <= 5; k++)
        {
            var (q, r) = DivMod(phis[k], s10);
            var bq = ToSBasis(q, s).Where(kv => kv.Value != 0).ToDictionary(kv => kv.Key, kv => kv.Value);
            var want = expectedQ[k].Where(kv => kv.Value != 0).ToDictionary(kv => kv.Key, kv => kv.Value);
            if (bq.Count != want.Count || bq.Any(kv => !want.TryGetValue(kv.Key, out long w) || w != kv.Value))
                quotOk = false;
            if (Deg(r) > 4 + k) remOk = false;
            if (k == 0 && Deg(r) != -1) k0Zero = false;
        }

        // the fence: l = 1 small, l = 2 overflow
        bool l1Ok = true, l2Ok = true;
        foreach (var nu in new[] { new[] { 2, 1, 0, 0, 0 }, new[] { 3, 1, 0, 0, 0 } })
        {
            var (_, r) = DivMod(PhiPoly(nu, eCos, eWin), s10);
            if (Deg(r) > 4) l1Ok = false;
        }
        foreach (var nu in new[] { new[] { 2, 2, 0, 0, 0 }, new[] { 3, 2, 0, 0, 0 }, new[] { 4, 2, 0, 0, 0 } })
        {
            var (_, r) = DivMod(PhiPoly(nu, eCos, eWin), s10);
            if (Deg(r) < 8) l2Ok = false;
        }

        // corruption control: bump the PARITY-LIVE ψ_{0,4} in Φ₀ (t even is the live parity at
        // k = 0, i = 0; the bump feeds y¹³ − y¹¹, which the k = 0 row reads); must break
        var bumped = ToSBasis(PhiPoly(new[] { 0, 0, 0, 0, 0 }, eCos, eWin, (0, 4)), s);
        bool broke = false;
        for (int lam1 = 0; lam1 <= 10 && !broke; lam1 += 2)
        {
            long b = bumped.TryGetValue(lam1 + 5, out long bv) ? bv : 0;
            long want = tableMap.TryGetValue($"{lam1},0", out long tv) ? tv : 0;
            if (b != want) broke = true;
        }

        return new SeamIdentityReport(skewOk, degOk, cells, mism, quotOk, remOk, k0Zero,
                                      l1Ok, l2Ok, broke);
    }

    /// <summary>The read-off-derived nonzero coefficients (λ, n_λ) over the window: a code path
    /// disjoint from the embedded table (the GF(p) certificate consumes THESE, not the table).</summary>
    public static List<(int[] Lam, long N)> DeriveCoefficients(Halves h)
    {
        var result = new List<(int[], long)>();
        foreach (var lam in WindowLambdas())
        {
            long raw = NRaw(lam, h);
            if (raw != 0) result.Add((lam, raw / 2));
        }
        return result;
    }

    // ------------------------------------------------------------------- the GF(p) certificate ----

    /// <summary>Deterministic xorshift64* stream (never System.Random, so the witness reproduces
    /// bit-for-bit), matching the <see cref="CrossFormCertificate"/> idiom.</summary>
    private sealed class Stream(ulong seed)
    {
        private ulong _s = seed == 0 ? 0x9E3779B97F4A7C15UL : seed;

        public long Next(long loInclusive, long hiExclusive)
        {
            _s ^= _s >> 12; _s ^= _s << 25; _s ^= _s >> 27;
            ulong r = _s * 0x2545F4914F6CDD1DUL;
            return loInclusive + (long)(r % (ulong)(hiExclusive - loInclusive));
        }
    }

    /// <summary>What <see cref="CertifyGfpSlice"/> recomputes at random GF(p) points: the identity
    /// C₆·SP == 2⁻⁵·(2i)⁻⁴⁶·Σ_λ n_λ·A_{λ+ρ} (mismatches must be 0), and the corruption control
    /// (one n_λ bumped) which must BREAK it at nearly every point.</summary>
    public sealed record GfpSliceReport(int Points, int Mismatches, int ControlChecks, int ControlMismatches);

    private static long DetMod(long[][] mat, long p)
    {
        int n = mat.Length;
        var m = new long[n][];
        for (int i = 0; i < n; i++) m[i] = (long[])mat[i].Clone();
        long det = 1;
        for (int col = 0; col < n; col++)
        {
            int piv = -1;
            for (int k = col; k < n; k++) if (m[k][col] % p != 0) { piv = k; break; }
            if (piv < 0) return 0;
            if (piv != col) { (m[col], m[piv]) = (m[piv], m[col]); det = det == 0 ? 0 : p - det; }
            det = det * m[col][col] % p;
            long inv = CrossFormCertificate.Inv(m[col][col], p);
            for (int k = col + 1; k < n; k++)
            {
                long f = m[k][col] * inv % p;
                if (f == 0) continue;
                for (int j = col; j < n; j++)
                    m[k][j] = ((m[k][j] - f * m[col][j]) % p + p) % p;
            }
        }
        return det % p;
    }

    /// <summary>One prime's slice of the closed-form identity C₆·SP == 2⁻⁵·(2i)⁻⁴⁶·Σ_λ n_λ·A_{λ+ρ},
    /// evaluated at <paramref name="samples"/> random points t_u ∈ GF(p). LHS = SP · (1/64) Σ_ε sgn(ε)
    /// Δ̂(ε∘x)/sin s(ε∘x) (the literal 64-flip sum, SP flip-invariant so C₆·SP = 𝒪[X]); RHS uses the
    /// symplectic alternants A_μ = det[t_u^{2μ_j} − t_u^{−2μ_j}] weighted by <paramref name="coeffs"/>
    /// (the read-off-derived n_λ). Points with sin s(ε∘x) = 0 on any flip are skipped. The corruption
    /// control bumps one n_λ by +1; the identity must then fail. Requires p ≡ 1 mod 4.</summary>
    public static GfpSliceReport CertifyGfpSlice(long p, int samples, IReadOnlyList<(int[] Lam, long N)> coeffs)
    {
        long iRoot = CrossFormCertificate.SqrtMinusOne(p);
        if (iRoot < 0 || iRoot * iRoot % p != p - 1)
            throw new InvalidOperationException($"p = {p} is not split (need p = 1 mod 4)");
        long Inv(long a) => CrossFormCertificate.Inv(a, p);
        long inv2i = Inv(2 * iRoot % p);
        long inv64 = Inv(64);
        long SinM(long m) => (((m - Inv(m)) % p + p) % p) * inv2i % p;

        // RHS constant 2⁻⁵·(2i)⁻⁴⁶ and the constant symplectic alternant exponents μ = λ + ρ (2μ used).
        long rhsConst = Inv(32) * Inv(CrossFormCertificate.PowMod(2 * iRoot % p, 46, p)) % p;
        var mus = new int[coeffs.Count][];
        for (int c = 0; c < coeffs.Count; c++)
        {
            var lam = coeffs[c].Lam;
            var mu = new int[6];
            for (int j = 0; j < 6; j++) mu[j] = (j < lam.Length ? lam[j] : 0) + Rho[j];
            mus[c] = mu;
        }

        var sheets = new int[32][];
        for (int mask = 0; mask < 32; mask++)
        {
            var L = new int[6]; L[0] = 1;
            for (int j = 0; j < 5; j++) L[1 + j] = ((mask >> (4 - j)) & 1) == 1 ? -1 : 1;
            sheets[mask] = L;
        }

        var stream = new Stream((ulong)p * 0xF1357AEA2E62A9C5UL);
        int points = 0, mismatches = 0, ctrlChecks = 0, ctrlMismatches = 0;
        var t = new long[6];
        var tinv = new long[6];
        int guard = 0;
        while (points < samples && guard++ < samples * 50)
        {
            for (int u = 0; u < 6; u++) { t[u] = stream.Next(2, p); tinv[u] = Inv(t[u]); }

            // SP = Π over the 32 canonical sheets of sin(s_L), T_L = Π t_u^{L_u}.
            long sp = 1;
            for (int s = 0; s < 32; s++)
            {
                long tl = 1;
                for (int u = 0; u < 6; u++) tl = tl * (sheets[s][u] == 1 ? t[u] : tinv[u]) % p;
                sp = sp * SinM(tl) % p;
            }

            // C₆ = (1/64) Σ_ε sgn(ε) Δ̂(ε∘x) / sin s(ε∘x), τ_u = t_u^{ε_u}.
            long c6 = 0;
            bool skip = false;
            for (int mask = 0; mask < 64 && !skip; mask++)
            {
                int sgn = 1;
                long te = 1;
                var tau = new long[6];
                for (int u = 0; u < 6; u++)
                {
                    if ((mask >> u & 1) == 0) tau[u] = t[u];
                    else { tau[u] = tinv[u]; sgn = -sgn; }
                    te = te * tau[u] % p;
                }
                if (te * te % p == 1) { skip = true; break; }   // sin s = 0 on this flip: skip the point
                long dh = 1;
                for (int u = 0; u < 6; u++)
                    for (int v = u + 1; v < 6; v++)
                        dh = dh * SinM(tau[u] * Inv(tau[v]) % p) % p;
                long term = sgn == 1 ? dh * Inv(SinM(te)) % p : (p - dh * Inv(SinM(te)) % p) % p;
                c6 = (c6 + term) % p;
            }
            if (skip) continue;
            c6 = c6 * inv64 % p;
            long lhs = c6 * sp % p;

            // RHS = 2⁻⁵·(2i)⁻⁴⁶ Σ_λ n_λ A_{λ+ρ}, A_μ = det[t_u^{2μ_j} − t_u^{−2μ_j}].
            long sumA = 0;
            var alt = new long[coeffs.Count];
            for (int c = 0; c < coeffs.Count; c++)
            {
                var mu = mus[c];
                var m = new long[6][];
                for (int u = 0; u < 6; u++)
                {
                    m[u] = new long[6];
                    for (int j = 0; j < 6; j++)
                    {
                        long up = CrossFormCertificate.PowMod(t[u], 2 * mu[j], p);
                        m[u][j] = ((up - Inv(up)) % p + p) % p;
                    }
                }
                long a = DetMod(m, p);
                alt[c] = a;
                sumA = (sumA + ((coeffs[c].N % p + p) % p) * a) % p;
            }
            long rhs = rhsConst * sumA % p;
            points++;
            if (lhs != rhs) mismatches++;

            // corruption control: bump one n_λ by +1 (choose one with a nonzero alternant so it bites).
            int bump = -1;
            for (int c = 0; c < coeffs.Count; c++) if (alt[c] != 0) { bump = c; break; }
            if (bump >= 0)
            {
                long sumBad = (sumA + alt[bump]) % p;
                long rhsBad = rhsConst * sumBad % p;
                ctrlChecks++;
                if (lhs != rhsBad) ctrlMismatches++;
            }
        }
        return new GfpSliceReport(points, mismatches, ctrlChecks, ctrlMismatches);
    }

    // ------------------------------------------------------- the committed table (embedded) --------

    /// <summary>The 143 committed χ^{C₆} coefficients n_λ from
    /// <c>simulations/results/f133_w_closed_form/chiC_coeffs.txt</c>, embedded so the witness is
    /// self-contained. Lines are "λ_parts;n"; the empty partition is the leading ";−1". The aggregate
    /// checksums (count 143, max 8, sum 359, n_() = −1) are asserted in the tests, so a transcription
    /// slip is caught.</summary>
    public static List<(int[] Lam, long N)> TableCoefficients()
    {
        var result = new List<(int[], long)>();
        foreach (var raw in TableText.Split('\n'))
        {
            var line = raw.Trim();
            if (line.Length == 0 || line.StartsWith("#")) continue;
            int semi = line.IndexOf(';');
            var lamPart = line.Substring(0, semi);
            long n = long.Parse(line.Substring(semi + 1));
            int[] lam = lamPart.Length == 0
                ? Array.Empty<int>()
                : lamPart.Split(',').Select(int.Parse).ToArray();
            result.Add((lam, n));
        }
        return result;
    }

    private const string TableText = @";-1
1,1;1
2;-1
1,1,1,1;-1
1,1,1,1,1,1;1
2,1,1,1,1;1
2,2,1,1;-1
2,2,2;1
3,1,1,1;1
3,3;-1
4,1,1;-1
4,2;3
5,1;-2
2,2,1,1,1,1;1
2,2,2,2;2
3,1,1,1,1,1;-3
3,2,2,1;-2
3,3,1,1;5
3,3,2;-3
4,1,1,1,1;2
4,2,1,1;-5
4,2,2;4
4,3,1;-1
4,4;1
5,1,1,1;3
5,3;-3
6,1,1;-1
6,2;3
8;-1
2,2,2,2,1,1;-4
2,2,2,2,2;1
3,2,2,1,1,1;3
3,3,1,1,1,1;-7
3,3,2,1,1;2
3,3,2,2;-7
3,3,3,1;5
4,2,1,1,1,1;4
4,2,2,1,1;-1
4,2,2,2;7
4,3,3;-3
4,4,1,1;-4
4,4,2;4
5,1,1,1,1,1;-4
5,2,2,1;-4
5,3,1,1;6
5,3,2;-1
6,1,1,1,1;2
6,2,1,1;-3
6,2,2;3
6,3,1;-1
6,4;1
7,3;-1
9,1;1
10;-1
2,2,2,2,2,2;4
3,3,2,2,1,1;8
3,3,2,2,2;-4
3,3,3,1,1,1;-6
3,3,3,3;6
4,2,2,2,1,1;-8
4,2,2,2,2;2
4,3,2,2,1;1
4,3,3,2;-5
4,4,1,1,1,1;7
4,4,2,1,1;-3
4,4,2,2;7
4,4,3,1;-2
4,4,4;1
5,2,2,1,1,1;4
5,3,1,1,1,1;-5
5,3,2,2;-4
5,3,3,1;6
5,4,2,1;-2
5,5,1,1;1
5,5,2;-1
6,2,1,1,1,1;1
6,2,2,2;3
6,3,3;-3
6,4,1,1;-1
6,4,2;1
7,1,1,1,1,1;-2
7,2,2,1;-1
7,3,1,1;1
8,1,1,1,1;2
9,1,1,1;-1
3,3,2,2,2,2;-4
3,3,3,2,2,1;1
3,3,3,3,1,1;-7
3,3,3,3,2;3
4,2,2,2,2,2;3
4,3,3,2,1,1;4
4,3,3,2,2;-1
4,4,2,2,1,1;-8
4,4,2,2,2;4
4,4,3,1,1,1;4
4,4,3,3;-3
4,4,4,2;3
5,2,2,2,2,1;1
5,3,2,2,1,1;4
5,3,2,2,2;-2
5,3,3,1,1,1;-4
5,3,3,3;3
5,4,2,1,1,1;1
5,4,4,1;-2
5,5,1,1,1,1;-2
5,5,2,1,1;1
5,5,2,2;-1
5,5,3,1;1
6,2,2,2,1,1;-4
6,3,2,2,1;1
6,3,3,2;-2
6,4,1,1,1,1;1
6,4,2,2;1
6,4,4;1
7,2,2,1,1,1;2
7,3,2,1,1;-1
7,3,3,1;1
8,2,1,1,1,1;-1
3,3,3,3,2,2;2
3,3,3,3,3,1;-1
4,4,2,2,2,2;2
4,4,3,2,2,1;-1
4,4,3,3,1,1;4
4,4,3,3,2;-2
4,4,4,2,1,1;-3
4,4,4,2,2;1
5,3,2,2,2,2;-3
5,3,3,3,1,1;-3
5,3,3,3,2;1
5,4,2,2,2,1;1
5,4,4,1,1,1;1
5,5,2,2,1,1;1
5,5,2,2,2;-1
5,5,3,1,1,1;-1
6,2,2,2,2,2;2
6,3,3,2,1,1;1
6,4,2,2,1,1;-1
3,3,3,3,3,3;1
4,3,3,3,3,2;-1
4,4,3,3,2,2;-1
4,4,3,3,3,1;1
4,4,4,2,2,2;1
5,3,3,3,2,2;1";
}
