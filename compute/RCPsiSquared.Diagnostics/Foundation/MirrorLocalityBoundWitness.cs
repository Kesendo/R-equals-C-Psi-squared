using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the mirror's locality bound under a temperature gradient
/// (<c>experiments/THE_CUBE_SITE_BY_SITE.md</c>, finding 1, the paragraph on the Z corner; claim
/// <see cref="Core.Symmetry.MirrorLocalityBoundClaim"/>). The system: an open chain of N sites, a number-conserving
/// Hamiltonian H = Σ (J/2)(X_lX_{l+1} + Y_lY_{l+1}) + (Δ/2) Z_lZ_{l+1}, and on every site a thermal bath, cooling σ⁻ at
/// γ↓_l and heating σ⁺ at γ↑_l. The argument the witness checks the ingredients of:
///
/// <list type="number">
///   <item>Z⊗N is an eigenvector of L at −c, c = Σ_l (γ↓_l + γ↑_l): each site charges Z the height −(γ↓ + γ↑) and
///   [H, Z⊗N] = 0. Checked here entry for entry in integers (the generator is built with integer entries).</item>
///   <item>Any element E with E·L·E⁻¹ = −L − c carries the −c mode onto the 0 mode: L(E v) = (−λ − c)·E v for
///   L v = λ v, so E(Z⊗N) lies in the kernel of L, which is the steady state's line when the steady state is
///   unique (checked: the kernel has dimension 1).</item>
///   <item>A map of operator Schmidt rank r across a cut turns a product across that cut into an operator of rank at
///   most r, and Z⊗N is a product across every cut. Hence Schmidt-rank(E) ≥ Schmidt-rank(ρ_ss) across every cut,
///   and an invertible product map (rank 1) cannot exist where the steady state is not a product.</item>
/// </list>
///
/// <para>What is recomputed live: the integer generator in the Pauli-string basis (real there, since L preserves
/// Hermiticity), its kernel modulo primes (two at N = 2, more where the denominators are larger, up to
/// <see cref="MaxPrimes"/>), the rational steady state recovered from the residues by the Chinese remainder theorem
/// and rational reconstruction and then VERIFIED exactly (L·ρ_ss = 0 in BigInteger arithmetic and c_{I⊗N} = 2^{−N}
/// in ℚ, so the letters below are exact, not probabilistic), and the operator Schmidt rank of ρ_ss across the cut as
/// an exact rank over ℚ (BigRational elimination on the 4^k × 4^(N−k) coefficient matrix). Where the reconstruction
/// does not close, the witness falls back to the largest rank modulo a prime that did not divide the trace
/// coefficient, a lower bound on the rational rank (the reduction of ρ_ss modulo such a prime is the residue vector),
/// hence still a bound. No eigensolver anywhere.</para>
///
/// <para>Two controls run inside the same object: one temperature (every ratio equal), where the steady state is
/// the product Gibbs state and the rank is 1 (the site-local M_β exists there), and no hopping (J = 0), where
/// the gradient alone leaves a product. Guard: N ≤ <see cref="MaxN"/> (= 4, a 256-dimensional operator space; at
/// N = 5 the steady state's entries reach 469 bits, past what 16 primes reconstruct, and the elimination runs for
/// minutes).</para>
///
/// <para>Second route: <c>simulations/cube_site_by_site_findings.py</c> block G (sympy over ℚ(i), the matrix-unit
/// basis, a different reshaping), whose N = 2 letters ¼, 5/56, 27/224, 67/1568, ±1/112 the tests pin against
/// this witness.</para></summary>
public sealed class MirrorLocalityBoundWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The largest chain the live build admits (a 4^N-dimensional operator space; N = 4 is 256).</summary>
    public const int MaxN = 4;

    /// <summary>The largest number of primes the reconstruction will use (each below 2^30, so products of residues
    /// stay inside a long); the steady state's denominators grow fast with N (about 10^43 at N = 4), and the
    /// reconstruction needs a modulus beyond twice the square of the largest numerator or denominator.</summary>
    public const int MaxPrimes = 16;

    /// <summary>The primes, found by trial division below 2^30 at first use (nothing recalled from memory).</summary>
    public static readonly long[] Primes = FindPrimes(MaxPrimes, (1L << 30) - 1);

    private static long[] FindPrimes(int count, long below)
    {
        var found = new List<long>();
        for (long c = below; found.Count < count; c--) if (IsPrime(c)) found.Add(c);
        return found.ToArray();
    }

    private static bool IsPrime(long n)
    {
        if (n < 2) return false;
        if (n % 2 == 0) return n == 2;
        for (long d = 3; d * d <= n; d += 2) if (n % d == 0) return false;
        return true;
    }

    public int N { get; }
    public int Cut { get; }
    public IReadOnlyList<long> Down { get; }
    public IReadOnlyList<long> Up { get; }
    public long J { get; }
    public long Delta { get; }

    /// <summary>c = Σ_l (γ↓_l + γ↑_l), the Z corner's depth.</summary>
    public long C { get; }

    /// <summary>The generator ×2 in the Pauli-string basis: column s holds the coefficients of 2·L(P_s).</summary>
    public long[,] TwiceGenerator { get; }

    /// <summary>Whether 2L·Z⊗N = −2c·Z⊗N holds entry for entry in integers (the Z corner at the deepest rate).</summary>
    public bool ZCornerExact { get; }

    /// <summary>How many primes the reconstruction consumed (2 at N = 2; more where the denominators are larger).</summary>
    public int PrimesUsed { get; }

    /// <summary>How many primes were skipped on the way, because the kernel modulo p was not one-dimensional or the
    /// trace coefficient vanished modulo p (a bad prime for this generator, not a statement about ℚ).</summary>
    public int PrimesSkipped { get; }

    /// <summary>Kernel dimension modulo each prime used, in the order of <see cref="Primes"/>.</summary>
    public IReadOnlyList<int> KernelDimensionModP { get; }

    /// <summary>The steady state's coefficients on the 4^N Pauli strings, exact over ℚ and verified by L·ρ = 0 in
    /// BigInteger arithmetic, normalised to trace 1 (c_{I⊗N} = 2^{−N}); null where the reconstruction did not close.</summary>
    public BigRational[]? SteadyState { get; }

    /// <summary>The operator Schmidt rank of the steady state across the cut, exact over ℚ when
    /// <see cref="SteadyState"/> is set, else the largest of the ranks modulo the primes used (a lower bound).</summary>
    public int SchmidtRank { get; }

    /// <summary>Whether <see cref="SchmidtRank"/> is the exact rational rank (true) or only a lower bound (false).</summary>
    public bool RankIsExact { get; }

    /// <summary>The ranks modulo each prime used, in the order of <see cref="Primes"/>.</summary>
    public IReadOnlyList<int> SchmidtRankModP { get; }

    /// <summary>4^min(k, N−k), the largest rank an operator can have across the cut.</summary>
    public long MaxRank => IntPow(4, Math.Min(Cut, N - Cut));

    public MirrorLocalityBoundWitness(int n = 2, IReadOnlyList<long>? down = null, IReadOnlyList<long>? up = null,
                                      long j = 1, long delta = 0, int cut = 1)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n), $"N must be ≥ 2 (a cut needs two sides); got {n}");
        if (n > MaxN) throw new ArgumentOutOfRangeException(nameof(n), $"N = {n} exceeds the live-build guard {MaxN}");
        if (cut < 1 || cut >= n) throw new ArgumentOutOfRangeException(nameof(cut), $"cut must lie in 1..N−1; got {cut}");
        // the documented chains as defaults: N = 2 (3,2 | 1,1), N = 3 (3,2,1 | 1,1,1), N = 4 (6,5,4,3 | 2,2,2,2)
        if (down is null && up is null && n == 4) { down = new long[] { 6, 5, 4, 3 }; up = new long[] { 2, 2, 2, 2 }; }
        down ??= Enumerable.Range(0, n).Select(l => (long)(n + 1 - l)).ToArray();
        up ??= Enumerable.Repeat(1L, n).ToArray();
        if (down.Count != n || up.Count != n)
            throw new ArgumentOutOfRangeException(nameof(down), $"one cooling and one heating rate per site: got {down.Count} and {up.Count} for N = {n}");
        for (int l = 0; l < n; l++)
        {
            if (down[l] < 0 || up[l] < 0) throw new ArgumentOutOfRangeException(nameof(down), $"rates must be ≥ 0; site {l} has ({down[l]}, {up[l]})");
            if (down[l] + up[l] == 0) throw new ArgumentOutOfRangeException(nameof(down), $"site {l} has no bath (γ↓ + γ↑ = 0); the steady state would not be unique");
        }

        N = n; Cut = cut; Down = down.ToArray(); Up = up.ToArray(); J = j; Delta = delta;
        C = Down.Zip(Up, (d, u) => d + u).Sum();
        TwiceGenerator = BuildTwiceGenerator(n, Down, Up, j, delta);
        ZCornerExact = CheckZCorner(TwiceGenerator, n, C);

        // the kernel modulo primes, skipping a prime where the kernel is not one-dimensional or the trace
        // coefficient vanishes (a bad prime for this generator, not a statement about ℚ); CRT + reconstruction
        // as soon as two primes are in hand, more until it closes or the primes run out
        var dims = new List<int>();
        var ranksModP = new List<int>();
        var residues = new List<long[]>();
        var primesUsed = new List<long>();
        BigRational[]? exact = null;
        int skipped = 0;
        foreach (long p in Primes)
        {
            var ker = NullspaceModP(TwiceGenerator, p);
            if (ker.Count != 1 || ker[0][0] == 0) { skipped++; continue; }
            dims.Add(ker.Count);
            primesUsed.Add(p);
            var r = NormaliseTraceModP(ker[0], n, p);
            residues.Add(r);
            ranksModP.Add(RankModP(Reshape(r, n, cut), p));
            if (residues.Count >= 2)
            {
                exact = Reconstruct(residues, primesUsed, TwiceGenerator, n);
                if (exact is not null) break;
            }
        }
        if (primesUsed.Count == 0) throw new InvalidOperationException($"no usable prime among {Primes.Length}: the kernel is not one-dimensional (or every prime was bad)");
        PrimesUsed = primesUsed.Count;
        PrimesSkipped = skipped;
        KernelDimensionModP = dims;
        SchmidtRankModP = ranksModP;
        SteadyState = exact;
        if (SteadyState is not null)
        {
            SchmidtRank = RankOverQ(ReshapeQ(SteadyState, n, cut));
            RankIsExact = true;
        }
        else
        {
            SchmidtRank = ranksModP.Max();
            RankIsExact = false;
        }
    }

    // ----- the generator, integer entries in the Pauli-string basis -----

    // letter digits: 0 = I, 1 = X, 2 = Y, 3 = Z; (x, z) bits of each letter
    private static readonly int[] LetterX = { 0, 1, 1, 0 };
    private static readonly int[] LetterZ = { 0, 0, 1, 1 };
    private static int LetterOf(int x, int z) => x == 0 ? (z == 0 ? 0 : 3) : (z == 0 ? 1 : 2);
    // P_a·P_b = i^Phase[a,b]·P_{a⊕b}: XY = iZ, YZ = iX, ZX = iY and their reverses with −i
    private static readonly int[,] Phase =
    {
        { 0, 0, 0, 0 },
        { 0, 0, 1, 3 },
        { 0, 3, 0, 1 },
        { 0, 1, 3, 0 },
    };

    public static long IntPow(long b, int e) { long r = 1; for (int i = 0; i < e; i++) r *= b; return r; }

    private static int Digit(long s, int n, int site) => (int)((s / IntPow(4, n - 1 - site)) % 4);

    /// <summary>The coefficients of 2·L(P_s) for every string s, as the columns of a 4^N × 4^N integer matrix.
    /// Thermal site l on the letter at l: I → 2(γ↓ − γ↑)·Z, X → −(γ↓ + γ↑)·X, Y → −(γ↓ + γ↑)·Y, Z → −2(γ↓ + γ↑)·Z;
    /// a bond term (J/2)·Q with Q ∈ {XX, YY} (and (Δ/2)·ZZ) contributes 2·(−i)[(J/2)Q, P] = −2iJ·QP on the strings P
    /// that anticommute with Q, a real integer multiple of the product string.</summary>
    public static long[,] BuildTwiceGenerator(int n, IReadOnlyList<long> down, IReadOnlyList<long> up, long j, long delta)
    {
        long dim = IntPow(4, n);
        var M = new long[dim, dim];
        for (long s = 0; s < dim; s++)
        {
            // dissipator, site by site
            for (int l = 0; l < n; l++)
            {
                int a = Digit(s, n, l);
                long gammaSum = down[l] + up[l];
                long stride = IntPow(4, n - 1 - l);
                switch (a)
                {
                    case 0: M[s + 3 * stride, s] += 2 * (down[l] - up[l]); break;   // I → Z, the arrow
                    case 1: case 2: M[s, s] -= gammaSum; break;                       // X, Y: height −Γ/2
                    case 3: M[s, s] -= 2 * gammaSum; break;                           // Z: height −Γ
                }
            }
            // Hamiltonian, bond by bond: (J/2)(XX + YY) + (Δ/2) ZZ on (l, l+1)
            for (int l = 0; l + 1 < n; l++)
            {
                foreach (var (q, coeff) in new[] { (1, j), (2, j), (3, delta) })
                {
                    if (coeff == 0) continue;
                    int a0 = Digit(s, n, l), a1 = Digit(s, n, l + 1);
                    bool anti0 = a0 != 0 && a0 != q, anti1 = a1 != 0 && a1 != q;
                    if (anti0 == anti1) continue;                     // commutes (0 or 2 anticommuting sites)
                    int phase = Phase[q, a0] + Phase[q, a1];           // QP = i^phase · R
                    int r0 = LetterOf(LetterX[q] ^ LetterX[a0], LetterZ[q] ^ LetterZ[a0]);
                    int r1 = LetterOf(LetterX[q] ^ LetterX[a1], LetterZ[q] ^ LetterZ[a1]);
                    long stride0 = IntPow(4, n - 1 - l), stride1 = IntPow(4, n - 2 - l);
                    long t = s + (r0 - a0) * stride0 + (r1 - a1) * stride1;
                    // −2i·coeff·i^phase: real exactly when phase is odd (it is, one anticommuting site)
                    int total = (phase + 3) % 4;                       // −i = i^3
                    if (total % 2 != 0) throw new InvalidOperationException("a bond term produced a non-real coefficient");
                    M[t, s] += (total == 0 ? 2 : -2) * coeff;
                }
            }
        }
        return M;
    }

    private static bool CheckZCorner(long[,] M2, int n, long c)
    {
        long dim = IntPow(4, n), corner = dim - 1;
        for (long t = 0; t < dim; t++)
        {
            long expected = t == corner ? -2 * c : 0;
            if (M2[t, corner] != expected) return false;
        }
        return true;
    }

    // ----- modular linear algebra -----

    private static long Mod(long a, long p) { long r = a % p; return r < 0 ? r + p : r; }
    private static long MulMod(long a, long b, long p) => (long)((Int128)a * b % p);
    private static long PowMod(long b, long e, long p) { long r = 1; b = Mod(b, p); for (; e > 0; e >>= 1) { if ((e & 1) == 1) r = MulMod(r, b, p); b = MulMod(b, b, p); } return r; }
    private static long InvMod(long a, long p) => PowMod(a, p - 2, p);

    /// <summary>Row-reduced echelon form modulo p; returns the pivot column per row and the reduced rows.</summary>
    private static (long[][] Rows, List<int> Pivots) ReduceModP(long[][] rows, int cols, long p)
    {
        var pivots = new List<int>();
        int r = 0;
        for (int c = 0; c < cols && r < rows.Length; c++)
        {
            int pr = -1;
            for (int i = r; i < rows.Length; i++) if (rows[i][c] != 0) { pr = i; break; }
            if (pr < 0) continue;
            (rows[r], rows[pr]) = (rows[pr], rows[r]);
            long inv = InvMod(rows[r][c], p);
            for (int k = 0; k < cols; k++) rows[r][k] = MulMod(rows[r][k], inv, p);
            for (int i = 0; i < rows.Length; i++)
            {
                if (i == r || rows[i][c] == 0) continue;
                long f = rows[i][c];
                for (int k = 0; k < cols; k++) rows[i][k] = Mod(rows[i][k] - MulMod(f, rows[r][k], p), p);
            }
            pivots.Add(c);
            r++;
        }
        return (rows, pivots);
    }

    private static long[][] ToRowsModP(long[,] M, long p)
    {
        int n = M.GetLength(0), m = M.GetLength(1);
        var rows = new long[n][];
        for (int i = 0; i < n; i++) { rows[i] = new long[m]; for (int k = 0; k < m; k++) rows[i][k] = Mod(M[i, k], p); }
        return rows;
    }

    public static int RankModP(long[,] M, long p) => ReduceModP(ToRowsModP(M, p), M.GetLength(1), p).Pivots.Count;

    /// <summary>A basis of {v : M v = 0 (mod p)}.</summary>
    public static List<long[]> NullspaceModP(long[,] M, long p)
    {
        int cols = M.GetLength(1);
        var (rows, pivots) = ReduceModP(ToRowsModP(M, p), cols, p);
        var free = Enumerable.Range(0, cols).Where(c => !pivots.Contains(c)).ToList();
        var basis = new List<long[]>();
        foreach (int f in free)
        {
            var v = new long[cols];
            v[f] = 1;
            for (int i = 0; i < pivots.Count; i++) v[pivots[i]] = Mod(-rows[i][f], p);
            basis.Add(v);
        }
        return basis;
    }

    /// <summary>Scale the kernel vector so that c_{I⊗N} = 2^{−N} (trace one) modulo p.</summary>
    private static long[] NormaliseTraceModP(long[] v, int n, long p)
    {
        if (v[0] == 0) throw new InvalidOperationException("the steady state has zero trace coefficient modulo p");
        long f = MulMod(InvMod(v[0], p), InvMod(IntPow(2, n), p), p);
        return v.Select(x => MulMod(x, f, p)).ToArray();
    }

    /// <summary>The coefficient vector as the 4^k × 4^(N−k) matrix c[left string, right string].</summary>
    public static long[,] Reshape(long[] v, int n, int cut)
    {
        long rows = IntPow(4, cut), cols = IntPow(4, n - cut);
        var S = new long[rows, cols];
        for (long i = 0; i < rows; i++) for (long k = 0; k < cols; k++) S[i, k] = v[i * cols + k];
        return S;
    }

    private static BigRational[,] ReshapeQ(BigRational[] v, int n, int cut)
    {
        long rows = IntPow(4, cut), cols = IntPow(4, n - cut);
        var S = new BigRational[rows, cols];
        for (long i = 0; i < rows; i++) for (long k = 0; k < cols; k++) S[i, k] = v[i * cols + k];
        return S;
    }

    /// <summary>Exact rank over ℚ by Gaussian elimination in BigRational.</summary>
    public static int RankOverQ(BigRational[,] A)
    {
        int n = A.GetLength(0), m = A.GetLength(1);
        var rows = new BigRational[n][];
        for (int i = 0; i < n; i++) { rows[i] = new BigRational[m]; for (int k = 0; k < m; k++) rows[i][k] = A[i, k]; }
        int r = 0;
        for (int c = 0; c < m && r < n; c++)
        {
            int pr = -1;
            for (int i = r; i < n; i++) if (!rows[i][c].IsZero) { pr = i; break; }
            if (pr < 0) continue;
            (rows[r], rows[pr]) = (rows[pr], rows[r]);
            for (int i = r + 1; i < n; i++)
            {
                if (rows[i][c].IsZero) continue;
                var f = rows[i][c] / rows[r][c];
                for (int k = c; k < m; k++) rows[i][k] = rows[i][k] - f * rows[r][k];
            }
            r++;
        }
        return r;
    }

    // ----- the exact steady state: CRT over the primes used, rational reconstruction, exact verification -----

    private static BigRational[]? Reconstruct(IReadOnlyList<long[]> residues, IReadOnlyList<long> primes, long[,] M2, int n)
    {
        int dim = residues[0].Length;
        // CRT over the primes used, entry by entry
        BigInteger modulus = BigInteger.One;
        var x = new BigInteger[dim];
        for (int i = 0; i < primes.Count; i++)
        {
            BigInteger p = primes[i];
            if (i == 0) { for (int k = 0; k < dim; k++) x[k] = residues[0][k]; modulus = p; continue; }
            BigInteger inv = BigInteger.ModPow(modulus % p, p - 2, p);   // modulus^{-1} mod p (p prime)
            for (int k = 0; k < dim; k++)
            {
                BigInteger t = (residues[i][k] - x[k] % p) % p; if (t < 0) t += p;
                t = t * inv % p;
                x[k] += t * modulus;
            }
            modulus *= p;
        }
        var result = new BigRational[dim];
        for (int k = 0; k < dim; k++)
        {
            var q = RationalReconstruct(x[k], modulus);
            if (q is null) return null;
            result[k] = q.Value;
        }
        // the scale fixed in ℚ, not only modulo p: trace one means c_{I⊗N} = 2^{−N}
        if (result[0] != new BigRational(1, IntPow(2, n))) return null;
        // exact verification: M2 * result = 0 in BigRational arithmetic (so the letters are exact, not probable)
        for (int t = 0; t < dim; t++)
        {
            var acc = BigRational.Zero;
            for (int k = 0; k < dim; k++)
            {
                long m = M2[t, k];
                if (m != 0) acc = acc + new BigRational(m) * result[k];
            }
            if (!acc.IsZero) return null;
        }
        return result;
    }

    /// <summary>Wang's rational reconstruction: the unique n/d with |n|, d ≤ √(m/2) and n ≡ d·x (mod m), if any.</summary>
    public static BigRational? RationalReconstruct(BigInteger x, BigInteger m)
    {
        BigInteger bound = Isqrt(m / 2);
        BigInteger r0 = m, r1 = x % m, s0 = 0, s1 = 1;
        if (r1 < 0) r1 += m;
        while (r1 > bound)
        {
            BigInteger q = r0 / r1;
            (r0, r1) = (r1, r0 - q * r1);
            (s0, s1) = (s1, s0 - q * s1);
        }
        if (s1.IsZero || BigInteger.Abs(s1) > bound) return null;
        if (s1.Sign < 0) { s1 = -s1; r1 = -r1; }
        if (BigInteger.GreatestCommonDivisor(r1, s1) != 1) return null;
        return new BigRational(r1, s1);
    }

    /// <summary>Floor of the square root by Newton's iteration on BigInteger (a double seed would leave 10^56 steps
    /// of linear correction at a 480-bit modulus).</summary>
    private static BigInteger Isqrt(BigInteger v)
    {
        if (v < 2) return v;
        BigInteger x = BigInteger.One << ((int)(v.GetBitLength() + 1) / 2);
        while (true)
        {
            BigInteger y = (x + v / x) >> 1;
            if (y >= x) return x;
            x = y;
        }
    }

    // ----- reading -----

    /// <summary>The string s as letters, site 0 first.</summary>
    public static string Label(long s, int n)
    {
        const string letters = "IXYZ";
        var chars = new char[n];
        for (int l = 0; l < n; l++) chars[l] = letters[Digit(s, n, l)];
        return new string(chars);
    }

    /// <summary>The nonzero exact coefficients of the steady state, (label, value), or empty when not exact.</summary>
    public IEnumerable<(string Label, BigRational Value)> Letters()
    {
        if (SteadyState is null) yield break;
        for (long s = 0; s < SteadyState.Length; s++)
            if (!SteadyState[s].IsZero) yield return (Label(s, N), SteadyState[s]);
    }

    private string Rates => "γ↓=(" + string.Join(",", Down) + "), γ↑=(" + string.Join(",", Up) + ")";

    public string DisplayName =>
        $"the mirror's locality bound: N={N}, cut {Cut}|{N - Cut}, {Rates}, J={J}, Δ={Delta}";

    public string Summary =>
        $"Z corner Z⊗N at −c = −{C}: {(ZCornerExact ? "exact" : "NOT an eigenvector (construction broken)")}; " +
        $"steady state unique mod p: yes ({PrimesUsed} prime(s)); " +
        $"Schmidt rank of ρ_ss across the cut = {SchmidtRank} of {MaxRank} ({(RankIsExact ? "exact over ℚ" : "lower bound, mod p")}); " +
        $"so every element E with E·L·E⁻¹ = −L − c has Schmidt rank ≥ {SchmidtRank} across this cut" +
        (SchmidtRank > 1 ? ", and no product map is one." : "; a product map is not excluded by this cut.");

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode(
                displayName: "1. the Z corner",
                summary: $"2L·Z⊗N = −2c·Z⊗N checked entry for entry in the integer generator (c = Σ(γ↓ + γ↑) = {C}): " +
                         (ZCornerExact ? "exact. Each site charges Z its height −(γ↓ + γ↑), and the hopping and ZZ terms commute with Z⊗N."
                                         : "FAILS."),
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "2. the steady mode",
                summary: $"kernel dimension modulo {PrimesUsed} prime(s) below 2^30 ({PrimesSkipped} skipped of {PrimesUsed + PrimesSkipped} tried): {string.Join(", ", KernelDimensionModP)}. " +
                         (SteadyState is null
                             ? (PrimesUsed < 2 ? "Reconstruction needs two residues; only one prime was usable, the ranks below are modular lower bounds."
                                               : $"Rational reconstruction did not close with {PrimesUsed} primes; the ranks below are modular lower bounds.")
                             : $"Exact over ℚ (CRT over {PrimesUsed} primes + rational reconstruction, then L·ρ_ss = 0 and the trace verified in BigInteger arithmetic), trace 1: " +
                               string.Join(" + ", Letters().Select(l => $"{l.Value}·{l.Label}")) + "."),
                provenance: NodeProvenance.Live);

            yield return new InspectableNode(
                displayName: "3. the bound",
                summary: $"operator Schmidt rank of ρ_ss across the cut {Cut}|{N - Cut}: {SchmidtRank} of at most {MaxRank} " +
                         $"({(RankIsExact ? "exact over ℚ" : "max of the ranks mod p " + string.Join("/", SchmidtRankModP) + ", a lower bound")}). " +
                         "E(Z⊗N) ∝ ρ_ss for every element of the shape −L − c, and a map of Schmidt rank r turns a product into " +
                         $"rank ≤ r: Schmidt-rank(E) ≥ {SchmidtRank}. " +
                         (SchmidtRank > 1 ? "No invertible product map across this cut, the site-local NONE of finding 1 on a second route."
                                          : "A product is not excluded here."),
                provenance: NodeProvenance.Live);

            // controls, built from the same code on the same chain
            var oneTemperature = new MirrorLocalityBoundWitness(N, Enumerable.Repeat(Down[0], N).ToArray(), Enumerable.Repeat(Up[0], N).ToArray(), J, Delta, Cut);
            yield return new InspectableNode(
                displayName: "4. control: one temperature",
                summary: $"every site at ({Down[0]}, {Up[0]}): Schmidt rank {oneTemperature.SchmidtRank} " +
                         $"({(oneTemperature.RankIsExact ? "exact" : "mod p")}); the product of the one-site steady states" +
                         (Down[0] > 0 && Up[0] > 0 ? ", the Gibbs state, where the site-local M_β exists." : "; at ratio zero or infinity M_β is singular (finding 1)."),
                provenance: NodeProvenance.Live);

            var noHop = new MirrorLocalityBoundWitness(N, Down, Up, 0, 0, Cut);
            yield return new InspectableNode(
                displayName: "5. control: no hopping",
                summary: $"J = 0, Δ = 0, the same rates: Schmidt rank {noHop.SchmidtRank} " +
                         $"({(noHop.RankIsExact ? "exact" : "mod p")}); without the hop the gradient leaves a product and the bound is empty.",
                provenance: NodeProvenance.Live);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
