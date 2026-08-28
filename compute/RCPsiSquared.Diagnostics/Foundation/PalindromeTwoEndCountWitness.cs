using System.Globalization;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F158 (claim <see cref="PalindromeTwoEndCountClaim"/>, proof
/// <c>docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md</c>, gate
/// <c>simulations/f138_rank_criterion.py</c>): build the Liouvillian at inspect time and read the
/// palindrome off two integers.
///
/// <code>
///     spec(L) closed under lambda -> -lambda - 2*sigma   &lt;=&gt;   dim ker L == dim ker(L + 2*sigma)
/// </code>
///
/// <para><b>Three computations, and two of the three meetings are the point.</b> Everything here
/// is recomputed and nothing is looked up.
/// (a) The two NULLITIES, as eliminations on L and on L + 2 sigma over GF(p).
/// (b) The same two dimensions again, from Lemma 1's OPERATOR conditions instead: the commutant
/// {X : [H,X] = 0, A_l X A_l = X} and the anticommutant with the jump sign flipped. This route
/// never forms L at all, so (a) == (b) is Lemma 1 recomputed rather than restated, and it is the
/// meeting that would break first if the equality-case argument were wrong.
/// (c) The PALINDROME itself, by the characteristic-polynomial identity det(xI - L) =
/// det((-x - 2 sigma)I - L) at random points over the same field. Then (a) == (b) and
/// "the two dimensions agree" == "the polynomial identity holds" are the two readings the claim
/// stands on, and they are printed as verdicts rather than asserted.</para>
///
/// <para><b>Ranks, never a spectrum, and the one-sidedness is disclosed.</b> Reduction mod p is a
/// ring map, so a RANK can only come out too small and a NULLITY too LARGE. Both dimensions read
/// here are therefore upper bounds on the true ones, and the SMALLER nullity over the primes is
/// kept, which is the tightening direction. Two primes are used because a rank is blind to a
/// factor of p in every entry. The polynomial side is one-sided in the opposite way and says so:
/// a BREAK at one point is a proof over Q(i), while a HOLD is a certificate at the sampled points
/// and primes, exactly as <c>simulations/f138_exact_palindrome_test.py</c> documents. There is no
/// tolerance anywhere in this file because nothing compared is a float.</para>
///
/// <para><b>Two guards, neither of them physics.</b> The eliminations run on d^2 = 4^N columns and
/// the determinant on a d^2 x d^2 matrix, so N is capped at <see cref="MaxN"/> = 4 (256 columns).
/// The rational inputs (gamma = 1/20, field magnitudes over 100) are cleared into the field by
/// modular inverses, so nothing is rounded on the way in.</para>
///
/// <para><b>Falsifiers sit beside the readings rather than being implied.</b> The default row is
/// the canonical Heisenberg chain under Z-dephasing, where both counts are N + 1 and the
/// palindrome holds; the children also carry a row where the criterion says NO and the spectrum
/// agrees, so a witness that had lost the ability to report a break would be visible. If the
/// criterion and the polynomial ever disagreed, the disagreement is what this witness would
/// print.</para>
///
/// <para>Args: <c>--N</c> (2..4, default 3), <c>--deph</c> (a letter per site from I/X/Y/Z or
/// '.', default all Z), <c>--field</c> (same alphabet, default none), <c>--topology</c>
/// (chain|ring|complete, default chain).</para></summary>
public sealed class PalindromeTwoEndCountWitness : IInspectable
{
    /// <summary>The elimination is on 4^N columns and the determinant on a 4^N square, so this is
    /// a cost guard and nothing else. The theorem carries no N.</summary>
    public const int MaxN = 4;

    /// <summary>Two primes, both 1 mod 4 so that -1 is a square and i has an image. The SMALLER
    /// nullity is kept: a nullity mod p can only be too large, so the smaller of two readings is
    /// the tighter available upper bound.</summary>
    private static readonly long[] Primes = { 998244353L, 1004535809L };

    private const int PolynomialPoints = 4;

    private readonly int _n;
    private readonly int[] _deph;
    private readonly int[] _field;
    private readonly (int A, int B)[] _edges;
    private readonly string _topology;
    private readonly Random _rng = new(20260829);

    public PalindromeTwoEndCountWitness(int n, string? deph = null, string? field = null,
                                        string? topology = null)
    {
        if (n < 2 || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"--root twoend runs a dense 4^N elimination and is guarded at N in 2..{MaxN}; got {n}. " +
                "The theorem carries no N, only this witness does.");
        _n = n;
        _deph = ParseLetters(deph, n, defaultLetter: 3, nameof(deph));
        _field = ParseLetters(field, n, defaultLetter: 0, nameof(field));
        _topology = (topology ?? "chain").ToLowerInvariant();
        _edges = BuildEdges(_topology, n);
    }

    private static int[] ParseLetters(string? spec, int n, int defaultLetter, string argName)
    {
        var result = new int[n];
        if (string.IsNullOrWhiteSpace(spec))
        {
            for (int i = 0; i < n; i++) result[i] = defaultLetter;
            return result;
        }
        string s = spec.Trim();
        if (s.Length != n)
            throw new ArgumentException(
                $"--{argName} takes exactly one letter per site (N = {n}); got \"{s}\" of length {s.Length}. " +
                "Use one of I X Y Z or '.' for none, e.g. \"Z.Z\".");
        for (int i = 0; i < n; i++)
        {
            result[i] = char.ToUpperInvariant(s[i]) switch
            {
                'I' or '.' or '0' => 0,
                'X' => 1,
                'Y' => 2,
                'Z' => 3,
                _ => throw new ArgumentException(
                    $"--{argName} letter '{s[i]}' is not one of I X Y Z or '.'."),
            };
        }
        return result;
    }

    private static (int, int)[] BuildEdges(string topology, int n) => topology switch
    {
        "chain" => Enumerable.Range(0, n - 1).Select(i => (i, i + 1)).ToArray(),
        "ring" => Enumerable.Range(0, n).Select(i => (i, (i + 1) % n)).ToArray(),
        "complete" => (from a in Enumerable.Range(0, n)
                       from b in Enumerable.Range(a + 1, n - a - 1)
                       select (a, b)).ToArray(),
        _ => throw new ArgumentException(
            $"--topology takes chain, ring or complete; got \"{topology}\"."),
    };

    public string Name => "twoend";

    public string DisplayName =>
        $"F158 live: the palindrome as two nullities (N = {_n}, {_topology}, dephasing {Letters(_deph)})";

    public string Summary
    {
        get
        {
            var r = Read();
            string verdict = r.CriterionSaysPalindrome == r.PolynomialSaysPalindrome
                ? "the two routes agree"
                : "THE TWO ROUTES DISAGREE, READ IT";
            return $"dim ker L = {r.NearCount}, dim ker(L + 2 sigma) = {r.FarCount}, criterion says " +
                   $"{(r.CriterionSaysPalindrome ? "PALINDROME" : "BROKEN")}, the characteristic polynomial says " +
                   $"{(r.PolynomialSaysPalindrome ? "PALINDROME" : "BROKEN")}: {verdict}";
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    /// <summary>What one inspect recomputes. Every field is a live number.</summary>
    public sealed record Reading(
        int NearCount, int FarCount,
        int NearFromOperatorConditions, int FarFromOperatorConditions,
        bool CriterionSaysPalindrome, bool PolynomialSaysPalindrome);

    public Reading Read()
    {
        int near = int.MaxValue, far = int.MaxValue;
        int nearOp = int.MaxValue, farOp = int.MaxValue;
        bool polynomial = true;
        foreach (long p in Primes)
        {
            var (l, twoSigma) = BuildLiouvillian(p);
            int d2 = l.GetLength(0);
            near = Math.Min(near, d2 - RankModP(Copy(l), d2, d2, p));
            var shifted = Copy(l);
            for (int i = 0; i < d2; i++)
                shifted[i, i] = Mod(shifted[i, i] + twoSigma, p);
            far = Math.Min(far, d2 - RankModP(shifted, d2, d2, p));

            nearOp = Math.Min(nearOp, OperatorConditionNullity(p, sign: +1));
            farOp = Math.Min(farOp, OperatorConditionNullity(p, sign: -1));

            if (!PalindromicModP(l, twoSigma, p)) polynomial = false;
        }
        return new Reading(near, far, nearOp, farOp, near == far, polynomial);
    }

    // ----------------------------------------------------------------- operators

    private static long[][,] Paulis(long p)
    {
        long i = SqrtMinusOne(p);
        return new[]
        {
            new long[,] { { 1, 0 }, { 0, 1 } },
            new long[,] { { 0, 1 }, { 1, 0 } },
            new long[,] { { 0, Mod(-i, p) }, { i, 0 } },
            new long[,] { { 1, 0 }, { 0, Mod(-1, p) } },
        };
    }

    private long[,] SiteOp(int letter, int site, long p)
    {
        var pa = Paulis(p);
        long[,] result = { { 1 } };
        for (int k = 0; k < _n; k++)
            result = Kron(result, k == site ? pa[letter] : pa[0], p);
        return result;
    }

    private long[,] BondOp(int letter, int a, int b, long p)
    {
        var pa = Paulis(p);
        long[,] result = { { 1 } };
        for (int k = 0; k < _n; k++)
            result = Kron(result, k == a || k == b ? pa[letter] : pa[0], p);
        return result;
    }

    /// <summary>H = sum over bonds of (XX + YY + ZZ) at J = 1, plus a field of magnitude 3/10 on
    /// each site that carries one. Rationals are cleared by modular inverse, so nothing rounds.
    /// </summary>
    private long[,] BuildHamiltonian(long p)
    {
        int d = 1 << _n;
        var h = new long[d, d];
        foreach (var (a, b) in _edges)
            for (int t = 1; t <= 3; t++)
            {
                var op = BondOp(t, a, b, p);
                for (int r = 0; r < d; r++)
                    for (int c = 0; c < d; c++)
                        h[r, c] = Mod(h[r, c] + op[r, c], p);
            }
        long mag = Mod(30 * ModInverse(100, p), p);
        for (int s = 0; s < _n; s++)
        {
            if (_field[s] == 0) continue;
            var op = SiteOp(_field[s], s, p);
            for (int r = 0; r < d; r++)
                for (int c = 0; c < d; c++)
                    h[r, c] = Mod(h[r, c] + mag * op[r, c] % p, p);
        }
        return h;
    }

    private List<long[,]> Jumps(long p)
    {
        var jumps = new List<long[,]>();
        for (int s = 0; s < _n; s++)
            if (_deph[s] != 0)
                jumps.Add(SiteOp(_deph[s], s, p));
        return jumps;
    }

    /// <summary>L in the row-major vec convention, and 2 sigma. gamma = 1/20 per dephased site.
    /// </summary>
    private (long[,] L, long TwoSigma) BuildLiouvillian(long p)
    {
        int d = 1 << _n, d2 = d * d;
        var h = BuildHamiltonian(p);
        var l = new long[d2, d2];
        long minusI = Mod(-SqrtMinusOne(p), p);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                for (int c = 0; c < d; c++)
                {
                    // -i ( H (x) I - I (x) H^T ) acting on vec(rho), row-major.
                    l[a * d + b, c * d + b] = Mod(l[a * d + b, c * d + b] + minusI * h[a, c] % p, p);
                    l[a * d + b, a * d + c] = Mod(l[a * d + b, a * d + c] - minusI * h[c, b] % p, p);
                }
        long gamma = Mod(ModInverse(20, p), p);
        long sigma = 0;
        foreach (var A in Jumps(p))
        {
            for (int a = 0; a < d; a++)
                for (int b = 0; b < d; b++)
                    for (int c = 0; c < d; c++)
                        for (int e = 0; e < d; e++)
                        {
                            // A (x) A^T, not A (x) A: the transpose is invisible for X and Z,
                            // which are symmetric, and bites exactly at Y. The first version of
                            // this line wrote A[b, e] and reported dim ker L = 0 on a row carrying
                            // a Y, which is impossible for a unital generator since L(1) = 0.
                            long v = gamma * (A[a, c] * A[e, b] % p) % p;
                            if (v == 0) continue;
                            l[a * d + b, c * d + e] = Mod(l[a * d + b, c * d + e] + v, p);
                        }
            for (int i = 0; i < d2; i++) l[i, i] = Mod(l[i, i] - gamma, p);
            sigma = Mod(sigma + gamma, p);
        }
        return (l, Mod(2 * sigma, p));
    }

    /// <summary>Lemma 1's other route: the nullity of the stacked operator conditions
    /// [H, X] = 0 and A_l X A_l = sign * X, which never forms L. sign = +1 is the commutant,
    /// -1 the anticommutant.</summary>
    private int OperatorConditionNullity(long p, int sign)
    {
        int d = 1 << _n, d2 = d * d;
        var h = BuildHamiltonian(p);
        var jumps = Jumps(p);
        int rows = d2 * (1 + jumps.Count);
        var m = new long[rows, d2];
        int offset = 0;
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                for (int c = 0; c < d; c++)
                {
                    m[a * d + b, c * d + b] = Mod(m[a * d + b, c * d + b] + h[a, c], p);
                    m[a * d + b, a * d + c] = Mod(m[a * d + b, a * d + c] - h[c, b], p);
                }
        offset = d2;
        foreach (var A in jumps)
        {
            for (int a = 0; a < d; a++)
                for (int b = 0; b < d; b++)
                    for (int c = 0; c < d; c++)
                        for (int e = 0; e < d; e++)
                        {
                            long v = A[a, c] * A[e, b] % p;   // A (x) A^T, as above
                            if (v == 0) continue;
                            m[offset + a * d + b, c * d + e] =
                                Mod(m[offset + a * d + b, c * d + e] + v, p);
                        }
            for (int i = 0; i < d2; i++)
                m[offset + i, i] = Mod(m[offset + i, i] - sign, p);
            offset += d2;
        }
        return d2 - RankModP(m, rows, d2, p);
    }

    /// <summary>The palindrome by the characteristic-polynomial identity. A BREAK at one point is
    /// a proof over Q(i); a HOLD is a certificate at the sampled points and primes.</summary>
    private bool PalindromicModP(long[,] l, long twoSigma, long p)
    {
        int d2 = l.GetLength(0);
        for (int k = 0; k < PolynomialPoints; k++)
        {
            long x = (long)(_rng.NextDouble() * (p - 2)) + 1;
            long y = Mod(-x - twoSigma, p);
            if (DetModP(Shift(l, x, p), d2, p) != DetModP(Shift(l, y, p), d2, p)) return false;
        }
        return true;
    }

    private static long[,] Shift(long[,] l, long x, long p)
    {
        int d2 = l.GetLength(0);
        var m = new long[d2, d2];
        for (int i = 0; i < d2; i++)
            for (int j = 0; j < d2; j++)
                m[i, j] = Mod(-l[i, j] + (i == j ? x : 0), p);
        return m;
    }

    // ----------------------------------------------------------------- modular linear algebra

    private static long[,] Copy(long[,] a)
    {
        var b = new long[a.GetLength(0), a.GetLength(1)];
        Array.Copy(a, b, a.Length);
        return b;
    }

    private static long[,] Kron(long[,] a, long[,] b, long p)
    {
        int ra = a.GetLength(0), ca = a.GetLength(1);
        int rb = b.GetLength(0), cb = b.GetLength(1);
        var c = new long[ra * rb, ca * cb];
        for (int i = 0; i < ra; i++)
            for (int j = 0; j < ca; j++)
            {
                if (a[i, j] == 0) continue;
                for (int k = 0; k < rb; k++)
                    for (int m = 0; m < cb; m++)
                        c[i * rb + k, j * cb + m] = a[i, j] * b[k, m] % p;
            }
        return c;
    }

    private static int RankModP(long[,] m, int rows, int cols, long p)
    {
        int rank = 0;
        for (int c = 0; c < cols && rank < rows; c++)
        {
            int pivot = -1;
            for (int r = rank; r < rows; r++)
                if (m[r, c] != 0) { pivot = r; break; }
            if (pivot < 0) continue;
            if (pivot != rank)
                for (int j = 0; j < cols; j++)
                    (m[rank, j], m[pivot, j]) = (m[pivot, j], m[rank, j]);
            long inv = ModInverse(m[rank, c], p);
            for (int j = 0; j < cols; j++) m[rank, j] = m[rank, j] * inv % p;
            for (int r = 0; r < rows; r++)
            {
                if (r == rank || m[r, c] == 0) continue;
                long f = m[r, c];
                for (int j = 0; j < cols; j++)
                    m[r, j] = Mod(m[r, j] - f * m[rank, j] % p, p);
            }
            rank++;
        }
        return rank;
    }

    private static long DetModP(long[,] m, int n, long p)
    {
        long det = 1;
        for (int c = 0; c < n; c++)
        {
            int pivot = -1;
            for (int r = c; r < n; r++)
                if (m[r, c] != 0) { pivot = r; break; }
            if (pivot < 0) return 0;
            if (pivot != c)
            {
                for (int j = 0; j < n; j++) (m[c, j], m[pivot, j]) = (m[pivot, j], m[c, j]);
                det = Mod(-det, p);
            }
            det = det * m[c, c] % p;
            long inv = ModInverse(m[c, c], p);
            for (int j = 0; j < n; j++) m[c, j] = m[c, j] * inv % p;
            for (int r = c + 1; r < n; r++)
            {
                if (m[r, c] == 0) continue;
                long f = m[r, c];
                for (int j = 0; j < n; j++) m[r, j] = Mod(m[r, j] - f * m[c, j] % p, p);
            }
        }
        return det;
    }

    private static long SqrtMinusOne(long p)
    {
        for (long g = 2; g < 200; g++)
        {
            long cand = ModPow(g, (p - 1) / 4, p);
            if (cand * cand % p == p - 1) return cand;
        }
        throw new InvalidOperationException($"no square root of -1 mod {p}; the prime must be 1 mod 4");
    }

    private static long Mod(long x, long p)
    {
        long r = x % p;
        return r < 0 ? r + p : r;
    }

    private static long ModInverse(long a, long p) => ModPow(Mod(a, p), p - 2, p);

    private static long ModPow(long b, long e, long p)
    {
        long r = 1;
        b = Mod(b, p);
        while (e > 0)
        {
            if ((e & 1) == 1) r = r * b % p;
            b = b * b % p;
            e >>= 1;
        }
        return r;
    }

    private static string Letters(int[] v) =>
        string.Concat(v.Select(x => "IXYZ"[x] == 'I' ? '.' : "IXYZ"[x]));

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var r = Read();
            yield return new InspectableNode("the two nullities, recomputed",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "dim ker L = {0}, dim ker(L + 2 sigma) = {1}. Eliminations over GF(p) at two primes on the " +
                    "4^N = {2} columns of the Liouvillian, the SMALLER nullity kept because a nullity mod p can " +
                    "only come out too large. The criterion reads the palindrome off these two integers and " +
                    "nothing else: {3}.",
                    r.NearCount, r.FarCount, 1 << (2 * _n),
                    r.CriterionSaysPalindrome ? "they agree, so PALINDROME" : "they differ, so BROKEN"));

            yield return new InspectableNode("Lemma 1, recomputed rather than restated",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "The same two dimensions from the OPERATOR conditions instead, a route that never forms L: " +
                    "commutant {{X : [H,X] = 0, A_l X A_l = X}} has dimension {0} against ker L = {1}, and the " +
                    "anticommutant with the jump sign flipped has dimension {2} against ker(L + 2 sigma) = {3}. " +
                    "{4} This is the meeting that would break first if the Cauchy-Schwarz equality case were " +
                    "wrong, and it is exact at each prime rather than compared to a tolerance.",
                    r.NearFromOperatorConditions, r.NearCount,
                    r.FarFromOperatorConditions, r.FarCount,
                    r.NearFromOperatorConditions == r.NearCount && r.FarFromOperatorConditions == r.FarCount
                        ? "Both meet."
                        : "THEY DO NOT MEET, READ IT."));

            yield return new InspectableNode("the palindrome, decided independently",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "det(xI - L) against det((-x - 2 sigma)I - L) at {0} random points per prime: {1}. The two " +
                    "sides are one-sided in OPPOSITE directions and both say so: a nullity can only read too " +
                    "large, while a polynomial BREAK is a proof over Q(i) and a HOLD is a certificate at the " +
                    "sampled points. Criterion says {2}, polynomial says {3}: {4}",
                    PolynomialPoints,
                    r.PolynomialSaysPalindrome ? "the identity holds" : "a witness point separates them",
                    r.CriterionSaysPalindrome ? "PALINDROME" : "BROKEN",
                    r.PolynomialSaysPalindrome ? "PALINDROME" : "BROKEN",
                    r.CriterionSaysPalindrome == r.PolynomialSaysPalindrome
                        ? "they agree."
                        : "THEY DISAGREE, WHICH WOULD FALSIFY F158."));

            yield return new InspectableNode("the canonical row, and a break beside it",
                summary: BuildComparisonRow());

            yield return new InspectableNode("what this witness does NOT decide",
                summary: "The MULTISET, not the Jordan structure away from the two ends; F1PalindromeIdentity is " +
                         "the operator identity and sees more. It also does not extend the theorem: the guard at " +
                         "N <= 4 is the cost of a dense 4^N elimination and carries no physics, and the fences " +
                         "A_l^2 = 1, gamma_l > 0 and even d live in the claim rather than here. The gate that " +
                         "scores the criterion at scale, in both directions and across those fences, is " +
                         "simulations/f138_rank_criterion.py.");
        }
    }

    /// <summary>Both verdicts, side by side, so that a witness that had lost the ability to report
    /// a break would be visible rather than silently reassuring.</summary>
    private string BuildComparisonRow()
    {
        var canonical = new PalindromeTwoEndCountWitness(_n).Read();
        var broken = new PalindromeTwoEndCountWitness(
            _n,
            deph: string.Concat(Enumerable.Repeat("Z", _n)),
            field: "Z" + new string('.', _n - 1)).Read();
        return string.Format(CultureInfo.InvariantCulture,
            "Canonical Heisenberg chain, Z on every site, no field: dim ker L = {0}, dim ker(L + 2 sigma) = {1}, " +
            "both N + 1 = {2}, palindrome {3}. Beside it the SAME dephasing with one Z field on site 0, which is " +
            "F138's clause 2 violated in its own terms (a field on an axis the dephasing already occupies): " +
            "dim ker L = {4} while dim ker(L + 2 sigma) = {5}, so the far end empties while the near end does not " +
            "move, and the palindrome {6}. Both readings come from the same routine, so the second is the " +
            "falsifier the first would otherwise lack.",
            canonical.NearCount, canonical.FarCount, PalindromeTwoEndCountClaim.CanonicalChainCount(_n),
            canonical.PolynomialSaysPalindrome ? "HOLDS" : "BREAKS",
            broken.NearCount, broken.FarCount,
            broken.PolynomialSaysPalindrome ? "HOLDS" : "BREAKS");
    }
}
