using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The triangle's palindromic locus at N = 3, live (claim <see cref="PalindromeComplementConnectionClaim"/>,
/// sections "The triangle's palindromic locus at N = 3" and "The triangle with planar fields" of
/// <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>). Heisenberg triangle on (d, u, v) = (0, 1, 2) with
/// couplings J_du, J_dv, J_uv (at most one zero), the jump Z on d alone. Two settings, two rules.
///
/// <para><b>Letter fields</b>, a single-letter field h_l·σ_{w_l} on every site. A row pairs exactly when it lies in</para>
/// <code>
///     A  a colouring: every nonzero field along one letter p in {X, Y};
///     B  J_du = J_dv and a pi rotation about an axis n ⊥ z with R_n f_d = f_d, R_n f_u = f_v;
///     C  e₂ = J_duJ_dv + J_duJ_uv + J_dvJ_uv = 0, f_d = 0, f_u and f_v ⊥ z, |f_u| = |f_v|.
/// </code>
///
/// <para><b>Planar fields</b>, no field on d and fields f_u, f_v of any direction perpendicular to z. A row pairs
/// exactly when it lies in</para>
/// <code>
///     A  f_u ∥ f_v (one of them zero included);
///     B  J_du = J_dv and |f_u| = |f_v|;
///     C  e₂ = 0 and |f_u| = |f_v|, at any angle;
///     E  R1 = J_uv(J_du − J_dv)² − J_duJ_dv(J_du + J_dv) = 0, R2 = f_u·f_v + 2e₂ = 0,
///        R3 = (J_du + J_uv)|f_u|² + (J_dv + J_uv)|f_v|² + 2(J_du + J_dv + J_uv)(f_u·f_v) = 0.
/// </code>
///
/// <para><b>A road of its own.</b> The proofs are a Python elimination with a stored certificate
/// (<c>simulations/n3_palindrome_covering.py</c>, letter fields) and a symbolic sympy derivation
/// (<c>simulations/n3_planar_field_locus_gate.py</c>, planar fields). This witness shares neither: it builds the two
/// ends of F158 itself, as the real kernels of W ↦ [H, W] on the Pauli strings whose d-letter is in {I, Z} (near)
/// or {X, Y} (far), with <see cref="BigRational"/> entries and an exact nullspace, and evaluates the classes on the
/// point. A field with two components is two Pauli terms on its site, so the same kernels read both settings. It
/// recomputes the rule at the row asked for and on two fixed grids: special loci over all 27 letter words, and
/// planar rows on the swap, reciprocal and R1 couplings with parallel, equal-length, generic and E fields and their
/// perturbations, and a row with R2 = R3 = 0 but R1 ≠ 0 that must break. Rational rows of E with J_uv ≠ 0 are scarce;
/// the grid knows one, w = 1, ξ = 8 in the proof's chart, J = (1, 2, 6), and carries E through that row's images under the mirror η₂ → −η₂, a rational turn about z, the swap u ↔ v and a common scale, and
/// through the rows of its corner J_uv = 0, f_u·f_v = 2J_du², which are rational at every rational κ. It is
/// evidence that the rules hold where they read, not the exhaustion, which stays the certificate's and the
/// derivation's.</para>
///
/// <para><b>A field on d.</b> The planar setting also takes a field f_d perpendicular to z on d. The kernels read it as
/// they read every field, two Pauli terms on its site; no rule is stated there (the classes above are proven for
/// f_d = 0), so such a row reports its two ends and its far element and no class.</para>
///
/// <para><b>The far element's locality.</b> Where the row pairs through a one-dimensional far end the witness keeps
/// its element W and reads its operator Schmidt ranks across d | uv, u | dv and v | du (<see cref="OperatorSchmidtRank"/>).
/// A mirror, M·L·M⁻¹ = −L − 2σ or the reflection onto −L† − 2σ (the jump is Hermitian, so ker(L† + 2σ) = ker(L + 2σ)),
/// sends ker L onto that far end, so it sends the identity to a multiple of W, and across every cut its own Schmidt
/// rank is at least W's. A mirror that factorises over the sites, linear or antilinear, after a frame change per site,
/// a transposition or a relabelling of the sites by conjugation, sends the identity to a product. So where W is not a
/// product no such mirror carries the row: the lemma of F116 §4 and, in rank form, of finding 1 of
/// <c>experiments/THE_CUBE_SITE_BY_SITE.md</c>. A mirror conditioned on the bits of a site is not covered (it sends
/// the identity to a sum over sectors), and a one-sided multiplication by a site permutation, such as B's
/// SWAP_uv·(m·σ)^⊗3, is not a relabelling.</para>
///
/// <para>Args: <c>--word</c> (three letters X/Y/Z, one per site d, u, v; default ZXY), <c>--J</c> (J_du,J_dv,J_uv,
/// integers or p/q; default 3,6,-2), <c>--h</c> (h_d,h_u,h_v; default 0,1,1); or, for planar fields, <c>--fu</c> and
/// <c>--fv</c> (x,y each, integers or p/q) with <c>--J</c>, optionally <c>--fd</c> (x,y) for a field on d, and no
/// <c>--word</c> or <c>--h</c>.</para></summary>
public sealed class TrianglePalindromeLocusWitness : IInspectable
{
    private readonly string _word;
    private readonly BigRational[] _j, _h;
    private readonly BigRational[]? _fu, _fv, _fd;

    public TrianglePalindromeLocusWitness(string? word = null, string? couplings = null, string? fields = null,
        string? fieldU = null, string? fieldV = null, string? fieldD = null)
    {
        _j = ParseTuple(couplings ?? "3,6,-2", "J", 3);
        if (_j.Count(x => x.IsZero) > 1)
            throw new ArgumentException("--J leaves the triangle disconnected (two couplings zero); the locus is stated for connected rows.");
        if (fieldU is not null || fieldV is not null || fieldD is not null)
        {
            if (word is not null || fields is not null)
                throw new ArgumentException("--fu/--fv/--fd read the planar setting (fields perpendicular to z); " +
                                            "--word/--h the letter one. Give one of the two.");
            _fu = ParseTuple(fieldU ?? "0,0", "fu", 2);
            _fv = ParseTuple(fieldV ?? "0,0", "fv", 2);
            _fd = ParseTuple(fieldD ?? "0,0", "fd", 2);
            _word = "";
            _h = new[] { BigRational.Zero, BigRational.Zero, BigRational.Zero };
            return;
        }
        _word = (word ?? "ZXY").Trim().ToUpperInvariant();
        if (_word.Length != 3 || _word.Any(c => c is not ('X' or 'Y' or 'Z')))
            throw new ArgumentException($"--word takes three letters from X, Y, Z (sites d, u, v); got \"{word}\".");
        _h = ParseTuple(fields ?? "0,1,1", "h", 3);
    }

    private static BigRational[] ParseTuple(string spec, string arg, int count)
    {
        var parts = spec.Split(',', StringSplitOptions.TrimEntries);
        if (parts.Length != count) throw new ArgumentException($"--{arg} takes {count} comma-separated values; got \"{spec}\".");
        return parts.Select(p =>
        {
            var q = p.Split('/');
            var num = BigInteger.Parse(q[0], CultureInfo.InvariantCulture);
            var den = q.Length > 1 ? BigInteger.Parse(q[1], CultureInfo.InvariantCulture) : BigInteger.One;
            return new BigRational(num, den);
        }).ToArray();
    }

    // ---- the two ends, exactly --------------------------------------------------------------------------

    private static (int Phase, int Letter) Mul1(int a, int b)
    {
        if (a == 0) return (0, b);
        if (b == 0) return (0, a);
        if (a == b) return (0, 0);
        bool cyclic = (a, b) is (1, 2) or (2, 3) or (3, 1);
        return (cyclic ? 1 : 3, 6 - a - b);
    }

    private static bool Anticommute(int[] s, int[] t)
    {
        int n = 0;
        for (int i = 0; i < 3; i++) if (s[i] != 0 && t[i] != 0 && s[i] != t[i]) n++;
        return (n & 1) == 1;
    }

    /// <summary>The terms of H: (coefficient, Pauli string over 0 = I, 1 = X, 2 = Y, 3 = Z); the field of site l is
    /// the vector f[l] = (x, y, z), one term per component.</summary>
    private static List<(BigRational C, int[] S)> Terms(BigRational[] j, BigRational[][] f)
    {
        var terms = new List<(BigRational, int[])>();
        var bonds = new[] { (0, 1), (0, 2), (1, 2) };
        for (int b = 0; b < 3; b++)
            for (int c = 1; c <= 3; c++)
            {
                var s = new int[3]; s[bonds[b].Item1] = c; s[bonds[b].Item2] = c;
                terms.Add((j[b], s));
            }
        for (int l = 0; l < 3; l++)
            for (int c = 1; c <= 3; c++)
            {
                var s = new int[3]; s[l] = c;
                terms.Add((f[l][c - 1], s));
            }
        return terms;
    }

    /// <summary>dim of {W real on the strings with d-letter in {I, Z} (near) or {X, Y} (far) : [H, W] = 0}.</summary>
    public static int End(BigRational[] j, BigRational[][] f, bool far) => EndBasis(j, f, far).Count;

    /// <summary>A basis of that end, each W as its 64 real coefficients on the Pauli strings (index 16·a + 4·b + c
    /// for the letters a, b, c of d, u, v; 0 = I, 1 = X, 2 = Y, 3 = Z). [t, s] = 2·i^k·u for anticommuting strings
    /// with k odd; divided by 2i every entry is ± a coefficient.</summary>
    public static List<BigRational[]> EndBasis(BigRational[] j, BigRational[][] f, bool far)
    {
        var cols = new List<int[]>();
        for (int a = 0; a < 4; a++)
            for (int b = 0; b < 4; b++)
                for (int c = 0; c < 4; c++)
                    if ((a is 1 or 2) == far) cols.Add(new[] { a, b, c });
        var rows = new Dictionary<int, BigRational[]>();
        foreach (var (coef, t) in Terms(j, f))
        {
            if (coef.IsZero) continue;
            for (int col = 0; col < cols.Count; col++)
            {
                var s = cols[col];
                if (!Anticommute(t, s)) continue;
                int phase = 0, key = 0;
                for (int i = 0; i < 3; i++)
                {
                    var (p, letter) = Mul1(t[i], s[i]);
                    phase += p; key = key * 4 + letter;
                }
                if (!rows.TryGetValue(key, out var row)) rows[key] = row = Enumerable.Repeat(BigRational.Zero, cols.Count).ToArray();
                row[col] = row[col] + ((phase & 3) == 1 ? coef : -coef);
            }
        }
        List<BigRational[]> kernel;
        if (rows.Count == 0)
            kernel = Enumerable.Range(0, cols.Count).Select(c =>
            {
                var e = Enumerable.Repeat(BigRational.Zero, cols.Count).ToArray(); e[c] = BigRational.One; return e;
            }).ToList();
        else
        {
            var m = new BigRational[rows.Count, cols.Count];
            int r = 0;
            foreach (var row in rows.Values) { for (int c = 0; c < cols.Count; c++) m[r, c] = row[c]; r++; }
            kernel = BigRationalLinearAlgebra.Nullspace(m);
        }
        return kernel.Select(k =>
        {
            var w = Enumerable.Repeat(BigRational.Zero, 64).ToArray();
            for (int c = 0; c < cols.Count; c++) w[16 * cols[c][0] + 4 * cols[c][1] + cols[c][2]] = k[c];
            return w;
        }).ToList();
    }

    /// <summary>The far element's operator Schmidt ranks across d | uv, u | dv and v | du, where the row pairs
    /// through a one-dimensional far end (else null: a row that does not pair has no mirror at all, and a larger far
    /// end leaves the identity's image open).</summary>
    public sealed record FarLocality(int DUv, int UDv, int VDu)
    {
        public bool Product => DUv == 1 && UDv == 1 && VDu == 1;
    }

    /// <summary>The kind patterns of the far element where the row pairs through one far element (else null): per
    /// string the bits (d, u, v), 0 for a letter in {I, Z} and 1 for one in {X, Y}, the per-site share of the cube's
    /// k_Z; a product's patterns form a product set. Sorted, distinct.</summary>
    public static IReadOnlyList<string>? ReadFarPatterns(BigRational[] j, BigRational[][] f)
    {
        var far = EndBasis(j, f, true);
        if (far.Count != 1 || End(j, f, false) != 1) return null;
        static char Kind(int letter) => letter is 1 or 2 ? '1' : '0';
        return Enumerable.Range(0, 64).Where(s => !far[0][s].IsZero)
            .Select(s => $"{Kind(s / 16)}{Kind(s / 4 % 4)}{Kind(s % 4)}").Distinct().OrderBy(p => p, StringComparer.Ordinal).ToList();
    }

    public static FarLocality? ReadFarLocality(BigRational[] j, BigRational[][] f)
    {
        var far = EndBasis(j, f, true);
        if (far.Count != 1 || End(j, f, false) != 1) return null; // near = far = 1: the row pairs, one far element
        var w = far[0];
        return new FarLocality(OperatorSchmidtRank.Across(w, 3, new[] { 0 }), OperatorSchmidtRank.Across(w, 3, new[] { 1 }),
            OperatorSchmidtRank.Across(w, 3, new[] { 2 }));
    }

    /// <summary>The letter setting's end: site l carries h_l along the letter word[l].</summary>
    public static int End(string word, BigRational[] j, BigRational[] h, bool far) =>
        End(j, Enumerable.Range(0, 3).Select(l => Vec(word[l], h[l])).ToArray(), far);

    /// <summary>The planar setting's fields as site vectors: none on d, (x, y, 0) on u and v.</summary>
    public static BigRational[][] PlanarFields(BigRational[] fu, BigRational[] fv) =>
        PlanarFields(new[] { BigRational.Zero, BigRational.Zero }, fu, fv);

    /// <summary>The planar setting with a field on d as well: (x, y, 0) on every site.</summary>
    public static BigRational[][] PlanarFields(BigRational[] fd, BigRational[] fu, BigRational[] fv) => new[]
    {
        new[] { fd[0], fd[1], BigRational.Zero },
        new[] { fu[0], fu[1], BigRational.Zero },
        new[] { fv[0], fv[1], BigRational.Zero },
    };

    // ---- the classes ------------------------------------------------------------------------------------

    /// <summary>The letter class a row lies in (A, B or C, the first that holds), or null. <paramref name="without"/>
    /// leaves one class out, for the planted mutation in the tests.</summary>
    public static char? Class(string word, BigRational[] j, BigRational[] h, char? without = null)
    {
        var f = Enumerable.Range(0, 3).Select(l => Vec(word[l], h[l])).ToArray();
        foreach (int p in new[] { 0, 1 })
            if (without != 'A' && f.All(v => Enumerable.Range(0, 3).All(q => q == p || v[q].IsZero))) return 'A';
        var rotations = new Func<BigRational[], BigRational[]>[]
        {
            v => new[] { v[0], -v[1], -v[2] }, v => new[] { -v[0], v[1], -v[2] },
            v => new[] { v[1], v[0], -v[2] }, v => new[] { -v[1], -v[0], -v[2] },
        };
        if (without != 'B' && j[0] == j[1] && rotations.Any(R => Same(R(f[0]), f[0]) && Same(R(f[1]), f[2]))) return 'B';
        var e2 = j[0] * j[1] + j[0] * j[2] + j[1] * j[2];
        if (without != 'C' && e2.IsZero && f[0].All(x => x.IsZero) && f[1][2].IsZero && f[2][2].IsZero
            && Norm(f[1]) == Norm(f[2])) return 'C';
        return null;
    }

    /// <summary>The planar class a row lies in (A, B, C or E, the first that holds), or null; <paramref name="without"/>
    /// as in <see cref="Class"/>.</summary>
    public static char? PlanarClass(BigRational[] j, BigRational[] fu, BigRational[] fv, char? without = null)
    {
        BigRational s = j[0], t = j[1], r = j[2];
        var nu = fu[0] * fu[0] + fu[1] * fu[1];
        var nv = fv[0] * fv[0] + fv[1] * fv[1];
        var dot = fu[0] * fv[0] + fu[1] * fv[1];
        var e2 = s * t + s * r + t * r;
        if (without != 'A' && (fu[0] * fv[1] - fu[1] * fv[0]).IsZero) return 'A';
        if (without != 'B' && s == t && nu == nv) return 'B';
        if (without != 'C' && e2.IsZero && nu == nv) return 'C';
        var r1 = r * (s - t) * (s - t) - s * t * (s + t);
        var r2 = dot + 2 * e2;
        var r3 = (s + r) * nu + (t + r) * nv + 2 * (s + t + r) * dot;
        if (without != 'E' && r1.IsZero && r2.IsZero && r3.IsZero) return 'E';
        return null;
    }

    private static BigRational[] Vec(char letter, BigRational h)
    {
        var v = new[] { BigRational.Zero, BigRational.Zero, BigRational.Zero };
        v["XYZ".IndexOf(letter)] = h;
        return v;
    }

    private static bool Same(BigRational[] a, BigRational[] b) => a.Zip(b).All(p => p.First == p.Second);
    private static BigRational Norm(BigRational[] v) => v.Aggregate(BigRational.Zero, (s, x) => s + x * x);

    // ---- one reading ------------------------------------------------------------------------------------

    /// <summary>One row read: its two ends, the class it lies in, whether the rule and the ends agree. With a field
    /// on d no rule is stated: <c>RuleStated</c> is false, <c>InClass</c> null, and <c>Agrees</c> carries nothing (a
    /// grid refuses such rows).</summary>
    public sealed record Reading(int Near, int Far, bool Palindrome, char? InClass, bool Agrees, bool RuleStated = true);

    public static Reading ReadRow(string word, BigRational[] j, BigRational[] h, char? without = null)
    {
        int near = End(word, j, h, false), far = End(word, j, h, true);
        bool pal = far >= 1 && near == far;
        var cls = Class(word, j, h, without);
        return new Reading(near, far, pal, cls, pal == (cls is not null));
    }

    /// <summary>A planar row with a field on d: its two ends and no rule (see <see cref="Reading"/>); with f_d = 0
    /// the planar rule as in the overload without it.</summary>
    public static Reading ReadPlanarRow(BigRational[] j, BigRational[] fd, BigRational[] fu, BigRational[] fv)
    {
        if (fd.All(x => x.IsZero)) return ReadPlanarRow(j, fu, fv);
        var f = PlanarFields(fd, fu, fv);
        int near = End(j, f, false), far = End(j, f, true);
        return new Reading(near, far, far >= 1 && near == far, null, true, RuleStated: false);
    }

    public static Reading ReadPlanarRow(BigRational[] j, BigRational[] fu, BigRational[] fv, char? without = null)
    {
        var f = PlanarFields(fu, fv);
        int near = End(j, f, false), far = End(j, f, true);
        bool pal = far >= 1 && near == far;
        var cls = PlanarClass(j, fu, fv, without);
        return new Reading(near, far, pal, cls, pal == (cls is not null));
    }

    private Reading? _reading;
    public Reading Read() => _reading ??= _fu is null ? ReadRow(_word, _j, _h) : ReadPlanarRow(_j, _fd!, _fu, _fv!);

    private BigRational[][] Fields() => _fu is null
        ? Enumerable.Range(0, 3).Select(l => Vec(_word[l], _h[l])).ToArray()
        : PlanarFields(_fd!, _fu, _fv!);

    private FarLocality? _locality;
    private bool _localityRead;
    public FarLocality? ReadFarLocality()
    {
        if (!_localityRead) { _locality = ReadFarLocality(_j, Fields()); _localityRead = true; }
        return _locality;
    }

    // ---- the grids --------------------------------------------------------------------------------------

    private static BigRational Q(long p, long q = 1) => new(p, q);

    /// <summary>Every connected row of the letter grid: the 27 words, couplings equal, opposite, reciprocal (e₂ = 0,
    /// alone and on the swap branch), one zero and generic, fields equal, opposite, zero and generic.</summary>
    public static IEnumerable<(string Word, BigRational[] J, BigRational[] H)> Grid()
    {
        var couplings = new[]
        {
            new[] { Q(1), Q(2), Q(3) }, new[] { Q(2), Q(2), Q(3) }, new[] { Q(2), Q(-2), Q(3) },
            new[] { Q(3), Q(6), Q(-2) }, new[] { Q(1), Q(-2), Q(-2) }, new[] { Q(2), Q(2), Q(-1) },
            new[] { Q(0), Q(2), Q(3) }, new[] { Q(2), Q(5), Q(-10, 7) },
        };
        var fields = new[]
        {
            new[] { Q(0), Q(1), Q(1) }, new[] { Q(0), Q(1), Q(-1) }, new[] { Q(1), Q(1), Q(1) },
            new[] { Q(0), Q(2), Q(3) }, new[] { Q(2), Q(0), Q(0) }, new[] { Q(1), Q(2), Q(2) },
        };
        foreach (var w in from a in "XYZ" from b in "XYZ" from c in "XYZ" select $"{a}{b}{c}")
            foreach (var j in couplings)
                foreach (var h in fields)
                    yield return (w, j, h);
    }

    /// <summary>The planar grid's couplings: generic, swap, two reciprocal ones, four on R1 (w = 1, its swap, w = ½,
    /// the w = 1 row scaled by two), the corner J_uv = 0 at two scales, d at an end, one row just off R1, and one off R1
    /// where R2 and R3 hold for the last field pair (R1 = 11), so E's first relation is tested on its own.</summary>
    public static readonly BigRational[][] PlanarCouplings =
    {
        new[] { Q(1), Q(2), Q(3) }, new[] { Q(2), Q(2), Q(3) }, new[] { Q(3), Q(6), Q(-2) }, new[] { Q(2), Q(5), Q(-10, 7) },
        new[] { Q(1), Q(2), Q(6) }, new[] { Q(2), Q(1), Q(6) }, new[] { Q(1), Q(3), Q(3) }, new[] { Q(2), Q(4), Q(12) },
        new[] { Q(2), Q(-2), Q(0) }, new[] { Q(4), Q(-4), Q(0) }, new[] { Q(0), Q(2), Q(3) }, new[] { Q(1), Q(2), Q(7) },
        new[] { Q(-4), Q(-1), Q(-1) },
    };

    /// <summary>The planar grid's field pairs (f_u, f_v): parallel, one zero, equal lengths at two angles, generic;
    /// the E row of J = (1, 2, 6) with its mirror, its turn by (3/5, 4/5) about z, its swap (for J = (2, 1, 6)), its
    /// double (for J = (2, 4, 12)) and a perturbation; the corner at κ = 2 in the symmetric gauge and turned, at κ = 4
    /// scaled to J = (4, −4, 0), and a perturbation; and the pair with R2 = R3 = 0 on J = (−4, −1, −1).</summary>
    public static readonly (BigRational[] Fu, BigRational[] Fv)[] PlanarFieldPairs =
    {
        (new[] { Q(1), Q(2) }, new[] { Q(2), Q(4) }), (new[] { Q(0), Q(0) }, new[] { Q(3), Q(1) }),
        (new[] { Q(3), Q(4) }, new[] { Q(5), Q(0) }), (new[] { Q(3), Q(4) }, new[] { Q(4), Q(3) }),
        (new[] { Q(1), Q(2) }, new[] { Q(3), Q(1) }),
        (new[] { Q(8), Q(0) }, new[] { Q(-5), Q(3) }), (new[] { Q(8), Q(0) }, new[] { Q(-5), Q(-3) }),
        (new[] { Q(24, 5), Q(32, 5) }, new[] { Q(-27, 5), Q(-11, 5) }), (new[] { Q(-5), Q(3) }, new[] { Q(8), Q(0) }),
        (new[] { Q(16), Q(0) }, new[] { Q(-10), Q(6) }), (new[] { Q(8), Q(0) }, new[] { Q(-5), Q(4) }),
        (new[] { Q(3), Q(1) }, new[] { Q(3), Q(-1) }), (new[] { Q(1), Q(4) }, new[] { Q(4), Q(1) }),
        (new[] { Q(9), Q(7) }, new[] { Q(9), Q(-7) }), (new[] { Q(3), Q(1) }, new[] { Q(3), Q(-2) }),
        (new[] { Q(6), Q(0) }, new[] { Q(-3), Q(3) }),
    };

    public static IEnumerable<(BigRational[] J, BigRational[] Fu, BigRational[] Fv)> PlanarGrid() =>
        from j in PlanarCouplings from p in PlanarFieldPairs select (j, p.Fu, p.Fv);

    public sealed record GridReading(int Rows, int Palindromic, int A, int B, int C, int E, int Disagreements);

    public static GridReading ReadGrid(char? without = null) =>
        Tally(Grid().Select(g => ReadRow(g.Word, g.J, g.H, without)));

    public static GridReading ReadPlanarGrid(char? without = null) =>
        Tally(PlanarGrid().Select(g => ReadPlanarRow(g.J, g.Fu, g.Fv, without)));

    private static GridReading Tally(IEnumerable<Reading> readings)
    {
        int rows = 0, pal = 0, a = 0, b = 0, c = 0, e = 0, bad = 0;
        foreach (var r in readings)
        {
            if (!r.RuleStated) throw new ArgumentException("a grid reads a rule; a row with a field on d states none.");
            rows++;
            if (r.Palindrome) pal++;
            switch (r.InClass) { case 'A': a++; break; case 'B': b++; break; case 'C': c++; break; case 'E': e++; break; }
            if (!r.Agrees) bad++;
        }
        return new GridReading(rows, pal, a, b, c, e, bad);
    }

    // ---- IInspectable -----------------------------------------------------------------------------------

    public string DisplayName => _fu is null
        ? $"the triangle's palindromic locus at N = 3 (word {_word}, J = {string.Join(", ", _j)}, h = {string.Join(", ", _h)})"
        : $"the triangle's palindromic locus at N = 3 with planar fields (J = {string.Join(", ", _j)}, " +
          $"f_u = ({string.Join(", ", _fu)}), f_v = ({string.Join(", ", _fv!)}), " +
          (_fd!.All(x => x.IsZero) ? "no field on d)" : $"f_d = ({string.Join(", ", _fd!)}))");

    public string Summary
    {
        get
        {
            var r = Read();
            var loc = ReadFarLocality();
            string locality = loc is null ? "" : string.Format(CultureInfo.InvariantCulture,
                "; far element of Schmidt ranks {0}, {1}, {2} across d|uv, u|dv, v|du{3}", loc.DUv, loc.UDv, loc.VDu,
                loc.Product ? ", a product" : ", no product, so no site-local mirror (conditioned ones aside)");
            if (!r.RuleStated)
                return string.Format(CultureInfo.InvariantCulture, "near {0}, far {1}: {2}; no rule is stated with a field on d{3}",
                    r.Near, r.Far, r.Palindrome ? "palindrome" : "broken", locality);
            return string.Format(CultureInfo.InvariantCulture, "near {0}, far {1}: {2}; {3}; the rule {4}{5}",
                r.Near, r.Far, r.Palindrome ? "palindrome" : "broken",
                r.InClass is { } k ? $"in class {k}" : _fu is null ? "in none of A, B, C" : "in none of A, B, C, E",
                r.Agrees ? "holds" : "FAILS", locality);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var loc = ReadFarLocality();
            yield return new InspectableNode("the far element's locality",
                summary: loc is null
                    ? "The row does not pair, so no mirror carries it, or its far end is not one-dimensional, so the " +
                      "identity's image is not fixed; the locality reading is not made."
                    : string.Format(CultureInfo.InvariantCulture,
                        "Operator Schmidt ranks of the far element: {0} across d|uv, {1} across u|dv, {2} across v|du. {3} " +
                        "A mirror, onto −L − 2σ or onto −L† − 2σ, sends the identity onto the far end, so across every " +
                        "cut its own Schmidt rank is at least the element's; a site-factorised one sends it onto a product " +
                        "(F116 §4; in rank form finding 1 of experiments/THE_CUBE_SITE_BY_SITE.md). Mirrors conditioned on " +
                        "the bits of a site are not covered.",
                        loc.DUv, loc.UDv, loc.VDu,
                        loc.Product ? "It is a product, so a site-local mirror is not ruled out."
                                    : "It is no product, so no site-local mirror carries this row, in any per-site frame."),
                provenance: NodeProvenance.Live);
            var g = ReadGrid();
            yield return new InspectableNode("the letter grid of special loci, recomputed",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "{0} rows over the 27 words: {1} palindromic ({2} in A, {3} in B, {4} in C, first class that " +
                    "holds), {5} rows where the rule and the two ends disagree", g.Rows, g.Palindromic, g.A, g.B, g.C,
                    g.Disagreements), provenance: NodeProvenance.Live);
            var pg = ReadPlanarGrid();
            yield return new InspectableNode("the planar grid, recomputed",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "{0} rows with no field on d and fields perpendicular to z: {1} palindromic ({2} in A, {3} in B, " +
                    "{4} in C, {5} in E, first class that holds), {6} rows where the rule and the two ends disagree",
                    pg.Rows, pg.Palindromic, pg.A, pg.B, pg.C, pg.E, pg.Disagreements), provenance: NodeProvenance.Live);
            yield return new InspectableNode("where the exhaustion lives",
                summary: "This witness reads the rules at points. That every palindromic connected row with letter " +
                         "fields lies in A, B or C is the elimination's, certificate " +
                         "simulations/results/n3_palindrome_covering.json (1,305 far leaves, 0 open) and its independent " +
                         "gate simulations/n3_palindrome_covering_gate.py; that with planar fields and none on d it pairs " +
                         "exactly in A, B, C or E is the derivation of the section 'The triangle with planar fields', " +
                         "gate simulations/n3_planar_field_locus_gate.py.");
        }
    }
}
