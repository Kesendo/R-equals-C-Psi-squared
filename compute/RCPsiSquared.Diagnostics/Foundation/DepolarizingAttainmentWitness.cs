using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for F5 (registry <c>docs/ANALYTICAL_FORMULAS.md</c> F5, claim
/// <see cref="RCPsiSquared.Core.Symmetry.F5DepolarizingErrorPi2Inheritance"/>, gate
/// <c>simulations/f5_depolarizing_attainment.py</c>): the depolarizing error is a lower bound for
/// every H and exact when ad_H has an eigenvector among the operators traceless on every site.
///
/// <code>
///     error = 2 sigma - r_max  >=  (2/3) sigma,     equality  iff  ad_H has an eigenvector in E
/// </code>
///
/// with sigma the total rate, r_max the fastest decay rate present, and E the span of the Pauli
/// strings with no identity letter. Off the isotropic class, a Pauli channel with its own rate per
/// letter and site, the same reading holds with the bound 2 sum_l min_P gamma_P^l and E the span of
/// the strings whose letter at each site is one of the site's fastest letters.
///
/// <para><b>Three computations, two of them exact.</b> Everything is recomputed and nothing is
/// looked up. (a) THE CERTIFICATE, exactly: the Hamiltonian is built with Gaussian-integer entries,
/// every candidate string of the fast span is a signed permutation, and [H, P] = 0 is decided with
/// == 0 on the exact products; for a commuting string the dissipator is applied exactly in units
/// of 1/30 of a rate (isotropic rates in tenths split three ways, Pauli-channel rates in tenths
/// times three), and 30 L(P) = -(2 sigma - bound) P is checked entry by entry. (b) THE FASTEST
/// RATE, exactly: every string's decay rate is read off the dissipator's diagonal as the sum of
/// 2 gamma over the (site, letter) pairs that anticommute with it, and 2 sigma minus the maximum
/// over all 4^N strings is the bound, in integers. (c) THE SPECTRUM, in floats: the Liouvillian's
/// eigenvalues, r_max, the shortfall against the bound, with the deviation read against the
/// eigensolver's rounding model eps * ||L||_F rather than against a threshold (on an attained row
/// the edge eigenvalue lies on the boundary of L's numerical range, where L acts normally, so it
/// is semisimple with condition number 1 and eps ||L|| is its scale; on a strict row the deviation
/// is the physics, 1e12 times that scale and more); and the general criterion, the largest
/// dimension of an eigenspace of ad_H intersected with E, by singular values, which is what
/// decides the rows where no single string commutes (the ladder). The criterion carries readings
/// beside its two cut-offs, so that neither cut-off stands alone. The differences of H's
/// eigenvalues are grouped to 1e3 eps ||H||_F, membership decided against the group's first
/// member, and the smallest gap between the groups' means over the largest spread inside a group
/// is printed in decades: an over-split (one difference cut in two) closes that gap, and a small
/// spread with a wide gap shows that no two distinct differences sit near the tolerance, which is
/// the merge the reading cannot see directly, since a merged group has a tiny spread by
/// construction. The singular values are counted at 1 within 1e-6; on attained rows the nearest
/// uncounted value's distance from 1 is printed in decades above the worst counted one (floored
/// at eps) and the worst counted one in units of eps, on rows where nothing is counted the largest
/// singular value.</para>
///
/// <para><b>Falsifiers sit beside the readings.</b> The default row is the Heisenberg chain, where
/// the three uniform strings commute and the bound is attained; the children carry the generic
/// two-local chain (bonds and fields) at the same N and rates, where no string commutes, the
/// criterion reads 0 and the shortfall exceeds the bound, so a witness that had lost the ability
/// to report "not attained" would be visible.</para>
///
/// <para>Args: <c>--N</c> (2..4, default 3), <c>--model</c> (heisenberg | xy | xx | ising | xxz2 |
/// dm | heisdm | isingx | ladder | generic, default heisenberg), <c>--rates</c> (per-site
/// isotropic totals in tenths, e.g. 3,7,11, a longer list cut to N; default 3,7,11,5), <c>--pauli</c> (per-site
/// X:Y:Z rates in tenths, e.g. 1:5:7,2:9:4,1:3:8, a longer list cut to N, replacing --rates with
/// an anisotropic Pauli channel; zero rates allowed, a site with every rate zero is dark, and at
/// least one rate must be positive: the registry's dark-site clause, "no rate needs to be
/// positive", is read here per site, and a profile with no jump at all is outside the
/// question).</para></summary>
public sealed class DepolarizingAttainmentWitness : IInspectable
{
    /// <summary>A cost guard only: the spectrum is a dense 4^N eigenproblem.</summary>
    public const int MaxN = 4;

    public static readonly string[] Models =
        { "heisenberg", "xy", "xx", "ising", "xxz2", "dm", "heisdm", "isingx", "ladder", "generic" };

    private const double Eps = 2.220446049250313e-16;

    private readonly int _n;
    private readonly string _model;
    /// <summary>Per site and letter (X, Y, Z), the rate in THIRTIETHS; isotropic rows carry the
    /// site's tenths at every letter (a third of the site's total, three times).</summary>
    private readonly int[][] _c;
    private readonly bool _isotropic;
    private readonly string _ratesText;
    /// <summary>H with Gaussian-integer entries; the physical H is this divided by <see cref="_scale"/>.</summary>
    private readonly Complex[,] _h;
    private readonly int _scale;

    public DepolarizingAttainmentWitness(int n, string? model = null, string? rates = null, string? pauliRates = null)
    {
        if (n < 2 || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), n,
                $"--root depol reads a dense 4^N spectrum and is guarded at N in 2..{MaxN}; got {n}. " +
                "The theorem carries no N, only this witness does.");
        _n = n;
        _model = (model ?? "heisenberg").ToLowerInvariant();
        if (Array.IndexOf(Models, _model) < 0)
            throw new ArgumentException($"unknown model '{model}'; one of {string.Join(", ", Models)}", nameof(model));
        if (pauliRates is not null)
        {
            _isotropic = false;
            _c = ParsePauli(pauliRates, n);
            _ratesText = "Pauli channel, X:Y:Z per site in tenths: " + pauliRates;
        }
        else
        {
            _isotropic = true;
            var tenths = ParseIsotropic(rates, n);
            _c = tenths.Select(t => new[] { t, t, t }).ToArray();
            _ratesText = "isotropic, per-site totals in tenths: " + string.Join(",", tenths);
        }
        (_h, _scale) = BuildHamiltonian(_model, n);
    }

    // ------------------------------------------------------------------ parsing

    private static int[] ParseIsotropic(string? spec, int n)
    {
        var text = spec ?? "3,7,11,5";
        var parts = text.Split(',', StringSplitOptions.TrimEntries | StringSplitOptions.RemoveEmptyEntries);
        if (parts.Length < n)
            throw new ArgumentException($"--rates needs at least one entry per site ({n}); got '{text}'", nameof(spec));
        parts = parts.Take(n).ToArray();                          // a longer list is cut to N sites
        var r = parts.Select(p => ParseInt(p, nameof(spec))).ToArray();
        if (r.Any(t => t <= 0))
            throw new ArgumentException("F5's class needs a positive rate at every site; use --pauli for zero rates", nameof(spec));
        return r;
    }

    private static int ParseInt(string text, string argName)
    {
        if (!int.TryParse(text, NumberStyles.Integer, CultureInfo.InvariantCulture, out int v))
            throw new ArgumentException($"'{text}' is not an integer", argName);
        return v;
    }

    private static int[][] ParsePauli(string spec, int n)
    {
        var sites = spec.Split(',', StringSplitOptions.TrimEntries | StringSplitOptions.RemoveEmptyEntries);
        if (sites.Length < n)
            throw new ArgumentException($"--pauli needs at least one X:Y:Z triple per site ({n}); got '{spec}'", nameof(spec));
        var c = sites.Take(n).Select(s =>
        {
            var t = s.Split(':').Select(v => ParseInt(v, nameof(spec))).ToArray();
            if (t.Length != 3 || t.Any(v => v < 0))
                throw new ArgumentException($"a site's rates are three non-negative tenths X:Y:Z; got '{s}'", nameof(spec));
            return new[] { 3 * t[0], 3 * t[1], 3 * t[2] };          // tenths -> thirtieths
        }).ToArray();
        if (c.All(site => site.All(v => v == 0)))
            throw new ArgumentException("every rate is zero: no jump at all, outside F5's question", nameof(spec));
        return c;
    }

    // ------------------------------------------------------------------ exact operators

    private static readonly Complex[][,] Pauli =
    {
        new Complex[,] { { 1, 0 }, { 0, 1 } },
        new Complex[,] { { 0, 1 }, { 1, 0 } },
        new Complex[,] { { 0, new Complex(0, -1) }, { new Complex(0, 1), 0 } },
        new Complex[,] { { 1, 0 }, { 0, -1 } },
    };

    private static Complex[,] Kron(Complex[,] a, Complex[,] b)
    {
        int ra = a.GetLength(0), ca = a.GetLength(1), rb = b.GetLength(0), cb = b.GetLength(1);
        var o = new Complex[ra * rb, ca * cb];
        for (int i = 0; i < ra; i++)
            for (int j = 0; j < ca; j++)
            {
                if (a[i, j] == Complex.Zero) continue;
                for (int k = 0; k < rb; k++)
                    for (int l = 0; l < cb; l++)
                        o[i * rb + k, j * cb + l] = a[i, j] * b[k, l];
            }
        return o;
    }

    /// <summary>The Pauli string with letter[l] at site l (0 = I), site 0 the most significant.</summary>
    private static Complex[,] StringOp(int[] letters)
    {
        var o = new Complex[,] { { 1 } };
        foreach (int t in letters) o = Kron(o, Pauli[t]);
        return o;
    }

    private static Complex[,] SiteOp(int letter, int site, int n) =>
        StringOp(Enumerable.Range(0, n).Select(k => k == site ? letter : 0).ToArray());

    private static Complex[,] BondOp(int a, int b, int site, int n) =>
        StringOp(Enumerable.Range(0, n).Select(k => k == site ? a : k == site + 1 ? b : 0).ToArray());

    private static Complex[,] Mul(Complex[,] a, Complex[,] b)
    {
        int d = a.GetLength(0);
        var o = new Complex[d, d];
        for (int i = 0; i < d; i++)
            for (int k = 0; k < d; k++)
            {
                if (a[i, k] == Complex.Zero) continue;
                for (int j = 0; j < d; j++)
                    if (b[k, j] != Complex.Zero) o[i, j] += a[i, k] * b[k, j];
            }
        return o;
    }

    private static void AddScaled(Complex[,] acc, Complex[,] m, int c)
    {
        int d = acc.GetLength(0);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                acc[i, j] += c * m[i, j];
    }

    private static bool IsZero(Complex[,] m)
    {
        foreach (var z in m) if (z != Complex.Zero) return false;
        return true;
    }

    private static Complex[,] Commutator(Complex[,] a, Complex[,] b)
    {
        var o = Mul(a, b);
        var ba = Mul(b, a);
        int d = o.GetLength(0);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                o[i, j] -= ba[i, j];
        return o;
    }

    /// <summary>A small deterministic integer stream for the ladder and the generic chain.</summary>
    private static IEnumerable<int> Lcg(int seed)
    {
        uint s = (uint)seed;
        while (true)
        {
            s = unchecked(s * 1664525u + 1013904223u);
            yield return (int)((s >> 16) % 5) - 2;                   // -2 .. 2
        }
    }

    private static (Complex[,] H, int Scale) BuildHamiltonian(string model, int n)
    {
        int d = 1 << n;
        var h = new Complex[d, d];
        void Bond(int a, int b, int c) { for (int s = 0; s < n - 1; s++) AddScaled(h, BondOp(a, b, s, n), c); }
        void Field(int a, int c) { for (int s = 0; s < n; s++) AddScaled(h, SiteOp(a, s, n), c); }
        switch (model)
        {
            case "heisenberg": Bond(1, 1, 1); Bond(2, 2, 1); Bond(3, 3, 1); return (h, 1);
            case "xy": Bond(1, 1, 1); Bond(2, 2, 1); return (h, 1);
            case "xx": Bond(1, 1, 1); return (h, 1);
            case "ising": Bond(3, 3, 1); return (h, 1);
            case "xxz2": Bond(1, 1, 1); Bond(2, 2, 1); Bond(3, 3, 2); return (h, 1);
            case "dm": Bond(1, 2, 1); Bond(2, 1, -1); return (h, 1);
            case "heisdm": Bond(1, 1, 10); Bond(2, 2, 10); Bond(3, 3, 10); Bond(1, 2, 3); Bond(2, 1, -3); return (h, 10);
            case "isingx": Bond(3, 3, 10); Field(1, 7); return (h, 10);
            case "generic":
            {
                // every two-letter word on every bond and every letter on every site with a
                // coefficient in -2..2, the same stream at every N so the rows are reproducible.
                // The fields matter: a bond-only H at N = 2 is itself a sum of strings with no
                // identity letter, so it lies in E and is its own ad_H eigenvector there.
                var g = Lcg(20260928).GetEnumerator();
                for (int s = 0; s < n - 1; s++)
                    for (int a = 1; a <= 3; a++)
                        for (int b = 1; b <= 3; b++)
                        {
                            g.MoveNext();
                            if (g.Current != 0) AddScaled(h, BondOp(a, b, s, n), g.Current);
                        }
                for (int s = 0; s < n; s++)
                    for (int a = 1; a <= 3; a++)
                    {
                        g.MoveNext();
                        if (g.Current != 0) AddScaled(h, SiteOp(a, s, n), g.Current);
                    }
                return (h, 1);
            }
            case "ladder":
            {
                // |0..0> and |1..1> eigenvectors at the distinct energies 7 and -5, a Hermitian
                // Gaussian-integer block with no structure on the rest
                var g = Lcg(20260929).GetEnumerator();
                h[0, 0] = 7;
                h[d - 1, d - 1] = -5;
                for (int i = 1; i < d - 1; i++)
                    for (int j = i; j < d - 1; j++)
                    {
                        g.MoveNext(); int re = g.Current;
                        if (i == j) { h[i, i] = re; continue; }
                        g.MoveNext(); int im = g.Current;
                        h[i, j] = new Complex(re, im);
                        h[j, i] = new Complex(re, -im);
                    }
                return (h, 1);
            }
            default: throw new ArgumentException(model);
        }
    }

    // ------------------------------------------------------------------ the exact readings

    /// <summary>The site's fastest letters (0 = I included at a dark site) and its fastest rate in
    /// thirtieths: a letter P decays at 2 (sum of the site's rates) - 2 gamma_P, the identity at 0.</summary>
    private (int[] Letters, int Rate) FastestLetters(int site)
    {
        int tot = _c[site].Sum();
        var rate = new[] { 0, 2 * (tot - _c[site][0]), 2 * (tot - _c[site][1]), 2 * (tot - _c[site][2]) };
        int top = rate.Max();
        return (Enumerable.Range(0, 4).Where(t => rate[t] == top).ToArray(), top);
    }

    /// <summary>The decay rate of a string under the dissipator, in thirtieths, read off the
    /// dissipator itself: 2 gamma_Q^l over the (l, Q) at which Q_l P Q_l = -P (decided exactly).</summary>
    private int RateOfString(int[] letters)
    {
        var p = StringOp(letters);
        int rate = 0;
        for (int l = 0; l < _n; l++)
            for (int q = 1; q <= 3; q++)
            {
                var ql = SiteOp(q, l, _n);
                var conj = Mul(Mul(ql, p), ql);
                bool plus = true, minus = true;
                for (int i = 0; i < p.GetLength(0) && (plus || minus); i++)
                    for (int j = 0; j < p.GetLength(1); j++)
                    {
                        if (conj[i, j] != p[i, j]) plus = false;
                        if (conj[i, j] != -p[i, j]) minus = false;
                    }
                if (!plus && !minus) throw new InvalidOperationException("a Pauli string conjugated by a Pauli is +- itself");
                if (minus) rate += 2 * _c[l][q - 1];
            }
        return rate;
    }

    public int SigmaThirtieths => _c.Sum(site => site.Sum());

    public int BoundThirtieths => 2 * _c.Sum(site => site.Min());

    private IEnumerable<int[]> FastSpan()
    {
        var choices = Enumerable.Range(0, _n).Select(l => FastestLetters(l).Letters).ToArray();
        var cur = new int[_n];
        IEnumerable<int[]> Walk(int site)
        {
            if (site == _n) { yield return (int[])cur.Clone(); yield break; }
            foreach (int t in choices[site])
            {
                cur[site] = t;
                foreach (var s in Walk(site + 1)) yield return s;
            }
        }
        return Walk(0);
    }

    // ------------------------------------------------------------------ the float readings

    private ComplexMatrix Liouvillian()
    {
        int d = 1 << _n;
        var hTrue = ComplexMatrix.Build.Dense(d, d, (i, j) => _h[i, j] / _scale);
        var id = ComplexMatrix.Build.DenseIdentity(d);
        var id2 = ComplexMatrix.Build.DenseIdentity(d * d);
        var l = (hTrue.KroneckerProduct(id) - id.KroneckerProduct(hTrue.Transpose())) * new Complex(0, -1);
        for (int s = 0; s < _n; s++)
            for (int q = 1; q <= 3; q++)
            {
                double g = _c[s][q - 1] / 30.0;
                if (g == 0) continue;
                var op = SiteOp(q, s, _n);
                var a = ComplexMatrix.Build.Dense(d, d, (i, j) => op[i, j]);
                l += (a.KroneckerProduct(a.Conjugate()) - id2) * g;
            }
        return l;
    }

    /// <summary>The largest dimension of (an eigenspace of ad_H) intersected with the fast span,
    /// by singular values; beside it the separation between the counted and the uncounted
    /// values, the largest singular value seen, and the grouping's own quality (the smallest gap
    /// between distinct eigenvalue differences against the largest spread inside a group).</summary>
    private (int Dim, double SeparationDecades, double WorstCountedOverEps, double MaxSingular, double GroupingDecades) Criterion()
    {
        int d = 1 << _n;
        var hTrue = ComplexMatrix.Build.Dense(d, d, (i, j) => _h[i, j] / _scale);
        var evd = hTrue.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian);
        var w = evd.EigenValues.Select(z => z.Real).ToArray();
        var v = evd.EigenVectors;
        double hNorm = hTrue.FrobeniusNorm();
        var span = FastSpan().ToArray();
        var e = ComplexMatrix.Build.Dense(d * d, span.Length);
        for (int k = 0; k < span.Length; k++)
        {
            var p = StringOp(span[k]);
            for (int i = 0; i < d; i++)
                for (int j = 0; j < d; j++)
                    e[i * d + j, k] = p[i, j] / Math.Sqrt(d);
        }
        // group the differences w_a - w_b to the eigensolver's resolution
        var groups = new List<(double Mu, List<(int, int)> Pairs)>();
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
            {
                double mu = w[a] - w[b];
                int gi = groups.FindIndex(g => Math.Abs(g.Mu - mu) <= 1e3 * Eps * Math.Max(1.0, hNorm));
                if (gi < 0) groups.Add((mu, new List<(int, int)> { (a, b) }));
                else groups[gi].Pairs.Add((a, b));
            }
        double spread = 0.0, gap = double.PositiveInfinity;
        for (int gi = 0; gi < groups.Count; gi++)
        {
            double mean = groups[gi].Pairs.Average(ab => w[ab.Item1] - w[ab.Item2]);
            spread = Math.Max(spread, groups[gi].Pairs.Max(ab => Math.Abs(w[ab.Item1] - w[ab.Item2] - mean)));
            for (int gj = 0; gj < gi; gj++)
                gap = Math.Min(gap, Math.Abs(mean - groups[gj].Pairs.Average(ab => w[ab.Item1] - w[ab.Item2])));
        }
        double groupingDecades = Math.Log10(gap / Math.Max(spread, Eps * Math.Max(1.0, hNorm)));
        int best = 0;
        double worstCounted = 0.0, nearestUncounted = double.PositiveInfinity, maxSingular = 0.0;
        foreach (var (_, pairs) in groups)
        {
            var k = ComplexMatrix.Build.Dense(d * d, pairs.Count);
            for (int c = 0; c < pairs.Count; c++)
            {
                var (a, b) = pairs[c];
                for (int i = 0; i < d; i++)
                    for (int j = 0; j < d; j++)
                        k[i * d + j, c] = v[i, a] * Complex.Conjugate(v[j, b]);
            }
            var sv = (k.ConjugateTranspose() * e).Svd(false).S;
            int count = 0;
            foreach (var sz in sv)
            {
                double s = sz.Real;
                maxSingular = Math.Max(maxSingular, s);
                if (Math.Abs(1 - s) < 1e-6) { count++; worstCounted = Math.Max(worstCounted, Math.Abs(1 - s)); }
                else nearestUncounted = Math.Min(nearestUncounted, Math.Abs(1 - s));
            }
            best = Math.Max(best, count);
        }
        // with nothing counted (not attained) or nothing uncounted there is no separation to read
        double sep = best == 0 || double.IsPositiveInfinity(nearestUncounted)
            ? double.NaN
            : Math.Log10(nearestUncounted / Math.Max(worstCounted, Eps));
        return (best, sep, worstCounted / Eps, maxSingular, groupingDecades);
    }

    // ------------------------------------------------------------------ the reading

    /// <param name="Certificates">Every string of the fast span commuting with H, exactly.</param>
    /// <param name="Certificate">The first of them, or null.</param>
    /// <param name="CertificateEigenvalueExact">30 L(P) = -(2 sigma - bound) P entry by entry, for every certificate.</param>
    /// <param name="FastestRateThirtieths">The largest decay rate among all 4^N strings, off the dissipator.</param>
    /// <param name="CriterionDimension">The largest dim of an ad_H eigenspace meeting the fast span (floats).</param>
    /// <param name="CriterionSeparationDecades">Nearest uncounted |1 - s| over the worst counted (floored at eps), in decades; NaN with nothing counted.</param>
    /// <param name="CriterionWorstCountedOverEps">The worst counted |1 - s| in units of eps, the rounding read on the count.</param>
    /// <param name="CriterionMaxSingular">The largest singular value seen, the reading on rows where nothing is counted.</param>
    /// <param name="GroupingDecades">The smallest gap between distinct eigenvalue differences over the largest spread inside a group.</param>
    /// <param name="ShortfallOverBound">(2 sigma - r_max) / bound from the spectrum; NaN when the bound is 0.</param>
    /// <param name="ShortfallDeviationOverModel">|shortfall - bound| / (eps ||L||_F), the rounding read.</param>
    public sealed record Reading(
        IReadOnlyList<string> Certificates, string? Certificate, bool CertificateEigenvalueExact,
        int SigmaThirtieths, int BoundThirtieths, int FastestRateThirtieths,
        int CriterionDimension, double CriterionSeparationDecades, double CriterionWorstCountedOverEps,
        double CriterionMaxSingular, double GroupingDecades,
        double RMax, double Shortfall, double ShortfallOverBound, double ShortfallDeviationOverModel);

    private Reading? _reading;

    /// <summary>Computed once per witness; every node reads the same numbers.</summary>
    public Reading Read() => _reading ??= Compute();

    private Reading Compute()
    {
        // (a) the certificates, exactly
        var certs = new List<string>();
        bool exact = true;
        foreach (var letters in FastSpan())
        {
            if (letters.All(t => t == 0)) continue;              // only an all-dark profile reaches this, and those are rejected
            var p = StringOp(letters);
            if (!IsZero(Commutator(_h, p))) continue;
            certs.Add(string.Concat(letters.Select(t => "IXYZ"[t])));
            // 30 L(P): the Hamiltonian part is exactly zero, the dissipator part is
            // sum_l sum_Q c_Q^l (Q P Q - P), compared with -(2 sigma - bound) P
            int d = 1 << _n;
            var acc = new Complex[d, d];
            for (int l = 0; l < _n; l++)
                for (int q = 1; q <= 3; q++)
                {
                    if (_c[l][q - 1] == 0) continue;
                    var ql = SiteOp(q, l, _n);
                    var conj = Mul(Mul(ql, p), ql);
                    for (int i = 0; i < d; i++)
                        for (int j = 0; j < d; j++)
                            acc[i, j] += _c[l][q - 1] * (conj[i, j] - p[i, j]);
                }
            int rate = 2 * SigmaThirtieths - BoundThirtieths;
            for (int i = 0; i < d && exact; i++)
                for (int j = 0; j < d; j++)
                    if (acc[i, j] != -rate * p[i, j]) { exact = false; break; }
        }
        string? cert = certs.Count > 0 ? certs[0] : null;
        if (certs.Count == 0) exact = false;
        // (b) the fastest rate, exactly, over every string
        int fastest = 0;
        var all = new int[_n];
        for (int code = 0; code < (1 << (2 * _n)); code++)
        {
            for (int l = 0; l < _n; l++) all[l] = (code >> (2 * (_n - 1 - l))) & 3;
            fastest = Math.Max(fastest, RateOfString(all));
        }
        // (c) the spectrum and the criterion, in floats
        var lv = Liouvillian();
        var ev = lv.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Asymmetric).EigenValues;
        double rMax = ev.Max(z => -z.Real);
        double sigma = SigmaThirtieths / 30.0, bound = BoundThirtieths / 30.0;
        double shortfall = 2 * sigma - rMax;
        double model = Eps * lv.FrobeniusNorm();
        var (dim, sep, worstEps, maxS, grouping) = Criterion();
        return new Reading(certs, cert, exact, SigmaThirtieths, BoundThirtieths, fastest, dim, sep, worstEps, maxS, grouping,
            rMax, shortfall, bound > 0 ? shortfall / bound : double.NaN, Math.Abs(shortfall - bound) / model);
    }

    // ------------------------------------------------------------------ the inspectable

    public string Name => "depol";

    public string DisplayName =>
        $"F5 the depolarizing error as a bound and its attainment, live (N = {_n}, {_model}, {_ratesText})";

    public string Summary
    {
        get
        {
            var r = Read();
            return string.Format(CultureInfo.InvariantCulture,
                "bound {0}/30 = {1:F4}, shortfall from the spectrum {2:F6} ({3}); certificate {4}; criterion dim {5}",
                r.BoundThirtieths, r.BoundThirtieths / 30.0, r.Shortfall,
                r.CriterionDimension > 0 ? "ATTAINED" : "NOT ATTAINED, the bound is strict",
                r.Certificate ?? "none (no string of the fast span commutes with H)", r.CriterionDimension);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var r = Read();
            yield return new InspectableNode("the certificate, exactly",
                summary: r.Certificate is null
                    ? string.Format(CultureInfo.InvariantCulture,
                        "No string of the fast span ({0} strings tried, each a signed permutation, [H, P] " +
                        "compared with == 0 on Gaussian-integer entries) commutes with H. The bound may still be " +
                        "attained through an eigenvector of ad_H that is no string; the criterion node decides.",
                        FastSpan().Count())
                    : string.Format(CultureInfo.InvariantCulture,
                        "{0} commute(s) with H exactly, and for each 30 L(P) = -{1} P entry by entry, {2}: so " +
                        "-{1}/30 = -{3:F4} is an eigenvalue of L, the fastest rate the dissipator allows, and " +
                        "the shortfall 2 sigma - r_max equals the bound {4}/30 = {5:F4}.",
                        string.Join(", ", r.Certificates), 2 * r.SigmaThirtieths - r.BoundThirtieths,
                        r.CertificateEigenvalueExact ? "exact" : "NOT EXACT, READ IT",
                        (2 * r.SigmaThirtieths - r.BoundThirtieths) / 30.0, r.BoundThirtieths, r.BoundThirtieths / 30.0));

            yield return new InspectableNode("the fastest rate, off the dissipator",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "Over all 4^N strings the largest decay rate is {0}/30, read exactly as 2 gamma over the " +
                    "(site, letter) pairs that anticommute with the string; 2 sigma - that = {1}/30 {2} the bound " +
                    "2 sum_l min_P gamma_P^l = {3}/30. The isotropic case is (2/3) sigma.",
                    r.FastestRateThirtieths, 2 * r.SigmaThirtieths - r.FastestRateThirtieths,
                    2 * r.SigmaThirtieths - r.FastestRateThirtieths == r.BoundThirtieths ? "==" : "!=, READ IT",
                    r.BoundThirtieths));

            yield return new InspectableNode("the spectrum, in floats",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "r_max = {0:F6}, shortfall 2 sigma - r_max = {1:F6} against the bound {2:F6}: " +
                    "|shortfall - bound| is {3:E1} times eps ||L||_F (the eigensolver's rounding model, read " +
                    "rather than gated){4}.",
                    r.RMax, r.Shortfall, r.BoundThirtieths / 30.0, r.ShortfallDeviationOverModel,
                    double.IsNaN(r.ShortfallOverBound) ? "" :
                        string.Format(CultureInfo.InvariantCulture, "; shortfall / bound = {0:F6}", r.ShortfallOverBound)));

            yield return new InspectableNode("the criterion: an eigenvector of ad_H in the fast span",
                summary: string.Format(CultureInfo.InvariantCulture,
                    "Largest dimension of an eigenspace of ad_H meeting the span: {0} ({1}); {2}; the eigenvalue " +
                    "differences group with {3:F1} decades between the groups' gap and their inner spread. This is " +
                    "the general reading, the one that decides the ladder row at N >= 3, where |0..0><1..1| is the " +
                    "eigenvector and no string commutes.",
                    r.CriterionDimension, r.CriterionDimension > 0 ? "attained" : "not attained",
                    double.IsNaN(r.CriterionSeparationDecades)
                        ? string.Format(CultureInfo.InvariantCulture,
                            "nothing counted at 1, the largest singular value being {0:F4}, 1 - s = {1:E1}, " +
                            "{2:F1} decades above the 1e-6 cut",
                            r.CriterionMaxSingular, 1 - r.CriterionMaxSingular, Math.Log10((1 - r.CriterionMaxSingular) / 1e-6))
                        : string.Format(CultureInfo.InvariantCulture,
                            "singular values counted at 1 against the rest, {0:F1} decades apart, the worst counted " +
                            "{1:F1} eps from 1", r.CriterionSeparationDecades, r.CriterionWorstCountedOverEps),
                    r.GroupingDecades));

            if (_model != "generic")
            {
                var g = new DepolarizingAttainmentWitness(_n, "generic",
                    _isotropic ? string.Join(",", _c.Select(s => s[0])) : null,
                    _isotropic ? null : string.Join(",", _c.Select(s => $"{s[0] / 3}:{s[1] / 3}:{s[2] / 3}"))).Read();
                yield return new InspectableNode("the bound off the locus: a generic chain beside",
                    summary: string.Format(CultureInfo.InvariantCulture,
                        "The same N and rates on a generic two-local chain: certificate {0}, criterion dim {1}, " +
                        "shortfall {2:F6} against the bound {3:F6}{4}. The falsifier the attained row would " +
                        "otherwise lack.",
                        g.Certificate ?? "none", g.CriterionDimension, g.Shortfall, g.BoundThirtieths / 30.0,
                        double.IsNaN(g.ShortfallOverBound) ? "" :
                            string.Format(CultureInfo.InvariantCulture, ", ratio {0:F4}", g.ShortfallOverBound)));
            }

            yield return new InspectableNode("what this witness does NOT decide",
                summary: "The complex pairing distance and the best-pairing error of the depolarizing experiments, " +
                         "which are different numbers (registry F5, 'Which number'); the Jordan structure at the " +
                         "far end; and nothing beyond N = 4, a cost guard. The window-edge lemma behind the " +
                         "criterion is PROOF_CODIM1_BY_ADDITIVITY section 6; the gate that scores the criterion " +
                         "on random chains in three classes and on the Pauli channel is " +
                         "simulations/f5_depolarizing_attainment.py.");
        }
    }
}
