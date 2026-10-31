using System.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>A polynomial over GF(p) in <see cref="PackedMonomial"/> terms, stored largest first in grevlex,
/// coefficients in 0..p−1, none zero.</summary>
public sealed class ModPPolynomial
{
    public PackedMonomial[] Monomials { get; }
    public long[] Coefficients { get; }

    public ModPPolynomial(PackedMonomial[] monomials, long[] coefficients)
    {
        if (monomials.Length != coefficients.Length) throw new ArgumentException("one coefficient per monomial");
        Monomials = monomials; Coefficients = coefficients;
    }

    public static readonly ModPPolynomial Zero = new(Array.Empty<PackedMonomial>(), Array.Empty<long>());

    public bool IsZero => Monomials.Length == 0;
    public bool IsConstant => Monomials.Length == 1 && Monomials[0].Degree == 0;
    public PackedMonomial Lead => Monomials[0];
    public int Count => Monomials.Length;

    /// <summary>The polynomial from (monomial, coefficient) pairs in any order, like terms summed modulo p.</summary>
    public static ModPPolynomial FromTerms(IEnumerable<(PackedMonomial Monomial, long Coefficient)> terms, long p)
    {
        var acc = new Dictionary<PackedMonomial, long>();
        foreach (var (m, c) in terms)
        {
            long r = ((c % p) + p) % p;
            acc[m] = acc.TryGetValue(m, out long old) ? (old + r) % p : r;
        }
        var kept = acc.Where(kv => kv.Value != 0).OrderBy(kv => kv.Key, PackedMonomial.Descending).ToArray();
        return new ModPPolynomial(kept.Select(kv => kv.Key).ToArray(), kept.Select(kv => kv.Value).ToArray());
    }
}

/// <summary>Arithmetic of ℤ[i] and ℚ(i) modulo a prime p ≡ 1 (mod 4), where −1 has a square root r and i ↦ r is
/// a ring map ℤ[i] → GF(p). The root is found by search, a^((p−1)/4) for the first a that gives one, and checked
/// (r² ≡ −1), nothing recalled.</summary>
public sealed class GaussianModP
{
    public long P { get; }

    /// <summary>The image of i: r with r² ≡ −1 (mod p).</summary>
    public long RootOfMinusOne { get; }

    public GaussianModP(long p)
    {
        if (p < 5 || p % 4 != 1) throw new ArgumentOutOfRangeException(nameof(p), $"p must be a prime ≡ 1 (mod 4) above 4; got {p}");
        if (p >= 1L << 31) throw new ArgumentOutOfRangeException(nameof(p), "p must stay below 2^31 so products fit in a long");
        for (long d = 3; d * d <= p; d += 2)
            if (p % d == 0) throw new ArgumentOutOfRangeException(nameof(p), $"{p} is not prime ({d} divides it)");
        P = p;
        for (long a = 2; a < p; a++)
        {
            long r = Pow(a, (p - 1) / 4);
            if (Mul(r, r) == p - 1) { RootOfMinusOne = r; return; }
        }
        throw new ArgumentException($"no square root of −1 modulo {p}; is it prime?");
    }

    public long Mod(long a) { long r = a % P; return r < 0 ? r + P : r; }
    public long Mul(long a, long b) => a * b % P;
    public long Pow(long b, long e) { long r = 1; b = Mod(b); for (; e > 0; e >>= 1) { if ((e & 1) == 1) r = Mul(r, b); b = Mul(b, b); } return r; }
    public long Inverse(long a) => Pow(a, P - 2);

    public long Map(BigInteger n) => (long)(((n % P) + P) % P);

    public long Map(GaussianInteger z) => Mod(Map(z.Re) + Mul(Map(z.Im), RootOfMinusOne));

    /// <summary>The image of a Gaussian rational, null where a denominator vanishes modulo p.</summary>
    public long? Map(GaussianRational z)
    {
        long dr = Map(z.Re.Denominator), di = Map(z.Im.Denominator);
        if (dr == 0 || di == 0) return null;
        long re = Mul(Map(z.Re.Numerator), Inverse(dr)), im = Mul(Map(z.Im.Numerator), Inverse(di));
        return Mod(re + Mul(im, RootOfMinusOne));
    }

    /// <summary>The primes p ≡ 1 (mod 4) just below the bound, by trial division, largest first.</summary>
    public static long[] PrimesOneModFour(int count, long below = (1L << 31) - 1)
    {
        var found = new List<long>();
        for (long c = below; found.Count < count && c > 4; c--)
        {
            if (c % 4 != 1) continue;
            bool prime = true;
            for (long d = 3; d * d <= c; d += 2) if (c % d == 0) { prime = false; break; }
            if (prime) found.Add(c);
        }
        return found.ToArray();
    }
}

/// <summary>Buchberger's algorithm over GF(p) in grevlex, keeping a TRACE: for every polynomial it adds, where it
/// came from (an input, an earlier trace polynomial reduced again, or the S-polynomial of two earlier ones) and the
/// reduction steps that followed, each a reductor and the monomial it was multiplied by. The strategy is modelled on
/// sympy's <c>groebner</c> (Becker–Weispfenning GROEBNERNEWS2) up to the breaking of ties, which sympy leaves to the
/// iteration order of its sets and this run to list order: the inputs inter-reduced until stable (each by the ones
/// before it, in list order), entered by increasing leading monomial through Gebauer–Möller's update, the pair with
/// the smallest lcm first, every S-polynomial reduced by the basis sorted by increasing leading monomial, and at the
/// end the reduced basis. The strategy
/// matters for the exact half: <see cref="GroebnerExactReplay"/> repeats the trace over ℚ(i), and the size of the
/// numbers it meets is the size along this route.
///
/// <para>The run stops at the first nonzero constant. Over GF(p) that constant says only that 1 lies in the ideal
/// modulo p; when its exact replay also reaches a constant, that proves 1 ∈ I outright, whatever happened at p.</para></summary>
public sealed class GroebnerTraceModP
{
    public readonly record struct Step(int Reductor, PackedMonomial Multiplier);

    /// <summary>Where a trace polynomial came from: Input ≥ 0 for that input; otherwise, with Right &lt; 0, the trace
    /// polynomial Left reduced again; otherwise the S-polynomial of the trace polynomials Left and Right. Then the
    /// reduction steps in order.</summary>
    public sealed record Origin(int Input, int Left, int Right, IReadOnlyList<Step> Steps);

    public long P { get; }
    public int VariableCount { get; }

    /// <summary>Every polynomial the run added, monic, in order.</summary>
    public IReadOnlyList<ModPPolynomial> Polynomials => _polys;

    public IReadOnlyList<Origin> Origins => _origins;

    /// <summary>The trace index of the nonzero constant, or −1 if the run closed without one.</summary>
    public int UnitIndex { get; private set; } = -1;

    /// <summary>Whether 1 lies in the ideal modulo p (the run met a nonzero constant).</summary>
    public bool ContainsOneModP => UnitIndex >= 0;

    /// <summary>The reduced Gröbner basis modulo p at the end, trace indices, largest leading monomial first; empty
    /// when the run reached a constant (and when every input was zero).</summary>
    public IReadOnlyList<int> Basis => _reduced;

    public int PairsTreated { get; private set; }

    private readonly List<ModPPolynomial> _polys = new();
    private readonly List<Origin> _origins = new();
    private List<int> _basis = new();
    private List<int> _reduced = new();
    private List<(int I, int J, PackedMonomial Lcm)> _pairs = new();

    private static readonly IComparer<PackedMonomial> Ascending = Comparer<PackedMonomial>.Create(PackedMonomial.CompareGrevlex);

    private GroebnerTraceModP(long p, int nvars) { P = p; VariableCount = nvars; }

    /// <summary>Run Buchberger on the inputs (already reduced modulo p). maxPairs bounds the work; past it the run
    /// throws rather than report a basis it did not finish.</summary>
    public static GroebnerTraceModP Compute(IReadOnlyList<ModPPolynomial> inputs, long p, int variableCount, int maxPairs = 2_000_000)
    {
        if (variableCount < 1 || variableCount > PackedMonomial.MaxVariables)
            throw new ArgumentOutOfRangeException(nameof(variableCount), $"1..{PackedMonomial.MaxVariables} variables; got {variableCount}");
        var run = new GroebnerTraceModP(p, variableCount);
        // the inputs, monic, then inter-reduced until stable: each by the ones before it, in list order
        var cur = new List<int>();
        for (int k = 0; k < inputs.Count; k++)
        {
            if (inputs[k].IsZero) continue;
            int idx = run.Append(inputs[k], new Origin(k, -1, -1, Array.Empty<Step>()));
            if (run.UnitIndex >= 0) return run;
            cur.Add(idx);
        }
        while (true)
        {
            var next = new List<int>();
            bool changed = false;
            for (int i = 0; i < cur.Count; i++)
            {
                var steps = new List<Step>();
                var r = run.Reduce(run._polys[cur[i]], cur.GetRange(0, i), steps);
                if (r.IsZero) { changed = true; continue; }
                if (steps.Count == 0) { next.Add(cur[i]); continue; }
                changed = true;
                int idx = run.Append(r, new Origin(-1, cur[i], -1, steps));
                if (run.UnitIndex >= 0) return run;
                next.Add(idx);
            }
            cur = next;
            if (!changed) break;
        }
        foreach (int idx in cur.OrderBy(k => run._polys[k].Lead, Ascending))
            run.Update(idx);
        while (run._pairs.Count > 0)
        {
            if (run.PairsTreated >= maxPairs) throw new InvalidOperationException($"Buchberger passed {maxPairs} pairs without closing");
            int best = 0;
            for (int q = 1; q < run._pairs.Count; q++)
                if (PackedMonomial.CompareGrevlex(run._pairs[q].Lcm, run._pairs[best].Lcm) < 0) best = q;
            var pair = run._pairs[best];
            run._pairs.RemoveAt(best);
            run.PairsTreated++;
            var s = run.SPolynomial(pair.I, pair.J, pair.Lcm);
            var steps = new List<Step>();
            var h = run.Reduce(s, run.SortedBasis(), steps);
            if (h.IsZero) continue;
            int idx = run.Append(h, new Origin(-1, pair.I, pair.J, steps));
            if (run.UnitIndex >= 0) return run;
            run.Update(idx);
        }
        // the reduced basis: each element's tail reduced by the others
        var reduced = new List<int>();
        var sorted = run.SortedBasis();
        foreach (int g in run._basis)
        {
            var others = sorted.Where(o => o != g).ToList();
            var steps = new List<Step>();
            var r = run.Reduce(run._polys[g], others, steps);
            if (r.IsZero) continue;                           // a lead an earlier-entered element divides
            reduced.Add(steps.Count == 0 ? g : run.Append(r, new Origin(-1, g, -1, steps)));
        }
        run._reduced = reduced.OrderBy(k => run._polys[k].Lead, PackedMonomial.Descending).ToList();
        return run;
    }

    private List<int> SortedBasis() => _basis.OrderBy(k => _polys[k].Lead, Ascending).ToList();

    private long Mul(long a, long b) => a * b % P;

    private long Inverse(long a)
    {
        long r = 1, e = P - 2, b = a;
        for (; e > 0; e >>= 1) { if ((e & 1) == 1) r = Mul(r, b); b = Mul(b, b); }
        return r;
    }

    private ModPPolynomial Monic(ModPPolynomial f)
    {
        long inv = Inverse(f.Coefficients[0]);
        if (inv == 1) return f;
        return new ModPPolynomial(f.Monomials, f.Coefficients.Select(c => Mul(c, inv)).ToArray());
    }

    private int Append(ModPPolynomial h, Origin origin)
    {
        int index = _polys.Count;
        _polys.Add(Monic(h));
        _origins.Add(origin);
        if (_polys[index].IsConstant) UnitIndex = index;
        return index;
    }

    /// <summary>S(f, g) = (L/lt f)·f − (L/lt g)·g for the monic f, g.</summary>
    private ModPPolynomial SPolynomial(int i, int j, PackedMonomial lcm)
    {
        var f = _polys[i]; var g = _polys[j];
        var uf = f.Lead.DivideInto(lcm); var ug = g.Lead.DivideInto(lcm);
        var terms = new List<(PackedMonomial, long)>(f.Count + g.Count);
        for (int k = 1; k < f.Count; k++) terms.Add((uf * f.Monomials[k], f.Coefficients[k]));
        for (int k = 1; k < g.Count; k++) terms.Add((ug * g.Monomials[k], P - g.Coefficients[k]));
        return ModPPolynomial.FromTerms(terms, P);
    }

    /// <summary>Full reduction of f by the divisors (the first one in list order whose leading monomial divides a
    /// term reduces it), every step recorded.</summary>
    private ModPPolynomial Reduce(ModPPolynomial f, IReadOnlyList<int> divisors, List<Step> steps)
    {
        var acc = new Dictionary<PackedMonomial, long>(f.Count * 2);
        var heap = new PriorityQueue<PackedMonomial, PackedMonomial>(PackedMonomial.Descending);
        for (int k = 0; k < f.Count; k++) { acc[f.Monomials[k]] = f.Coefficients[k]; heap.Enqueue(f.Monomials[k], f.Monomials[k]); }
        var outM = new List<PackedMonomial>();
        var outC = new List<long>();
        while (heap.Count > 0)
        {
            var m = heap.Dequeue();
            if (!acc.Remove(m, out long c) || c == 0) continue;
            int reductor = -1;
            foreach (int b in divisors)
                if (_polys[b].Lead.Divides(m)) { reductor = b; break; }
            if (reductor < 0) { outM.Add(m); outC.Add(c); continue; }
            var g = _polys[reductor];
            var u = g.Lead.DivideInto(m);
            steps.Add(new Step(reductor, u));
            for (int k = 1; k < g.Count; k++)
            {
                var t = u * g.Monomials[k];
                long sub = Mul(c, g.Coefficients[k]);
                if (acc.TryGetValue(t, out long old)) acc[t] = (old - sub + P) % P;
                else { acc[t] = (P - sub) % P; heap.Enqueue(t, t); }
            }
        }
        return new ModPPolynomial(outM.ToArray(), outC.ToArray());
    }

    /// <summary>Gebauer–Möller's update (Becker–Weispfenning, UPDATE) for the new basis element h.</summary>
    private void Update(int h)
    {
        var lh = _polys[h].Lead;
        var c = _basis.Select(g => (G: g, Lcm: PackedMonomial.Lcm(lh, _polys[g].Lead))).ToList();
        var d = new List<(int G, PackedMonomial Lcm)>();
        for (int q = 0; q < c.Count; q++)
        {
            var (g1, l1) = c[q];
            bool keep = PackedMonomial.AreCoprime(lh, _polys[g1].Lead);
            if (!keep)
            {
                keep = true;
                for (int r = q + 1; r < c.Count && keep; r++) if (c[r].Lcm.Divides(l1)) keep = false;
                for (int r = 0; r < d.Count && keep; r++) if (d[r].Lcm.Divides(l1)) keep = false;
            }
            if (keep) d.Add((g1, l1));
        }
        var e = d.Where(x => !PackedMonomial.AreCoprime(lh, _polys[x.G].Lead)).ToList();
        var kept = new List<(int I, int J, PackedMonomial Lcm)>(_pairs.Count + e.Count);
        foreach (var pr in _pairs)
        {
            bool drop = lh.Divides(pr.Lcm)
                        && PackedMonomial.Lcm(_polys[pr.I].Lead, lh) != pr.Lcm
                        && PackedMonomial.Lcm(lh, _polys[pr.J].Lead) != pr.Lcm;
            if (!drop) kept.Add(pr);
        }
        foreach (var (g, l) in e) kept.Add((g, h, l));
        _pairs = kept;
        var nb = _basis.Where(g => !lh.Divides(_polys[g].Lead)).ToList();
        nb.Add(h);
        _basis = nb;
    }

    /// <summary>The trace indices the unit depends on, the unit included, ascending.</summary>
    public IReadOnlyList<int> AncestorsOfUnit() => UnitIndex < 0 ? Array.Empty<int>() : Ancestors(new[] { UnitIndex });

    /// <summary>The trace indices the targets depend on, the targets included, ascending.</summary>
    public IReadOnlyList<int> Ancestors(IEnumerable<int> targets)
    {
        var seen = new HashSet<int>();
        var stack = new Stack<int>(targets);
        while (stack.Count > 0)
        {
            int k = stack.Pop();
            if (!seen.Add(k)) continue;
            var o = _origins[k];
            if (o.Input < 0) { stack.Push(o.Left); if (o.Right >= 0) stack.Push(o.Right); }
            foreach (var s in o.Steps) stack.Push(s.Reductor);
        }
        return seen.OrderBy(k => k).ToArray();
    }
}
