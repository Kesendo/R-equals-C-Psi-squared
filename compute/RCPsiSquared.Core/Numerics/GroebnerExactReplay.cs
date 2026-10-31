using System.Numerics;

namespace RCPsiSquared.Core.Numerics;

/// <summary>The exact half of a modular Gröbner run: the steps a <see cref="GroebnerTraceModP"/> took, repeated over
/// ℚ(i) on the exact inputs, every replayed polynomial made monic. Every step is an ideal operation (an
/// S-polynomial, a reduction f ↦ f − a·u·g, a nonzero scaling), so every polynomial the replay produces lies in the
/// ideal the exact inputs generate, whatever the prime did; when the replayed unit comes out a nonzero constant,
/// 1 ∈ I over ℚ(i) and the variety is empty over ℂ. That is a proof by construction and needs no Gröbner criterion.
///
/// <para>What the prime can spoil is only the route: where p divides a coefficient that is not zero exactly, the
/// exact polynomial has a term the modular one lacks, and later steps aimed at the modular leading terms would
/// miss; and where p divides the denominator of a real or imaginary part, the coefficient has no image here (the
/// map is taken part by part, so this can also set aside a prime at which the coefficient is integral). The replay therefore checks after every polynomial that its support equals the modular one and that it
/// reduces to a nonzero multiple of it modulo p; a mismatch is reported (a bad prime), never papered over, and the
/// caller tries another prime. Only the ancestors of the targets are replayed. The coefficients met are those of
/// the polynomials along the trace's route, so their size is the route's (see the strategy note on the trace).</para></summary>
public static class GroebnerExactReplay
{
    /// <summary>The outcome: whether the replay proved 1 ∈ I over ℚ(i) (the unit came out a nonzero constant, monic
    /// so 1), how many trace polynomials were replayed, the largest numerator or denominator of the replayed monic
    /// polynomials in bits (not of the intermediate sums inside a reduction), and where it stopped if it did.</summary>
    public sealed record Result(bool ProvesUnit, int Replayed, long MaxBits, string? Mismatch);

    public static Result ReplayUnit(GroebnerTraceModP trace, IReadOnlyList<Dictionary<PackedMonomial, GaussianInteger>> exactInputs, GaussianModP field)
    {
        if (trace.UnitIndex < 0) return new Result(false, 0, 0, "the modular run reached no constant");
        if (field.P != trace.P) throw new ArgumentException("the field's prime is not the trace's");
        var (exact, replayed, maxBits, mismatch0) = Replay(trace, new[] { trace.UnitIndex }, exactInputs, field);
        if (mismatch0 is not null) return new Result(false, replayed, maxBits, mismatch0);
        var unit = exact[trace.UnitIndex];
        if (unit.Count != 1 || !unit.ContainsKey(PackedMonomial.One))
            return new Result(false, replayed, maxBits, "the replayed unit is not a constant");
        return new Result(true, replayed, maxBits, null);
    }

    /// <summary>The exact monic polynomials of the trace indices the targets depend on (the targets included),
    /// each checked against its modular image; a mismatch stops the replay and is returned.</summary>
    private static (Dictionary<int, Dictionary<PackedMonomial, GaussianRational>> Exact, int Replayed, long MaxBits, string? Mismatch) Replay(
        GroebnerTraceModP trace, IReadOnlyList<int> targets, IReadOnlyList<Dictionary<PackedMonomial, GaussianInteger>> exactInputs,
        GaussianModP field)
    {
        var exact = new Dictionary<int, Dictionary<PackedMonomial, GaussianRational>>();
        long maxBits = 0;
        int replayed = 0;
        foreach (int k in trace.Ancestors(targets))
        {
            var o = trace.Origins[k];
            Dictionary<PackedMonomial, GaussianRational> f;
            if (o.Input >= 0) f = exactInputs[o.Input].ToDictionary(kv => kv.Key, kv => (GaussianRational)kv.Value);
            else if (o.Right < 0) f = new Dictionary<PackedMonomial, GaussianRational>(exact[o.Left]);
            else f = SPolynomial(exact[o.Left], trace.Polynomials[o.Left].Lead, exact[o.Right], trace.Polynomials[o.Right].Lead);
            foreach (var step in o.Steps)
            {
                var g = exact[step.Reductor];                   // monic, lead at the modular lead
                var target = step.Multiplier * trace.Polynomials[step.Reductor].Lead;
                if (!f.TryGetValue(target, out var a) || a.IsZero) continue;
                foreach (var (m, c) in g)
                {
                    var t = step.Multiplier * m;
                    var v = (f.TryGetValue(t, out var old) ? old : GaussianRational.Zero) - a * c;
                    if (v.IsZero) f.Remove(t); else f[t] = v;
                }
            }
            string? mismatch = Compare(f, trace.Polynomials[k], field);
            if (mismatch is not null) return (exact, replayed + 1, maxBits, $"trace polynomial {k}: {mismatch}");
            var inv = GaussianRational.One / f[trace.Polynomials[k].Lead];
            var monic = new Dictionary<PackedMonomial, GaussianRational>(f.Count);
            foreach (var (m, c) in f)
            {
                var v = c * inv;
                monic[m] = v;
                maxBits = Math.Max(maxBits, Math.Max(Bits(v.Re), Bits(v.Im)));
            }
            replayed++;
            exact[k] = monic;
        }
        return (exact, replayed, maxBits, null);
    }

    private static long Bits(BigRational q) => Math.Max((long)BigInteger.Abs(q.Numerator).GetBitLength(), (long)q.Denominator.GetBitLength());

    private static Dictionary<PackedMonomial, GaussianRational> SPolynomial(
        Dictionary<PackedMonomial, GaussianRational> f, PackedMonomial lf,
        Dictionary<PackedMonomial, GaussianRational> g, PackedMonomial lg)
    {
        // both monic: S = (L/lf)·f − (L/lg)·g
        var lcm = PackedMonomial.Lcm(lf, lg);
        var uf = lf.DivideInto(lcm); var ug = lg.DivideInto(lcm);
        var s = new Dictionary<PackedMonomial, GaussianRational>(f.Count + g.Count);
        foreach (var (m, c) in f) s[uf * m] = c;
        foreach (var (m, c) in g)
        {
            var t = ug * m;
            var v = (s.TryGetValue(t, out var old) ? old : GaussianRational.Zero) - c;
            if (v.IsZero) s.Remove(t); else s[t] = v;
        }
        return s;
    }

    /// <summary>Whether the exact f has the modular polynomial's support and reduces to a nonzero multiple of it.</summary>
    private static string? Compare(Dictionary<PackedMonomial, GaussianRational> f, ModPPolynomial modular, GaussianModP field)
    {
        if (f.Count != modular.Count) return $"{f.Count} exact terms against {modular.Count} modulo p";
        long lambda = -1;
        for (int k = 0; k < modular.Count; k++)
        {
            if (!f.TryGetValue(modular.Monomials[k], out var c)) return $"monomial {modular.Monomials[k]} missing in the exact polynomial";
            if (field.Map(c) is not long r) return "a denominator vanishes modulo p";
            if (k == 0)
            {
                if (r == 0) return "the exact leading coefficient vanishes modulo p";
                lambda = field.Mul(r, field.Inverse(modular.Coefficients[0]));
            }
            if (r != field.Mul(lambda, modular.Coefficients[k])) return "the exact polynomial is not a multiple of the modular one modulo p";
        }
        return null;
    }
}
