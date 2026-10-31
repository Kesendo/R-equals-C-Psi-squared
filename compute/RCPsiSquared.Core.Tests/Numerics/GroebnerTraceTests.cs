using RCPsiSquared.Core.Numerics;
using Xunit;

namespace RCPsiSquared.Core.Tests.Numerics;

/// <summary>The modular Gröbner engine with its exact replay: grevlex on packed monomials, Buchberger over GF(p)
/// with a trace, and the replay over ℚ(i) that turns a modular unit into a proof. Systems whose answer is known by
/// hand, a bad prime the replay must refuse, and systems that are not empty.</summary>
public class GroebnerTraceTests
{
    private static readonly long Prime = GaussianModP.PrimesOneModFour(1)[0];
    private static readonly GaussianModP Field = new(Prime);

    private static PackedMonomial Mono(params int[] e)
    {
        var m = PackedMonomial.One;
        for (int v = 0; v < e.Length; v++) if (e[v] > 0) m *= PackedMonomial.Variable(v, e[v]);
        return m;
    }

    private static Dictionary<PackedMonomial, GaussianInteger> Poly(params (GaussianInteger C, int[] E)[] terms)
    {
        var f = new Dictionary<PackedMonomial, GaussianInteger>();
        foreach (var (c, e) in terms)
        {
            var m = Mono(e);
            f[m] = f.TryGetValue(m, out var old) ? old + c : c;
        }
        return f;
    }

    private static (GroebnerTraceModP Trace, GroebnerExactReplay.Result Replay) Run(int nvars, params Dictionary<PackedMonomial, GaussianInteger>[] exact)
    {
        var modular = exact.Select(f => ModPPolynomial.FromTerms(f.Select(kv => (kv.Key, Field.Map(kv.Value))), Prime)).ToArray();
        var trace = GroebnerTraceModP.Compute(modular, Prime, nvars);
        return (trace, GroebnerExactReplay.ReplayUnit(trace, exact, Field));
    }

    [Fact]
    public void TheRootOfMinusOne_SquaresToMinusOne()
    {
        Assert.Equal(1, Prime % 4);
        Assert.Equal(Prime - 1, Field.Mul(Field.RootOfMinusOne, Field.RootOfMinusOne));
        Assert.Equal(Field.RootOfMinusOne, Field.Map(GaussianInteger.I));
    }

    [Fact]
    public void Grevlex_OrdersTheDegreeTwoMonomialsOfThreeVariables()
    {
        // x > y > z: x² > xy > y² > xz > yz > z²
        var expected = new[] { Mono(2, 0, 0), Mono(1, 1, 0), Mono(0, 2, 0), Mono(1, 0, 1), Mono(0, 1, 1), Mono(0, 0, 2) };
        var sorted = expected.Reverse().OrderBy(m => m, PackedMonomial.Descending).ToArray();
        Assert.Equal(expected, sorted);
        Assert.True(PackedMonomial.CompareGrevlex(Mono(0, 0, 3), Mono(2, 0, 0)) > 0);   // degree first
    }

    [Fact]
    public void Monomials_DivideMultiplyAndTakeTheLcm()
    {
        var a = Mono(2, 0, 1, 0, 0, 0, 0, 0, 0, 3);      // across a word boundary (variable 9)
        var b = Mono(1, 4, 0, 0, 0, 0, 0, 0, 0, 1);
        Assert.Equal(Mono(2, 4, 1, 0, 0, 0, 0, 0, 0, 3), PackedMonomial.Lcm(a, b));
        Assert.Equal(10, PackedMonomial.Lcm(a, b).Degree);
        Assert.False(a.Divides(b));
        Assert.True(Mono(1, 0, 1).Divides(a));
        Assert.Equal(Mono(1, 0, 0, 0, 0, 0, 0, 0, 0, 3), Mono(1, 0, 1).DivideInto(a));
        Assert.Equal(Mono(3, 4, 1, 0, 0, 0, 0, 0, 0, 4), a * b);
        Assert.True(PackedMonomial.AreCoprime(Mono(2, 0, 0), Mono(0, 1, 5)));
        Assert.False(PackedMonomial.AreCoprime(a, b));
    }

    [Fact]
    public void Monomials_InTheLastWord_OrderAndOverflowAsTheyShould()
    {
        var hi = PackedMonomial.Variable(31);                 // the last variable, byte 7 of word 3
        var lo = PackedMonomial.Variable(24);                 // byte 0 of word 3
        Assert.True(PackedMonomial.CompareGrevlex(lo, hi) > 0);   // x_24 > x_31 at equal degree
        Assert.True(PackedMonomial.CompareGrevlex(PackedMonomial.Variable(0) * hi, lo * lo) < 0);
        Assert.Equal(1, (lo * hi)[31]);
        Assert.Throws<ArgumentOutOfRangeException>(() => hi[32]);
        var top = PackedMonomial.Variable(31, PackedMonomial.MaxExponent);
        Assert.Throws<OverflowException>(() => top * hi);
        Assert.Equal(PackedMonomial.MaxExponent, (top * PackedMonomial.Variable(30))[31]);
    }

    [Fact]
    public void TwoValuesForOneVariable_IsEmpty_AndTheReplayProvesIt()
    {
        var (trace, replay) = Run(1, Poly((1, new[] { 1 }), (-1, new[] { 0 })), Poly((1, new[] { 1 }), (-2, new[] { 0 })));
        Assert.True(trace.ContainsOneModP);
        Assert.True(replay.ProvesUnit, replay.Mismatch);
    }

    [Fact]
    public void XEqualsI_IsARootOfXSquaredPlusOne_ButNotOfXSquaredMinusOne()
    {
        var xMinusI = Poly((1, new[] { 1 }), (-GaussianInteger.I, new[] { 0 }));
        var (consistent, _) = Run(1, xMinusI, Poly((1, new[] { 2 }), (1, new[] { 0 })));
        Assert.False(consistent.ContainsOneModP);
        var (empty, replay) = Run(1, xMinusI, Poly((1, new[] { 2 }), (-1, new[] { 0 })));
        Assert.True(empty.ContainsOneModP);
        Assert.True(replay.ProvesUnit, replay.Mismatch);
    }

    [Fact]
    public void Rabinowitsch_TwoInvertibleFactorsWithAZeroProduct_IsEmpty()
    {
        // x·y = 0, x·z = 1, y·w = 1
        var (trace, replay) = Run(4, Poly((1, new[] { 1, 1 })), Poly((1, new[] { 1, 0, 1 }), (-1, new[] { 0 })),
                                  Poly((1, new[] { 0, 1, 0, 1 }), (-1, new[] { 0 })));
        Assert.True(trace.ContainsOneModP);
        Assert.True(replay.ProvesUnit, replay.Mismatch);
    }

    [Fact]
    public void CyclicThree_HasSolutions_AndPinningXToTwoEmptiesIt()
    {
        // x + y + z = 0, xy + yz + zx = 0, xyz = 1: the solutions are the cube roots of unity, permuted
        var e1 = Poly((1, new[] { 1, 0, 0 }), (1, new[] { 0, 1, 0 }), (1, new[] { 0, 0, 1 }));
        var e2 = Poly((1, new[] { 1, 1, 0 }), (1, new[] { 0, 1, 1 }), (1, new[] { 1, 0, 1 }));
        var e3 = Poly((1, new[] { 1, 1, 1 }), (-1, new[] { 0, 0, 0 }));
        var (open, _) = Run(3, e1, e2, e3);
        Assert.False(open.ContainsOneModP);
        // the basis contains z³ − 1 up to scale, as by hand
        Assert.Contains(open.Basis, k => open.Polynomials[k].Count == 2 && open.Polynomials[k].Lead == Mono(0, 0, 3));
        var (empty, replay) = Run(3, e1, e2, e3, Poly((1, new[] { 1, 0, 0 }), (-2, new[] { 0, 0, 0 })));
        Assert.True(empty.ContainsOneModP);
        Assert.True(replay.ProvesUnit, replay.Mismatch);
        Assert.True(replay.Replayed >= 2);
    }

    [Fact]
    public void ABadPrime_InTheMiddleOfTheRoute_IsCaughtAfterTheReduction()
    {
        // x − 1 and x·y − (1 + p)·y + 1: exactly y = 1/p, a point; modulo p the second reduces by the first to 1,
        // and the exact reduction leaves −p·y + 1, a term the modular one lacks
        var (trace, replay) = Run(2, Poly((1, new[] { 1, 0 }), (-1, new[] { 0, 0 })),
                                  Poly((1, new[] { 1, 1 }), (new GaussianInteger(-1 - Prime, 0), new[] { 0, 1 }), (1, new[] { 0, 0 })));
        Assert.True(trace.ContainsOneModP);
        Assert.False(replay.ProvesUnit);
        Assert.Contains("exact terms against", replay.Mismatch);
    }

    [Fact]
    public void ABadPrime_GivesAModularUnit_ThatTheReplayRefuses()
    {
        // p·x − 1 has the root 1/p, but modulo p it is the constant −1
        var (trace, replay) = Run(1, Poly((new GaussianInteger(Prime, 0), new[] { 1 }), (-1, new[] { 0 })));
        Assert.True(trace.ContainsOneModP);
        Assert.False(replay.ProvesUnit);
        Assert.NotNull(replay.Mismatch);
    }
}
