using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the binomial machinery the theodolite reads with: the exact Gaussian square root against
// squaring, the Gaussian common divisor against both of its multiples, a long elimination whose constants must stay
// small, and a rational point that needs a rational square root.
public class BinomialsTests
{
    static Binomial B(string monomial, GaussBig num, GaussBig den)
    {
        var e = new Dictionary<string, int>();
        foreach (var part in monomial.Split('·'))
        {
            var bits = part.Split('^');
            e[bits[0]] = bits.Length > 1 ? int.Parse(bits[1]) : 1;
        }
        return new Binomial(e, num, den, "test");
    }

    [Fact]
    public void The_Square_Root_Is_Exact_Or_Absent()
    {
        foreach (var (re, im) in new (long, long)[] { (-1, 0), (3, 4), (-3, 4), (0, 2), (5, -12), (9, 0), (-7, 24) })
        {
            var g = GaussBig.Of(re, im);
            var root = Binomials.SquareRoot(g);
            Assert.NotNull(root);
            Assert.Equal(g, root!.Value * root.Value);
        }
        foreach (var (re, im) in new (long, long)[] { (2, 0), (3, 0), (1, 1), (5, 0), (-3, 0) })
            Assert.Null(Binomials.SquareRoot(GaussBig.Of(re, im)));
    }

    [Fact]
    public void The_Common_Divisor_Divides_Both()
    {
        var two = GaussBig.Of(1, 1) * GaussBig.Of(1, -1);                     // 2 = (1 + i)(1 - i)
        var a = GaussBig.Of(2, 1) * GaussBig.Of(2, 1) * two;                 // (2 + i)^2 · 2
        var b = GaussBig.Of(2, 1) * GaussBig.Of(3, 0) * GaussBig.Of(1, 1);   // (2 + i) · 3 · (1 + i)
        var g = GaussBig.Gcd(a, b);
        Assert.Equal(new BigInteger(10), g.Norm);                           // (2 + i)(1 + i), norm 5·2, up to a unit
        Assert.Equal(a, GaussBig.ExactQuotient(a, g) * g);
        Assert.Equal(b, GaussBig.ExactQuotient(b, g) * g);
    }

    [Fact]
    public void A_Long_Elimination_Keeps_Its_Constants_Small()
    {
        // a chain of forty ratios u_k / u_(k+1) = (3 + 4i)/5 and its closure u_40 / u_0 = ((3 + 4i)/5)^-40: consistent,
        // and every reduced constant stays small because each step cancels the common divisor
        int n = 41;
        var unknowns = Enumerable.Range(0, n).Select(k => $"u{k}").ToList();
        var step = (GaussBig.Of(3, 4), GaussBig.Of(5));
        var eqs = Enumerable.Range(0, n - 1).Select(k => B($"u{k}·u{k + 1}^-1", step.Item1, step.Item2)).ToList();
        var closure = B($"u{n - 1}·u0^-1", GaussBig.Of(5).Pow(n - 1), GaussBig.Of(3, 4).Pow(n - 1));
        var red = Binomials.Solve(eqs.Append(closure), unknowns);
        Assert.NotNull(red);
        Assert.All(red!.Constants, c => Assert.True(c.Num.Norm < BigInteger.Pow(10, 4) && c.Den.Norm < BigInteger.Pow(10, 4)));
        var wrong = B($"u{n - 1}·u0^-1", GaussBig.Of(1), GaussBig.Of(1));
        Assert.Null(Binomials.Solve(eqs.Append(wrong), unknowns));
    }

    [Fact]
    public void A_Rational_Point_Takes_A_Rational_Square_Root()
    {
        // x/y = -4 and x·y = -1: y^2 = 1/4, so y = ±1/2 and x = -4y, no irrational root needed
        var red = Binomials.Solve(new[] { B("x·y^-1", GaussBig.Of(-4), GaussBig.One), B("x·y", GaussBig.Of(-1), GaussBig.One) },
            new[] { "x", "y" });
        var sol = Binomials.RationalSolution(red!);
        Assert.NotNull(sol);
        var (xn, xd) = sol!["x"];
        var (yn, yd) = sol["y"];
        Assert.Equal(GaussBig.Of(-4) * yn * xd, xn * yd);                    // x = -4y
        Assert.Equal(-(xd * yd), xn * yn);                                   // x·y = -1
    }
}
