using BigInteger = System.Numerics.BigInteger;
using System.Text;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Quadrupolar;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Quadrupolar;

/// <summary>The planes of coinciding population rates, J₂/J₁ = t with J₀ free, where the sodium page's parity route
/// leaves a centre possible at spins I ≥ 2, read by the numerical-range bound at 2I = 4 … 12, together with the
/// lines J₁ = 0 (every spin) and J₂ = 0 (half-integer spins). The verifier
/// (simulations/sodium_two_mirrors.py) finds every plane with t &gt; 0 exactly and checks that each strip in
/// <see cref="Strips"/> holds exactly one. This test checks by a second route that each strip holds a coincidence: the
/// order of the population rates even and odd under m ↦ −m, read exactly by inertia (Sylvester), differs at its two
/// ends. Then it reads the bound on the strip, J = (x, 1, t) with x = J₀/J₁ ≥ 0 and t in the strip: decided on a
/// window of x, certified by its four corners, since the region where D + 2s·1 is positive definite is a convex cone
/// in the densities; or silent at every x ≥ 0, certified by rational test vectors y whose quadratic forms
/// yᵀ(D + 2s·1)y, affine in (x, t), are ≤ 0 at the corners of pieces that cover x ≥ 0 across the strip. Floating
/// point only proposes the corners and the vectors; every verdict is exact.</summary>
public class NumericalRangePlanesTests
{
    private readonly ITestOutputHelper _out;
    public NumericalRangePlanesTests(ITestOutputHelper output) => _out = output;

    /// <summary>(2I, t_lo, t_hi): one plane J₂/J₁ = t with t &gt; 0 in each open strip, every such plane at
    /// 2I = 4 … 12 listed once. The verifier reads this table from this file between the two markers.</summary>
    // PLANES-BEGIN
    private static readonly (int TwoI, string Lo, string Hi)[] Strips =
    {
        (4, "28469/2000000", "71173/5000000"),
        (4, "177429/100000", "17743/10000"),
        (6, "123497/50000000", "49399/20000000"),
        (6, "25947/2500000", "103789/10000000"),
        (6, "389163/1000000", "97291/250000"),
        (6, "12398/3125", "396737/100000"),
        (6, "288017/50000", "115207/20000"),
        (7, "53639/500000", "107279/1000000"),
        (7, "15467/100000", "154671/1000000"),
        (7, "999999/1000000", "1000001/1000000"),
        (7, "32721/2000", "81803/5000"),
        (8, "152367/200000000", "190459/250000000"),
        (8, "46187/25000000", "184749/100000000"),
        (8, "818267/100000000", "204567/25000000"),
        (8, "771631/10000000", "48227/625000"),
        (8, "28153/4000", "351913/50000"),
        (8, "170181/5000", "340363/10000"),
        (9, "441327/10000000", "27583/625000"),
        (9, "81239/500000", "162479/1000000"),
        (9, "999999/1000000", "1000001/1000000"),
        (9, "5463/5000", "109261/100000"),
        (9, "70723/25000", "282893/100000"),
        (9, "173773/10000", "86887/5000"),
        (9, "305629/5000", "611259/10000"),
        (10, "77689/250000000", "310757/1000000000"),
        (10, "599969/1000000000", "599971/1000000000"),
        (10, "72481/50000000", "144963/100000000"),
        (10, "136091/20000000", "85057/12500000"),
        (10, "269483/10000000", "67371/2500000"),
        (10, "700027/10000000", "175007/2500000"),
        (10, "388049/1000000", "7761/20000"),
        (10, "155673/50000", "311347/100000"),
        (10, "54967/5000", "21987/2000"),
        (10, "195843/5000", "391687/10000"),
        (10, "100211/1000", "25053/250"),
        (11, "178171/10000000", "44543/2500000"),
        (11, "384359/10000000", "9609/250000"),
        (11, "109893/500000", "219787/1000000"),
        (11, "9167/5000", "183341/100000"),
        (11, "261429/50000", "522859/100000"),
        (11, "56883/10000", "568831/100000"),
        (11, "33899/2000", "21187/1250"),
        (11, "358471/5000", "716943/10000"),
        (11, "154159/1000", "3854/25"),
        (12, "1501089/10000000000", "1501091/10000000000"),
        (12, "127183/500000000", "254367/1000000000"),
        (12, "487789/1000000000", "48779/100000000"),
        (12, "2943/2500000", "117721/100000000"),
        (12, "589057/100000000", "294529/50000000"),
        (12, "60253/5000000", "120507/10000000"),
        (12, "232367/10000000", "14523/625000"),
        (12, "133073/2000000", "332683/5000000"),
        (12, "137403/1000000", "34351/250000"),
        (12, "31679/25000", "126717/100000"),
        (12, "214523/100000", "53631/25000"),
        (12, "92437/10000", "924371/100000"),
        (12, "79189/5000", "158379/10000"),
        (12, "436357/10000", "218179/5000"),
        (12, "117791/1000", "14724/125"),
        (12, "226127/1000", "28266/125"),
    };
    // PLANES-END

    /// <summary>The table's shape per spin, 2I = 4 … 12 (sixty planes in all; that these are all of them is the
    /// verifier's claim, not this test's), and the planes on which the bound decides on a window (32, computed here).
    /// The silent ones are those with t &lt; 0.071 from 2I = 6 on, and with t &gt; 70 at 2I = 11 and t &gt; 15 at 2I = 12:
    /// a description of these sixty, pinned, not a law.</summary>
    private static readonly int[] PlaneCounts = { 2, 0, 5, 4, 6, 7, 11, 9, 16 };
    private static readonly int[] DecidedCounts = { 2, 0, 3, 4, 3, 6, 5, 5, 4 };

    [Fact]
    public void EveryStrip_HoldsACoincidenceOfPopulationRates_ReadByInertiaAtItsEnds()
    {
        foreach (var (twoI, lo, hi) in Strips)
        {
            var parts = Parts.Of(SpinSpace.Of(twoI));
            string a = parts.ParityOrder(Parse(lo)), b = parts.ParityOrder(Parse(hi));
            Assert.True(a != b, $"2I = {twoI}, strip ({lo}, {hi}): the even/odd order {a} is the same at both ends");
        }
    }

    [Fact]
    public void TheBound_OnEveryPlane_DecidesOnACertifiedWindow_OrIsCertifiedSilent()
    {
        Assert.Equal(PlaneCounts.Sum(), Strips.Length);
        var sb = new StringBuilder();
        var decided = new int[PlaneCounts.Length];
        for (int twoI = 4; twoI <= 12; twoI++)
        {
            var parts = Parts.Of(SpinSpace.Of(twoI));
            int planes = 0;
            foreach (var (_, lo, hi) in Strips.Where(s => s.TwoI == twoI))
            {
                planes++;
                var family = new Family(Base: (0, 1, 0), Direction: (0, 0, 1), Parse(lo), Parse(hi));
                var verdict = parts.Read(family);
                if (verdict.Window is not null) decided[twoI - 4]++;
                // where the silent planes lie, as the page says it: J₂ far below J₁ from 2I = 6 on, far above at 2I = 11, 12
                double t = D(family.TMid);
                bool silentAsSaid = (twoI >= 6 && t < 0.071) || (twoI == 11 && t > 70) || (twoI == 12 && t > 15);
                Assert.True(silentAsSaid == (verdict.Window is null), $"2I = {twoI}, t ≈ {t}: the silent planes are not where the page says");
                sb.AppendLine($"2I = {twoI}, t in ({lo}, {hi}): {verdict}");
            }
            Assert.Equal(PlaneCounts[twoI - 4], planes);
        }
        _out.WriteLine(sb.ToString());
        Assert.Equal(DecidedCounts, decided);
    }

    /// <summary>The two lines, J₁ = 0 at every spin (J = (x, 0, 1)) and J₂ = 0 at half-integer spins (J = (x, 1, 0)):
    /// the verdicts pinned are the ones the page quotes.</summary>
    [Fact]
    public void TheBound_OnTheTwoLines()
    {
        var sb = new StringBuilder();
        var j1Zero = new List<int>();
        var j2Zero = new List<int>();
        for (int twoI = 4; twoI <= 12; twoI++)
        {
            var parts = Parts.Of(SpinSpace.Of(twoI));
            var v = parts.Read(new Family(Base: (0, 0, 1), Direction: (0, 0, 0), BigRational.Zero, BigRational.Zero));
            sb.AppendLine($"2I = {twoI}, J₁ = 0, x = J₀/J₂: {v}");
            if (v.Window is not null) j1Zero.Add(twoI);
            if (twoI % 2 == 1)
            {
                var u = parts.Read(new Family(Base: (0, 1, 0), Direction: (0, 0, 0), BigRational.Zero, BigRational.Zero));
                sb.AppendLine($"2I = {twoI}, J₂ = 0, x = J₀/J₁: {u}");
                if (u.Window is not null) j2Zero.Add(twoI);
            }
        }
        _out.WriteLine(sb.ToString());
        Assert.Equal(new[] { 4, 5, 6, 7, 8, 9 }, j1Zero);
        Assert.Equal(new[] { 5 }, j2Zero);
    }

    /// <summary>The blocks the certificates use are the generator's: D + 2s·1 is linear in the densities, and the
    /// combination of the three unit-density blocks equals the generator's own block at a sample point.</summary>
    [Theory]
    [InlineData(4)]
    [InlineData(9)]
    [InlineData(12)]
    public void TheCombinedBlocks_AreTheGenerators(int twoI)
    {
        var s = SpinSpace.Of(twoI);
        var parts = Parts.Of(s);
        var j = (new BigRational(7, 3), new BigRational(2), new BigRational(5, 11));
        var g = new QuadrupolarGenerator(s, new QuadrupolarRates(j.Item1, j.Item2, j.Item3), SecularHamiltonian.None);
        var twoS = 2 * g.MeanRate();
        for (int p = -(s.Levels - 1); p <= s.Levels - 1; p++)
        {
            var direct = g.DissipatorBlock(p);
            for (int i = 0; i < direct.GetLength(0); i++) direct[i, i] += twoS;
            var combined = parts.Shifted(p, j);
            for (int r = 0; r < direct.GetLength(0); r++)
                for (int c = 0; c < direct.GetLength(0); c++)
                    Assert.True(direct[r, c] == combined[r, c], $"2I = {twoI}, p = {p}, [{r},{c}]");
        }
    }

    /// <summary>A control for the silence certificate: at extreme narrowing the bound decides, so no family through
    /// it can be certified silent, and the reader must find a window there (2I = 7, where extreme narrowing lies on
    /// the plane t = 1, the window includes x = 1).</summary>
    [Fact]
    public void Control_AtExtremeNarrowing_TheReaderFindsAWindowContainingIt()
    {
        var parts = Parts.Of(SpinSpace.Of(7));
        var v = parts.Read(new Family(Base: (0, 1, 0), Direction: (0, 0, 1), new BigRational(999999, 1000000), new BigRational(1000001, 1000000)));
        Assert.NotNull(v.Window);
        Assert.True((v.Window!.Value.Lo - 1).Sign <= 0 && (v.Window!.Value.Hi - 1).Sign >= 0, v.ToString());
    }

    /// <summary>Controls for the coincidence check: across a strip that holds no plane (2I = 4, t in (1/2, 1), between
    /// the two planes of the table) the order of the even and odd rates stays the same.</summary>
    [Fact]
    public void Control_AStripWithoutAPlane_KeepsTheOrder()
    {
        var parts = Parts.Of(SpinSpace.Of(4));
        Assert.Equal(parts.ParityOrder(new BigRational(1, 2)), parts.ParityOrder(BigRational.One));
    }

    /// <summary>Control for the silence certificate: on a plane where the bound decides (2I = 7, extreme narrowing on
    /// t = 1) it cannot be built.</summary>
    [Fact]
    public void Control_SilenceCannotBeCertified_WhereTheBoundDecides()
    {
        var parts = Parts.Of(SpinSpace.Of(7));
        var f = new Family(Base: (0, 1, 0), Direction: (0, 0, 1), new BigRational(999999, 1000000), new BigRational(1000001, 1000000));
        var e = Assert.Throws<Xunit.Sdk.TrueException>(() => parts.CertifySilence(f));
        Assert.Contains("the vector does not hold", e.Message);
    }

    /// <summary>The quadratic form the silence certificate reads, yᵀΦ⁻²·B·y with B the gauged block, is vᵀSv for
    /// v = Φ⁻¹y and S the physical block of D + 2s·1: checked against the physical generator built in floating point
    /// from the float spin matrices (an independent route through no gauge), on every order p ≠ 0, where the gauge
    /// weights are not all 1. Error model: a float quadratic form carries at most a few n·eps of Σ|v_r S_rc v_c|.</summary>
    [Theory]
    [InlineData(8)]
    [InlineData(12)]
    public void TheSilenceForm_IsThePhysicalQuadraticForm(int twoI)
    {
        var s = SpinSpace.Of(twoI);
        var parts = Parts.Of(s);
        var j = (new BigRational(3, 2), BigRational.One, new BigRational(1, 50));
        var g = new QuadrupolarGenerator(s, new QuadrupolarRates(j.Item1, j.Item2, j.Item3), SecularHamiltonian.None);
        var lf = g.Float();
        double twoS = 2 * D(g.MeanRate());
        int d = s.Levels;
        for (int p = -(d - 1); p <= d - 1; p++)
        {
            if (p == 0) continue;
            var (y, coef, phi2) = parts.FormFromLowestVector(p, (D(j.Item1), D(j.Item2), D(j.Item3)));
            var exact = j.Item1 * coef[0] + j.Item2 * coef[1] + j.Item3 * coef[2];
            var cells = QuadrupolarGenerator.OrderCells(d, p);
            int n = cells.Length;
            var v = new double[n];
            for (int r = 0; r < n; r++) v[r] = D(y[r]) / Math.Sqrt(D(phi2[r]));
            double q = 0, scale = 0;
            for (int r = 0; r < n; r++)
                for (int c = 0; c < n; c++)
                {
                    double sv = lf[cells[r], cells[c]].Real + (r == c ? twoS : 0);
                    q += v[r] * sv * v[c];
                    scale += Math.Abs(v[r] * sv * v[c]);
                }
            double diff = Math.Abs(D(exact) - q);
            Assert.True(diff <= 16 * n * 2.220446049250313e-16 * scale, $"2I = {twoI}, p = {p}: exact {D(exact)}, float {q}, |diff| {diff:E2}, scale {scale:E2}");
        }
    }

    // ------------------------------------------------------------------ machinery

    private static BigRational Parse(string q)
    {
        var parts = q.Split('/');
        return parts.Length == 1 ? new BigRational(BigInteger.Parse(parts[0])) : new BigRational(BigInteger.Parse(parts[0]), BigInteger.Parse(parts[1]));
    }

    /// <summary>A double as the exact dyadic rational it is.</summary>
    private static BigRational Exact(double x)
    {
        if (x == 0) return BigRational.Zero;
        long bits = BitConverter.DoubleToInt64Bits(x);
        bool negative = bits < 0;
        int exponent = (int)((bits >> 52) & 0x7FF);
        long mantissa = bits & 0xFFFFFFFFFFFFFL;
        if (exponent == 0) exponent = 1; else mantissa |= 1L << 52;
        exponent -= 1075;
        BigInteger num = mantissa, den = BigInteger.One;
        if (exponent > 0) num <<= exponent; else den <<= -exponent;
        return new BigRational(negative ? -num : num, den);
    }

    /// <summary>x rounded to four significant digits, toward the inside of a window (down for the upper end).</summary>
    private static BigRational Inward(double x, bool down)
    {
        if (x == 0) return BigRational.Zero;
        int e = (int)Math.Floor(Math.Log10(Math.Abs(x)));
        double scale = Math.Pow(10, 3 - e);
        double m = down ? Math.Floor(x * scale) : Math.Ceiling(x * scale);
        return e <= 3 ? new BigRational((BigInteger)m, BigInteger.Pow(10, 3 - e)) : new BigRational((BigInteger)m * BigInteger.Pow(10, e - 3));
    }

    private static double D(BigRational r) => SpinSpace.ToDouble(r);

    /// <summary>The densities J(x, t) = x·(1, 0, 0) + Base + t·Direction with t in [TLo, THi].</summary>
    private sealed record Family((int, int, int) Base, (int, int, int) Direction, BigRational TLo, BigRational THi)
    {
        public (BigRational, BigRational, BigRational) At(BigRational x, BigRational t) =>
            (x + Base.Item1 + t * Direction.Item1, Base.Item2 + t * Direction.Item2, Base.Item3 + t * Direction.Item3);
        public BigRational TMid => (TLo + THi) / 2;
    }

    private readonly record struct Verdict((BigRational Lo, BigRational Hi)? Window, int Pieces)
    {
        public override string ToString() => Window is { } w
            ? $"decided for x in [{w.Lo}, {w.Hi}] (≈ [{D(w.Lo):G4}, {D(w.Hi):G4}]), four corners exact"
            : $"silent at every x ≥ 0, {Pieces} pieces exact";
    }

    /// <summary>The three unit-density blocks of D + 2s·1 by coherence order, in the gauge, with the cells' gauge
    /// weights φ² (Φ scales |i⟩⟨j| by σ_i/σ_j; the physical block is Φ⁻¹·B·Φ, real symmetric).</summary>
    private sealed class Parts
    {
        private readonly SpinSpace _s;
        private readonly Dictionary<int, BigRational[][,]> _blocks = new();   // p → [J₀, J₁, J₂] parts, each with its 2s
        private readonly Dictionary<int, BigRational[]> _phi2 = new();
        private readonly Dictionary<int, Matrix<double>[]> _physicalFloat = new();   // p → the three unit parts, physical, symmetric
        private Parts(SpinSpace s) => _s = s;

        public static Parts Of(SpinSpace s)
        {
            var parts = new Parts(s);
            var units = new[] { (1, 0, 0), (0, 1, 0), (0, 0, 1) }
                .Select(u => new QuadrupolarGenerator(s, new QuadrupolarRates(u.Item1, u.Item2, u.Item3), SecularHamiltonian.None)).ToArray();
            int d = s.Levels;
            for (int p = -(d - 1); p <= d - 1; p++)
            {
                var three = new BigRational[3][,];
                for (int k = 0; k < 3; k++)
                {
                    var b = units[k].DissipatorBlock(p);
                    var twoS = 2 * units[k].MeanRate();
                    for (int i = 0; i < b.GetLength(0); i++) b[i, i] += twoS;
                    three[k] = b;
                }
                parts._blocks[p] = three;
                var cells = QuadrupolarGenerator.OrderCells(d, p);
                parts._phi2[p] = cells.Select(c => s.GaugeWeightSquared(c / d) / s.GaugeWeightSquared(c % d)).ToArray();
                var phi2 = parts._phi2[p];
                parts._physicalFloat[p] = three.Select(bk =>
                {
                    int n = phi2.Length;
                    var m = Matrix<double>.Build.Dense(n, n);
                    for (int r = 0; r < n; r++)
                        for (int c = 0; c < n; c++)
                            m[r, c] = D(bk[r, c]) * Math.Sqrt(D(phi2[c]) / D(phi2[r]));
                    return (m + m.Transpose()) / 2;
                }).ToArray();
            }
            return parts;
        }

        public BigRational[,] Shifted(int p, (BigRational, BigRational, BigRational) j)
        {
            var three = _blocks[p];
            int n = three[0].GetLength(0);
            var m = new BigRational[n, n];
            for (int r = 0; r < n; r++)
                for (int c = 0; c < n; c++)
                    m[r, c] = j.Item1 * three[0][r, c] + j.Item2 * three[1][r, c] + j.Item3 * three[2][r, c];
            return m;
        }

        private Matrix<double> PhysicalFloat(int p, (double, double, double) j)
        {
            var f = _physicalFloat[p];
            return j.Item1 * f[0] + j.Item2 * f[1] + j.Item3 * f[2];
        }

        /// <summary>min over orders of the smallest eigenvalue of D + 2s·1, in floating point, with its order and vector.</summary>
        private (double Value, int P, Vector<double> V) Smallest((double, double, double) j)
        {
            (double, int, Vector<double>) best = (double.PositiveInfinity, 0, Vector<double>.Build.Dense(1));
            foreach (int p in _blocks.Keys)
            {
                var evd = PhysicalFloat(p, j).Evd(Symmetricity.Symmetric);
                double v = evd.EigenValues[0].Real;
                int arg = 0;
                for (int i = 1; i < evd.EigenValues.Count; i++)
                    if (evd.EigenValues[i].Real < v) { v = evd.EigenValues[i].Real; arg = i; }
                if (v < best.Item1) best = (v, p, evd.EigenVectors.Column(arg));
            }
            return best;
        }

        private bool DecidedExactly((BigRational, BigRational, BigRational) j) =>
            NumericalRangeBound.ExcludesEveryCentre(new QuadrupolarGenerator(_s, new QuadrupolarRates(j.Item1, j.Item2, j.Item3), SecularHamiltonian.None));

        public Verdict Read(Family f)
        {
            double tm = D(f.TMid);
            (double, double, double) J(double x) { var j = f.At(Exact(x), f.TMid); return (D(j.Item1), D(j.Item2), D(j.Item3)); }
            double G(double x) => Smallest(J(x)).Value;
            // proposal: the largest smallest eigenvalue over x ≥ 0 (concave in x): a scan, then a golden section
            double span = 4 * (1 + Math.Abs(tm)) + 20;
            const int steps = 400;
            int bestI = 0;
            double bestG = G(0);
            for (int i = 1; i <= steps; i++)
            {
                double g = G(span * i / steps);
                if (g > bestG) { bestG = g; bestI = i; }
            }
            double a = span * Math.Max(0, bestI - 1) / steps, b = span * Math.Min(steps, bestI + 1) / steps;
            for (int i = 0; i < 80; i++)
            {
                double m1 = a + 0.382 * (b - a), m2 = a + 0.618 * (b - a);
                if (G(m1) < G(m2)) a = m1; else b = m2;
            }
            double bestX = (a + b) / 2;
            bestG = Math.Max(bestG, G(bestX));
            return bestG > 0 ? CertifyWindow(f, G, bestX) : CertifySilence(f);
        }

        private Verdict CertifyWindow(Family f, Func<double, double> g, double xBest)
        {
            double Edge(double inside, double outside)
            {
                for (int i = 0; i < 100; i++)
                {
                    double c = (inside + outside) / 2;
                    if (g(c) > 0) inside = c; else outside = c;
                }
                return inside;
            }
            double lo = g(0) > 0 ? 0 : Edge(xBest, 0);
            double hi = Edge(xBest, 1e6);
            Assert.True(hi < 1e6, "a window without an upper end: the bound would decide at J₀ alone");
            double margin = 0.01 * (hi - lo);
            var xa = lo == 0 ? BigRational.Zero : Inward(lo + margin, down: false);
            var xb = Inward(hi - margin, down: true);
            Assert.True((xb - xa).Sign > 0, $"window too narrow to certify: [{lo}, {hi}]");
            foreach (var x in new[] { xa, xb })
                foreach (var t in new[] { f.TLo, f.THi })
                    Assert.True(DecidedExactly(f.At(x, t)), $"corner (x = {x}, t = {t}) not decided");
            // the same exact decision says no just outside the window, so the corners are not a check that cannot fail
            var outside = new List<BigRational> { Inward(1.1 * hi + 0.1, down: false) };
            if (lo > 0) outside.Add(Inward(0.9 * lo, down: true));
            foreach (var x in outside)
                foreach (var t in new[] { f.TLo, f.THi })
                    Assert.False(DecidedExactly(f.At(x, t)), $"decided outside the window at x = {x}, t = {t}");
            return new Verdict((xa, xb), 0);
        }

        /// <summary>For a float physical vector v of order p: y = Φv as exact dyadic rationals, and the three
        /// coefficients of q = yᵀΦ⁻²·B·y = v′ᵀSv′ (v′ = Φ⁻¹y, the vector the certificate actually uses) on the unit
        /// densities, so q(J) = J₀·a + J₁·b + J₂·c.</summary>
        public (BigRational[] Y, BigRational[] Coef) Form(int p, Vector<double> v)
        {
            var phi2 = _phi2[p];
            var y = new BigRational[v.Count];
            for (int i = 0; i < v.Count; i++) y[i] = Exact(v[i] * Math.Sqrt(D(phi2[i])));
            var coef = new BigRational[3];
            for (int k = 0; k < 3; k++)
            {
                var bk = _blocks[p][k];
                var acc = BigRational.Zero;
                for (int r = 0; r < y.Length; r++)
                {
                    if (y[r].IsZero) continue;
                    var row = BigRational.Zero;
                    for (int c = 0; c < y.Length; c++) if (!y[c].IsZero && !bk[r, c].IsZero) row += bk[r, c] * y[c];
                    acc += y[r] / phi2[r] * row;
                }
                coef[k] = acc;
            }
            return (y, coef);
        }

        /// <summary>The form of the lowest float eigenvector of order p at the densities j, with the cells' φ².</summary>
        public (BigRational[] Y, BigRational[] Coef, BigRational[] Phi2) FormFromLowestVector(int p, (double, double, double) j)
        {
            var evd = PhysicalFloat(p, j).Evd(Symmetricity.Symmetric);
            int arg = 0;
            for (int i = 1; i < evd.EigenValues.Count; i++)
                if (evd.EigenValues[i].Real < evd.EigenValues[arg].Real) arg = i;
            var (y, coef) = Form(p, evd.EigenVectors.Column(arg));
            return (y, coef, _phi2[p]);
        }

        /// <summary>Pieces [x_k, x_{k+1}] (the last one unbounded) with a rational vector y in one order block each:
        /// q(x, t) = yᵀΦ⁻²·B(x, t)·y = vᵀ S v for v = Φ⁻¹y, S the physical block of D + 2s·1, affine in (x, t), checked
        /// ≤ 0 at the piece's corners (and of non-positive slope in x on the last piece). Then S is not positive definite
        /// anywhere on x ≥ 0 across the strip, so the bound decides nowhere there.</summary>
        public Verdict CertifySilence(Family f)
        {
            var x = BigRational.Zero;
            for (int piece = 1; piece <= 500; piece++)
            {
                var (_, p, v) = Smallest(JOf(f, D(x)));
                var (_, coef) = Form(p, v);
                BigRational Q(BigRational xx, BigRational t)
                {
                    var j = f.At(xx, t);
                    return j.Item1 * coef[0] + j.Item2 * coef[1] + j.Item3 * coef[2];
                }
                Assert.True(Q(x, f.TLo).Sign <= 0 && Q(x, f.THi).Sign <= 0, $"piece {piece} at x = {D(x)}: the vector does not hold");
                var slope = coef[0];   // ∂q/∂x: x enters J₀ only
                if (slope.Sign <= 0) return new Verdict(null, piece);
                // q stays ≤ 0 up to the smaller root over the two ends of the strip
                var next = new[] { f.TLo, f.THi }.Select(t => x - Q(x, t) / slope).Aggregate((a, b) => (a - b).Sign <= 0 ? a : b);
                Assert.True((next - x).Sign > 0, $"no progress at x = {D(x)}");
                x = next;
            }
            throw new InvalidOperationException("silence certificate: too many pieces");
        }

        private static (double, double, double) JOf(Family f, double x)
        {
            var j = f.At(Exact(x), f.TMid);
            return (D(j.Item1), D(j.Item2), D(j.Item3));
        }

        /// <summary>The order of the population rates at J = (0, 1, t), each labelled e or o by its parity under
        /// m ↦ −m, read exactly: rational separators proposed between the float rates, the count of even and odd
        /// rates below each separator read by the signs of the pivots of W_± − μ·1 (W_± the parity blocks in the
        /// unnormalized basis e_k ± e_{d−1−k}, a diagonal similarity of a real symmetric matrix, so its pivots are the
        /// symmetric one's and Sylvester's law counts its negative eigenvalues).</summary>
        public string ParityOrder(BigRational t)
        {
            int d = _s.Levels;
            var pop = new QuadrupolarGenerator(_s, new QuadrupolarRates(BigRational.Zero, BigRational.One, t), SecularHamiltonian.None).DissipatorBlock(0);
            for (int r = 0; r < d; r++)
                for (int c = 0; c < d; c++)
                    Assert.True(pop[r, c] == pop[c, r] && pop[r, c] == pop[d - 1 - r, d - 1 - c],
                        "the populations block is not symmetric, or not invariant under m ↦ −m");
            var even = ParityBlock(pop, +1);
            var odd = ParityBlock(pop, -1);
            var rates = new List<(double Value, char Label)>();
            foreach (var (blk, label) in new[] { (even, 'e'), (odd, 'o') })
            {
                if (blk.GetLength(0) == 0) continue;
                foreach (var ev in SymmetricFloat(blk, label == 'e').Evd(Symmetricity.Symmetric).EigenValues) rates.Add((ev.Real, label));
            }
            rates.Sort((a, b) => a.Value.CompareTo(b.Value));
            var order = new StringBuilder();
            int ePrev = 0, oPrev = 0;
            for (int k = 1; k <= rates.Count; k++)
            {
                int eNow, oNow;
                if (k == rates.Count) { eNow = even.GetLength(0); oNow = odd.GetLength(0); }
                else
                {
                    // the float only proposes the separator; the exact counts below decide whether it separates
                    Assert.True(rates[k].Value > rates[k - 1].Value, $"two float rates coincide at t = {t}");
                    var mu = Exact((rates[k].Value + rates[k - 1].Value) / 2);
                    eNow = NegativeCount(even, mu);
                    oNow = NegativeCount(odd, mu);
                    Assert.Equal(k, eNow + oNow);
                }
                if (eNow - ePrev == 1 && oNow == oPrev) order.Append('e');
                else if (oNow - oPrev == 1 && eNow == ePrev) order.Append('o');
                else throw new InvalidOperationException("one rate per gap expected");
                (ePrev, oPrev) = (eNow, oNow);
            }
            return order.ToString();
        }

        private BigRational[,] ParityBlock(BigRational[,] w, int sign)
        {
            int d = w.GetLength(0);
            var basis = new List<int>();
            for (int k = 0; k < d; k++)
            {
                int partner = d - 1 - k;
                if (k < partner || (k == partner && sign == 1)) basis.Add(k);
            }
            int n = basis.Count;
            var b = new BigRational[n, n];
            for (int r = 0; r < n; r++)
                for (int c = 0; c < n; c++)
                {
                    int kc = basis[c], partner = d - 1 - kc;
                    b[r, c] = kc == partner ? w[basis[r], kc] : w[basis[r], kc] + sign * w[basis[r], partner];
                }
            return b;
        }

        /// <summary>The parity block in the orthonormal basis: the pair vectors carry 1/√2, the middle level 1.</summary>
        private Matrix<double> SymmetricFloat(BigRational[,] b, bool even)
        {
            int n = b.GetLength(0);
            int d = _s.Levels;
            bool hasMiddle = even && d % 2 == 1;
            double Norm(int i) => hasMiddle && i == n - 1 ? 1.0 : Math.Sqrt(2.0);
            var m = Matrix<double>.Build.Dense(n, n);
            for (int r = 0; r < n; r++)
                for (int c = 0; c < n; c++)
                    m[r, c] = D(b[r, c]) * Norm(r) / Norm(c);
            return (m + m.Transpose()) / 2;
        }

        private static int NegativeCount(BigRational[,] a, BigRational mu)
        {
            int n = a.GetLength(0);
            var m = (BigRational[,])a.Clone();
            for (int i = 0; i < n; i++) m[i, i] -= mu;
            int negative = 0;
            for (int k = 0; k < n; k++)
            {
                if (m[k, k].IsZero) throw new InvalidOperationException("a zero pivot: choose another separator");
                if (m[k, k].Sign < 0) negative++;
                for (int r = k + 1; r < n; r++)
                {
                    if (m[r, k].IsZero) continue;
                    var factor = m[r, k] / m[k, k];
                    for (int c = k; c < n; c++) m[r, c] -= factor * m[k, c];
                }
            }
            return negative;
        }
    }
}
