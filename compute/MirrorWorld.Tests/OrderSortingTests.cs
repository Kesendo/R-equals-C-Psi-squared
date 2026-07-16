using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the mirror's order-sorting law (F131, adopted 2026-07-16 from
// docs/proofs/PROOF_MIRROR_ORDER_SORTING.md, Theorem A only -- the unitary column). The law:
// mirror conjugation reflects a parameter scan, (R x R) L(base + t*dir) (R x R) = L(base - t*dir)
// for an R-odd direction (sigma_eff = -1), and for an OPERATOR-R-even preparation the response
// orders sort by the readout parity q: q = +1 reads EVEN in t, q = -1 reads ODD in t; with an
// R-even direction (sigma_eff = +1) the q = -1 readout is IDENTICALLY ZERO and the q = +1 readout
// is generic. Trajectory face only, twin RK4 on the gamma axis -- no eigensolver anywhere; the
// spectral evenness corollary and Theorem B (antiunitary, Floquet) stay in the main repo.
public class OrderSortingTests
{
    static readonly World W = new();

    static double[] Uniform(int len, double v) => Enumerable.Repeat(v, len).ToArray();

    static double[] Scan(double[] baseP, double[] dir, double t)
    {
        var r = new double[baseP.Length];
        for (int l = 0; l < r.Length; l++) r[l] = baseP[l] + t * dir[l];
        return r;
    }

    // T0, the pencil face: the conjugation identity holds entry-wise on every joint-popcount block,
    // on all three parameter axes (gamma per site, J per bond, h per site), machine zero; and the
    // fence is real: a direction that does not negate EVERY component is rejected at O(1).
    [Fact]
    public void Conjugation_Identity_Holds_On_All_Three_Axes_And_The_Fence_Rejects()
    {
        var klein = new ParameterKlein(W, 5);
        var jU = Uniform(4, 1.0);
        var hZ = Uniform(5, 0.0);
        var gBase = OrderSorting.SymmetricBase(5);
        var gDir = OrderSorting.OddDirection(5);

        // gamma axis (F91's axis)
        Assert.True(klein.MirrorConjugationResidual(
            Scan(gBase, gDir, +0.37), jU, hZ,
            Scan(gBase, gDir, -0.37), jU, hZ) < 1e-12);

        // J axis: R-even base + R-odd bond direction (bond b mirrors to N-2-b)
        var jBase = new[] { 1.0, 0.85, 0.85, 1.0 };
        var jDir = new[] { -0.3, -0.1, 0.1, 0.3 };
        var gU = Uniform(5, 0.45);
        Assert.True(klein.MirrorConjugationResidual(
            gU, Scan(jBase, jDir, +0.5), hZ,
            gU, Scan(jBase, jDir, -0.5), hZ) < 1e-12);

        // h axis: R-odd longitudinal detuning
        var hDir = new[] { -0.2, -0.1, 0.0, 0.1, 0.2 };
        Assert.True(klein.MirrorConjugationResidual(
            gU, jU, Scan(Uniform(5, 0.3), hDir, +0.4),
            gU, jU, Scan(Uniform(5, 0.3), hDir, -0.4)) < 1e-12);

        // the fence: zero one component of the odd gamma direction -> not R-odd -> rejected O(1)
        var broken = (double[])gDir.Clone();
        broken[0] = 0.0;
        Assert.True(klein.MirrorConjugationResidual(
            Scan(gBase, broken, +0.37), jU, hZ,
            Scan(gBase, broken, -0.37), jU, hZ) > 0.01);
    }

    // the hypotheses as operator identities: the preparation is OPERATOR-R-even (R rho R = rho,
    // entry-wise rho[rev i, rev j] = rho[i, j]), and the leak's admixture is operator-R-odd.
    [Fact]
    public void The_Preparation_Parities_Are_Operator_Identities()
    {
        int n = 5, dim = 1 << n;
        var even = new Restless(W, n, 1.0, 0.05);
        OrderSorting.SeedEvenPrep(even, n);
        var leak = new Restless(W, n, 1.0, 0.05);
        OrderSorting.SeedLeakPrep(leak, n, 0.04);

        double worstEven = 0, worstOddPart = 0;
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                int ri = OrderSorting.Reverse(i, n), rj = OrderSorting.Reverse(j, n);
                worstEven = Math.Max(worstEven, (even[ri, rj] - even[i, j]).Magnitude);
                // leak = even + eps*odd: the difference to the even prep must be R-ODD exactly
                var d = leak[i, j] - even[i, j];
                var dRev = leak[ri, rj] - even[ri, rj];
                worstOddPart = Math.Max(worstOddPart, (dRev + d).Magnitude);
            }
        Assert.True(worstEven < 1e-14, $"even prep must be operator-R-even, worst {worstEven:E1}");
        Assert.True(worstOddPart < 1e-14, $"leak admixture must be operator-R-odd, worst {worstOddPart:E1}");
    }

    // the two response cells (sigma_eff = -1): the q = +1 readout is EVEN in the scan parameter and
    // the q = -1 readout is ODD, at every probed evolution time, machine precision from twin RK4
    // runs; the odd cell is non-vacuous (the forbidden orders are absent, the allowed ones are O(1)).
    [Fact]
    public void The_Even_And_Odd_Cells_Sort_By_Readout_Parity()
    {
        foreach (int n in new[] { 4, 5 })
        {
            var rep = new OrderSorting(W, n).Run();
            Assert.True(rep.EvenCellResidual < 1e-12, $"N={n}: even cell {rep.EvenCellResidual:E1}");
            Assert.True(rep.OddCellResidual < 1e-12, $"N={n}: odd cell {rep.OddCellResidual:E1}");
            Assert.True(rep.OddMagnitude > 1e-3, $"N={n}: odd cell vacuous, |<O_odd>| max {rep.OddMagnitude:E1}");
        }
    }

    // the sigma_eff = +1 column: with an R-even generator (no scan, or an R-even direction) the
    // R-odd readout vanishes IDENTICALLY (every RK4 tick, not just probe times) -- the zero cell is
    // a pure selection rule; and the (+,+) cell is genuinely generic (no even/odd constraint links
    // +t to -t), asserted O(1) from below.
    [Fact]
    public void The_Zero_Cell_Vanishes_Identically_And_The_Generic_Cell_Does_Not()
    {
        var rep = new OrderSorting(W, 5).Run();
        Assert.True(rep.ZeroCellWorst < 1e-12, $"zero cell {rep.ZeroCellWorst:E1}");
        Assert.True(rep.GenericGap > 1e-3, $"generic cell must respond, gap {rep.GenericGap:E1}");
    }

    // the leak (hypothesis violation): prep = R-even + eps * (R-odd) opens the forbidden odd-in-t
    // channel; the opening is EXACTLY affine in eps (the master equation and RK4 are linear in
    // rho0), so halving eps halves the slope to machine precision.
    [Fact]
    public void The_Leak_Opens_Exactly_Affine_In_Eps()
    {
        var rep = new OrderSorting(W, 5).Run();
        Assert.True(Math.Abs(rep.LeakSlope) > 1e-6, $"the leak must exist, slope {rep.LeakSlope:E1}");
        Assert.True(Math.Abs(rep.LeakRatio - 2.0) < 1e-6, $"eps-halving ratio {rep.LeakRatio:F7}");
    }

    // T0 rides inside the report too (the gamma-axis conjugation residual the twin runs rest on).
    [Fact]
    public void The_Report_Carries_The_Conjugation_Residual()
    {
        var rep = new OrderSorting(W, 4).Run();
        Assert.True(rep.ConjugationResidual < 1e-12, $"T0 residual {rep.ConjugationResidual:E1}");
    }
}
