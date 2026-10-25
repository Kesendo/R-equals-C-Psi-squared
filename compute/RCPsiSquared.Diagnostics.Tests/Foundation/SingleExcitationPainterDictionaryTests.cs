using System;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.States;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The α-dictionary in the single-excitation block against the painters' 4^N pipeline, and what it
/// reads past N = 6. Two error models: (1) the α fit stops at a golden-section width of 1e-7 in α, so a letter
/// read twice from trajectories that agree to rounding differs by at most 1e-7/δJ = 5e-6 in f; (2) the fit's
/// 512-point grid seed picks a basin by argmin, and a near-tie between seeds can flip the basin and move f by
/// about 1. The handshake gate 1e-4 sits between the two, so it passes the first and catches the second.</summary>
public class SingleExcitationPainterDictionaryTests
{
    private const double BasinFlipGate = 1e-4;
    // two fits of mirror-related inputs each stop within 1e-7 of their minimiser in α, so their letters
    // differ by at most 2·1e-7/δJ = 1e-5 in f where the misfit has a well-defined minimum
    private const double TwoStops = 2e-7 / SingleExcitationPainterDictionary.DefaultDeltaJ;

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    public void Letters_AgreeWithTheFullPipeline(int n)
    {
        var block = new SingleExcitationPainterDictionary(n);
        for (int b = 0; b < n - 1; b++)
        {
            var s = new Symphony(n: n, j: 1.0, gamma: 0.05, hType: HamiltonianType.XY,
                                 initialState: InitialStateKind.BondingMode, defectBond: b, deltaJ: 0.02);
            var pm = ((IInspectable)s).Children.OfType<PaintersMovement>().Single();
            for (int i = 0; i < n; i++)
            {
                Assert.True(Math.Abs(pm.F[i] - block.F[b][i]) < BasinFlipGate,
                    $"N={n} bond {b} site {i}: pipeline f = {pm.F[i]}, block f = {block.F[b][i]}");
                Assert.Equal(pm.Reliable[i], block.Reliable[b][i]);
            }
        }
    }

    /// <summary>The premise of the chain-mirror identity: the reliability flags are mirror-consistent exactly,
    /// and on the reliable sites the letter of bond N − 2 − b is R applied to the letter of bond b to two
    /// stopping widths (on a rejected site the misfit is flat and rounding moves its minimiser, which is why
    /// the painters reject it). At even N the centre bond is its own mirror, so its letter is mirror-even to the
    /// same bound.</summary>
    [Theory]
    [InlineData(6)]
    [InlineData(7)]
    [InlineData(8)]
    [InlineData(9)]
    [InlineData(10)]
    [InlineData(11)]
    [InlineData(12)]
    public void MirrorBondsLetter_IsTheMirrorImage(int n)
    {
        var d = new SingleExcitationPainterDictionary(n);
        for (int b = 0; b < n - 1; b++)
            for (int i = 0; i < n; i++)
            {
                Assert.Equal(d.Reliable[b][n - 1 - i], d.Reliable[n - 2 - b][i]);
                if (!d.Reliable[n - 2 - b][i]) continue;
                Assert.True(Math.Abs(d.F[n - 2 - b][i] - d.F[b][n - 1 - i]) <= TwoStops,
                    $"N={n} bond {b} site {i}: {d.F[n - 2 - b][i]} against {d.F[b][n - 1 - i]}");
            }
    }

    /// <summary>The identity itself is algebra: for any f and its mirror image Rf, cos(f, Rf) = 1 − 2w(f).
    /// Checked on each letter and its exact mirror image, to machine rounding.</summary>
    [Theory]
    [InlineData(7)]
    [InlineData(12)]
    public void CosOfALetterWithItsMirrorImage_IsOneMinusTwiceTheOddWeight(int n)
    {
        var d = new SingleExcitationPainterDictionary(n);
        for (int b = 0; b < n - 1; b++)
        {
            var f = d.F[b];
            var rf = Enumerable.Range(0, n).Select(i => f[n - 1 - i]).ToArray();
            double cos = SingleExcitationPainterDictionary.Cos(f, rf);
            Assert.True(Math.Abs(cos - (1 - 2 * d.MirrorOddWeight(b))) < 1e-12,
                $"N={n} bond {b}: cos = {cos}, 1 − 2w = {1 - 2 * d.MirrorOddWeight(b)}");
        }
    }

    /// <summary>The reading over the sites reliable for every bond, pinned to the printed three decimals
    /// (the fit's width moves a cosine by far less than the 5e-4 of rounding), with the class of the worst
    /// pair: a chain-mirror pair (a + b = N − 2) or a pair across the centre (a &lt; (N − 2)/2 &lt; b).</summary>
    [Theory]
    [InlineData(3, -0.976, true)]
    [InlineData(5, -0.965, false)]
    [InlineData(6, -0.962, true)]
    [InlineData(7, -0.940, false)]
    [InlineData(8, -0.927, false)]
    [InlineData(9, -0.992, true)]
    [InlineData(10, -0.950, false)]
    [InlineData(11, -0.914, true)]
    [InlineData(12, -0.961, true)]
    public void OverReliableSites_TheWorstPair(int n, double expected, bool mirrorPair)
    {
        var w = new SingleExcitationPainterDictionary(n).WorstAntiCollinearPair(reliableOnly: true);
        Assert.True(Math.Abs(w.Cos - expected) < 5e-4, $"N={n}: {w.Cos:F4} at ({w.A},{w.B})");
        if (mirrorPair) Assert.Equal(n - 2, w.A + w.B);
        else Assert.True(w.A < (n - 2) / 2.0 && (n - 2) / 2.0 < w.B, $"N={n}: ({w.A},{w.B}) is not across the centre");
    }

    /// <summary>N = 4: every site reliable, the one chain-mirror pair (0,2) at w = 0.707, cos = −0.414, and the
    /// worst pair (0,1) at −0.541, an edge bond with the mirror-even centre (neither class above).</summary>
    [Fact]
    public void AtNFour_TheThreeBondInventory()
    {
        var d = new SingleExcitationPainterDictionary(4);
        Assert.All(d.SitesReliableForEveryBond(), Assert.True);
        Assert.True(Math.Abs(d.MirrorOddWeight(0) - 0.707) < 5e-4, $"w(0) = {d.MirrorOddWeight(0)}");
        Assert.True(Math.Abs(SingleExcitationPainterDictionary.Cos(d.F[0], d.F[2]) - (-0.414)) < 5e-4);
        var w = d.WorstAntiCollinearPair(reliableOnly: true);
        Assert.True(Math.Abs(w.Cos - (-0.541)) < 5e-4 && (w.A, w.B) == (0, 1), $"{w.Cos:F4} at ({w.A},{w.B})");
    }

    /// <summary>Over all sites N = 6 reads −0.378 and N = 10 −0.701; the sites the painters' rule rejects carry
    /// those values (two of six at N = 6, with |f| up to 216; two of ten at N = 10, |f| about 42).</summary>
    [Theory]
    [InlineData(6, -0.378, 4, 215.6)]
    [InlineData(10, -0.701, 8, 41.8)]
    public void OverAllSites_TheWeakReadingIsCarriedByRejectedSites(int n, double allSites, int reliableSites, double maxF)
    {
        var d = new SingleExcitationPainterDictionary(n);
        Assert.True(Math.Abs(d.WorstAntiCollinearPair(reliableOnly: false).Cos - allSites) < 5e-4);
        var keep = d.SitesReliableForEveryBond();
        Assert.Equal(reliableSites, keep.Count(x => x));
        double maxRejected = Enumerable.Range(0, n).Where(i => !keep[i])
            .Max(i => Enumerable.Range(0, n - 1).Max(b => Math.Abs(d.F[b][i])));
        Assert.True(Math.Abs(maxRejected - maxF) < 0.05, $"N={n}: max |f| on rejected sites {maxRejected:F2}");
    }

    /// <summary>N = 5: the painted chain-mirror pairs read −0.783 and −0.894; the worst pair, −0.965, is
    /// (0,2)/(1,3), across the centre.</summary>
    [Fact]
    public void AtNFive_TheMirrorPairsAreNotTheWorstPair()
    {
        var d = new SingleExcitationPainterDictionary(5);
        Assert.True(Math.Abs(SingleExcitationPainterDictionary.Cos(d.F[0], d.F[3]) - (-0.783)) < 5e-4);
        Assert.True(Math.Abs(SingleExcitationPainterDictionary.Cos(d.F[1], d.F[2]) - (-0.894)) < 5e-4);
        Assert.True(Math.Abs(SingleExcitationPainterDictionary.Cos(d.F[1], d.F[3]) - (-0.965)) < 5e-4);
    }

    [Fact]
    public void AtNTwelve_BondThreeIsNearlyMirrorOdd() =>
        Assert.True(Math.Abs(new SingleExcitationPainterDictionary(12).MirrorOddWeight(3) - 0.980) < 5e-4);

    /// <summary>The α that stretches the EXACT clean purity onto P_B: the root of the misfit's stationarity
    /// h(α) = Σ_t t·Ṗ_A(αt)·(P_A(αt) − P_B(t)), bisected to the rounding of h (a root, unlike a minimum, is
    /// located to ε rather than √ε).</summary>
    private static double ExactStretch(SingleExcitationPainterDictionary d, int site, double[] pB, double lo, double hi)
    {
        double H(double alpha)
        {
            var times = d.TimeGrid.Select(t => alpha * t).ToArray();
            var (p, rate) = d.CleanPopulations(times);
            double h = 0;
            for (int s = 0; s < times.Length; s++)
            {
                double x = 1 - 2 * p[site][s];
                h += d.TimeGrid[s] * (-2 * x * rate[site][s]) * (0.5 * (1 + x * x) - pB[s]);
            }
            return h;
        }
        double hl = H(lo);
        Assert.True(hl * H(hi) < 0, $"site {site}: no sign change on [{lo}, {hi}]");
        for (int it = 0; it < 200 && hi - lo > 1e-15; it++)
        {
            double mid = 0.5 * (lo + hi), hm = H(mid);
            if (hm * hl > 0) { lo = mid; hl = hm; } else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    /// <summary>The linear letter is the δJ → 0 limit of the exactly stretched fit: f_fit = f_lin + a₂·δJ + O(δJ²)
    /// at every site, so e(δJ) = (f_fit − f_lin)/δJ converges and its successive differences halve as δJ halves
    /// (gated at ≤ 0.6 to leave room for the O(δJ) correction to the ratio). Two routes: the eigenbasis
    /// derivative on one side, three defected runs and a root on the other. The floor is the rounding of the
    /// trajectories: η, read off the populations' sum, which the exact dynamics holds at 1, taken as N·η per
    /// population and 2N·η per purity; a purity error η_P moves the root by at most η_P·Σ|g|/Σg², and e divides that
    /// by δJ², twice for the two subtractions.</summary>
    [Theory]
    [InlineData(5, 1)]
    [InlineData(8, 2)]
    [InlineData(12, 3)]
    public void LinearLetter_IsTheLimitOfTheExactlyStretchedFit(int n, int bond)
    {
        var d = new SingleExcitationPainterDictionary(n);
        var lin = d.LinearF[bond];
        double[] deltas = { 0.004, 0.002, 0.001 };
        double eta = 0;
        double TraceDefect(double[][] pop) =>
            Enumerable.Range(0, pop[0].Length).Max(s => Math.Abs(pop.Sum(row => row[s]) - 1));
        var (pClean, rate) = d.CleanPopulations(d.TimeGrid);
        eta = Math.Max(eta, TraceDefect(pClean));
        var e = deltas.Select(dj =>
        {
            var pop = d.DefectPopulations(bond, dj);
            eta = Math.Max(eta, TraceDefect(pop));
            var pB = pop.Select(row => row.Select(x => 0.5 * (1 + (1 - 2 * x) * (1 - 2 * x))).ToArray()).ToArray();
            return Enumerable.Range(0, n).Select(i =>
                (ExactStretch(d, i, pB[i], 1 - 20 * dj, 1 + 20 * dj) - 1) / dj - lin[i]).Select(x => x / dj).ToArray();
        }).ToArray();
        for (int i = 0; i < n; i++)
        {
            var g = Enumerable.Range(0, d.TimeGrid.Count)
                .Select(s => d.TimeGrid[s] * -2 * (1 - 2 * pClean[i][s]) * rate[i][s]).ToArray();
            double floor = 2 * (2 * n * eta) * g.Sum(Math.Abs) / g.Sum(x => x * x) / (deltas[^1] * deltas[^1]);
            double first = Math.Abs(e[0][i] - e[1][i]), second = Math.Abs(e[1][i] - e[2][i]);
            Assert.True(second <= 0.6 * first + floor,
                $"N={n} bond {bond} site {i}: e = {e[0][i]}, {e[1][i]}, {e[2][i]}, floor {floor}");
            // the floor may not carry the gate: with floor ≤ 0.3·first the gate still demands second ≤ 0.9·first
            Assert.True(floor <= 0.3 * first, $"N={n} site {i}: floor {floor} against {first}");
        }
    }

    /// <summary>The carrier coupling's closed form against the mode matrix element ψ_k(b)ψ₁(b+1) + ψ_k(b+1)ψ₁(b):
    /// both reduce their integer angles into [0, 2π) before the sine or cosine, so each side is a few products of
    /// values ≤ 1 rounded to a few ε at every N (8ε gated). At k = N the closed form vanishes, the K-partner rule.</summary>
    [Fact]
    public void CarrierCoupling_ClosedFormIsTheModeMatrixElement()
    {
        for (int n = 3; n <= 16; n++)
        {
            var psi1 = SingleExcitationPainterDictionary.Mode(n, 1);
            for (int k = 1; k <= n; k++)
            {
                var psik = SingleExcitationPainterDictionary.Mode(n, k);
                for (int b = 0; b < n - 1; b++)
                {
                    double element = psik[b] * psi1[b + 1] + psik[b + 1] * psi1[b];
                    double closed = SingleExcitationPainterDictionary.CarrierCoupling(n, k, b);
                    Assert.True(Math.Abs(element - closed) <= 8 * Math.Pow(2, -52),
                        $"N={n} k={k} b={b}: {element} against {closed}");
                    if (k == n) Assert.True(Math.Abs(closed) <= 8 * Math.Pow(2, -52), $"N={n} b={b}: c_N = {closed}");
                }
            }
        }
    }

    private static double[] CarrierModel(SingleExcitationPainterDictionary d, int bond)
    {
        var f = new double[d.N];
        for (int k = 1; k <= d.N; k++)
        {
            double c = SingleExcitationPainterDictionary.CarrierCoupling(d.N, k, bond);
            var u = d.CarrierChannelLetter(k);
            for (int i = 0; i < d.N; i++) f[i] += c * u[i];
        }
        return f;
    }

    private static double OddWeight(double[] f)
    {
        double odd = 0, all = 0;
        for (int i = 0; i < f.Length; i++) { double o = (f[i] - f[^(i + 1)]) / 2; odd += o * o; all += f[i] * f[i]; }
        return odd / all;
    }

    private static double Norm(double[] f) => Math.Sqrt(f.Sum(x => x * x));

    /// <summary>The carrier channels (1, k) carry most of each letter: the remainder ‖F − Σ_k c_k·u_k‖/‖F‖ over the
    /// whole dictionary, pinned to six decimals, 8 % to 17 % at N = 4..12 (the rest is the mode pairs that do not
    /// touch ψ₁).</summary>
    [Theory]
    [InlineData(4, 0.168539)]
    [InlineData(5, 0.108456)]
    [InlineData(6, 0.125317)]
    [InlineData(7, 0.084931)]
    [InlineData(8, 0.080397)]
    [InlineData(9, 0.078617)]
    [InlineData(10, 0.089745)]
    [InlineData(11, 0.088618)]
    [InlineData(12, 0.098031)]
    public void CarrierChannels_CarryTheLetters(int n, double residual)
    {
        var d = new SingleExcitationPainterDictionary(n);
        double r2 = 0, f2 = 0;
        for (int b = 0; b < n - 1; b++)
        {
            var m = CarrierModel(d, b);
            for (int i = 0; i < n; i++) { r2 += Math.Pow(d.LinearF[b][i] - m[i], 2); f2 += Math.Pow(d.LinearF[b][i], 2); }
        }
        Assert.True(Math.Abs(Math.Sqrt(r2 / f2) - residual) < 5e-6, $"N={n}: residual {Math.Sqrt(r2 / f2):F6}");
    }

    /// <summary>The carrier model's mirror-odd weight against the linear letter's, the largest departure over the
    /// bonds: within 0.06 at N = 5..10, within 0.1 at N = 4 and N = 11..13 (pinned to four decimals).</summary>
    [Theory]
    [InlineData(4, 0.0922)]
    [InlineData(5, 0.0226)]
    [InlineData(6, 0.0426)]
    [InlineData(7, 0.0266)]
    [InlineData(8, 0.0282)]
    [InlineData(9, 0.0238)]
    [InlineData(10, 0.0512)]
    [InlineData(11, 0.0684)]
    [InlineData(12, 0.0950)]
    [InlineData(13, 0.0922)]
    public void CarrierModel_ReadsTheMirrorOddWeight(int n, double maxDeparture)
    {
        var d = new SingleExcitationPainterDictionary(n);
        double worst = Enumerable.Range(0, n - 1)
            .Max(b => Math.Abs(OddWeight(d.LinearF[b]) - OddWeight(CarrierModel(d, b))));
        Assert.True(Math.Abs(worst - maxDeparture) < 5e-5, $"N={n}: {worst:F5}");
    }

    /// <summary>The halves V± = (V_b ± V_{N−2−b})/2 of nearly pure seesaws. N = 6 bond 1: c₃ near its node, the even
    /// half 0.339 against 2.883. N = 12 bond 3: c₃ near its node, the even half 2.650 against
    /// 12.548. N = 8 bonds 2 and 1 share |c₃| = 1/6; bond 1's even half is the larger, 3.198 against 2.246, and its w
    /// 0.866 against 0.913. Pinned to five decimals, with the carrier model's w beside the linear one.</summary>
    [Theory]
    [InlineData(6, 1, -0.05379, 0.33893, 2.88313, 0.98637, 0.98916)]
    [InlineData(8, 2, -0.16667, 2.24587, 7.25454, 0.91254, 0.89941)]
    [InlineData(8, 1, 0.16667, 3.19768, 8.14439, 0.86644, 0.88524)]
    [InlineData(12, 3, -0.02441, 2.65035, 12.54802, 0.95729, 0.96534)]
    public void NearlyPureSeesaw_TheHalves(int n, int bond, double c3, double even, double odd, double w, double wCarrier)
    {
        var d = new SingleExcitationPainterDictionary(n);
        Assert.True(Math.Abs(SingleExcitationPainterDictionary.CarrierCoupling(n, 3, bond) - c3) < 5e-5);
        double[,] Half(int sign)
        {
            var v = SingleExcitationPainterDictionary.BondOperator(n, bond);
            var m = SingleExcitationPainterDictionary.BondOperator(n, n - 2 - bond);
            var h = new double[n, n];
            for (int x = 0; x < n; x++) for (int y = 0; y < n; y++) h[x, y] = (v[x, y] + sign * m[x, y]) / 2;
            return h;
        }
        Assert.True(Math.Abs(Norm(d.LinearLetter(Half(+1))) - even) < 5e-5);
        Assert.True(Math.Abs(Norm(d.LinearLetter(Half(-1))) - odd) < 5e-5);
        Assert.True(Math.Abs(OddWeight(d.LinearF[bond]) - w) < 5e-5);
        Assert.True(Math.Abs(OddWeight(CarrierModel(d, bond)) - wCarrier) < 5e-5);
    }

    /// <summary>N = 6: the two sites the painters reject at δJ = 0.02 (|f| 216 and 42, a far minimum of the α fit
    /// at that δJ, which the exactly stretched purity has too) carry ordinary linear letters, 1.0471 and −0.9726
    /// on bond 1.</summary>
    [Fact]
    public void AtNSix_TheRejectedSitesHaveOrdinaryLinearLetters()
    {
        var d = new SingleExcitationPainterDictionary(6);
        Assert.False(d.Reliable[1][4]);
        Assert.True(Math.Abs(d.F[1][4] - 215.64) < 0.005 && Math.Abs(d.F[1][1] - (-42.27)) < 0.005);
        Assert.True(Math.Abs(d.LinearF[1][4] - 1.04712) < 5e-5 && Math.Abs(d.LinearF[1][1] - (-0.97263)) < 5e-5);
    }

    /// <summary>The linear dictionary over all sites, no reliability mask needed: the worst signed cos is
    /// anti-collinear at every N = 3, 5..12 and −0.536 at N = 4.</summary>
    [Theory]
    [InlineData(3, -0.9646)]
    [InlineData(4, -0.5364)]
    [InlineData(5, -0.9826)]
    [InlineData(6, -0.9727)]
    [InlineData(7, -0.9281)]
    [InlineData(8, -0.9595)]
    [InlineData(9, -0.9675)]
    [InlineData(10, -0.9306)]
    [InlineData(11, -0.9049)]
    [InlineData(12, -0.9146)]
    public void LinearDictionary_WorstPairOverAllSites(int n, double expected)
    {
        var f = new SingleExcitationPainterDictionary(n).LinearF;
        double worst = double.PositiveInfinity;
        for (int a = 0; a < n - 1; a++)
            for (int b = a + 1; b < n - 1; b++) worst = Math.Min(worst, SingleExcitationPainterDictionary.Cos(f[a], f[b]));
        Assert.True(Math.Abs(worst - expected) < 5e-5, $"N={n}: {worst:F5}");
    }

    [Fact]
    public void Guard_RejectsOutOfRange()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new SingleExcitationPainterDictionary(2));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SingleExcitationPainterDictionary(25));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SingleExcitationPainterDictionary(5, deltaJ: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SingleExcitationPainterDictionary(5, timePoints: 1));
    }
}
