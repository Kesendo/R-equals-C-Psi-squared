using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the atmosphere (adopted 2026-10-11 from
// experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md, typed in the main repo as
// CorrelatedZDephasingPriceClaim): correlated Z dephasing as the Mirror's price turned into an
// operator. Integer inputs; every residual is compared with == 0.0 and every count exactly.
public class AtmosphereTests
{
    static readonly World W = new();

    static Atmosphere Make(int n, long gamma, long[,] cov, double j = 1.0)
        => new(new Mirror(W, n, j, gamma), cov);

    // a deterministic integer pattern, every covariance nonzero, entries in {-3..3}
    static long[,] Pattern(int n, int seed)
    {
        var rng = new Random(seed);
        var m = new long[n, n];
        for (int k = 0; k < n; k++)
            for (int l = k + 1; l < n; l++)
            {
                long v = 0;
                while (v == 0) v = rng.Next(-3, 4);
                m[k, l] = m[l, k] = v;
            }
        return m;
    }

    // a fixed integer pattern with every covariance positive, C_kl = 1 + (k + l) mod 3 (row sums unequal at N = 3 and 4)
    static long[,] Formula(int n)
    {
        var m = new long[n, n];
        for (int k = 0; k < n; k++)
            for (int l = 0; l < n; l++)
                if (k != l) m[k, l] = 1 + (k + l) % 3;
        return m;
    }

    static IEnumerable<(int P, int Q)> Blocks(int n)
    {
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
                yield return (p, q);
    }

    static long Binomial(int n, int k)
    {
        long r = 1;
        for (int i = 1; i <= k; i++) r = r * (n - k + i) / i;
        return r;
    }

    // ---- the door ----

    [Fact]
    public void The_Mirrors_Gamma_Must_Be_An_Integer()
        => Assert.Throws<ArgumentException>(() => new Atmosphere(new Mirror(W, 3, 1.0, 0.5), new long[3, 3]));

    [Fact]
    public void The_Mirrors_J_Must_Be_An_Integer()
        => Assert.Throws<ArgumentException>(() => new Atmosphere(new Mirror(W, 3, 0.9, 10), new long[3, 3]));

    [Fact]
    public void The_Covariance_Must_Be_Symmetric()
    {
        var m = new long[3, 3];
        m[0, 1] = 1;
        Assert.Throws<ArgumentException>(() => Make(3, 10, m));
    }

    [Fact]
    public void The_Covariance_Has_No_Diagonal()
    {
        var m = new long[3, 3];
        m[1, 1] = 2;
        Assert.Throws<ArgumentException>(() => Make(3, 10, m));
    }

    [Fact]
    public void The_Covariance_Has_The_Mirrors_Size()
        => Assert.Throws<ArgumentException>(() => Make(3, 10, new long[4, 4]));

    [Fact]
    public void The_Covariance_Is_Bounded()
    {
        var m = new long[3, 3];
        m[0, 1] = m[1, 0] = 2_000_000_000L;
        Assert.Throws<ArgumentOutOfRangeException>(() => Make(3, 10, m));
    }

    [Fact]
    public void The_Cell_Sums_Stop_At_Eight_Sites()
        => Assert.Throws<ArgumentOutOfRangeException>(() => Make(9, 10, new long[9, 9]));

    // ---- the ontology: the variance is the mirror's, the covariance the atmosphere's ----

    [Fact]
    public void The_Variance_Is_Inherited_The_Covariance_Owned()
    {
        var atm = Make(4, 10, Atmosphere.Uniform(4, 3));
        Assert.Equal(new[] { "rate", "ising", "footprint", "centres", "dark" }, atm.Own);
        Assert.IsType<Mirror>(atm.Parent);
        Assert.Equal(atm.Mirror.Gamma, (double)atm.Gamma);   // the variance is read from the Mirror, not owned
        foreach (var name in new[] { "legs", "price", "orbit", "x", "y", "z" })
            Assert.Contains(name, atm.Inherited);
    }

    // ---- a featureless atmosphere is the Mirror, entry by entry ----

    [Theory, InlineData(2), InlineData(3), InlineData(4), InlineData(5)]
    public void A_Featureless_Atmosphere_Is_The_Mirror(int n)
    {
        var atm = Make(n, 10, new long[n, n]);
        foreach (var (p, q) in Blocks(n))
        {
            var mine = atm.BuildBlock(p, q, out _, out _);
            var mirrors = atm.Mirror.BuildBlock(p, q, out _, out _);
            for (int i = 0; i < mine.GetLength(0); i++)
                for (int j = 0; j < mine.GetLength(1); j++)
                    Assert.True(mine[i, j] == mirrors[i, j], $"N={n} block ({p},{q}) entry ({i},{j})");
        }
        for (int i = 0; i < (1 << n); i++)
            for (int j = 0; j < (1 << n); j++)
                Assert.Equal(atm.Mirror.Price, (double)atm.FoldPrice(i, j));
        Assert.Equal(0, atm.FootprintByCells());
        Assert.Equal(Enumerable.Range(0, n + 1).ToArray(), atm.QualifyingPopcounts());
        Assert.Equal(1 << n, atm.DarkCount());
    }

    // ---- rate: what a cell pays ----

    [Theory, InlineData(2, 1), InlineData(3, 2), InlineData(4, 3), InlineData(5, 4), InlineData(6, 5)]
    public void The_Rate_Read_Two_Ways_Agrees(int n, int seed)
    {
        var atm = Make(n, 10, Pattern(n, seed));
        for (int i = 0; i < (1 << n); i++)
            for (int j = 0; j < (1 << n); j++)
                Assert.Equal(atm.RateFromKossakowski(i, j), atm.Rate(i, j));
    }

    [Theory, InlineData(3), InlineData(4), InlineData(5)]
    public void The_Single_Disagreement_Blocks_Run_As_Under_Local_Noise_And_The_Image_Of_Zero_One_Does_Not(int n)
    {
        var atm = Make(n, 10, Formula(n));
        foreach (var (p, q) in new[] { (0, 1), (1, 0), (n - 1, n), (n, n - 1) })
        {
            var mine = atm.BuildBlock(p, q, out _, out _);
            var mirrors = atm.Mirror.BuildBlock(p, q, out _, out _);
            for (int i = 0; i < mine.GetLength(0); i++)
                for (int j = 0; j < mine.GetLength(1); j++)
                    Assert.True(mine[i, j] == mirrors[i, j], $"N={n} single-disagreement block ({p},{q}) entry ({i},{j})");
        }
        var image = atm.BuildBlock(1, n, out _, out _);          // the Pi_Z image of (0, 1)
        var local = atm.Mirror.BuildBlock(1, n, out _, out _);
        bool differs = false;
        for (int i = 0; i < image.GetLength(0); i++)
            differs |= image[i, i] != local[i, i];
        Assert.True(differs, $"N={n}: the image (1, N) of a single-disagreement block carries the covariance");
    }

    // the dressed block's diagonal is the cell rate itself, read both ways (the legs alone would not catch a
    // dressing applied to the flipped bra: the fold identity holds for it too)
    [Theory, InlineData(3), InlineData(4), InlineData(5)]
    public void The_Dressed_Diagonal_Is_The_Cell_Rate(int n)
    {
        var atm = Make(n, 10, Formula(n));
        foreach (var (p, q) in Blocks(n))
        {
            var block = atm.BuildBlock(p, q, out var kets, out var bras);
            for (int ka = 0; ka < kets.Count; ka++)
                for (int bb = 0; bb < bras.Count; bb++)
                {
                    int cell = ka * bras.Count + bb;
                    long rate = atm.Rate(kets[ka], bras[bb]);
                    Assert.True(block[cell, cell] == new Complex(rate, 0), $"N={n} block ({p},{q}) cell {cell}");
                    Assert.Equal(atm.RateFromKossakowski(kets[ka], bras[bb]), rate);
                }
        }
    }

    // ---- ising: the folds at the cell's own price ----

    [Theory, InlineData(2, 21), InlineData(3, 22), InlineData(4, 23), InlineData(5, 24), InlineData(6, 25)]
    public void Every_Leg_Holds_With_The_Cell_Price(int n, int seed)
    {
        foreach (var cov in new[] { Pattern(n, seed), Atmosphere.Uniform(n, 3), Atmosphere.Uniform(n, -2) })
        {
            var atm = Make(n, 10, cov, j: 2.0);
            foreach (var (p, q) in Blocks(n))
            {
                double bra = atm.BraFoldResidual(p, q), ket = atm.KetFoldResidual(p, q);
                Assert.True(bra == 0.0, $"N={n} bra fold ({p},{q}): {bra}");
                Assert.True(ket == 0.0, $"N={n} ket fold ({p},{q}): {ket}");
                Assert.True(atm.TransposeResidual(p, q) == 0.0, $"N={n} transpose ({p},{q})");
                Assert.True(atm.KleinResidual(p, q) == 0.0, $"N={n} Klein ({p},{q})");
            }
        }
    }

    [Fact]
    public void The_Mirrors_Scalar_Price_Fails_Where_The_Covariance_Speaks()
    {
        var atm = Make(4, 10, Formula(4));                 // row sums 6, 5, 7, 6: K is not one number on popcount 1
        Assert.True(atm.BraFoldResidualAtScalarPrice(1, 1, (long)atm.Mirror.Price) > 0.0,
            "the anti-vacuity control: the one number 2*N*gamma must NOT fold a block the covariance reads");
    }

    [Theory, InlineData(2, 31), InlineData(3, 32), InlineData(4, 33), InlineData(5, 34)]
    public void The_Mirrors_Price_Is_The_Atmospheres_Mean(int n, int seed)
    {
        var atm = Make(n, 10, Pattern(n, seed));
        Assert.Equal((1L << n) * n * 10, atm.IsingSum());       // mean K = N*gamma, mean fold price 2*N*gamma
    }

    // ---- footprint: what the mirror no longer cancels ----

    [Theory, InlineData(2, 41), InlineData(3, 42), InlineData(4, 43), InlineData(5, 44), InlineData(6, 45)]
    public void The_Footprint_Read_Two_Ways_Agrees(int n, int seed)
    {
        var atm = Make(n, 10, Pattern(n, seed));
        Assert.Equal(atm.FootprintClosedForm(), atm.FootprintByCells());
        Assert.True(atm.FootprintByCells() > 0);
    }

    [Fact]
    public void The_Shared_Field_Footprint_At_N3_Is_The_Closed_Form()
        => Assert.Equal(13824, Make(3, 10, Atmosphere.Uniform(3, 3)).FootprintByCells());   // 8 * 4^3 * 3 * 3^2

    [Fact]
    public void A_Featureless_Atmosphere_Leaves_No_Footprint()
        => Assert.Equal(0, Make(4, 10, new long[4, 4]).FootprintByCells());

    // ---- centres: where the price is one number per block ----

    [Theory, InlineData(2), InlineData(3), InlineData(4), InlineData(5)]
    public void Under_One_Common_Covariance_Every_Block_Folds_At_One_Number(int n)
    {
        const long gamma = 10, kappa = 3;
        var atm = Make(n, gamma, Atmosphere.Uniform(n, kappa));
        Assert.Equal(Enumerable.Range(0, n + 1).ToArray(), atm.QualifyingPopcounts());
        foreach (var (p, q) in Blocks(n))
        {
            Assert.Equal(2 * Atmosphere.SharedBlockSigma(n, gamma, kappa, p, q), atm.BlockPrice(p, q));
            Assert.True(atm.BraFoldResidualAtScalarPrice(p, q, atm.BlockPrice(p, q)) == 0.0, $"N={n} block ({p},{q})");
        }
    }

    [Fact]
    public void The_Pages_Two_Centres()
    {
        Assert.Equal(28, Atmosphere.SharedBlockSigma(4, 10, 3, 2, 2));     // 2.8 gamma, the half-filled block
        Assert.Equal(76, Atmosphere.SharedBlockSigma(4, 10, 3, 0, 4));     // 7.6 gamma, the (0, 4) block
        Assert.Equal(40, Atmosphere.SharedBlockSigma(4, 10, 0, 1, 3));     // featureless: sigma on every block
    }

    [Fact]
    public void At_Two_Sites_Every_Pattern_Is_Shared()
        => Assert.Equal(new[] { 0, 1, 2 }, Make(2, 10, Pattern(2, 51)).QualifyingPopcounts());

    [Fact]
    public void Two_Patterns_With_Unequal_Row_Sums_Qualify_Only_The_Empty_And_The_Full_Configuration()
    {
        // N = 3: row sums 3, 4, 5. N = 4: row sums 6, 2, 2, 6, and popcount 2 sees two values of K.
        var c3 = new long[3, 3] { { 0, 1, 2 }, { 1, 0, 3 }, { 2, 3, 0 } };
        var c4 = new long[4, 4] { { 0, 1, 2, 3 }, { 1, 0, -1, 2 }, { 2, -1, 0, 1 }, { 3, 2, 1, 0 } };
        Assert.Equal(new[] { 0, 3 }, Make(3, 10, c3).QualifyingPopcounts());
        Assert.Equal(new[] { 0, 4 }, Make(4, 10, c4).QualifyingPopcounts());
        Assert.Throws<InvalidOperationException>(() => Make(3, 10, c3).BlockPrice(1, 1));
    }

    // the source page's sign pattern 3*s_k*s_l, s = (1, 1, 1, -1): row sums 3, 3, 3, -9, yet popcount 2 sees one
    // value of K, so row sums decide only popcounts 1 and N - 1
    [Fact]
    public void A_Sign_Pattern_Qualifies_Popcount_Two_Despite_Unequal_Row_Sums()
    {
        var s = new long[] { 1, 1, 1, -1 };
        var c = new long[4, 4];
        for (int k = 0; k < 4; k++)
            for (int l = 0; l < 4; l++)
                if (k != l) c[k, l] = 3 * s[k] * s[l];
        var atm = Make(4, 10, c);
        Assert.Equal(new[] { 0, 2, 4 }, atm.QualifyingPopcounts());
        AssertQualifyingBlocksFoldAtOneNumber(atm, 9);
    }

    // two disjoint pairs, C_01 = C_23 = 3: every row sum 3, so popcounts 1 and 3 qualify, popcount 2 does not
    [Fact]
    public void Two_Disjoint_Pairs_Qualify_Sixteen_Blocks()
    {
        var c = new long[4, 4];
        c[0, 1] = c[1, 0] = c[2, 3] = c[3, 2] = 3;
        var atm = Make(4, 10, c);
        Assert.Equal(new[] { 0, 1, 3, 4 }, atm.QualifyingPopcounts());
        AssertQualifyingBlocksFoldAtOneNumber(atm, 16);
    }

    static void AssertQualifyingBlocksFoldAtOneNumber(Atmosphere atm, int expectedBlocks)
    {
        var ok = atm.QualifyingPopcounts();
        int count = 0;
        foreach (int p in ok)
            foreach (int q in ok)
            {
                double r = atm.BraFoldResidualAtScalarPrice(p, q, atm.BlockPrice(p, q));
                Assert.True(r == 0.0, $"block ({p},{q}) at its one price: {r}");
                count++;
            }
        Assert.Equal(expectedBlocks, count);
    }

    // ---- dark: what the atmosphere does not see ----

    [Theory, InlineData(2), InlineData(3), InlineData(4), InlineData(5), InlineData(6)]
    public void What_The_Atmosphere_Does_Not_See(int n)
    {
        int d = 1 << n;
        Assert.Equal(d, Make(n, 10, new long[n, n]).DarkCount());                       // featureless: the populations
        Assert.Equal(d, Make(n, 10, Atmosphere.Uniform(n, 3)).DarkCount());             // eta = 3/10, positive definite
        Assert.Equal(Binomial(2 * n, n), Make(n, 10, Atmosphere.Uniform(n, 10)).DarkCount());   // eta = 1: C(2N, N)
        var antiEnd = Make(n, n - 1, Atmosphere.Uniform(n, -1));                       // eta = -1/(N - 1): rows sum to 0
        Assert.Equal(d + 2, antiEnd.DarkCount());                                       // the populations and the GHZ pair
        Assert.Equal(0, antiEnd.Rate(0, d - 1));
        Assert.Equal(0, antiEnd.Rate(d - 1, 0));
    }
}
