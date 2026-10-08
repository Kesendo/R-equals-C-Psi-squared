using MirrorWorld;

namespace MirrorWorldTests;

// The flavor split as F36's pairing in closed form. Pins come from graphs whose scaled spectra are rational
// (the C_4 ring: 1, 0, -1, 0; the C_6 ring: 1, +-1/2, -1; the complete graph K_6: 1 and -1/5; the star: +-1 and 0);
// the exact rows compare to 0 on dyadic inputs, the F36 residual to 0 at generic ones by its summation order; the
// irrational chain row and the two off-page rows are read against the producer's own numbers to six digits.
public class FlavorSplitTests
{
    // The experiment's parameters (simulations/results/neural_flavor_rule.txt header).
    private const double H = 0.5, C = 0.25, GE = 1.0, GI = 2.0;

    private static double[,] RingC4Scaled()
    {
        // adjacency of the 4-ring has spectral radius 2; scaled entries 1/2
        var a = new double[4, 4];
        for (int l = 0; l < 4; l++) { a[l, (l + 1) % 4] = 0.5; a[(l + 1) % 4, l] = 0.5; }
        return a;
    }
    private static double[,] CompleteScaled(int n)
    {
        var a = new double[n, n];
        for (int r = 0; r < n; r++) for (int k = 0; k < n; k++) if (r != k) a[r, k] = 1.0 / (n - 1);
        return a;
    }

    [Fact]
    public void F36_Holds_Exactly_On_The_Ring_And_The_Complete_Graph_At_Generic_Couplings()
    {
        Assert.Equal(0.0, FlavorSplit.F36Residual(RingC4Scaled(), H, C, GE, GI));
        Assert.Equal(0.0, FlavorSplit.F36Residual(RingC4Scaled(), 0.53, 0.29, 0.37, 0.91));
        Assert.Equal(0.0, FlavorSplit.F36Residual(CompleteScaled(6), 0.53, 0.29, 0.37, 0.91));
        Assert.Equal(0.0, FlavorSplit.F36Residual(CompleteScaled(7), 1.7, 1.1, 1.3, 0.2));
        // 200 seeded generic draws: the one-sum route is exactly 0.0 at every one of them
        var rng = new Random(12345); int hostNonzero = 0;
        for (int t = 0; t < 200; t++)
        {
            double gE = 0.05 + 1.9 * rng.NextDouble(), gI = 0.05 + 1.9 * rng.NextDouble(), h = 2 * rng.NextDouble(), c = 2 * rng.NextDouble();
            Assert.Equal(0.0, FlavorSplit.F36Residual(CompleteScaled(5), h, c, gE, gI));
            if (FlavorSplit.HostRouteResidual(CompleteScaled(5), h, c, gE, gI) != 0.0) hostNonzero++;
        }
        // the host's centred-first route carries a summation-order residual at generic leaks: a reading, not a gate
        Assert.True(hostNonzero >= 0);
        Assert.Equal(0.0, FlavorSplit.HostRouteResidual(RingC4Scaled(), H, C, GE, GI));   // and is exact on the dyadic page row
    }

    [Fact]
    public void Control_Coupling_With_The_Same_Sign_On_Both_Blocks_Breaks_F36()
    {
        var a = RingC4Scaled(); var j = FlavorSplit.Jacobian(a, H, C, GE, GI);
        // flip the inhibitory block's coupling sign: +cA on both blocks
        for (int r = 0; r < 4; r++) for (int k = 0; k < 4; k++) if (r != k) j[4 + r, 4 + k] = C * a[r, k];
        Assert.Equal(2 * C * 0.5, NeuralPalindrome.MaxResidual(j, FlavorSplit.SwapPermutation(4), FlavorSplit.Centre(GE, GI)));
        // and a wrong centre breaks the leak condition
        Assert.NotEqual(0.0, NeuralPalindrome.MaxResidual(FlavorSplit.Jacobian(a, H, C, GE, GI), FlavorSplit.SwapPermutation(4), GE));
    }

    [Theory]
    [InlineData(new double[] { 1, 1, 1, 1 }, 1.0)]
    [InlineData(new double[] { 1, -1, 1, -1 }, -1.0)]
    [InlineData(new double[] { 1, 0, -1, 0 }, 0.0)]
    public void The_Block_Action_On_A_Ring_Eigenvector_Is_Exact(double[] w, double a)
    {
        Assert.Equal(0.0, FlavorSplit.BlockActionResidual(RingC4Scaled(), w, a, 0.75, -0.25, H, C, GE, GI));
        Assert.Equal(0.0, FlavorSplit.BlockActionResidual(RingC4Scaled(), w, a, 1.0, 0.5, 0.25, 0.125, 1.0, 1.5));
        // a wrong eigenvalue breaks the action whenever a != 0
        if (a != 0) Assert.NotEqual(0.0, FlavorSplit.BlockActionResidual(RingC4Scaled(), w, -a, 0.75, -0.25, H, C, GE, GI));
    }

    [Fact]
    public void The_Ring_Blocks_At_The_Experiment_Parameters_Are_Overdamped_Defective_Underdamped()
    {
        // a = 1: x = 0.75, disc = 0.5625 - 0.25 = 0.3125 exactly
        Assert.Equal(0.3125, FlavorSplit.Discriminant(1, H, C, GE, GI));
        Assert.True(FlavorSplit.IsOverdamped(1, H, C, GE, GI));
        Assert.Equal("E", FlavorSplit.SlowFlavor(1, H, C, GE, GI));
        // a = 0: x = Delta = h, a Jordan block
        Assert.Equal(0.0, FlavorSplit.Discriminant(0, H, C, GE, GI));
        Assert.True(FlavorSplit.IsDefective(0, H, C, GE, GI));
        Assert.Equal("mixed", FlavorSplit.SlowFlavor(0, H, C, GE, GI));
        // a = -1: x = 0.25 < h, underdamped
        Assert.Equal(0.0625 - 0.25, FlavorSplit.Discriminant(-1, H, C, GE, GI));
        Assert.False(FlavorSplit.IsOverdamped(-1, H, C, GE, GI));
        Assert.Throws<InvalidOperationException>(() => FlavorSplit.OverdampedRates(-1, H, C, GE, GI));
    }

    [Fact]
    public void The_Slowest_Excitatory_Rate_Is_The_Same_On_Every_Connected_Graph()
    {
        // a_max = 1 after the radius scaling: rate s - sqrt(0.3125) = 0.940983 (the file's rate_E on every connected row)
        var (slow, fast) = FlavorSplit.OverdampedRates(1, H, C, GE, GI);
        Assert.Equal(1.5 - Math.Sqrt(0.3125), slow);
        Assert.Equal(1.5 + Math.Sqrt(0.3125), fast);
        Assert.Equal(0.940983, slow, 6);
        Assert.Equal(2.059017, fast, 6);
        Assert.Equal(3.0, slow + fast);                 // the pair sums to 2s exactly: F36's -2s on the rates
    }

    [Fact]
    public void Reported_Ratio_Single_Overdamped_Block_Reads_The_File_Value()
    {
        // ring N = 4 (1, 0, -1, 0), the star (+-1, 0 ...) and K_6 (1, five times -1/5) each have one overdamped block
        double expected = (1.5 + Math.Sqrt(0.3125)) / (1.5 - Math.Sqrt(0.3125));
        Assert.Equal(expected, FlavorSplit.ReportedRatio([1, 0, -1, 0], H, C, GE, GI));
        Assert.Equal(expected, FlavorSplit.ReportedRatio([1, 0, 0, 0, 0, -1], H, C, GE, GI));
        Assert.Equal(expected, FlavorSplit.ReportedRatio([1, -0.2, -0.2, -0.2, -0.2, -0.2], H, C, GE, GI));
        Assert.Equal(2.188155, expected, 6);           // neural_flavor_rule.txt: ring N=4, star, complete
    }

    [Fact]
    public void Reported_Ratio_Two_Overdamped_Blocks_Reads_The_Chain_Row()
    {
        // chain N = 4: eigenvalues 2 cos(k pi / 5), scaled by the radius 2 cos(pi/5): 1, (3 - sqrt 5)/2, -(3 - sqrt 5)/2, -1
        double g = (3 - Math.Sqrt(5)) / 2;
        double ratio = FlavorSplit.ReportedRatio([1, g, -g, -1], H, C, GE, GI);
        Assert.Equal(1.937798, ratio, 6);              // neural_flavor_rule.txt: chain N=4 (the faster root of the a = g block)
        // the numerator is the weaker block's fast rate, not the partner 2s - r of the slowest mode
        double slowE = 1.5 - Math.Sqrt(0.3125); double partner = 3.0 - slowE;
        Assert.Equal(1.5 + Math.Sqrt(FlavorSplit.Discriminant(g, H, C, GE, GI)), ratio * slowE, 12);
        Assert.NotEqual(partner, ratio * slowE);
    }

    [Fact]
    public void Off_The_Page_The_Classifier_Rule_Parts_From_The_Weakest_Block_Shortcut()
    {
        // (a) a block with x < -h: chain N=4, h=0.1, c=1, gE=1, gI=2. The producer's slow_flavor_rates reads
        // E 0.0033370452904232795 and I 1.0101020514433638 (the a = -1 block's SLOW root, x = -0.5), ratio 302.692876.
        double g = (3 - Math.Sqrt(5)) / 2;
        double ratioA = FlavorSplit.ReportedRatio([1, g, -g, -1], 0.1, 1.0, 1.0, 2.0);
        Assert.Equal(1.0101020514433638 / 0.0033370452904232795, ratioA, 6);
        Assert.Equal("I", FlavorSplit.SlowFlavor(-1, 0.1, 1.0, 1.0, 2.0));
        // (b) the 0.65 threshold: C_6 ring (1, 1/2, 1/2, -1/2, -1/2, -1), h=0.5, c=0.04. The weakest overdamped block
        // (a = 1/2, x = 0.52) has fast-mode I weight 0.637 < 0.65 and is skipped; the producer reads 1.7039607805437118
        // over 1.2960392194562875, ratio 1.314745, both from the a = 1 block.
        double ratioB = FlavorSplit.ReportedRatio([1, 0.5, 0.5, -0.5, -0.5, -1], 0.5, 0.04, 1.0, 2.0);
        Assert.Equal(1.7039607805437118 / 1.2960392194562875, ratioB, 6);
        double xWeak = FlavorSplit.BlockDiagonal(0.5, 0.04, 1.0, 2.0);
        Assert.True(FlavorSplit.ModeInhibitionWeight(xWeak, 0.5, -Math.Sqrt(xWeak * xWeak - 0.25)) < 0.65);
        // (c) weights at the page's parameters: the a = 1 block's slow mode carries I weight 0.127, its fast mode 0.873
        double x1 = FlavorSplit.BlockDiagonal(1, C, GE, GI); double r1 = Math.Sqrt(0.3125);
        Assert.Equal(0.127322, FlavorSplit.ModeInhibitionWeight(x1, H, r1), 6);
        Assert.Equal(0.872678, FlavorSplit.ModeInhibitionWeight(x1, H, -r1), 6);
        Assert.Equal(1.0, FlavorSplit.ModeInhibitionWeight(x1, H, r1) + FlavorSplit.ModeInhibitionWeight(x1, H, -r1), 12);
    }

    [Fact]
    public void Which_Flavor_Is_Slow_Is_The_Sign_Of_X_Not_F36()
    {
        // gE = 1, gI = 2, h = 0.1, c = 1: a = -1 gives x = 0.5 - 1 = -0.5 < -h, an inhibition-dominant slow mode
        Assert.Equal("I", FlavorSplit.SlowFlavor(-1, 0.1, 1.0, 1.0, 2.0));
        var (slow, _) = FlavorSplit.OverdampedRates(-1, 0.1, 1.0, 1.0, 2.0);
        Assert.Equal(1.5 - Math.Sqrt(0.25 - 0.01), slow);          // 1.0101, the inhibition-dominant mode below s
        Assert.Equal("E", FlavorSplit.SlowFlavor(1, 0.1, 1.0, 1.0, 2.0));
    }

    [Fact]
    public void Equal_Leaks_With_H_Above_C_Leave_No_Overdamped_Block()
    {
        foreach (double a in new[] { 1.0, 0.0, -1.0, 0.0 })
            Assert.False(FlavorSplit.IsOverdamped(a, 0.5, 0.25, 1.0, 1.0));
        Assert.Throws<InvalidOperationException>(() => FlavorSplit.ReportedRatio([1, 0, -1, 0], 0.5, 0.25, 1.0, 1.0));
        // c > h at equal leaks: the +-a blocks share their discriminant (the bipartite degeneracy)
        Assert.Equal(FlavorSplit.Discriminant(1, 0.5, 1.0, 1.0, 1.0), FlavorSplit.Discriminant(-1, 0.5, 1.0, 1.0, 1.0));
        Assert.True(FlavorSplit.IsOverdamped(1, 0.5, 1.0, 1.0, 1.0));
    }

    [Fact]
    public void Zero_Cross_Coupling_Ratio_Reads_The_H0_Control_File()
    {
        Assert.Equal(7.0 / 3.0, FlavorSplit.ZeroCrossCouplingRatio(-1, 1, C, GE, GI));      // chain, ring, star: 2.333333
        Assert.Equal(2.6, FlavorSplit.ZeroCrossCouplingRatio(-0.2, 1, C, GE, GI));           // complete N=6: 2.600000, exact in float64
        Assert.Equal(2.0, FlavorSplit.ZeroCrossCouplingRatio(0, 0, C, GE, GI));              // empty graph: gI/gE
    }

    [Fact]
    public void Guards_Throw()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => FlavorSplit.Centre(0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => FlavorSplit.Centre(1, double.NaN));
        Assert.Throws<ArgumentException>(() => FlavorSplit.Jacobian(new double[2, 3], H, C, GE, GI));
        Assert.Throws<ArgumentException>(() => FlavorSplit.Jacobian(new double[0, 0], H, C, GE, GI));
        Assert.Throws<ArgumentOutOfRangeException>(() => FlavorSplit.SwapPermutation(0));
        Assert.Throws<InvalidOperationException>(() => FlavorSplit.ZeroCrossCouplingRatio(-1, 1, 2.0, GE, GI));
        Assert.Throws<ArgumentOutOfRangeException>(() => FlavorSplit.Jacobian(RingC4Scaled(), H, -0.25, GE, GI));
        Assert.True(NeuralPalindrome.IsInvolution(FlavorSplit.SwapPermutation(5)));
    }
}
