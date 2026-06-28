using MirrorWorld;

namespace MirrorWorldTests;

// From-below smoke test: every number MirrorWorld adopted from Tier-1 is recomputed or checked here,
// so a regression or a wrong adopted constant fails loudly. Nothing interpreted, just the numbers.
public class SmokeTests
{
    const double G = 0.5;
    static readonly World W = new();

    // --- the base: own (left) vs inherited (right) ---
    [Fact]
    public void World_Owns_The_Frame()
    {
        Assert.Equal(new[] { "x", "y", "z" }, W.Own);
        Assert.Empty(W.Inherited);                              // the root, no parent
    }

    [Fact]
    public void Pair_Reads_Disagreement_And_Inherits_The_Frame()
    {
        var p = new Pair(W, 0b011, 0b101, G);                   // popcount(011 ^ 101) = popcount(110) = 2
        Assert.Equal(2, p.Disagreement);
        Assert.Equal(-2.0 * G * 2, p.Rate, 12);
        Assert.Equal(new[] { "x", "y", "z" }, p.Inherited);     // inherited from the World
    }

    // --- Grading A: bare sectors 2^N C(N,k); the two bases (pairs, Pauli strings) agree in count ---
    [Theory]
    [InlineData(2, new[] { 4, 8, 4 })]
    [InlineData(3, new[] { 8, 24, 24, 8 })]
    [InlineData(4, new[] { 16, 64, 96, 64, 16 })]
    public void BareSectors_Match_And_Agree_Across_Bases(int n, int[] expected)
    {
        Assert.Equal(expected, Redistribution.Bare(n));
        var modes = PauliMode.Enumerate(W, n, G).ToList();
        for (int k = 0; k <= n; k++)
            Assert.Equal(expected[k], modes.Count(m => m.K == k));
        for (int k = 0; k <= n; k++)
            Assert.Equal(expected[k], expected[n - k]);         // palindrome fold k <-> N-k
    }

    // --- the superposition: each inner sector splits into four equal Klein cells ---
    [Theory]
    [InlineData(2, 1)]
    [InlineData(3, 1)]
    [InlineData(4, 2)]
    public void KleinQuartering_Four_Equal_Cells(int n, int k)
    {
        var cells = PauliMode.Enumerate(W, n, G).Where(m => m.K == k).GroupBy(m => m.Klein).ToList();
        Assert.Equal(4, cells.Count);
        int each = (1 << n) * (int)Block.Binomial(n, k) / 4;
        Assert.All(cells, c => Assert.Equal(each, c.Count()));
    }

    // --- H on: the on-grid d_total folds (DEGENERACY_PALINDROME), palindromic, edges N+1 ---
    [Theory]
    [InlineData(2, new[] { 3, 10, 3 }, 16)]
    [InlineData(3, new[] { 4, 14, 14, 4 }, 36)]
    [InlineData(4, new[] { 5, 20, 152, 20, 5 }, 202)]
    public void OnGridFold_Palindrome_Edges_Sum(int n, int[] expected, int sum)
    {
        var g = Redistribution.OnGrid(n)!;
        Assert.Equal(expected, g);
        Assert.Equal(sum, g.Sum());
        Assert.Equal(n + 1, g[0]);                             // kernel
        Assert.Equal(n + 1, g[n]);                             // drain
        for (int k = 0; k <= n; k++) Assert.Equal(g[k], g[n - k]);
    }

    // --- Grading B: blocks C(N,p)C(N,q) sum to 4^N ---
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void Blocks_Sum_To_4N(int n)
    {
        long total = 0;
        for (int p = 0; p <= n; p++)
            for (int q = 0; q <= n; q++)
            {
                var b = new Block(W, n, p, q);
                Assert.Equal(Block.Binomial(n, p) * Block.Binomial(n, q), b.Size);
                total += b.Size;
            }
        Assert.Equal(1L << (2 * n), total);
    }

    // --- the coherence horizon Q*(N) and the coherence hand omega_mem ---
    [Fact]
    public void Qstar_Exact()
    {
        Assert.Equal(1.0, Formulas.Qstar(2), 4);
        Assert.Equal(Math.Sqrt(2), Formulas.Qstar(3), 4);
        Assert.Equal(1.8787, Formulas.Qstar(4), 4);
        Assert.Equal(2.3737, Formulas.Qstar(5), 4);
    }

    [Fact]
    public void OmegaMem_Sqrt2_Phi_Sqrt3()
    {
        Assert.Equal(Math.Sqrt(2), Formulas.OmegaMem(3, 1, 0), 6);
        Assert.Equal((1 + Math.Sqrt(5)) / 2, Formulas.OmegaMem(4, 1, 0), 6);   // phi
        Assert.Equal(Math.Sqrt(3), Formulas.OmegaMem(5, 1, 0), 6);
    }

    // --- the dispersions ---
    [Fact]
    public void F2_Dispersion_N3_Is_2_6()
    {
        var w = Formulas.F2_Dispersion(3, 1.0);
        Assert.Equal(2.0, w[0], 6);
        Assert.Equal(6.0, w[1], 6);
    }

    [Fact]
    public void F2b_Is_Chiral_Paired()                          // E_k = -E_{N+1-k}
    {
        var e = Formulas.F2b_SingleExcitation(5, 1.0);
        for (int k = 0; k < 5; k++) Assert.Equal(e[k], -e[4 - k], 10);
    }

    [Fact]
    public void F7_Spread_Is_Cot2()
    {
        var (_, _, _, spread) = Formulas.F7_QSpectrum(4, 1.0, 0.5);
        double cot = 1.0 / Math.Tan(Math.PI / 8);
        Assert.Equal(cot * cot, spread, 6);
    }

    [Fact]
    public void F33_N3_Rates_Are_1_4over3_5over3()
    {
        var r = Formulas.F33_N3Rates(0.5);
        Assert.Equal(1.0, r[0], 10);
        Assert.Equal(4.0 / 3, r[1], 10);
        Assert.Equal(5.0 / 3, r[2], 10);
    }

    [Fact]
    public void F25_Crossing_Satisfies_3half()                  // f*(1+f*^2) = 3/2
    {
        double f = Formulas.F25_CrossingF;
        Assert.Equal(1.5, f * (1 + f * f), 3);
    }

    // --- the simple closed forms in one place ---
    [Fact]
    public void SimpleClosedForms()
    {
        Assert.Equal(5, Formulas.F4_KernelDim(4));                          // N+1
        Assert.Equal(5.0 / 256, Formulas.F23_XorFraction(4), 10);          // (N+1)/4^N
        Assert.Equal(8, Formulas.F50_Weight1Degeneracy(4));                // 2N
        Assert.Equal(8, Formulas.F50_Weight1Degeneracy(3, triangleK3: true));
        Assert.Equal(new[] { 0, 2 }, Formulas.F34_QubitNecessity());       // d^2-2d=0
        var (unp, pm) = Formulas.F8_DecayLaw(4, 0.5);
        Assert.Equal(2.0, unp / pm, 10);                                   // the 2x ratio
        var (mn, mx, bw) = Formulas.F3_RateBounds(4, 0.5);
        Assert.Equal(1.0, mn, 10);
        Assert.Equal(3.0, mx, 10);
        Assert.Equal(2.0, bw, 10);
        Assert.Equal(0.25, Formulas.F16_FoldBoundary);                     // the 1/4 fold
        Assert.Equal(Math.Log(2) / 8, Formulas.F27_KX, 10);               // K_X = ln2/8
        Assert.Equal(+1, Formulas.F38_PiSquared(0, 0));                    // Pi^2 parity
        Assert.Equal(-1, Formulas.F38_PiSquared(1, 0));
        Assert.Equal(-1, Formulas.F61_PiSquaredX(1, 0));                   // the k-parity
        Assert.Equal(8.0, Formulas.D1_Bandwidth(1000, 1.0), 1);           // -> 8J at large N
        Assert.Equal(0.5, Formulas.D4_CrossingRhs(2), 10);
        Assert.Equal(1.5, Formulas.D4_CrossingRhs(4), 10);
        Assert.Equal(1.0, Formulas.D6_Gap(0.5), 10);                       // gap = 2g
    }

    [Fact]
    public void CrossTerm_Determinant_PalindromicTime_Family()
    {
        Assert.Equal(-1, Formulas.F39_DetPi(1));
        Assert.Equal(+1, Formulas.F39_DetPi(4));
        Assert.Equal(0.0, Formulas.F49_CrossTerm(2), 10);                  // exact Pythagorean at N=2
        Assert.Equal(1.0 / Math.Sqrt(48), Formulas.F49_CrossTerm(3), 10);
        Assert.Equal(1.0 / Math.Sqrt(128), Formulas.F49_CrossTerm(4), 10);
        Assert.Equal(0.25 * 256 * 4, Formulas.F49b_CenteredDissipatorNormSq(4, 0.5), 6);   // g^2 4^N N = 256
        Assert.Equal(Math.PI / (4 * Math.Pow(Math.Sin(Math.PI / 8), 2)), Formulas.F41_PalindromicTime(4, 1.0), 9);
        // F44: a palindromic pair d_fast + d_slow = 2 Sg; here Sg=1, d_fast=1.5, d_slow=0.5
        Assert.Equal(Math.Log(1.5 / 0.5), Formulas.F44_LogRatio(1.5, 0.5, 1.0), 10);
    }

    // --- the even/odd self-mirror: half-filling survivor exists only at even N ---
    [Theory]
    [InlineData(2, true)]
    [InlineData(3, false)]
    [InlineData(4, true)]
    [InlineData(5, false)]
    public void HalfFillingSurvivor_Iff_Even(int n, bool expected)
    {
        Assert.Equal(expected, new Survivor(W, n).HasHalfFillingSurvivor);
    }
}
