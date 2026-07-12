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
        // N=6..8: the exact SE-EP (coherence_horizon_se_block.py qstar_se), NOT the 2N/pi asymptote.
        Assert.Equal(2.889253, Formulas.Qstar(6), 4);
        Assert.Equal(3.419782, Formulas.Qstar(7), 4);
        Assert.Equal(3.961618, Formulas.Qstar(8), 4);
        // the finite-N EP sits strictly BELOW the asymptote it approaches (the bug returned the asymptote):
        for (int n = 4; n <= 8; n++)
            Assert.True(Formulas.Qstar(n) < 2.0 * n / Math.PI, $"Q*({n}) must be below the 2N/pi asymptote");
        // monotone increasing across the exact->asymptote splice at N=8->9 (used by any horizon sweep):
        for (int n = 2; n <= 8; n++)
            Assert.True(Formulas.Qstar(n) < Formulas.Qstar(n + 1), $"Q*(N) not monotone at N={n}");
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
        Assert.Equal(Math.Sqrt(2.0 / 48), Formulas.F49c_CrossTermCrossing(3), 10);          // shadow-crossing
        Assert.Equal(Math.Log(10), Formulas.F55_KDeath, 9);                                 // K_death
        Assert.Equal(5, Formulas.F55_ImmortalModes(4));                                     // N+1 immortal
    }

    [Fact]
    public void Dwell_And_InitialCPsi_Family()
    {
        Assert.Equal(1.0 / 3, Formulas.F60_GhzCPsi0(2), 10);
        Assert.Equal(1.0 / 7, Formulas.F60_GhzCPsi0(3), 10);
        Assert.Equal(1.0 / 15, Formulas.F60_GhzCPsi0(4), 10);                               // GHZ_4 below 1/4
        Assert.Equal(1.0 / 3, Formulas.F62_WstateCPsi0(2), 10);                             // W_2 = Bell+, 1/3
        Assert.Equal(10.0 / 81, Formulas.F62_WstateCPsi0(3), 10);
        Assert.Equal(1.080, Formulas.F59_DwellPrefactor(2, 0.5, 0.3709), 3);               // reduces to F57
    }

    [Fact]
    public void SingleExcitation_And_BlockStructure()
    {
        Assert.Equal(64L, Formulas.F63_BlockDim(4));                                        // 4^(N-1)
        Assert.Equal((3, 2), Formulas.F63_ConservedPerSector(4));
        Assert.Equal((2, 1), Formulas.F63_ConservedPerSector(2));
        var a3 = Formulas.F65_SingleExcitationRates(3);
        Assert.Equal(0.5, a3[0], 10); Assert.Equal(1.0, a3[1], 10); Assert.Equal(0.5, a3[2], 10);
        var a5 = Formulas.F65_SingleExcitationRates(5);
        Assert.Equal(1.0 / 6, a5[0], 10); Assert.Equal(0.5, a5[1], 10); Assert.Equal(2.0 / 3, a5[2], 10);
        Assert.True(Formulas.F65_RatesRational(3));
        Assert.False(Formulas.F65_RatesRational(4));                                        // golden-irrational
        Assert.Equal(5, Formulas.F66_PoleMultiplicity(4));
    }

    [Fact]
    public void Partner_Optimum_Kinematic_Family()
    {
        Assert.Equal(2 * 0.05 - 0.0138, Formulas.F68_PartnerRate(0.0138, 0.05), 10);        // alpha_p = 2 g0 - alpha_b
        Assert.True(Formulas.F69_N3Optimum > 0.25);                                          // above the fold
        Assert.Equal(1, Formulas.F70_MaxVisibleDeltaN(1));
        Assert.Equal(2, Formulas.F70_MaxVisibleDeltaN(2));
        Assert.Equal(2, Formulas.F71_C1IndependentComponents(5));
        Assert.Equal(1, Formulas.F71_C1IndependentComponents(2));
        Assert.Equal(+1, Formulas.F71_ReflectionParity(1));
        Assert.Equal(-1, Formulas.F71_ReflectionParity(2));
    }

    [Fact]
    public void Dicke_Qudit_StructuralCeiling_Family()
    {
        Assert.Equal(0.3, Formulas.F98_DickeAsymptote(4), 10);                              // (N+2)/(4(N+1))
        Assert.Equal(8.0 / 28, Formulas.F98_DickeAsymptote(6), 10);
        Assert.Equal(16, Formulas.F121_PairedCeiling(2, 2));                                // d=2: fully paired = 2^(2N)
        Assert.Equal(54, Formulas.F121_PairedCeiling(3, 2));                                // d=3: 54/81 partial
        Assert.Equal(9, Formulas.F121_CoherenceCount(3, 2, 0));
        Assert.Equal(36, Formulas.F121_CoherenceCount(3, 2, 1));
        Assert.Equal(0.8, Formulas.F122_CompleteCeiling(5), 10);                            // K_5 = 4/5
        Assert.Equal(2.0 / 3, Formulas.F122_StarCeiling(7), 10);                            // star_7 = 4/6
        Assert.Equal(2 - 2 / Math.Sqrt(3), Formulas.F122_K4Ceiling(), 10);
        Assert.Equal(1.6, Formulas.F122_RingCommutant(5), 10);                              // ring-5 (1,1)
        Assert.Equal(1.0, Formulas.F122_RingCommutant(4), 10);                              // ring-4 co-occupies the band edge
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

    // --- F85: the k-body residual trichotomy. c(truly)=0, c(Pi^2-odd)=1, c(Pi^2-even non-truly)=2;
    // truly iff #Y even AND #Z even; ||M||^2 per term = 4 c ||H_k||^2 2^N. For k>=3 only the
    // Pi^2-class matters (YYY has n_YZ=3 but c=1, the registry's own flag example). ---
    [Fact]
    public void F85_KBody_Trichotomy_Counts_And_Residual()
    {
        Assert.Equal(0, Formulas.F85_FrobeniusFactor("YY"));                    // truly: both parities even
        Assert.Equal(1, Formulas.F85_FrobeniusFactor("XY"));                    // Pi^2-odd
        Assert.Equal(2, Formulas.F85_FrobeniusFactor("YZ"));                    // Pi^2-even, not truly
        Assert.Equal(1, Formulas.F85_FrobeniusFactor("YYY"));                   // n_YZ=3 yet c=1

        // the trichotomy enumeration over {X,Y,Z}^k, pinned to the registry table
        var expected = new Dictionary<int, (int Truly, int Odd, int Even)>
        {
            [2] = (3, 4, 2), [3] = (7, 14, 6), [4] = (21, 40, 20),
        };
        foreach (var (k, (truly, odd, even)) in expected)
        {
            int cTruly = 0, cOdd = 0, cEven = 0;
            foreach (var tuple in Tuples("XYZ", k))
                switch (Formulas.F85_FrobeniusFactor(tuple))
                {
                    case 0: cTruly++; break;
                    case 1: cOdd++; break;
                    default: cEven++; break;
                }
            Assert.Equal((truly, odd, even), (cTruly, cOdd, cEven));
            Assert.Equal(odd, Formulas.F85_Pi2OddCount(k));                     // (3^k - (-1)^k)/2
        }

        Assert.Equal(4.0 * 32, Formulas.F85_ResidualNormSqPerTerm(5, 1.0, 1), 10);   // 4 c ||H||^2 2^N
        Assert.Equal(8.0 * 32, Formulas.F85_ResidualNormSqPerTerm(5, 1.0, 2), 10);
        Assert.Equal(0.0, Formulas.F85_ResidualNormSqPerTerm(5, 3.7, 0), 12);
    }

    static IEnumerable<string> Tuples(string alphabet, int k)
        => k == 0 ? new[] { "" } : Tuples(alphabet, k - 1).SelectMany(t => alphabet.Select(c => t + c));

    // --- F97: the Mandelbrot main cardioid at b = 1/2. c(phi) = z*(1-z*) with z* = (1/2) e^{i phi};
    // cusp c(0) = 1/4 (the fold boundary F16 re-seen), tail c(pi) = -3/4, top c(pi/2) = 1/4 + i/2. ---
    [Fact]
    public void F97_Cardioid_Anchors_And_Identity()
    {
        Assert.Equal(0.25, Formulas.F97_Cardioid(0.0).Real, 12);
        Assert.Equal(0.0, Formulas.F97_Cardioid(0.0).Imaginary, 12);
        Assert.Equal(-0.75, Formulas.F97_Cardioid(Math.PI).Real, 12);
        Assert.Equal(0.25, Formulas.F97_Cardioid(Math.PI / 2).Real, 12);
        Assert.Equal(0.5, Formulas.F97_Cardioid(Math.PI / 2).Imaginary, 12);
        for (int s = 0; s < 1000; s++)
        {
            double phi = 2.0 * Math.PI * s / 1000.0;
            var z = Formulas.F97_FixedPoint(phi);
            Assert.Equal(0.5, z.Magnitude, 12);                                 // pinned to b = 1/2
            Assert.True((Formulas.F97_Cardioid(phi) - z * (System.Numerics.Complex.One - z)).Magnitude < 1e-15,
                $"the algebraic identity c = z*(1-z*) broke at phi={phi}");
        }
    }

    // --- F74: chromaticity of the (n, n+1) popcount coherence block at J = 0 -- the count
    // c(n, N) = min(n, N-1-n) + 1 and the odd rate ladder 2*gamma0*{1,3,..,2c-1}, pinned by DIRECT
    // ENUMERATION of all popcount-(n, n+1) basis pairs (HD = the Pair disagreement, rate -2*gamma0*HD).
    [Fact]
    public void F74_Chromaticity_And_Ladder_Match_Enumeration()
    {
        for (int N = 4; N <= 9; N++)
            for (int n = 0; n < N; n++)
            {
                var hds = new SortedSet<int>();
                for (int x = 0; x < (1 << N); x++)
                {
                    if (System.Numerics.BitOperations.PopCount((uint)x) != n) continue;
                    for (int y = 0; y < (1 << N); y++)
                        if (System.Numerics.BitOperations.PopCount((uint)y) == n + 1)
                            hds.Add(System.Numerics.BitOperations.PopCount((uint)(x ^ y)));
                }
                int c = Formulas.F74_Chromaticity(n, N);
                Assert.Equal(c, hds.Count);
                Assert.Equal(Enumerable.Range(0, c).Select(i => 2 * i + 1), hds);
                var ladder = Formulas.F74_RateLadder(n, N, 0.05);
                Assert.Equal(c, ladder.Length);
                Assert.Equal(0.1, ladder[0], 12);
                Assert.Equal(2 * 0.05 * (2 * c - 1), ladder[c - 1], 12);
            }
    }

    // --- F72: the DD + CC block split of per-site purity -- pinned FROM BELOW on a deterministic
    // dense Hermitian unit-trace rho at N = 5 with ALL joint-popcount sectors populated (NOT
    // positive-semidefinite: the split is a purely kinematic identity, deliberately pinned beyond
    // physical states). The site marginals are computed by a GENERIC partial trace that scans every
    // (u, v) entry and knows nothing about popcount blocks; the block content enters only through
    // projected COPIES of rho. The one-site invisibility of |DeltaN| >= 2 blocks is STRUCTURAL to
    // the single-site trace (that IS F70's proof: |Dpc| >= 2 forces a difference at another site),
    // so its zero is asserted together with the DISCRIMINATING contrast the registry's own k = 2
    // line supplies: the SAME |DeltaN| = 2 content is NONZERO in a two-site marginal -- the zero is
    // the physics of k = 1, not an artifact of the harness. Load-bearing coefficient pin: the
    // assembled 1/2 + DD + CC equals the directly computed Tr(rho_i^2) with both sectors active.
    // No dynamics run, any rho.
    [Fact]
    public void F72_SitePurity_Splits_Blockwise_With_No_Cross_Term()
    {
        int N = 5, dim = 1 << N;
        var rho = new System.Numerics.Complex[dim, dim];        // deterministic pseudo-dense Hermitian
        for (int i = 0; i < dim; i++)
            for (int j = 0; j < dim; j++)
            {
                var a = new System.Numerics.Complex((i * 3 + j * 7) % 11 - 5, (i * 5 + j * 2) % 7 - 3);
                var b = new System.Numerics.Complex((j * 3 + i * 7) % 11 - 5, (j * 5 + i * 2) % 7 - 3);
                rho[i, j] = (a + System.Numerics.Complex.Conjugate(b)) / 2.0;
            }
        double tr = 0;
        for (int i = 0; i < dim; i++) tr += rho[i, i].Real;
        for (int i = 0; i < dim; i++) rho[i, i] += (1.0 - tr) / dim;   // unit trace, still Hermitian

        int Dpc(int u, int v) => System.Numerics.BitOperations.PopCount((uint)u)
                               - System.Numerics.BitOperations.PopCount((uint)v);

        // generic partial trace over all sites but `site`: scans EVERY (u, v) pair of the projected
        // rho and keeps rho[u, v] iff u and v agree everywhere except possibly at `site`. No popcount
        // logic in here; `keep` projects rho onto a block class before tracing.
        System.Numerics.Complex[,] Marginal(int site, Func<int, int, bool> keep)
        {
            var m = new System.Numerics.Complex[2, 2];
            for (int u = 0; u < dim; u++)
                for (int v = 0; v < dim; v++)
                    if ((u & ~(1 << site)) == (v & ~(1 << site)) && keep(u, v))
                        m[(u >> site) & 1, (v >> site) & 1] += rho[u, v];
            return m;
        }

        int populatedHigh = 0;                                   // the |DeltaN| >= 2 sectors are genuinely there
        for (int u = 0; u < dim; u++)
            for (int v = 0; v < dim; v++)
                if (Math.Abs(Dpc(u, v)) >= 2 && rho[u, v].Magnitude > 0) populatedHigh++;
        Assert.True(populatedHigh > 100, "the pin needs a rho with real |DeltaN| >= 2 content");

        // the discriminating contrast (the registry's k = 2 doorway): the |DeltaN| = 2 content the
        // single site cannot see IS visible to a two-site marginal (pairs may differ at both sites,
        // so |Dpc| = 2 entries survive the trace). Sites 0 and 1, generic pair partial trace.
        double pairHigh = 0;
        int pairMask = ~((1 << 0) | (1 << 1));
        for (int u = 0; u < dim; u++)
            for (int v = 0; v < dim; v++)
                if ((u & pairMask) == (v & pairMask) && Math.Abs(Dpc(u, v)) == 2)
                    pairHigh += rho[u, v].Magnitude;
        Assert.True(pairHigh > 0.1, "|DeltaN| = 2 content must be visible to a two-site marginal");

        for (int site = 0; site < N; site++)
        {
            var full = Marginal(site, (u, v) => true);
            var dd = Marginal(site, (u, v) => Dpc(u, v) == 0);
            var cc = Marginal(site, (u, v) => Math.Abs(Dpc(u, v)) == 1);
            var high = Marginal(site, (u, v) => Math.Abs(Dpc(u, v)) >= 2);

            for (int a = 0; a < 2; a++)
                for (int b = 0; b < 2; b++)
                {
                    Assert.Equal(0.0, high[a, b].Magnitude, 15);            // F70: >= 2 steps invisible to one site
                    Assert.Equal(0.0, (full[a, b] - dd[a, b] - cc[a, b] - high[a, b]).Magnitude, 12);
                }
            Assert.Equal(0.0, dd[0, 1].Magnitude, 15);                      // DD carries only the diagonal (z)
            Assert.Equal(0.0, cc[0, 0].Magnitude + cc[1, 1].Magnitude, 15); // CC carries only the off-diagonal (x, y)

            double x = 2.0 * full[0, 1].Real, y = -2.0 * full[0, 1].Imaginary;
            double z = (full[0, 0] - full[1, 1]).Real;
            double direct = full[0, 0].Real * full[0, 0].Real + full[1, 1].Real * full[1, 1].Real
                          + 2.0 * full[0, 1].Magnitude * full[0, 1].Magnitude;
            Assert.Equal(direct, Formulas.F72_SitePurity(x, y, z), 9);      // 1/2 + DD + CC, no cross term
            Assert.True(Math.Abs(x) + Math.Abs(y) > 1e-3 && Math.Abs(z) > 1e-3,
                "the pin must exercise both sectors");
        }
    }

    // --- F73: the spatial-sum coherence closure -- sum_i 2*|(rho_i)_{01}(t)|^2 = (1/2) e^{-4 gamma0 t}
    // on the vac-SE coherent probe, pinned FROM BELOW by running the living world (Restless RK4) with a
    // deliberately asymmetric |alpha> under three different U(1) Hamiltonians (XY chain, XXZ chain
    // at zz = 1 i.e. Delta = 2, XY ring): the closure must be blind to all of them. The site-i marginal's
    // off-diagonal is computed inline by partial trace (sum over basis states with bit i = 0). Cross-dock:
    // the (vac, SE) block is F74's n = 0 mono-chromatic block, c(0, N) = 1, whose single Pair rate
    // -2*gamma0 the squared magnitude pays twice.
    [Fact]
    public void F73_SpatialSum_Closure_HBlind_Under_Restless()
    {
        int N = 5, dim = 1 << N;
        double gamma = 0.05;
        Assert.Equal(1, Formulas.F74_Chromaticity(0, N));       // the mono-chromatic end F73 rides
        double[] alpha = { 0.6, -0.3, 0.5, 0.2, -0.51 };        // asymmetric, no XY eigenmode
        double norm = Math.Sqrt(alpha.Sum(a => a * a));
        foreach (var (bonds, zz) in new[] { (Topology.Chain(N), 0.0), (Topology.Chain(N), 1.0), (Topology.Ring(N), 0.0) })
        {
            var w = new Restless(W, N, 1.0, gamma, bonds, zz: zz);
            for (int l = 0; l < N; l++)
                w.SeedCoherence(0, 1 << l, alpha[l] / norm / 2);    // rho0 = (|vac><alpha| + h.c.)/2
            for (int tick = 0; tick <= 1600; tick++)
            {
                if (tick % 400 == 0)
                {
                    double sum = 0, moved = 0;
                    for (int site = 0; site < N; site++)
                    {
                        var m01 = System.Numerics.Complex.Zero; // (rho_site)_{01} = sum_{s: bit site of s = 0} rho[s, s | bit]
                        for (int s = 0; s < dim; s++)
                            if (((s >> site) & 1) == 0) m01 += w[s, s | (1 << site)];
                        sum += 2.0 * m01.Magnitude * m01.Magnitude;
                        double atJ0 = alpha[site] * alpha[site] / (norm * norm)
                                    * Formulas.F73_SpatialSumClosure(gamma, w.T);  // the J = 0 per-site prediction
                        moved += Math.Abs(2.0 * m01.Magnitude * m01.Magnitude - atJ0);
                    }
                    Assert.Equal(Formulas.F73_SpatialSumClosure(gamma, w.T), sum, 9);
                    // the H-rotation must actually happen for "H-blind" to mean anything: at the final
                    // time the per-site pattern has left the J = 0 shape even though the sum has not.
                    if (tick == 1600) Assert.True(moved > 1e-3, "H did not redistribute the per-site coherence");
                }
                if (tick < 1600) w.Step(0.0025);
            }
        }
    }

    // --- F77: the multi-drop MM(0) saturates at 1 bit -- the closed asymptotic pinned against the
    // EXACT F75 mirror-pair sum at full best-k search (the registry's own verified rows), and the
    // rescaled deviation (MM-1)(N+1) descending toward 3/(4 ln 2) = 1.0820.
    [Fact]
    public void F77_MM0_Saturation_Against_Exact_F75_BestK()
    {
        Assert.Equal(1.0820, Formulas.F77_RescaledDeviationLimit(), 4);
        double prev = double.MaxValue;
        foreach (var (N, mm0, resc) in new[] { (101, 1.01078, 1.100), (201, 1.00540, 1.091), (1001, 1.00108, 1.084) })
        {
            double best = 0.0;
            for (int k = 1; k <= N; k++) best = Math.Max(best, Formulas.F75_MirrorPairSum(N, k));
            Assert.Equal(mm0, best, 5);
            Assert.Equal(resc, (best - 1.0) * (N + 1), 3);
            Assert.Equal(Formulas.F77_MMSaturation(N), best, 3);
            double dev = (best - 1.0) * (N + 1);
            Assert.True(dev > Formulas.F77_RescaledDeviationLimit() && dev < prev,
                "the rescaled deviation must descend toward 3/(4 ln 2) from above");
            prev = dev;
        }
    }

    // --- F75: mirror-pair mutual information for single-excitation mirror-symmetric states,
    // MI = 2 h(p) - h(2p) with p the site population; independent of the mirror sign; saturates at
    // 2 bits (a Bell pair) at p = 1/2. Bonding:k populations are the F65 amplitudes squared, and
    // MM(0) = the sum over mirror pairs is O(N), no propagation. Table pins from the registry. ---
    [Fact]
    public void F75_MirrorPair_MI_And_The_Bonding_Sum()
    {
        Assert.Equal(2.0, Formulas.F75_MirrorPairMI(0.5), 12);              // the Bell saturation
        Assert.Equal(0.0, Formulas.F75_MirrorPairMI(0.0), 12);
        Assert.Equal(0.25, Formulas.F75_BondingSitePopulation(5, 2, 0), 12); // (1/3) sin^2(pi/3) = 1/4
        Assert.Equal(0.800, Formulas.F75_MirrorPairSum(5, 1), 3);
        Assert.Equal(1.245, Formulas.F75_MirrorPairSum(5, 2), 3);           // the N=5 maximiser (even k: node at the centre)
        Assert.Equal(0.918, Formulas.F75_MirrorPairSum(5, 3), 3);
        Assert.Equal(1.245, Formulas.F75_MirrorPairSum(7, 4), 3);
        Assert.Equal(1.145, Formulas.F75_MirrorPairSum(11, 6), 3);
        Assert.Equal(0.961, Formulas.F75_MirrorPairSum(13, 7), 3);
    }

    // --- F76: pure-dephasing decay of the mirror-pair MI. The pair coherence decays at 4 gamma0
    // (lambda = e^{-4 gamma0 t}), populations stay; at lambda = 1 the pair entropy is h(2p) (F75
    // recovered), at lambda = 0 it is h(1-2p) + 2p. The 0.93 envelope at gamma0 = 0.05, t = 0.1 is
    // the gamma0 signature, not a hidden constant. ---
    [Fact]
    public void F76_Dephasing_Envelope_Recovers_F75_And_Explains_The_093()
    {
        foreach (double p in new[] { 0.1, 0.25, 0.4 })
        {
            Assert.Equal(Formulas.F75_MirrorPairMI(p), Formulas.F76_MirrorPairMI(p, 1.0), 12);
            Assert.Equal(2 * H2(p) - H2(1 - 2 * p) - 2 * p, Formulas.F76_MirrorPairMI(p, 0.0), 12);
        }
        // the registry table column, N=9/11/13 rows corrected 2026-07-04 against envelope_study.py
        // (the adoption's from-below pin caught the three stale cells; the exact column was right)
        Assert.Equal(0.936, Formulas.F76_Envelope(5, 2, 0.05, 0.1), 3);
        Assert.Equal(0.932, Formulas.F76_Envelope(7, 2, 0.05, 0.1), 3);
        Assert.Equal(0.929, Formulas.F76_Envelope(9, 4, 0.05, 0.1), 3);
        Assert.Equal(0.927, Formulas.F76_Envelope(11, 4, 0.05, 0.1), 3);
        Assert.Equal(0.926, Formulas.F76_Envelope(13, 4, 0.05, 0.1), 3);
        // the envelope is the gamma0 signature, monotone in the watching: gentler, higher
        Assert.Equal(0.964, Formulas.F76_Envelope(5, 2, 0.025, 0.1), 3);
        Assert.Equal(0.888, Formulas.F76_Envelope(5, 2, 0.10, 0.1), 3);
        Assert.True(Formulas.F76_Envelope(5, 2, 0.025, 0.1) > Formulas.F76_Envelope(5, 2, 0.05, 0.1));
        Assert.True(Formulas.F76_Envelope(5, 2, 0.05, 0.1) > Formulas.F76_Envelope(5, 2, 0.10, 0.1));
    }

    static double H2(double x) => x <= 0 || x >= 1 ? 0.0 : -x * Math.Log2(x) - (1 - x) * Math.Log2(1 - x);

    // --- F95: the theta-compass at the quadratic discriminant zero. For z^2 - 2bz + c = 0 the
    // complex-root angle is theta = arctan(sqrt(c/b^2 - 1)) above c = b^2; at b = 1/2 the threshold
    // is 1/4 and the Februar compass arctan(sqrt(4c-1)) is recovered (= the adopted F15); the
    // Lindblad specialization is theta = arctan(Q), the adopted Clock's angle. ---
    [Fact]
    public void F95_Theta_Compass_Docks_Onto_F15_And_The_Clock()
    {
        Assert.Equal(0.0, Formulas.F95_Theta(0.25, 0.5), 12);                       // the degenerate double root
        Assert.True(double.IsNaN(Formulas.F95_Theta(0.2, 0.5)));                    // below threshold: no angle
        foreach (double c in new[] { 1.0 / 3, 0.308, 0.5, 1.0 })
        {
            Assert.Equal(Formulas.F95_ThetaHalf(c), Formulas.F95_Theta(c, 0.5), 12);
            Assert.Equal(Formulas.F15_ThetaDeg(c), Formulas.F95_ThetaHalf(c) * 180.0 / Math.PI, 10);
        }
        foreach (var (j, g) in new[] { (1.0, 0.5), (2.0, 0.5), (0.3, 0.6) })        // theta = arctan(Q)
            Assert.Equal(new Clock(W, j, g).ThetaDeg, Formulas.F95_Theta(g * g + j * j, g) * 180.0 / Math.PI, 10);
    }

    // --- F99: the five canonical trig anchors. alpha(theta) = sin^2(theta)/2 at {0,30,45,60,90}
    // degrees gives the five Pi2 dyadic anchors {0, 1/8, 1/4, 3/8, 1/2}; the non-uniform Dicke
    // weight is c^2 = cos(theta)/(2 sin^2(theta/2)) with the silver ratio 1+sqrt2 at 45 degrees,
    // 2 sqrt3 + 3 at 30, and the uniform Dicke c = 1 at 60. ---
    [Fact]
    public void F99_Five_Canonical_Angles_Give_The_Dyadic_Anchors()
    {
        var anchors = new[] { (0.0, 0.0), (30.0, 0.125), (45.0, 0.25), (60.0, 0.375), (90.0, 0.5) };
        foreach (var (deg, alpha) in anchors)
            Assert.Equal(alpha, Formulas.F99_Alpha(deg * Math.PI / 180.0), 12);
        Assert.Equal(2.0 * Math.Sqrt(3.0) + 3.0, Formulas.F99_DickeWeightSq(30.0 * Math.PI / 180.0), 10);
        Assert.Equal(1.0 + Math.Sqrt(2.0), Formulas.F99_DickeWeightSq(45.0 * Math.PI / 180.0), 10);   // the silver ratio
        Assert.Equal(1.0, Formulas.F99_DickeWeightSq(60.0 * Math.PI / 180.0), 10);                    // uniform Dicke
        Assert.Equal(0.0, Formulas.F99_DickeWeightSq(90.0 * Math.PI / 180.0), 12);
        foreach (double deg in new[] { 30.0, 45.0, 60.0, 90.0 })                    // gamma = c^2/(1+c^2) = cos(theta)
        {
            double c2 = Formulas.F99_DickeWeightSq(deg * Math.PI / 180.0);
            Assert.Equal(Math.Cos(deg * Math.PI / 180.0), c2 / (1.0 + c2), 10);
        }
    }

    // --- F88b: the popcount-coherence Pi^2-odd / memory closed form for (|p> + |q>)/sqrt2. Three
    // alpha anchors from one Krawtchouk identity (0 at popcount-mirror, the K-intermediate ratio,
    // 1/2 generic), the static fraction s from the sector sizes, and the HD = N anchor (GHZ:
    // Pi^2-classical, zero odd content). The adjacent K-intermediate alpha IS the adopted F98. ---
    [Fact]
    public void F88b_Popcount_Coherence_Anchors_And_The_F98_Dock()
    {
        // the three alpha anchors
        Assert.Equal(0.0, Formulas.F88b_Alpha(6, 2, 4), 12);                        // popcount-mirror n_p + n_q = N
        Assert.Equal(0.0, Formulas.F88b_Alpha(6, 3, 3), 12);                        // intra-mirror at N/2 (even N)
        Assert.Equal(0.5, Formulas.F88b_Alpha(7, 1, 3), 12);                        // generic
        foreach (int n in new[] { 4, 6, 8 })                                        // adjacent K-intermediate = F98
            Assert.Equal(Formulas.F98_DickeAsymptote(n), Formulas.F88b_Alpha(n, n / 2, n / 2 + 1), 12);
        // the static fraction: inter- and intra-sector
        Assert.Equal(1.0 / (4 * 15) + 1.0 / (4 * 20), Formulas.F88b_StaticFraction(6, 2, 3), 12);
        Assert.Equal(1.0 / 20, Formulas.F88b_StaticFraction(6, 3, 3), 12);
        // HD = N is Pi^2-classical (GHZ_N, Bell at N=2): zero odd content in memory
        Assert.Equal(0.0, Formulas.F88b_Pi2OddInMemory(4, 0, 4, 4), 12);
        // otherwise (1/2 - alpha s)/(1 - s); generic small case pinned by direct evaluation
        double s = Formulas.F88b_StaticFraction(7, 1, 3);
        Assert.Equal((0.5 - 0.5 * s) / (1.0 - s), Formulas.F88b_Pi2OddInMemory(7, 1, 3, 2), 12);
        // the multi-state Dicke extension: alpha_total = (1 - gamma^2)/2, anchors {1/2, 3/8, 0}
        Assert.Equal(0.5, Formulas.F88b_DickeAlphaTotal(0.0), 12);
        Assert.Equal(0.375, Formulas.F88b_DickeAlphaTotal(0.5), 12);
        Assert.Equal(0.0, Formulas.F88b_DickeAlphaTotal(1.0), 12);
    }

    // --- F124: the band-edge transition invariant ||M||_F^2 + lambda_min = 2 (the coordination
    // number), split as (2 - E) + E with E = (4/(N+1)) sin^2(pi/(N+1)) -- exactly the k=1 rung of
    // the already-adopted F65 ladder (the carrier's weight on the two free ends). ---
    [Fact]
    public void F124_BandEdge_Invariant_Splits_On_The_F65_Rung()
    {
        Assert.Equal(0.5, Formulas.F124_EndWeight(3), 12);                      // sin^2(pi/4) = 1/2 exactly
        foreach (int n in new[] { 3, 4, 5, 8, 20 })
        {
            Assert.Equal(Formulas.F65_SingleExcitationRates(n)[0], Formulas.F124_EndWeight(n), 12);
            Assert.Equal(2.0, Formulas.F124_FrobeniusNormSq(n) + Formulas.F124_SpectralFloor(n), 12);
        }
        // the floor vanishes as (N+1)^-3: E (N+1)^3 -> 4 pi^2
        Assert.Equal(4.0 * Math.PI * Math.PI, Formulas.F124_EndWeight(2000) * Math.Pow(2001, 3), 3);
    }
}
