using System;
using System.Numerics;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Numerics;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>Route A of the Sturm design (sideways_spin_ladder arc; design note, local and NOT
/// in the repo: docs/superpowers/plans/2026-08-11-disc13-sturm-real-root-count-design.md, §2 +
/// §5 step 3; its content is summarized in the doc section this class gates):
/// the exact integer landing of disc_M(G) over ℤ. The evaluation grid takes the negative node
/// half for FREE from the derived-and-gated q-evenness of F_res (the bipartite gauge
/// T·L(q)·T = L(−q) acting within each R-sector at (1,3)@N=6), which is what brings the R-odd
/// cost under the design's ~30 min/parity threshold (quiet-machine rate 2026-08-11: 0.87 /
/// 2.37 s per node, measured after the load-contaminated scout figures were re-taken quiet).
/// Gates, all exact or at the two committed primes: deg/v of the landed polynomial against the
/// certified-backed 768/438 and 1038/594, fresh verification nodes INCLUDING a negative one
/// (the only end-to-end check of the mirroring), the composition identity
/// disc_Λ(F_res) = (−4)^m·f·disc_M(G)² at the POLYNOMIAL level mod both committed primes (sign
/// included; closes the degree-bookkeeping gap the design's machinery map names), and
/// gcd(u, disc_M G) = 1 exact over ℚ (the u-gate promoted from mod-p). The synthetic M-descent
/// control is the one control that can catch a G-reindex bug (design §6): every number it pins
/// is sympy-computed and recorded as a literal, coefficient for coefficient.</summary>
public class Disc13SturmTests
{
    private readonly ITestOutputHelper _out;
    public Disc13SturmTests(ITestOutputHelper output) => _out = output;

    /// <summary>The synthetic M-descent control (design §6): F(Λ, q) = G₀((Λ+12)², q) with
    /// G₀ = (M − (q²+1)²)·(M − 2), all reference values sympy-computed 2026-08-11 and pinned as
    /// literals: f = G₀(0, q) = 2(q²+1)² (a perfect square with lc 2, ũ = q²+1),
    /// disc_M(G₀) = ((q²+1)²−2)² = (q⁴+2q²−1)², landed coefficients low-to-high
    /// [1, 0, −4, 0, 2, 0, 4, 0, 1], deg_q disc_Λ(F) = 20, and the composition identity
    /// disc_Λ(F) = (−4)²·f·disc_M(G₀)² verified symbolically. The pipeline must reproduce the
    /// coefficients EXACTLY; the degree gates alone cannot see a reindexing bug.</summary>
    [Fact]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void SyntheticMDescent_LandsSympyPinnedDiscM()
    {
        // F(Λ, q) = ((Λ+12)² − (q²+1)²)·((Λ+12)² − 2), rows Λ^0..Λ^4, q-coefficients
        // low-to-high (sympy-expanded, pinned):
        var fRes = new[]
        {
            Row(20306, 0, -284, 0, -142),
            Row(6840, 0, -48, 0, -24),
            Row(861, 0, -2, 0, -1),
            Row(48),
            Row(1),
        };
        var r = FoldResultantCertificate.DiscMExactLandingFromBivariate(fRes, -12, certifiedDegQD: 20);
        _out.WriteLine($"synthetic: mDeg={r.MDeg} discDeg={r.DiscDeg} vQ={r.DiscVq} " +
                       $"verify={r.VerificationNodesPass} compId={r.CompositionIdentityModBothPrimes} " +
                       $"uGcd1={r.UGcdOne} layers=[{string.Join(",", r.StrippedLayerDegreesFirstPrime)}]");
        Assert.Equal(2, r.MDeg);
        Assert.Equal(8, r.DiscDeg);
        Assert.Equal(0, r.DiscVq);
        var expected = new BigInteger[] { 1, 0, -4, 0, 2, 0, 4, 0, 1 };
        Assert.Equal(expected.Length, r.DiscCoefficients.Length);
        for (int i = 0; i < expected.Length; i++)
            Assert.True(expected[i] == r.DiscCoefficients[i],
                $"disc_M coefficient q^{i}: expected {expected[i]}, landed {r.DiscCoefficients[i]}");
        Assert.True(r.VerificationNodesPass);
        Assert.True(r.CompositionIdentityModBothPrimes);
        Assert.True(r.UGcdOne, "gcd(q²+1, (q⁴+2q²−1)²) = 1 (value −2 at the u-roots)");
        Assert.Equal(new[] { 0, 4 }, r.StrippedLayerDegreesFirstPrime);
        Assert.Equal(new[] { 0, 4 }, r.StrippedLayerDegreesSecondPrime);
    }

    /// <summary>Route A landed on the census block, both parities, then the Sturm count and the
    /// psc device on the landed objects (the long pole: ~10 / ~47 min per parity, the R-odd
    /// bulk being the deg-954 psc gcd). Pins: deg/v of disc_M(G) = 768/438 (R-even) and
    /// 1038/594 (R-odd), both certified-backed (they chain from the certified deg/v of D, the
    /// exact deg f with v_q(f) = 0, and the proven composition identity); the mod-p layer split
    /// of the landed stripped object [90, 120] / [108, 168] at both committed primes (this pair
    /// is the design's EXPECTATION, the number this pipeline pins: the mod-p read of the exact
    /// object is still a floor for its true split, which step 4's squarefree PRS delivers);
    /// verification nodes (one negative), the polynomial-level composition identity at both
    /// primes, and the exact u-gate.</summary>
    [Theory]
    [InlineData(false, 1572, 768, 438, 90, 120)]
    [InlineData(true, 2124, 1038, 594, 108, 168)]
    [Trait("Category", "SLOW_DISC13_STURM")]
    public void N6_Disc13_DiscM_ExactLanding_AndSturmCount(bool rOdd, int certDegD, int discDeg,
        int vQ, int simpleLayerDeg, int doubleLayerDeg)
    {
        var r = FoldResultantCertificate.WeightDiscMExactLanding(6, 1, 3, rOdd, certDegD);
        _out.WriteLine($"N6 (1,3) rOdd={rOdd}: mDeg={r.MDeg} discDeg={r.DiscDeg} vQ={r.DiscVq} " +
                       $"verify={r.VerificationNodesPass} compId={r.CompositionIdentityModBothPrimes} " +
                       $"uGcd1={r.UGcdOne} layers=[{string.Join(",", r.StrippedLayerDegreesFirstPrime)}] / " +
                       $"[{string.Join(",", r.StrippedLayerDegreesSecondPrime)}]");
        Assert.Equal(rOdd ? 24 : 21, r.MDeg);
        Assert.Equal(discDeg, r.DiscDeg);
        Assert.Equal(vQ, r.DiscVq);
        Assert.True(r.VerificationNodesPass, "fresh nodes beyond the grid, incl. q0 < 0 " +
            "(the end-to-end check of the evenness mirroring)");
        Assert.True(r.CompositionIdentityModBothPrimes,
            "stripped-free polynomial identity disc_Λ(F_res) ≡ (−4)^m·f·disc_M(G)² mod both " +
            "committed primes, sign included");
        Assert.True(r.UGcdOne, "gcd(u, disc_M G) = 1 exact: the μ = 0 diabolics stay out of A");
        Assert.Equal(new[] { simpleLayerDeg, doubleLayerDeg }, r.StrippedLayerDegreesFirstPrime);
        Assert.Equal(new[] { simpleLayerDeg, doubleLayerDeg }, r.StrippedLayerDegreesSecondPrime);

        // The Sturm count on the exact split (design §3 + §4): the census pins of THIS parity
        // (experiments/F89_PATH_K_DIABOLIC.md, the census table) PLUS the two loci this count
        // found first (2026-08-11, sympy-isolated from the landed A): the ninth R-even locus
        // q ≈ 0.6133162 (between the census pins 0.4630 and 0.6400, missed by both census
        // channels) and the R-odd close-pair twin q ≈ 0.4627729 beside the census 0.4635508
        // (separation 7.8·10⁻⁴, below channel A's 0.001 grid step; the committed N=4 pair
        // 0.854438 / 0.857458 is the spacing precedent, though that one is same-species while
        // this pair is MIXED, fold beside axis). Brackets ±1/4000 (see SturmCountFromLanding).
        var pins = rOdd
            ? new[] { 0.4627729, 0.4635508, 0.8029541, 1.2273677, 1.9567664, 3.2608458 }
            : new[] { 0.4630013, 0.6133162, 0.6400121, 0.7026159, 0.7450837, 0.9184132, 1.5109958, 3.0922209, 3.3536850 };
        var s = FoldResultantCertificate.SturmCountFromLanding(r, pins);
        _out.WriteLine($"N6 (1,3) rOdd={rOdd} COUNT: degA={s.DegA} degB={s.DegB} split={s.SplitIdentityOk} " +
                       $"aEven={s.AEvenInQ} A(q>0)={s.CountAPositive} A(all)={s.CountAAllReal} " +
                       $"B(q>0)={s.CountBPositive} B(all)={s.CountBAllReal} " +
                       $"u(q>0)={s.CountUPositive} u(all)={s.CountUAllReal} " +
                       $"pins=[{string.Join(",", s.PinCountsA)}]");
        Assert.True(s.SplitIdentityOk, "P = A·B·B exact over ℤ + pairwise gcd/squarefree gates");
        Assert.Equal(simpleLayerDeg, s.DegA);   // the TRUE simple-layer degree: k = 0 pinned
        Assert.Equal(doubleLayerDeg, s.DegB);
        Assert.True(s.AEvenInQ, "A inherits the q-evenness (± root pairing)");
        Assert.Equal(2 * s.CountAPositive, s.CountAAllReal);

        // THE COUNT (first landed 2026-08-11, sympy-oracle-verified the same afternoon: A, B, u
        // counts and the root isolation all agree): the census floors 8 / 5 were both one
        // short: 9 R-even and 6 R-odd distinct defective loci on q > 0. With the two new loci
        // pinned, every bracket carries exactly one root.
        Assert.Equal(rOdd ? 6 : 9, s.CountAPositive);
        Assert.Equal(pins.Length, s.CountAPositive);
        foreach (var pc in s.PinCountsA)
            Assert.Equal(1, pc);

        // The corollaries, pinned as first-landed: real μ = 0 diabolics (u) and the doubled
        // off-line content (B; nonzero, so §7's free closure does not fire and the psc device
        // is the closure).
        Assert.Equal(rOdd ? 3 : 4, s.CountUPositive);
        Assert.Equal(2 * s.CountUPositive, s.CountUAllReal);
        Assert.Equal(rOdd ? 9 : 7, s.CountBPositive);
        Assert.Equal(2 * s.CountBPositive, s.CountBAllReal);

        // Export for the once-per-parity sympy oracle (design §4: an external count_roots on the
        // landed integer polynomials; a wrong chain still returns a plausible integer).
        var export = System.IO.Path.Combine(System.IO.Path.GetTempPath(),
            $"disc13_sturm_{(rOdd ? "odd" : "even")}.txt");
        var sb = new System.Text.StringBuilder();
        sb.AppendLine("A");
        foreach (var cf in s.ACoefficients) sb.AppendLine(cf.ToString());
        sb.AppendLine("B");
        foreach (var cf in s.BCoefficients) sb.AppendLine(cf.ToString());
        sb.AppendLine("U");
        foreach (var cf in r.UCoefficients) sb.AppendLine(cf.ToString());
        System.IO.File.WriteAllText(export, sb.ToString());
        _out.WriteLine($"oracle export: {export}");

        // The psc device in M (design §7's leading option), exact over ℤ on the landed objects:
        // gcd(psc₁, stripped disc) = 1 closes the mult-4 hideouts (two-double-M coincidences
        // everywhere, 3×3 Jordan cubic branches off μ = 0), and the Tarski chain reads the
        // species EXACTLY at the roots (−1 fold pair on the line, +1 conj pair on the axis).
        var d = FoldResultantCertificate.WeightPscDevice(r, s, pins);
        _out.WriteLine($"N6 (1,3) rOdd={rOdd} PSC: pscDeg={d.PscDeg} s10Deg={d.S10Deg} " +
                       $"coprime={d.DeviceCoprimeExact} pinSigns=[{string.Join(",", d.PinSpeciesSigns)}] " +
                       $"axis={d.AxisTotalPositive} fold={d.FoldTotalPositive}");
        Assert.True(d.DeviceCoprimeExact,
            "gcd(psc₁, stripped disc_M(G)) = 1 exact over ℤ: no two-double/triple-M point exists");
        Assert.Equal(rOdd ? 954 : 704, d.PscDeg);
        Assert.Equal(rOdd ? 956 : 704, d.S10Deg);
        // Census species column at the thirteen census pins (all confirmed by the exact chain),
        // plus the two new loci's species, first read off this chain 2026-08-11 and pinned as
        // measured: the ninth R-even locus 0.6133162 is a FOLD pair on the line, and the R-odd
        // twin 0.4627729 is a FOLD pair sitting 7.8·10⁻⁴ beside the AXIS locus 0.4635508, a
        // mixed-species close pair, which is exactly why the census missed it (channel B sees
        // only the axis species, and channel A's N=4 validation carried no such neighbourhood).
        var expectedSigns = rOdd
            ? new[] { -1, +1, +1, +1, -1, -1 }
            : new[] { -1, -1, -1, +1, -1, -1, -1, -1, +1 };
        Assert.Equal(expectedSigns, d.PinSpeciesSigns);
        // The species split of the WHOLE simple layer (the news the census never had):
        // 7 fold + 2 axis (R-even), 3 fold + 3 axis (R-odd). The cross-check is against the
        // pin signs, which the brackets pin root-for-root (axis+fold = n is true by
        // construction and checks nothing; a review round caught the tautology).
        Assert.Equal(rOdd ? 3 : 2, d.AxisTotalPositive);
        Assert.Equal(rOdd ? 3 : 7, d.FoldTotalPositive);
        int signSum = 0;
        foreach (var sg in d.PinSpeciesSigns) signSum += sg;
        Assert.Equal(d.AxisTotalPositive - d.FoldTotalPositive, signSum);
    }

    /// <summary>The N=4 (1,2) R-even control (design §6), the committed story end to end: the
    /// pipeline must reproduce the four committed √-branch loci 0.460212 / 0.854438 / 0.857458 /
    /// 1.738181 as EXACTLY the q &gt; 0 real roots of the simple layer, with the diabolic 0.659
    /// in the double layer and NOT counted. Committed shape (experiments doc):
    /// disc(F₈) = c·q²⁴·(3q⁴+q²−1)²·P₁₀(q²), so deg D = 52, v = 24, A = P₁₀(q²) (deg 20),
    /// B = 3q⁴+q²−1 (deg 4, one positive real root, the diabolic). Closes the recorded repo gap
    /// "no exact Sturm run on P₁₀ exists". The landing also upgrades the census's float
    /// sector-disc-reality to exact (Im D ≡ 0 identically, node count over the a-priori bound).</summary>
    [Fact]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void N4_Control_FourCommittedLoci_DiabolicExcluded()
    {
        var r = FoldResultantCertificate.WeightSturmN4Control();
        _out.WriteLine($"N4 (1,2) R-even: degD={r.DegD} vD={r.VqD} degA={r.DegA} degB={r.DegB} " +
                       $"split={r.SplitIdentityOk} A(q>0)={r.CountAPositive} B(q>0)={r.CountBPositive} " +
                       $"pins=[{string.Join(",", r.PinCountsA)}] pinB659={r.PinCountB659}");
        Assert.Equal(52, r.DegD);
        Assert.Equal(24, r.VqD);
        Assert.Equal(20, r.DegA);
        Assert.Equal(4, r.DegB);
        Assert.True(r.SplitIdentityOk);
        Assert.Equal(4, r.CountAPositive);
        Assert.Equal(new[] { 1, 1, 1, 1 }, r.PinCountsA);
        Assert.Equal(1, r.CountBPositive);
        Assert.Equal(1, r.PinCountB659);
    }

    /// <summary>The (1,2)@N=6 realness-contingency control (design §6): the census measured that
    /// block's sectors NOT conjugation-closed (disc complex, O(1) reality residual), so the ℤ
    /// landing must DECLINE loudly rather than silently produce a wrong real count.</summary>
    [Fact]
    [Trait("Category", "SLOW_DISC13_STURM")]
    public void N6_OneTwo_Control_LandingDeclines()
    {
        var ex = Assert.Throws<InvalidOperationException>(
            () => FoldResultantCertificate.WeightDiscMExactLanding(6, 1, 2, rOdd: false, certifiedDegQD: 0));
        _out.WriteLine($"decline message: {ex.Message}");
        Assert.Contains("not real", ex.Message);   // the RIGHT reason: five other gates also throw
    }

    private static RCPsiSquared.Core.Numerics.GaussianInteger[] Row(params int[] coefs)
    {
        var r = new RCPsiSquared.Core.Numerics.GaussianInteger[coefs.Length];
        for (int i = 0; i < coefs.Length; i++)
            r[i] = new RCPsiSquared.Core.Numerics.GaussianInteger(coefs[i], 0);
        return r;
    }
}
