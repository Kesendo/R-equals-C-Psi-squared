using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for GammaFold's PER-SITE turns (built 2026-08-20). The object's two mirrors turn
// the whole profile at once; these tests hold the law for turning ONE site, which is the move the
// arc site_resolved_vacuum_block asked for and the object did not carry.
//
// On a cell |i><j| the rate is -2*sum_{l differs} gamma_l, so s_l leaves every cell where site l
// agrees untouched and moves every cell where it differs by exactly +4*gamma_l: exactly half the
// cells, all by the same step. The whole turn s is a different animal, moving EVERY cell by twice
// its own rate, which is what makes it the reflection r -> -r.
//
// Two answers are pinned here and they are the arc's two questions. ON THE PROFILE the single turn
// is unconditional, exact, an involution, and the N of them commute, so they generate
// (Z/2)^|support|. ON THE RATE AXIS the turn has a shadow only when the profile restricted to its
// NONZERO sites has distinct subset sums; uniform gamma is the extreme failure, where the rate sees
// only how many sites disagree. And where the shadow exists it is PIECEWISE the identity and a
// translation, so it is no reflection, and the dihedral <s, s0> does not grow to absorb it.
//
// Two disciplines the file keeps. The subset-sum criterion is decided over the INTEGERS, and for a
// reason narrower than the tempting one: not because 1/10 + 2/10 = 3/10 slips past a float test
// (the stored doubles for those decimals have distinct sums and both routes agree), but because
// distinct exact sums can ROUND TOGETHER, so a float test invents collisions that are not there.
// And the rate identities are compared to literal zero on DYADIC profiles, where they are exact,
// with the price of a non-representable profile read in its own test rather than hidden in a
// tolerance; that price tracks representability and not magnitude, so scaling switches it off
// rather than sweeping it.
public class SiteTurnTests
{
    static readonly World W = new();

    // dyadic AND with distinct subset sums: exact in binary, and the rate pins down the
    // disagreement set.
    static readonly double[] Dyadic = { 0.125, 0.25, 0.5, 1.0 };

    [Fact]
    public void One_Site_Turned_Moves_Exactly_Half_The_Cells_By_Exactly_Four_Gamma()
    {
        var fold = new GammaFold(W, 4, siteGammas: Dyadic);
        for (int l = 0; l < 4; l++)
        {
            var r = fold.SiteTurn(l);
            Assert.Equal(1 << (2 * 4 - 1), r.CellsWhereSiteDiffers);      // exactly half of 4^N
            Assert.Equal(1 << (2 * 4 - 1), r.CellsWhereSiteAgrees);
            Assert.True(r.WorstStepResidual == 0.0,
                $"site {l}: the step is exactly 4*gamma_l where the site differs and exactly zero "
                + $"where it agrees; residual {r.WorstStepResidual:E1} on a dyadic profile would be "
                + "a finding, not noise");
            Assert.Equal(4.0 * Dyadic[l], r.Step);
        }
    }

    [Fact]
    public void The_Composite_Of_All_Single_Turns_Is_The_Whole_Turn_Exactly()
    {
        foreach (var n in new[] { 2, 3, 4 })
        {
            var fold = new GammaFold(W, n, siteGammas: Dyadic.Take(n).ToArray());
            var r = fold.SiteTurnGroup();
            Assert.True(r.WorstCompositeResidual == 0.0,
                $"N={n}: composing the N single turns gives r -> -r exactly, "
                + $"residual {r.WorstCompositeResidual:E1}");
            Assert.Equal(1 << n, r.OrbitSize);
            Assert.Equal(n, r.Support);
            Assert.True(r.AllInvolutions);
            Assert.True(r.AllCommute);
        }
    }

    // an unwatched site is turned to itself, so it contributes no factor: the group is
    // (Z/2)^support and not (Z/2)^N. The negative-zero trap is what would hide this.
    [Fact]
    public void An_Unwatched_Site_Contributes_No_Factor_To_The_Group()
    {
        var fold = new GammaFold(W, 4, siteGammas: new[] { 0.0, 0.25, 0.5, 1.0 });
        var r = fold.SiteTurnGroup();
        Assert.Equal(3, r.Support);
        Assert.Equal(8, r.OrbitSize);
    }

    [Fact]
    public void Sigma_Moves_By_Twice_The_Turned_Rates_Exactly()
    {
        var fold = new GammaFold(W, 4, siteGammas: Dyadic);
        var r = fold.SiteTurnGroup();
        Assert.True(r.WorstSigmaResidual == 0.0,
            $"sigma(s_S gamma) = sigma - 2*sum_(l in S) gamma_l, residual {r.WorstSigmaResidual:E1}");
    }

    // THE CRITERION, and the two rows that separate it from the obvious guess. Dissociation of the
    // WHOLE profile is sufficient and not necessary: a zero rate may collide freely, because the
    // turn moves it by 4*0 = 0. The criterion is the support.
    [Theory]
    [InlineData(new[] { 0.125, 0.25, 0.5, 1.0 }, true, true)]     // dyadic, all sums distinct
    [InlineData(new[] { 1.0, 2.0, 4.0, 8.0 }, true, true)]        // powers of two
    [InlineData(new[] { 0.0, 1.0, 2.0, 4.0 }, false, true)]       // a zero rate: whole fails, support does not
    [InlineData(new[] { 0.0, 0.25, 0.5, 1.0 }, false, true)]      // the same, dyadic
    [InlineData(new[] { 0.25, 0.25, 0.25, 0.25 }, false, false)]  // UNIFORM: the extreme failure
    [InlineData(new[] { 0.25, 0.5, 0.75, 2.0 }, false, false)]    // 1/4 + 1/2 = 3/4
    [InlineData(new[] { 0.25, 0.25, 0.5, 1.0 }, false, false)]    // two equal sites
    public void The_Shadow_On_The_Rate_Axis_Follows_The_SUPPORT_And_Not_The_Whole_Profile(
        double[] profile, bool wholeDistinct, bool expectedDescends)
    {
        var fold = new GammaFold(W, profile.Length, siteGammas: profile);
        var r = fold.SiteTurnGroup();
        Assert.Equal(wholeDistinct, r.WholeProfileSumsDistinct);
        Assert.Equal(expectedDescends, r.SupportSumsDistinct);
        Assert.Equal(expectedDescends, r.DescendsToRateAxis);
    }

    // the criterion is decided over the integers, and this is the profile that shows why it must
    // be. The tempting justification is the wrong one: 1/10 + 2/10 = 3/10 is a statement about
    // rationals a double array cannot hold, and the values actually stored for 0.1, 0.2 and 0.3
    // have distinct sums, so a float test gets THAT case right. The real gap is the false
    // collision: distinct exact sums rounding together. {1.0, 1e-20} is the smallest witness.
    [Fact]
    public void The_Criterion_Is_Decided_Exactly_Where_A_Float_Test_Would_Invent_A_Collision()
    {
        Assert.True(1.0 + 1e-20 == 1.0,
            "the premise: in doubles the two sums are the same number");
        var fold = new GammaFold(W, 2, siteGammas: new[] { 1.0, 1e-20 });
        var r = fold.SiteTurnGroup();
        Assert.True(r.WholeProfileSumsDistinct,
            "over the integers the four subset sums are distinct, so the profile qualifies; a "
            + "float HashSet would have collapsed 1 + 1e-20 onto 1 and answered no");
        Assert.True(r.DescendsToRateAxis);

        // and the case the tempting justification points at, where both routes agree
        var ten = new GammaFold(W, 3, siteGammas: new[] { 0.1, 0.2, 0.3 });
        Assert.True(ten.SiteTurnGroup().WholeProfileSumsDistinct,
            "the stored doubles for 0.1, 0.2, 0.3 do NOT collide; only the decimals they are "
            + "named after do");
    }

    // where the shadow exists it is PIECEWISE the identity and a translation. Two pieces, and the
    // moving one shifts by 4*gamma_l. That is why it is no reflection.
    [Fact]
    public void Where_It_Descends_The_Shadow_Is_Piecewise_Identity_And_Translation()
    {
        var fold = new GammaFold(W, 4, siteGammas: Dyadic);
        var r = fold.SiteTurnGroup();
        Assert.True(r.DescendsToRateAxis);
        Assert.True(r.DescendedIsInvolution);
        Assert.Equal(2, r.PieceCount);
        Assert.Equal(4.0 * Dyadic.Max(), r.PieceShift);
    }

    // the arc's second question: does the dihedral grow? It does not, and the measurement is that
    // (s_l o s0)^2 is not even total on the spectrum: the orbit leaves the set.
    [Fact]
    public void The_Dihedral_Does_Not_Absorb_The_Single_Turn()
    {
        var fold = new GammaFold(W, 4, siteGammas: Dyadic);
        var r = fold.SiteTurnGroup();
        Assert.True(r.DescendsToRateAxis);
        Assert.False(r.SquaredWithAntiWatchIsTotal);
    }

    // uniform gamma is where the question disappears: the rate sees only HOW MANY sites disagree,
    // so no single site can be addressed. The whole turn survives it, which is why the object never
    // noticed the difference.
    [Fact]
    public void Under_Uniform_Watching_There_Is_No_Single_Site_Shadow()
    {
        var fold = new GammaFold(W, 4, siteGammas: new[] { 0.25, 0.25, 0.25, 0.25 });
        var r = fold.SiteTurnGroup();
        Assert.False(r.SupportSumsDistinct);
        Assert.False(r.DescendsToRateAxis);
        Assert.Equal(5, r.DistinctRateValues);          // |D| = 0..N, so exactly N+1 values
        Assert.True(r.WorstCompositeResidual == 0.0);
    }

    // the price of a non-representable profile, read rather than gated, because it is a
    // deterministic function of an input the value does not contain. The input is REPRESENTABILITY,
    // not magnitude, and the measurement says so: scaling the same profile by powers of a hundred
    // does not sweep the residual smoothly, it switches it off, because the scaled entries land
    // back on representable values. An earlier version of this test gated a ratio across decades
    // as though the residual tracked the scale; three of its four points were exactly zero.
    [Fact]
    public void The_Price_Is_A_Function_Of_Representability_And_Not_Of_Scale()
    {
        const double eps = 2.220446049250313e-16;
        var seen = new List<(double Scale, double Worst)>();
        foreach (var scaleUp in new[] { 1.0, 1e2, 1e4, 1e6 })
        {
            var profile = new[] { 0.2, 0.35, 0.5, 0.9 }.Select(x => x * scaleUp).ToArray();
            var fold = new GammaFold(W, 4, siteGammas: profile);
            double worst = 0;
            for (int l = 0; l < 4; l++) worst = Math.Max(worst, fold.SiteTurn(l).WorstStepResidual);
            seen.Add((scaleUp, worst));
            Assert.True(worst <= eps * 4.0 * profile.Sum(),
                $"at scale {scaleUp:0e+0} the price is at most one rounding of the summed scale: "
                + $"{worst:E1} against {eps * 4.0 * profile.Sum():E1}");
        }
        Assert.Contains(seen, x => x.Worst > 0.0);
        Assert.Contains(seen, x => x.Worst == 0.0);
        Assert.True(seen.First().Worst > 0.0,
            "the unscaled profile is the non-representable one and must pay: "
            + string.Join(", ", seen.Select(x => $"{x.Scale:0e+0}->{x.Worst:E1}")));
    }
}
