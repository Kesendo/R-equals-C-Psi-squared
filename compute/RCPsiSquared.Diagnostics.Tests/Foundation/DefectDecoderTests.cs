using System;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The defect decoder (M2b): given an observed per-site α-profile, return the bond where the
/// defect sits and how strong it is, with the least-squares residual as confidence. The dictionary is
/// the painters' per-bond f-profiles, calibrated once at δJ_cal through the SAME PaintersMovement
/// pipeline the Symphony uses (one source of α semantics). These tests pin the spec's feasibility
/// table: exact-bond identification across the linear window, strength to ≤10% inside it, honest
/// degradation at the edge (bond still exact, strength may drift), the residual growing from window
/// center to edge (the confidence), and the N=5 mirror-pair degeneracy (site order distinguishes the
/// bonds, or it is reported as a real finding).</summary>
public class DefectDecoderTests
{
    // The decoder calibrates to the canonical PTF protocol the painters require: XY chain,
    // bonding-mode carrier, γ=0.05, J=1, δJ_cal=0.02 (the same construction PtfMovementTests.Canonical
    // uses). The decoder builds these planted-defect α-profiles through the same painters path.
    const double J = 1.0, Gamma = 0.05, DeltaJCal = 0.02;

    static DefectDecoder Calib(int n) => DefectDecoder.Calibrate(n, J, Gamma, DeltaJCal);

    /// <summary>Build the per-site α-profile a real run with a defect (bond b, δJ) would show, through
    /// the same painters pipeline — this is what the decoder must invert.</summary>
    static double[] AlphaProfile(int n, int bond, double deltaJ)
    {
        var s = new Symphony(n: n, j: J, gamma: Gamma, hType: HamiltonianType.XY,
                             initialState: InitialStateKind.BondingMode, defectBond: bond, deltaJ: deltaJ);
        var pm = ((IInspectable)s).Children.OfType<PaintersMovement>().Single();
        Assert.True(pm.HasLenses);
        return pm.Alphas.ToArray();
    }

    /// <summary>Build the per-site purity-DEVIATION profile a defect (bond b, δJ) shows, through the same
    /// painters pipeline — the de-lossed "observed reading" DecodeDeviation must invert.</summary>
    static double[] DeviationProfileObs(int n, int bond, double deltaJ)
    {
        var s = new Symphony(n: n, j: J, gamma: Gamma, hType: HamiltonianType.XY,
                             initialState: InitialStateKind.BondingMode, defectBond: bond, deltaJ: deltaJ);
        var pm = ((IInspectable)s).Children.OfType<PaintersMovement>().Single();
        Assert.True(pm.HasLenses);
        return pm.DeviationProfile.ToArray();
    }

    [Theory]
    [InlineData(0, 0.010)]
    [InlineData(1, 0.015)]
    [InlineData(2, 0.025)]
    [InlineData(2, -0.020)]
    public void N4_SpecTable_ExactBond_StrengthWithin10Percent(int bond, double deltaJ)
    {
        var dec = Calib(4);
        var alpha = AlphaProfile(4, bond, deltaJ);
        var r = dec.Decode(alpha);
        Assert.Equal(bond, r.Bond);
        double relErr = Math.Abs(r.DeltaJ - deltaJ) / Math.Abs(deltaJ);
        Assert.True(relErr <= 0.10, $"strength rel error {relErr:P2} should be ≤ 10% (decoded δĴ={r.DeltaJ}, truth {deltaJ})");
    }

    [Fact]
    public void N4_Residual_GrowsFromWindowCenterToEdge()
    {
        var dec = Calib(4);
        var center = dec.Decode(AlphaProfile(4, 1, 0.015));   // deep in the linear window
        var edge = dec.Decode(AlphaProfile(4, 0, 0.040));     // at the window edge
        Assert.True(center.Residual < edge.Residual,
            $"the residual is the confidence: it should grow toward the window edge " +
            $"(center {center.Residual:E3} < edge {edge.Residual:E3})");
    }

    [Fact]
    public void N4_WindowEdge_BondStillExact_StrengthMayDrift()
    {
        // (0, +0.04) sits at the edge of the linear window: the spec's table shows the bond is still
        // identified exactly but the strength error climbs to ~17% (the dictionary is linear, the true
        // response has curvature). Assert the bond only; document why the strength is not pinned here.
        var dec = Calib(4);
        var r = dec.Decode(AlphaProfile(4, 0, 0.040));
        Assert.Equal(0, r.Bond);
    }

    [Fact]
    public void N5_MirrorPair_Bond0_vs_Bond3_BothDecodeToCorrectBond()
    {
        // The symmetric bonding state makes bonds 0 and 3 a mirror pair; their f-profiles are
        // site-reversed, so the site order must distinguish them. If this fails it is NOT forced — the
        // degeneracy is a real result to report. Here we assert it succeeds (the order does resolve it).
        var dec = Calib(5);
        var r0 = dec.Decode(AlphaProfile(5, 0, 0.015));
        var r3 = dec.Decode(AlphaProfile(5, 3, 0.015));
        Assert.Equal(0, r0.Bond);
        Assert.Equal(3, r3.Bond);
    }

    [Fact]
    public void N5_NegativeEdgeDefect_ReportsAmbiguity_TruthAmongCandidates()
    {
        // Fix B (verified physics): at N=5 an edge bond weakened (bond 3, δJ = −0.02) has an f-profile
        // nearly ANTI-collinear with the complementary interior bond (bond 1 strengthened). The linear
        // dictionary cannot cleanly separate sign+location: bond 1 wins over the true bond 3 by only
        // ~1.5× in residual. The honest decoder must REPORT the ambiguity (not guess), and the truth
        // (bond 3) must be among the two candidates {Bond, RunnerUpBond}.
        var dec = Calib(5);
        var r = dec.Decode(AlphaProfile(5, 3, -0.02));
        Assert.True(r.IsAmbiguous,
            $"N=5 (3, −0.02) should be flagged ambiguous (residual ratio ≈ 1.5 < {DefectDecoder.AmbiguityFactor}); " +
            $"winner bond {r.Bond} res {r.Residual:E3}, runner-up bond {r.RunnerUpBond} res {r.RunnerUpResidual:E3}");
        Assert.True(r.Bond == 3 || r.RunnerUpBond == 3,
            $"the truth (bond 3) must be among the two candidates; got winner {r.Bond}, runner-up {r.RunnerUpBond}");
    }

    [Fact]
    public void N5_PositiveDefect_IsClean_NotAmbiguous_ExactBond()
    {
        // The mirror case: a positive defect is clean. Calibrating and decoding at the SAME point
        // (bond 3, +0.02) lands on the calibration profile, residual ~0, runner-up ratio ≫ 10 ⟹ not
        // ambiguous, and the bond is exact.
        var dec = Calib(5);
        var r = dec.Decode(AlphaProfile(5, 3, +0.02));
        Assert.False(r.IsAmbiguous,
            $"N=5 (3, +0.02) is a clean positive defect (ratio ≫ {DefectDecoder.AmbiguityFactor}); " +
            $"winner res {r.Residual:E3}, runner-up res {r.RunnerUpResidual:E3}");
        Assert.Equal(3, r.Bond);
    }

    [Theory]
    [InlineData(0, 0.010)]
    [InlineData(1, 0.015)]
    [InlineData(2, 0.025)]
    [InlineData(2, -0.020)]
    public void N4_SpecTable_AllUnambiguous_RatiosWellAboveThreshold(int bond, double deltaJ)
    {
        // At N=4 every spec-table case is unambiguous: the f-profiles are well-separated (residual
        // ratios ≫ 3), so the decoder never has to report a tie.
        var dec = Calib(4);
        var r = dec.Decode(AlphaProfile(4, bond, deltaJ));
        Assert.False(r.IsAmbiguous,
            $"N=4 ({bond}, {deltaJ}) should be unambiguous; winner bond {r.Bond} res {r.Residual:E3}, " +
            $"runner-up bond {r.RunnerUpBond} res {r.RunnerUpResidual:E3} (ratio {r.RunnerUpResidual / r.Residual:F2})");
    }

    [Fact]
    public void N5_DeLoss_DeviationPath_ResolvesMirrorPair_WithSign()
    {
        // The headline de-loss (spec §11 row [D]). The α path FLAGS the N=5 (bond 3 weakened) mirror pair
        // ambiguous (residual ratio ≈ 1.5); the SIGNED deviation path RESOLVES it (bond 3, δĴ < 0, ratio
        // ≫ AmbiguityFactor AND ≫ the α ratio). Bind QUALITATIVELY: the C# squared-residual ratio is ≈ 516
        // (the square of the ≈ 22.7 residual-norm ratio), do NOT pin the exact number. A ratio wildly off
        // ≈ 516 is a dense-vs-sector handshake break to investigate, not a number to loosen.
        var dec = Calib(5);

        var rAlpha = dec.Decode(AlphaProfile(5, 3, -0.02));
        Assert.True(rAlpha.IsAmbiguous, "the α path must still flag the N=5 mirror pair ambiguous");
        double ratioAlpha = rAlpha.RunnerUpResidual / rAlpha.Residual;

        var rDev = dec.DecodeDeviation(DeviationProfileObs(5, 3, -0.02));
        Assert.False(rDev.IsAmbiguous, "the signed deviation path must resolve the mirror pair (not ambiguous)");
        Assert.Equal(3, rDev.Bond);
        Assert.True(rDev.DeltaJ < 0, $"the deviation path must read the sign (weakened, δĴ < 0); got {rDev.DeltaJ}");
        double ratioDev = rDev.RunnerUpResidual / rDev.Residual;
        Assert.True(ratioDev > DefectDecoder.AmbiguityFactor,
            $"deviation ratio {ratioDev:F1} must clear the {DefectDecoder.AmbiguityFactor} threshold");
        Assert.True(ratioDev > ratioAlpha,
            $"the de-loss must beat the α path: deviation ratio {ratioDev:F1} > α ratio {ratioAlpha:F1}");
    }

    [Theory]
    [InlineData(4, 0, 0.010)]
    [InlineData(4, 1, 0.015)]
    [InlineData(4, 2, 0.025)]
    [InlineData(4, 2, -0.020)]
    [InlineData(5, 3, 0.020)]
    public void DeviationPath_RecoversBondAndSign_OnEasyCases(int n, int bond, double deltaJ)
    {
        // Parity: on the unambiguous cases the α path already handles, the deviation path recovers the
        // bond AND the sign.
        var dec = Calib(n);
        var r = dec.DecodeDeviation(DeviationProfileObs(n, bond, deltaJ));
        Assert.Equal(bond, r.Bond);
        Assert.Equal(Math.Sign(deltaJ), Math.Sign(r.DeltaJ));
    }

    [Theory]
    [InlineData(0, 0.010)]
    [InlineData(1, 0.015)]
    [InlineData(2, 0.025)]
    [InlineData(2, -0.020)]
    public void N4_DeviationPath_StrengthWithin10Percent(int bond, double deltaJ)
    {
        // Parity with N4_SpecTable_ExactBond_StrengthWithin10Percent (the α path): the deviation path's
        // signed δĴ recovers the strength MAGNITUDE too, not just bond + sign. This pins the
        // DeviationResponse = g/δJ normalization (a wrong divisor would leave bond + sign correct but δĴ
        // off by a factor of δJ_cal). Same deep-window N=4 cases and 10% tolerance as the α analogue.
        var dec = Calib(4);
        var r = dec.DecodeDeviation(DeviationProfileObs(4, bond, deltaJ));
        Assert.Equal(bond, r.Bond);
        double relErr = Math.Abs(r.DeltaJ - deltaJ) / Math.Abs(deltaJ);
        Assert.True(relErr <= 0.10, $"deviation strength rel error {relErr:P2} should be ≤ 10% (δĴ={r.DeltaJ}, truth {deltaJ})");
    }

    [Theory]
    [InlineData(2)]
    [InlineData(6)]
    public void Guard_RefusesNOutsideThreeToFive(int n)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => DefectDecoder.Calibrate(n, J, Gamma, DeltaJCal));
    }
}
