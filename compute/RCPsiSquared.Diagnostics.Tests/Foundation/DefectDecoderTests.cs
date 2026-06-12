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

    [Theory]
    [InlineData(2)]
    [InlineData(6)]
    public void Guard_RefusesNOutsideThreeToFive(int n)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => DefectDecoder.Calibrate(n, J, Gamma, DeltaJCal));
    }
}
