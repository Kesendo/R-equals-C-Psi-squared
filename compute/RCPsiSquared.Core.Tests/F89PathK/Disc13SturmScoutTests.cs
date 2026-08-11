using System;
using RCPsiSquared.Core.F89PathK;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>Scout steps 0-2 of the Sturm design (sideways_spin_ladder arc; design note, local:
/// docs/superpowers/plans/2026-08-11-disc13-sturm-real-root-count-design.md): §0 the exact
/// realness of F_res (expected real on the conjugation-closed (1,3) sectors; if this gate ever
/// fails it is the design's §0 contingency finding, not a code bug); §1 the exact G/f/u chain
/// (G a reindex of the certified-even F_res, f = G(0, q) with v_q(f) = 0 gated exactly, ũ the
/// monic square root with the LOAD-BEARING exact identity lc(f)·ũ² = f, upgrading the two-prime
/// perfect-square read to algebra); §2 the per-node disc_M(G) with the composition identity
/// disc_Λ(F_res) = (−4)^m·f·disc_M(G)² gated EXACTLY at fresh integer nodes (q0 = 2, 3, sign
/// included), plus the Route-A per-node timing the design's §5 decision leans on. The N=4 (1,2)
/// control exercises the declining branch (no evenness, f squarefree, the G-analysis void).</summary>
public class Disc13SturmScoutTests
{
    private readonly ITestOutputHelper _out;
    public Disc13SturmScoutTests(ITestOutputHelper output) => _out = output;

    /// <summary>The N=4 (1,2) control: evenness is void there (committed: expEven false, f-layers
    /// [8] squarefree, u-gate n/a), so the scout must DECLINE the G-analysis with sentinels, not
    /// silently produce a wrong G. f itself is still exact: deg f = 8 (the committed fDeg), and
    /// the square-root verification must return false on the squarefree f.</summary>
    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "EVEN_DISC13_SCOUT")]
    public void N4_Control_EvennessVoid_GAnalysisDeclines(bool rOdd)
    {
        var r = FoldResultantCertificate.WeightSturmScout(4, 1, 2, rOdd, timingNodes: 0);
        _out.WriteLine($"N4 (1,2) rOdd={rOdd}: real={r.ResidualIsReal} even={r.EvenInShiftedLambda} " +
                       $"mDeg={r.MDeg} fDeg={r.FDeg} vQf={r.VQf} square={r.FIsPerfectSquare} uDeg={r.UDeg}");
        Assert.False(r.EvenInShiftedLambda, "N=4 (1,2) has no (Λ+8)² evenness (committed control)");
        Assert.Equal(-1, r.MDeg);
        Assert.Equal(8, r.FDeg);
        Assert.False(r.FIsPerfectSquare, "the N=4 f is squarefree ([8]); the exact √ must refuse it");
        Assert.Equal(-1, r.UDeg);
        Assert.False(r.CompositionIdentityAtNodes, "no G, no identity: the sentinel is false/n-a");
    }

    /// <summary>Steps 0-2 on the census block, both parities. Pins from the design note: ResDeg
    /// 42/48, deg G in M = 21/24 (monic), deg f = 36/48 with v_q(f) = 0 (gated exactly for the
    /// first time here; the halving v(disc G) = v(D)/2 leans on it), lc(f)·ũ² = f EXACT (the
    /// load-bearing upgrade of the mod-p perfect-square read), deg ũ = 18/24, and the composition
    /// identity disc_Λ(F_res) = (−4)^m·f·disc_M(G)² exact at the fresh nodes q0 = 2, 3 with the
    /// (−1)^m sign (m = 21 R-even is odd: a genuine minus). The realness gate is §0 of the design;
    /// the per-node timing is logged for the §5 Route-A decision, gated only > 0 (quiet-machine
    /// numbers are re-measured before they are quoted anywhere).</summary>
    [Theory]
    [InlineData(false, 42, 21, 36, 18)]
    [InlineData(true, 48, 24, 48, 24)]
    [Trait("Category", "SLOW_EVEN_DISC13")]
    public void N6_Disc13_SturmScout_Steps0to2(bool rOdd, int resDeg, int mDeg, int fDeg, int uDeg)
    {
        var r = FoldResultantCertificate.WeightSturmScout(6, 1, 3, rOdd, timingNodes: 10);
        _out.WriteLine($"N6 (1,3) rOdd={rOdd}: real={r.ResidualIsReal} even={r.EvenInShiftedLambda} " +
                       $"mDeg={r.MDeg} fDeg={r.FDeg} vQf={r.VQf} square={r.FIsPerfectSquare} " +
                       $"uDeg={r.UDeg} uReal={r.UIsReal} compId={r.CompositionIdentityAtNodes} " +
                       $"ms/node(discM)={r.MsPerNodeDiscM:F1}");
        Assert.Equal(resDeg, r.ResDeg);
        Assert.True(r.ResidualIsReal, "§0: F_res expected real (conjugation-closed sector; " +
                                      "a failure here is the design's §0 contingency finding)");
        Assert.True(r.EvenInShiftedLambda, "the certified evenness, re-read on this path");
        Assert.Equal(mDeg, r.MDeg);
        Assert.Equal(fDeg, r.FDeg);
        Assert.Equal(0, r.VQf);
        Assert.True(r.FIsPerfectSquare, "lc(f)·ũ² = f must hold EXACTLY (load-bearing)");
        Assert.Equal(uDeg, r.UDeg);
        Assert.True(r.UIsReal);
        Assert.True(r.CompositionIdentityAtNodes,
            "disc_Λ(F_res)(q0) = (−4)^m·f(q0)·disc_M(G)(q0)² exact at q0 = 2, 3, sign included");
        Assert.True(r.MsPerNodeDiscM > 0.0, "timing nodes ran");
    }
}
