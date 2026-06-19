using System.Linq;
using RCPsiSquared.Core.Knowledge;
using Xunit;

namespace RCPsiSquared.Core.Tests.Knowledge;

/// <summary>Tom's marker idea (2026-05-28): the recognition itself is a marker, and
/// the markers are INDEPENDENT flags, not a single-valued classification. A thing can
/// be IsCore AND IsDeadEnd at once ("at least two sides"). The point of the test is the
/// "enttarnen": the two-flag combination IsCore ∧ IsDeadEnd unmasks a class that no
/// single flag gives — the load-bearing-unsolvable, the framework's living spine
/// (the carrier, g_eff, PTF), distinct from the closeable scaffold (the masks).</summary>
public class EpistemicFacetTests
{
    [Fact]
    public void Flags_AreIndependent_PtfIsCoreAndDeadEndAtOnce()
    {
        // PTF is not solvable (the mixing calc / alpha-closure) yet the framework
        // does not work without the perspectival structure. Two sides, both true.
        var ptf = EpistemicFacetMap.Facet("PTF");
        Assert.True(ptf.HasFlag(EpistemicFacet.IsCore));
        Assert.True(ptf.HasFlag(EpistemicFacet.IsDeadEnd));
    }

    [Fact]
    public void TwoFlags_Enttarnt_TheCoreDeadEndSpine()
    {
        // The carrier (gamma0), g_eff, and PTF each carry BOTH IsCore and IsDeadEnd:
        // essential and unsolvable. The two-flag query unmasks them as one class.
        var spine = EpistemicFacetMap.WithAll(EpistemicFacet.IsCore | EpistemicFacet.IsDeadEnd);
        Assert.Contains("gamma0", spine);
        Assert.Contains("g_eff", spine);
        Assert.Contains("PTF", spine);
    }

    [Fact]
    public void TwoFlags_Refine_WhatOneFlagCannot()
    {
        // The enttarnen, made precise: IsCore alone is too broad (it holds both the
        // closeable spine d^2-2d=0 and the unsolvable carrier). Adding IsDeadEnd
        // narrows to the load-bearing-unsolvable a single flag cannot isolate.
        var core = EpistemicFacetMap.WithAll(EpistemicFacet.IsCore);
        var spine = EpistemicFacetMap.WithAll(EpistemicFacet.IsCore | EpistemicFacet.IsDeadEnd);

        Assert.Contains("d2_minus_2d", core);          // core structure: in 'core'
        Assert.DoesNotContain("d2_minus_2d", spine);   // but NOT in the spine (it closes)
        Assert.Contains("gamma0", core);
        Assert.Contains("gamma0", spine);              // the carrier: in both
        Assert.True(spine.Count < core.Count);         // two flags narrow; the spine is unmasked
    }

    [Fact]
    public void PureStructure_IsCloseable_NotCore_NotDeadEnd()
    {
        // The masks (1/4 = 1/d^2): derivable scaffold, neither essential nor unsolvable.
        var quarter = EpistemicFacetMap.Facet("quarter");
        Assert.True(quarter.HasFlag(EpistemicFacet.IsStructure));
        Assert.False(quarter.HasFlag(EpistemicFacet.IsCore));
        Assert.False(quarter.HasFlag(EpistemicFacet.IsDeadEnd));

        // A pure mask is not part of the spine the two-flag query unmasks.
        var spine = EpistemicFacetMap.WithAll(EpistemicFacet.IsCore | EpistemicFacet.IsDeadEnd);
        Assert.DoesNotContain("quarter", spine);
    }

    [Fact]
    public void DeadEnd_WithoutCore_IsNotSpine_ButDeadEnd_WithRule_IsTheHeldExit()
    {
        // A pure dead-end (a number that does not close, nothing essential rides on it)
        // is IsDeadEnd without IsCore: not the spine.
        var xpeak = EpistemicFacetMap.Facet("x_peak");
        Assert.True(xpeak.HasFlag(EpistemicFacet.IsDeadEnd));
        Assert.False(xpeak.HasFlag(EpistemicFacet.IsCore));

        // The red button is a different two-flag object: a dead-end carrying a rule
        // (the self-erasing solution, held not pressed).
        var redButton = EpistemicFacetMap.Facet("red_button");
        Assert.True(redButton.HasFlag(EpistemicFacet.IsDeadEnd));
        Assert.True(redButton.HasFlag(EpistemicFacet.IsRule));
    }

    [Fact]
    public void IsAlive_IsIndependent_PtfIsCore_DeadEnd_AndAlive()
    {
        // PTF moves: it is the perspectival TIME field, the painters paint over time.
        // So besides core and dead-end it carries IsAlive. Three flags, all true.
        var ptf = EpistemicFacetMap.Facet("PTF");
        Assert.True(ptf.HasFlag(EpistemicFacet.IsCore));
        Assert.True(ptf.HasFlag(EpistemicFacet.IsDeadEnd));
        Assert.True(ptf.HasFlag(EpistemicFacet.IsAlive));
    }

    [Fact]
    public void IsAlive_Splits_LivingDeadEnds_FromStaticDeadEnds()
    {
        // The cut IsDeadEnd alone could not make: a dead-end that MOVES (PTF, the angle,
        // Im(lambda) != 0) versus one that just sits (g_eff, x_peak: a number that
        // neither closes nor stirs). IsAlive is the x,y,z + i / t axis.
        var deadEnds = EpistemicFacetMap.WithAll(EpistemicFacet.IsDeadEnd);
        var living = EpistemicFacetMap.WithAll(EpistemicFacet.IsDeadEnd | EpistemicFacet.IsAlive);

        Assert.Contains("PTF", living);
        Assert.Contains("the_angle", living);
        Assert.Contains("g_eff", deadEnds);          // g_eff IS a dead-end
        Assert.DoesNotContain("g_eff", living);      // but a static one, it does not move
        Assert.True(living.Count < deadEnds.Count);  // living dead-ends are a strict subset
    }

    [Fact]
    public void StaticThings_AreNotAlive()
    {
        // Masks, the contract, the carrier, the static dead-ends: real constants or pure
        // decay rates, no oscillation, no i over t. The metronome, not the dance.
        Assert.False(EpistemicFacetMap.Facet("quarter").HasFlag(EpistemicFacet.IsAlive));
        Assert.False(EpistemicFacetMap.Facet("d2_minus_2d").HasFlag(EpistemicFacet.IsAlive));
        Assert.False(EpistemicFacetMap.Facet("g_eff").HasFlag(EpistemicFacet.IsAlive));
        Assert.False(EpistemicFacetMap.Facet("gamma0").HasFlag(EpistemicFacet.IsAlive));
    }

    [Fact]
    public void HandshakeRk_IsAStaticDeadEnd_LikeXPeak()
    {
        // The bonding-mode per-bond defect response f(b) (handshake_decoder R_k): rate-protected
        // (U(1) + Π => re_shift(0,1) ~ 0) so it is eigenvector mixing, NOT a rate functional, and no
        // simple closed form fits it across N (EQ-021 power-law fails; the only-computable gate's
        // battery stays R^2 0.38-0.73 at N>=6). It is the PTF mixing calc instanced -- computable
        // (the painter, exact, carrier-independent to 0.999) but not derivable. A pure, static
        // dead-end like x_peak: nothing core rides on its value, and the per-N profile does not stir.
        var rk = EpistemicFacetMap.Facet("handshake_Rk");
        Assert.True(rk.HasFlag(EpistemicFacet.IsDeadEnd));
        Assert.False(rk.HasFlag(EpistemicFacet.IsCore));
        Assert.False(rk.HasFlag(EpistemicFacet.IsAlive));

        // It joins the static dead-ends (x_peak, g_eff), not the living ones (PTF, the_angle).
        var living = EpistemicFacetMap.WithAll(EpistemicFacet.IsDeadEnd | EpistemicFacet.IsAlive);
        Assert.DoesNotContain("handshake_Rk", living);
        Assert.Contains("handshake_Rk", EpistemicFacetMap.WithAll(EpistemicFacet.IsDeadEnd));
    }
}
