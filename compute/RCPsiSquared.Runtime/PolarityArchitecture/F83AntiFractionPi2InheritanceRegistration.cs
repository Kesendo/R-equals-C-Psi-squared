using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F83AntiFractionPi2Inheritance"/>:
/// F83's anti-fraction = 1/(2 + 4·r) closed form. F83 is the FIRST F-formula
/// whose primary anchor sits on <see cref="BilinearApexClaim"/> directly:
/// the maximum 1/2 at r=0 IS the bilinear-apex argmax. Three parent edges
/// (one explicit, two as registration discards):
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_2 = 1/2</c>
///         (max anti-fraction), <c>a_0 = 2</c>, <c>a_{-1} = 4</c>, <c>a_{-2} = 8</c>,
///         <c>a_{1-N} = 2^N</c> (Hilbert-space dimension).</item>
///   <item><see cref="BilinearApexClaim"/>: registration discard. F83 is the
///         first F-formula on the argmax-side of the bilinear-apex pair;
///         before F83 the BilinearApexClaim had 0 descendants. The
///         "max anti-fraction = 1/2 at r=0" reading IS the BilinearApex
///         apex value (ON_THE_HALF Coda 2026-05-07's argmax-maxval pair).</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: registration discard.
///         F83 crosses 1/4 at r=1/2 — same anchor as F57 + Dicke + F60 fold.
///         The crossover edge ties the BilinearApex argmax (1/2) to its
///         maxval shadow (1/4) through F83's continuous family.</item>
/// </list>
///
/// <para>Tier consistency: F83 is Tier 1 proven (PROOF_F83_PI_DECOMPOSITION_RATIO);
/// derived from F49's Frobenius identity. All four claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers BilinearApexClaim + QuarterAsBilinearMaxvalClaim) +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F83AntiFractionPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF83AntiFractionPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F83AntiFractionPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            _ = b.Get<BilinearApexClaim>();                  // FIRST direct edge: argmax 1/2
            _ = b.Get<QuarterAsBilinearMaxvalClaim>();       // crossover 1/4 at r=1/2
            return new F83AntiFractionPi2Inheritance(ladder);
        });
}
