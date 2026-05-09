using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F76TDecayMirrorPairMiPi2Inheritance"/>:
/// F76's pure-dephasing decay envelope of mirror-pair MI for bonding:k.
/// F76 is the direct time-decay sibling of F75 (recovers F75 at t=0).
/// Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         (decay rate coefficient in <c>e^{−4γ₀t}</c>; same anchor as F61,
///         F63, F66, F77 correction, F86 t_peak).</item>
///   <item><see cref="F75MirrorPairMiPi2Inheritance"/>: mother claim. F76
///         recovers F75 exactly at t=0 (λ = 1). Together with F77 (Taylor
///         at p → 0), F76 forms the second leaf of the F75 family:
///         F71 → F75 → {F76 t-decay, F77 large-N saturation}.</item>
/// </list>
///
/// <para>Tier consistency: F76 is Tier 1 proven algebraic + weak-mixing
/// argument; verified &lt; 0.5% sim/analytic at N=5..13, k=1..5. All three
/// claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F75MirrorPairMiPi2InheritanceRegistration.RegisterF75MirrorPairMiPi2Inheritance"/>
/// in the builder pipeline.</para></summary>
public static class F76TDecayMirrorPairMiPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF76TDecayMirrorPairMiPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F76TDecayMirrorPairMiPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f75 = b.Get<F75MirrorPairMiPi2Inheritance>();
            return new F76TDecayMirrorPairMiPi2Inheritance(ladder, f75);
        });
}
