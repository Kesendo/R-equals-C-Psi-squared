using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F88StaticDyadicAnchor"/>: F88's dyadic-N
/// singleton-mirror static-fraction inherits the Pi2 halving ladder. Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: the parent ladder; the dyadic anchors
///         a_(k+2) = 2^(−(k+1)) live there. F88StaticDyadicAnchor inherits the lattice
///         positions from it.</item>
///   <item><see cref="PopcountCoherenceClaim"/>: the F88 closed form that produces
///         StaticFraction = 1/(2·C(N, n_p)). The dyadic-N reduction lives inside the
///         binomial.</item>
/// </list>
///
/// <para>Tier consistency: F88StaticDyadicAnchor is Tier1Derived; both parents are
/// Tier1Derived. TierStrength inheritance trivially passes (5 ≥ 5). The wiring exposes
/// the inheritance reading explicitly: one algebra propagates from d=2 (Pi2 root)
/// through the binomial structure into F88's static fraction.</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>,
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>,
/// <see cref="F88PopcountCoherenceRegistration.RegisterF88PopcountCoherence"/> in the
/// builder pipeline before this extension. The builder errors with <c>MissingParent</c>
/// if any are absent.</para></summary>
public static class F88StaticDyadicAnchorRegistration
{
    public static ClaimRegistryBuilder RegisterF88StaticDyadicAnchor(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F88StaticDyadicAnchor>(b =>
        {
            _ = b.Get<Pi2DyadicLadderClaim>();      // parent: the Pi2 ladder
            _ = b.Get<PopcountCoherenceClaim>();     // parent: F88 closed form
            return new F88StaticDyadicAnchor();
        });
}
