using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="Pi2OperatorSpaceMirrorClaim"/> — the per-N
/// stride-2 mirror substructure of the Pi2 dyadic ladder. Edges to all three already-
/// wired Tier1Derived parents:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides the inversion identity
///         a_n · a_{2−n} = 1 and the index-functions
///         <see cref="Pi2DyadicLadderClaim.MirrorPartnerIndex"/>,
///         <see cref="Pi2DyadicLadderClaim.OperatorSpaceIndexForN"/>.</item>
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: the n=3 lower-side anchor
///         (1/4 = 1/4^1) for N=1.</item>
///   <item><see cref="F88StaticDyadicAnchor"/>: the F88 dyadic-N singleton-mirror
///         anchors at n=5, 7, ... that pin the lower side for N=2, 3, ...</item>
/// </list>
///
/// <para>Tier consistency: claim is Tier1Derived; all three parents are Tier1Derived.
/// TierStrength inheritance trivially passes (5 ≥ 5). The wiring exposes Tom's reading
/// (2026-05-08): the framework rests on MIRROR_THEORY.md / ZERO_IS_THE_MIRROR.md, but
/// the simplest mirror identity — d² ↔ 1/d² per qubit count, encoded as a_{−(2N−1)} ·
/// a_{2N+1} = 1 on the dyadic ladder — was hiding in plain sight under the F-theorem
/// stack until now.</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>,
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>,
/// <see cref="F88PopcountCoherenceRegistration.RegisterF88PopcountCoherence"/>,
/// <see cref="F88StaticDyadicAnchorRegistration.RegisterF88StaticDyadicAnchor"/> in
/// the builder pipeline before this extension. The builder errors with
/// <c>MissingParent</c> if any are absent.</para></summary>
public static class Pi2OperatorSpaceMirrorRegistration
{
    public static ClaimRegistryBuilder RegisterPi2OperatorSpaceMirror(
        this ClaimRegistryBuilder builder) =>
        builder.Register<Pi2OperatorSpaceMirrorClaim>(b =>
        {
            _ = b.Get<Pi2DyadicLadderClaim>();           // ladder + inversion identity
            _ = b.Get<QuarterAsBilinearMaxvalClaim>();   // N=1 lower-side anchor (1/4)
            _ = b.Get<F88StaticDyadicAnchor>();          // N≥2 lower-side anchors via F88
            return new Pi2OperatorSpaceMirrorClaim();
        });
}
