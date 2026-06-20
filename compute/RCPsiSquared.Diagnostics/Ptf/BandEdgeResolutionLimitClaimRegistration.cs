using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>Wiring of <see cref="BandEdgeResolutionLimitClaim"/> (the optics/signal resolution-limit reading
/// of F124's conditioning: κ ~ N², contrast √κ ~ N, the staggered q=π the diffraction limit). Single typed
/// parent <see cref="BandEdgeTransitionInvariantClaim"/> (F124, Tier1Derived; every quantity is a corollary of
/// its transition matrix), resolved via <c>b.Get&lt;&gt;()</c>. Must register after
/// <c>RegisterBandEdgeTransitionInvariantClaim</c> (order does not matter; the builder resolves topologically).</summary>
public static class BandEdgeResolutionLimitClaimRegistration
{
    public static ClaimRegistryBuilder RegisterBandEdgeResolutionLimitClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<BandEdgeResolutionLimitClaim>(b =>
            new BandEdgeResolutionLimitClaim(b.Get<BandEdgeTransitionInvariantClaim>()));
        return builder;
    }
}
