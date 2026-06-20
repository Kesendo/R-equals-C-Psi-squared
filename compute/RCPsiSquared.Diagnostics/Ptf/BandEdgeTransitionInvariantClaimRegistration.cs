using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>Wiring of <see cref="BandEdgeTransitionInvariantClaim"/> (F124, the band-edge transition
/// invariant ‖M‖_F² + λ_min(MMᵀ) = 2 with λ_min = E the Dirichlet-edge coupling). Two typed parents,
/// both Tier1Derived (so the child is Tier1Derived): <see cref="KPartnerSelectionRuleClaim"/> (the frame's
/// exact kernel IS the K-partner ψ_N; F124 is the same M completed with the strength column) and
/// <see cref="ClockHandLadderClaim"/> (the band edge E₁ = 2cos(π/(N+1)) the conserved envelope rides on,
/// which selects the carrier). The builder resolves both via <c>b.Get&lt;&gt;()</c> at Build() time, so
/// registration order does not matter.</summary>
public static class BandEdgeTransitionInvariantClaimRegistration
{
    public static ClaimRegistryBuilder RegisterBandEdgeTransitionInvariantClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<BandEdgeTransitionInvariantClaim>(b =>
            new BandEdgeTransitionInvariantClaim(
                b.Get<KPartnerSelectionRuleClaim>(),
                b.Get<ClockHandLadderClaim>()));
        return builder;
    }
}
