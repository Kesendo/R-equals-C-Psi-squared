using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wires the node-pair theorem to F157, its only typed premise. F157 owns the
/// zero-free Jacobi node lemma and blindness as a Krylov complement; the new claim adds the
/// reduced-resolvent zero and the sufficient non-incident-bond consequence.</summary>
public static class NodePairResolventClaimRegistration
{
    public static ClaimRegistryBuilder RegisterNodePairResolventClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<NodePairResolventClaim>(context =>
            new NodePairResolventClaim(context.Get<SeatCutBlindnessClaim>()));
}
