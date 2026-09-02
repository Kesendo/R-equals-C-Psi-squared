using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wiring of <see cref="CrackedRingExactCurveClaim"/> (F160; Tier1Derived). Two typed parents:
/// <see cref="TopologyBandEdgeClaim"/>, whose ring row is the curve's u = 1 end and whose §Scope fence
/// (PROOF_RING_GAP_DOMINANCE, the wrap bond detuned by 1e-4) is this curve's Perron root; and
/// <see cref="F2bXyChainSpectrumPi2Inheritance"/>, the open chain comb the curve reaches at u = 0. Both are
/// Tier1Derived, so the tier rule holds. Resolution is topological, so this registration may sit anywhere in
/// the chain after the two families that supply the parents.</summary>
public static class CrackedRingExactCurveClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCrackedRingExactCurveClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<CrackedRingExactCurveClaim>(b =>
            new CrackedRingExactCurveClaim(b.Get<TopologyBandEdgeClaim>(), b.Get<F2bXyChainSpectrumPi2Inheritance>()));
}
