using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="NivenRationalityRootClaim"/>: two typed parent edges, one per single-
/// excitation face whose arithmetic the Niven root reads — <see cref="TopologyBandEdgeClaim"/> (the IM-face
/// band edge 2cos(π/(N+1))) and <see cref="F65XxChainSpectrumPi2Inheritance"/> (the RE-face dissipator rates
/// α_k, the documented F65/F99 home of the rate-side Niven fact). Both are Tier1Derived, so the child stays
/// Tier1Derived (pure, sympy-proven number theory). Requires both registered first.</summary>
public static class NivenRationalityRootClaimRegistration
{
    public static ClaimRegistryBuilder RegisterNivenRationalityRootClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<NivenRationalityRootClaim>(b =>
            new NivenRationalityRootClaim(
                b.Get<TopologyBandEdgeClaim>(),
                b.Get<F65XxChainSpectrumPi2Inheritance>()));
}
