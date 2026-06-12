using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="CoherenceHorizonClaim"/>: typed parent edges to
/// <see cref="ClockHandLadderClaim"/> (the horizon IS the clock's Q-floor made exact, the Q
/// below which the coherence hand freezes) and <see cref="F2bXyChainSpectrumPi2Inheritance"/>
/// (the band edge 2cos(π/(N+1)) the horizon coincides with at N=2,3). Requires both to be
/// registered first.</summary>
public static class CoherenceHorizonClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCoherenceHorizonClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<CoherenceHorizonClaim>(b =>
            new CoherenceHorizonClaim(
                b.Get<ClockHandLadderClaim>(),
                b.Get<F2bXyChainSpectrumPi2Inheritance>()));
}
