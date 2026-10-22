using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="CoherenceHorizonClaim"/>: typed parent edges to
/// <see cref="ClockHandLadderClaim"/> (the horizon IS the clock's Q-floor made exact, the Q
/// below which the coherence hand freezes) and <see cref="F2bXyChainSpectrumPi2Inheritance"/>
/// (the band edge 2cos(π/(N+1)) the horizon coincides with at N=2,3). Requires both to be
/// registered first. A third edge, to <see cref="TwoByTwoExceptionalPointClaim"/>: at N = 2, 3 the horizon is the
/// clean closed 2 × 2, Q* = 2/√κ_N its meeting d = 0.</summary>
public static class CoherenceHorizonClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCoherenceHorizonClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<CoherenceHorizonClaim>(b =>
        {
            _ = b.Get<TwoByTwoExceptionalPointClaim>();
            return new CoherenceHorizonClaim(
                b.Get<ClockHandLadderClaim>(),
                b.Get<F2bXyChainSpectrumPi2Inheritance>());
        });
}
