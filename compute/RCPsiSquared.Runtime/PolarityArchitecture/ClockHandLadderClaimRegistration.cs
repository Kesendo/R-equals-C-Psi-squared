using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="ClockHandLadderClaim"/>: typed parent edges to
/// <see cref="F2bXyChainSpectrumPi2Inheritance"/> (the band edge the coherence hand is
/// γ-protected to), <see cref="AbsorptionTheoremClaim"/> (the −2γ that does the protecting),
/// and <see cref="UniversalCarrierClaim"/> (carrier-blindness, the γ-protection read from the
/// other side). Requires those three to be registered first.</summary>
public static class ClockHandLadderClaimRegistration
{
    public static ClaimRegistryBuilder RegisterClockHandLadderClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<ClockHandLadderClaim>(b =>
        {
            var bandEdge = b.Get<F2bXyChainSpectrumPi2Inheritance>();
            var absorption = b.Get<AbsorptionTheoremClaim>();
            var carrier = b.Get<UniversalCarrierClaim>();
            return new ClockHandLadderClaim(bandEdge, absorption, carrier);
        });
}
