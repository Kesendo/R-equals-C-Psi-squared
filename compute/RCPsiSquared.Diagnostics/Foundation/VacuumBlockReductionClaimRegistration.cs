using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Schicht-1 wiring of <see cref="VacuumBlockReductionClaim"/> (the N=5 birth-canal
/// boundary = the |1-exc&gt;&lt;vac| (0,1) Liouville sector block, Tier1Derived). Single typed
/// parent <see cref="AbsorptionTheoremClaim"/> (the rate = −2·Σ_l γ_l·⟨Δ_l⟩ law the block's
/// −2·diag(γ) decay and flat-γ blindness rest on; the (0,1)-block rate is that law restricted to a
/// single conserved sector). Must register after <c>RegisterAbsorptionTheoremClaim</c>.</summary>
public static class VacuumBlockReductionClaimRegistration
{
    public static ClaimRegistryBuilder RegisterVacuumBlockReductionClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<VacuumBlockReductionClaim>(b =>
        {
            var absorption = b.Get<AbsorptionTheoremClaim>();   // typed parent edge
            return new VacuumBlockReductionClaim(absorption);
        });
}
