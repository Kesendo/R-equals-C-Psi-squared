using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Wires the N = 7 missing-phase relaxation-scale germ to the three objects it
/// actually composes. The absorption theorem and F1 palindrome remain ancestors through the
/// composite <see cref="JDefectLightMigrationClaim"/>; this local result spends its absorption/light
/// half, while the F1 half is inherited graph context rather than a separate direct edge. The
/// vacuum-block reduction is not a premise.</summary>
public static class MissingPhaseRelaxationScaleClaimRegistration
{
    public static ClaimRegistryBuilder RegisterMissingPhaseRelaxationScaleClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<MissingPhaseRelaxationScaleClaim>(context =>
            new MissingPhaseRelaxationScaleClaim(
                context.Get<JDefectLightMigrationClaim>(),
                context.Get<SeatCutBlindnessClaim>(),
                context.Get<JointPopcountSectors>()));
}
