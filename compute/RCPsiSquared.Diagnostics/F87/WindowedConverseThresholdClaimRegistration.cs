using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="WindowedConverseThresholdClaim"/> (the Tier1Derived
/// two-reflection spine). Single typed parent <see cref="F87DiagonalCellBipartiteWitnessSet"/>
/// (Tier1Derived). Must register after the witness set and before
/// <see cref="WindowedConverseAllGammaClaimRegistration.RegisterWindowedConverseAllGammaClaim"/>
/// (which now takes this claim as a parent).</summary>
public static class WindowedConverseThresholdClaimRegistration
{
    public static ClaimRegistryBuilder RegisterWindowedConverseThresholdClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<WindowedConverseThresholdClaim>(b =>
        {
            _ = b.Get<F87DiagonalCellBipartiteWitnessSet>(); // typed parent edge (parent before child)
            return new WindowedConverseThresholdClaim();
        });
        return builder;
    }
}
