using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="WindowedConverseAllGammaClaim"/> (the F87 windowed-converse
/// all-γ residual lemma). Two Tier1Derived typed parents: <see cref="F87DiagonalCellBipartiteWitnessSet"/>
/// (the F103 §7 criterion + genericity carrier) and <see cref="WindowedConverseThresholdClaim"/> (the
/// proven Phase B two-reflection spine), which together leave this Tier1Candidate lemma proven modulo
/// the two residuals R-deg + R-sign.
///
/// <para>Requires both <see cref="F87DiagonalCellBipartiteWitnessSetRegistration.RegisterF87DiagonalCellBipartiteWitnessSet"/>
/// and <see cref="WindowedConverseThresholdClaimRegistration.RegisterWindowedConverseThresholdClaim"/>
/// to run first (parents before child).</para>
///
/// <para>Tier consistency: WindowedConverseAllGammaClaim is Tier1Candidate ← both parents Tier1Derived
/// (Derived ≥ Candidate, the registry's parent ≥ child rule).</para></summary>
public static class WindowedConverseAllGammaClaimRegistration
{
    public static ClaimRegistryBuilder RegisterWindowedConverseAllGammaClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<WindowedConverseAllGammaClaim>(b =>
        {
            _ = b.Get<F87DiagonalCellBipartiteWitnessSet>();   // genericity carrier (Tier1Derived)
            _ = b.Get<WindowedConverseThresholdClaim>();        // the proven Phase B spine (Tier1Derived)
            return new WindowedConverseAllGammaClaim();
        });
        return builder;
    }
}
