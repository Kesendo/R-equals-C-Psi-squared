using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="HardCellYInversionPattern"/>
/// (F110, 7th YParity-axis Claim). Standalone Claim: no ctor parents in this
/// registration extension.
///
/// <para><b>Layer-boundary note on the F87 dissipator-resonance law</b>:
/// `DissipatorResonanceLaw` (in `compute/RCPsiSquared.Diagnostics/F87/DissipatorResonanceLaw.cs`)
/// is the N = 4, k = 3 census of F110's Aspect A, not an input to it: Aspect A rests on
/// the colouring of the three non-diagonal Klein cells (PROOF_F110 §2), so no parent
/// edge to the census is called for. <c>RCPsiSquared.Runtime</c> does not reference
/// <c>RCPsiSquared.Diagnostics</c> in any case (per the PolarityCubeMap architectural
/// boundary).</para></summary>
public static class HardCellYInversionPatternRegistration
{
    public static ClaimRegistryBuilder RegisterHardCellYInversionPattern(
        this ClaimRegistryBuilder builder) =>
        builder.Register<HardCellYInversionPattern>(b =>
            new HardCellYInversionPattern(b.Get<KleinEightCellClaim>()));
}
