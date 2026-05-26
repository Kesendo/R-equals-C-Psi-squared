using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="HardCellYInversionPattern"/>
/// (F110, 7th YParity-axis Claim). Standalone Claim: no ctor parents in this
/// registration extension.
///
/// <para><b>Layer-boundary note on the F87 dissipator-resonance typed parent
/// edge</b>: F110's Aspect A derivation depends on `DissipatorResonanceLaw`
/// (Tier1Derived, in `compute/RCPsiSquared.Diagnostics/F87/DissipatorResonanceLaw.cs`)
/// as a structural input. The parent edge is NOT wired here because
/// `RCPsiSquared.Runtime` does not reference `RCPsiSquared.Diagnostics` (per the
/// PolarityCubeMap architectural boundary). Moving this Registration to
/// `RCPsiSquared.Diagnostics` would let us add `b.Get&lt;DissipatorResonanceLaw&gt;()`
/// as a typed parent edge, mirroring the
/// `DissipatorAxisSelectsPolarityClaim` pattern (see
/// `F87FamilyRegistration.cs:47-52`). Deferred as a separate Schicht-relocation
/// pass; the proof and Claim docstring cite DissipatorResonanceLaw explicitly so
/// the dependency is recorded even though the inheritance-graph edge is
/// untyped.</para></summary>
public static class HardCellYInversionPatternRegistration
{
    public static ClaimRegistryBuilder RegisterHardCellYInversionPattern(
        this ClaimRegistryBuilder builder) =>
        builder.Register<HardCellYInversionPattern>(b =>
            new HardCellYInversionPattern(b.Get<KleinEightCellClaim>()));
}
