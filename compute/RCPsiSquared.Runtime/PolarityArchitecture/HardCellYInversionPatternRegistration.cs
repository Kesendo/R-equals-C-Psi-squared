using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 registration for <see cref="HardCellYInversionPattern"/>
/// (F110, 7th YParity-axis Claim). Standalone Claim: no ctor parents. Closes
/// the third F87 trichotomy slot (after F107 truly + F109 mother soft) on the
/// YParity axis, completing the three-axis y_par classification of F87's
/// trichotomy classes.</summary>
public static class HardCellYInversionPatternRegistration
{
    public static ClaimRegistryBuilder RegisterHardCellYInversionPattern(
        this ClaimRegistryBuilder builder) =>
        builder.Register<HardCellYInversionPattern>(_ => new HardCellYInversionPattern());
}
