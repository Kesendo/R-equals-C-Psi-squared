using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>Schicht-1 wiring of <see cref="AntiFractionObstructionOrthogonalityClaim"/> (F81/F115
/// connector). Two typed parents: <see cref="F83AntiFractionPi2Inheritance"/> (Core, the bit_b
/// anti-fraction) and <see cref="WindowedHardnessClaim"/> (Diagnostics, the bit_a obstruction). Both
/// Tier1Derived, so the connector is Tier1Derived (5 >= 5). Requires both parents' registrations to
/// have run (the builder errors with MissingParent otherwise).</summary>
public static class AntiFractionObstructionOrthogonalityRegistration
{
    public static ClaimRegistryBuilder RegisterAntiFractionObstructionOrthogonalityClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<AntiFractionObstructionOrthogonalityClaim>(b =>
            new AntiFractionObstructionOrthogonalityClaim(
                b.Get<F83AntiFractionPi2Inheritance>(),
                b.Get<WindowedHardnessClaim>()));
}
