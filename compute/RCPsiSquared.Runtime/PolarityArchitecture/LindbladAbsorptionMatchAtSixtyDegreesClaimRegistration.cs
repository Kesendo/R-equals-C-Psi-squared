using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>
/// Wires the named two-level Q=sqrt(3) cross-reading to its sole genuine parent, F95.
/// The literal magnitude ratio 2 and the angle pi/3 do not create AbsorptionTheorem,
/// canonical-angle, or Pi2-ladder ancestry.
/// </summary>
public static class LindbladAbsorptionMatchAtSixtyDegreesClaimRegistration
{
    public static ClaimRegistryBuilder RegisterLindbladAbsorptionMatchAtSixtyDegreesClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<LindbladAbsorptionMatchAtSixtyDegreesClaim>(b =>
            new LindbladAbsorptionMatchAtSixtyDegreesClaim(
                b.Get<F95AngleAtQuadraticZeroPi2Inheritance>()));
}
