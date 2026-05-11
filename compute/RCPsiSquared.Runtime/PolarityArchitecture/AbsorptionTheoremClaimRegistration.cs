using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="AbsorptionTheoremClaim"/>. Single
/// parent edge to <see cref="Pi2DyadicLadderClaim"/> for the absorption-quantum
/// coefficient a_0 = 2. Children register their own back-edge via discard-Get.
///
/// <para>Requires:
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class AbsorptionTheoremClaimRegistration
{
    public static ClaimRegistryBuilder RegisterAbsorptionTheoremClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<AbsorptionTheoremClaim>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new AbsorptionTheoremClaim(ladder);
        });
}
