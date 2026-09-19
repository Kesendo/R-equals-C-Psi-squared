using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of the historically named <see cref="CpsiEnvelopeTheoremClaim"/> boundary with two
/// typed parent edges resolved from the registry — <see cref="F25CPsiBellPlusPi2Inheritance"/> (the
/// named Bell+/Z formula) and <see cref="QuarterAsBilinearMaxvalClaim"/> (the algebraic ¼ maxval).
/// Requires both parent registrations upstream.</summary>
public static class CpsiEnvelopeTheoremClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCpsiEnvelopeTheoremClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<CpsiEnvelopeTheoremClaim>(b =>
        {
            var f25 = b.Get<F25CPsiBellPlusPi2Inheritance>();
            var quarter = b.Get<QuarterAsBilinearMaxvalClaim>();
            return new CpsiEnvelopeTheoremClaim(f25, quarter);
        });
}
