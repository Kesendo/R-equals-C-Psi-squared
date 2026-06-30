using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>Schicht-1 wiring of <see cref="DecoderAntiCollinearityParityClaim"/> (Tier2Empirical, the
/// painted-dictionary odd-N parity refinement). Single typed parent
/// <see cref="DefectReadingEquivarianceClaim"/> (the structural mirror-pair equivariance this is the
/// empirical Q-instance of). Must register after <c>RegisterDefectReadingEquivarianceClaim</c>.</summary>
public static class DecoderAntiCollinearityParityClaimRegistration
{
    public static ClaimRegistryBuilder RegisterDecoderAntiCollinearityParityClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<DecoderAntiCollinearityParityClaim>(b =>
            new DecoderAntiCollinearityParityClaim(b.Get<DefectReadingEquivarianceClaim>()));
        return builder;
    }
}
