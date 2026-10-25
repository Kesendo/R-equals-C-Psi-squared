using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>Schicht-1 wiring of <see cref="DecoderAntiCollinearitySeesawClaim"/> (Tier2Empirical, the
/// painted dictionary's anti-collinearity read as the seesaw's). Single typed parent
/// <see cref="DefectReadingEquivarianceClaim"/> (the structural mirror-pair equivariance this is the
/// empirical Q-instance of). Must register after <c>RegisterDefectReadingEquivarianceClaim</c>.</summary>
public static class DecoderAntiCollinearitySeesawClaimRegistration
{
    public static ClaimRegistryBuilder RegisterDecoderAntiCollinearitySeesawClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<DecoderAntiCollinearitySeesawClaim>(b =>
            new DecoderAntiCollinearitySeesawClaim(b.Get<DefectReadingEquivarianceClaim>()));
        return builder;
    }
}
