using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>Schicht-1 wiring of <see cref="DefectReadingEquivarianceClaim"/> (M3, Tier1Derived,
/// the defect-reading spatial-reflection equivariance). Single load-bearing typed parent
/// <see cref="KPartnerSelectionRuleClaim"/> (it defines the location dictionary M[b,k], the
/// K-partner null column, and rank N−2). Must register after
/// <c>RegisterKPartnerSelectionRuleClaim</c>.</summary>
public static class DefectReadingEquivarianceClaimRegistration
{
    public static ClaimRegistryBuilder RegisterDefectReadingEquivarianceClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<DefectReadingEquivarianceClaim>(b =>
            new DefectReadingEquivarianceClaim(b.Get<KPartnerSelectionRuleClaim>()));
        return builder;
    }
}
