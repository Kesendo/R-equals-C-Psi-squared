using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>Schicht-1 wiring of <see cref="KPartnerSelectionRuleClaim"/> (the Tier1Derived
/// K-partner selection rule ⟨ψ_N|V_b|ψ_1⟩ = 0). Single typed parent
/// <see cref="ChiralMirrorTrajectoryClaim"/> (Tier1Derived: both ingredients of the two-line
/// derivation come from it). Must register after <c>RegisterChiralMirrorTrajectoryClaim</c>.</summary>
public static class KPartnerSelectionRuleClaimRegistration
{
    public static ClaimRegistryBuilder RegisterKPartnerSelectionRuleClaim(
        this ClaimRegistryBuilder builder)
    {
        builder.Register<KPartnerSelectionRuleClaim>(b =>
            new KPartnerSelectionRuleClaim(b.Get<ChiralMirrorTrajectoryClaim>()));
        return builder;
    }
}
