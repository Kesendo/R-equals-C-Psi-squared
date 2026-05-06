using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>Explicit strength ranking for <see cref="Tier"/>. The existing enum integer
/// order has Tier2Empirical (= 2) below Tier2Verified (= 3), which would make Tier2Verified
/// numerically "greater" but they are unordered as labels: Verified is hardware-confirmed
/// (stronger), Empirical is a numerical pattern without proof. This helper imposes the
/// strength order used for parent-child Tier inheritance validation.</summary>
public static class TierStrength
{
    public static int Of(Tier tier) => tier switch
    {
        Tier.Tier1Derived => 5,
        Tier.Tier1Candidate => 4,
        Tier.Tier2Verified => 3,
        Tier.Tier2Empirical => 2,
        Tier.OpenQuestion => 1,
        Tier.Retracted => 0,
        _ => throw new ArgumentOutOfRangeException(nameof(tier), $"unknown Tier {tier}"),
    };

    /// <summary>True when <paramref name="parent"/> is at least as strong as <paramref name="child"/>.
    /// Used by the registration-time Tier inheritance check: a child Claim cannot rest on a
    /// parent Claim of weaker strength.</summary>
    public static bool IsAtLeastAsStrong(Tier parent, Tier child) =>
        Of(parent) >= Of(child);
}
