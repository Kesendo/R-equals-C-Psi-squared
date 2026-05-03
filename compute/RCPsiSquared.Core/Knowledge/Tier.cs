namespace RCPsiSquared.Core.Knowledge;

/// <summary>Knowledge tier shared across all F-theorem typed knowledge bases. Mirrors the
/// framework-wide tier convention (see <c>hypotheses/</c> README and individual proof
/// docs). Used by <see cref="Claim"/> and every F-theorem KB.
///
/// <list type="bullet">
///   <item><see cref="Tier1Derived"/>: analytic proof, bit-exact verified.</item>
///   <item><see cref="Tier1Candidate"/>: strong numerical witness, no full algebraic
///         derivation yet; promotion to Tier 1 needs the missing piece.</item>
///   <item><see cref="Tier2Empirical"/>: empirical pattern across (c, N) without proof.</item>
///   <item><see cref="Tier2Verified"/>: hardware-confirmed prediction (Marrakesh / Kingston
///         entries in <c>ConfirmationsRegistry</c>).</item>
///   <item><see cref="OpenQuestion"/>: identified gap with a tractable approach.</item>
///   <item><see cref="Retracted"/>: claim was made, then refuted by extended-N data.</item>
/// </list>
/// </summary>
public enum Tier
{
    Tier1Derived,
    Tier1Candidate,
    Tier2Empirical,
    Tier2Verified,
    OpenQuestion,
    Retracted,
}

public static class TierExtensions
{
    public static string Label(this Tier tier) => tier switch
    {
        Tier.Tier1Derived => "Tier 1 (derived)",
        Tier.Tier1Candidate => "Tier 1 (candidate)",
        Tier.Tier2Empirical => "Tier 2 (empirical)",
        Tier.Tier2Verified => "Tier 2 (hardware-verified)",
        Tier.OpenQuestion => "Open question",
        Tier.Retracted => "Retracted",
        _ => "Unknown",
    };
}
