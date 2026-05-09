using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F69GhzWSexticAboveFoldPi2Inheritance"/>:
/// F69 GHZ+W mix lifts pair-CΨ(0) above the 1/4 fold at N=3 (sextic minimal
/// polynomial, irreducible over ℚ). Three typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_3 = 1/4</c>
///         (the bilinear-apex fold; same anchor F60 and F62 reference as
///         their FoldPosition).</item>
///   <item><see cref="F60GhzBornBelowFoldPi2Inheritance"/>: GHZ baseline
///         (CΨ(0) = 1/(2³−1) = 1/7 ≈ 0.143 at N=3, strictly below fold).</item>
///   <item><see cref="F62WStateBornBelowFoldPi2Inheritance"/>: W baseline
///         (CΨ(0) = 10/81 ≈ 0.1235 at N=3, strictly below fold).</item>
/// </list>
///
/// <para>F69's structural finding is that the mix of two below-fold pure states
/// can rise above the fold via interference (L1-off √6·α·√(1-α²) cross term).
/// The mechanism is N=3-specific; at N ≥ 4 the GHZ-purity scales as 1/(2^N − 1)
/// and collapses too fast for the W-contribution to lift above 1/4.</para>
///
/// <para>Tier consistency: F69 is Tier 1 — sextic minimal polynomial irreducible
/// over ℚ; numerical optimum at 25-digit precision; α²_opt is a degree-6
/// algebraic number with no radical form. All three parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F60GhzBornBelowFoldPi2InheritanceRegistration.RegisterF60GhzBornBelowFoldPi2Inheritance"/> +
/// <see cref="F62WStateBornBelowFoldPi2InheritanceRegistration.RegisterF62WStateBornBelowFoldPi2Inheritance"/>.</para></summary>
public static class F69GhzWSexticAboveFoldPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF69GhzWSexticAboveFoldPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F69GhzWSexticAboveFoldPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f60 = b.Get<F60GhzBornBelowFoldPi2Inheritance>();
            var f62 = b.Get<F62WStateBornBelowFoldPi2Inheritance>();
            return new F69GhzWSexticAboveFoldPi2Inheritance(ladder, f60, f62);
        });
}
