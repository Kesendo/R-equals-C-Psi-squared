using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89PathKHbMixedDegreesClaim"/>:
/// the path-k H_B-mixed sub-factor degree table {8, 18, 32, 53} for paths
/// {3, 4, 5, 6}. One typed parent edge:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: orbit-closure framework parent
///         (the H_B-mixed degrees are a combinatorial property of the path-k
///         orbit class).</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (combinatorial + verified). Parent
/// F89 also Tier 1 derived.</para>
///
/// <para>Requires: <see cref="F89TopologyOrbitClosureRegistration.RegisterF89TopologyOrbitClosure"/>.</para></summary>
public static class F89PathKHbMixedDegreesClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89PathKHbMixedDegreesClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89PathKHbMixedDegreesClaim>(b =>
        {
            var f89 = b.Get<F89TopologyOrbitClosure>();
            return new F89PathKHbMixedDegreesClaim(f89);
        });
}
