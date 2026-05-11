using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89AdditiveIdentityClaim"/>: the
/// mixed-topology additive identity from Lindbladian factorisation. One typed
/// parent edge:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: orbit-closure framework
///         parent (defines S_T(t) for ρ_cc + uniform-J multi-bond XY).</item>
/// </list>
///
/// <para>Tier consistency: F89AdditiveIdentity is Tier 1 derived (Lindbladian
/// factorisation + 27/27 N=7 CSV verification at 5.013·10⁻⁷ precision floor);
/// parent F89 also Tier 1 derived.</para>
///
/// <para>Requires: <see cref="F89TopologyOrbitClosureRegistration.RegisterF89TopologyOrbitClosure"/>.</para></summary>
public static class F89AdditiveIdentityClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89AdditiveIdentityClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89AdditiveIdentityClaim>(b =>
        {
            var f89 = b.Get<F89TopologyOrbitClosure>();
            return new F89AdditiveIdentityClaim(f89);
        });
}
