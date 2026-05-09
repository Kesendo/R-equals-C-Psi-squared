using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F84ThermalAmplitudeDampingPi2Inheritance"/>:
/// F84 thermal amplitude damping (cooling + heating). F84 → F82 mother-corollary.
/// Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         and <c>a_{2−N} = 2^(N−1)</c> via F82 delegation.</item>
///   <item><see cref="F82T1AmplitudeDampingPi2Inheritance"/>: mother claim.
///         F84 reduces to F82 at γ_↑ = 0 (vacuum bath); typed mother edge.</item>
/// </list>
///
/// <para>Tier consistency: F84 is Tier 1 proven (PROOF_F84_AMPLITUDE_DAMPING);
/// verified bit-exact at N=3 across 7 configurations + Pauli-channel
/// cancellation explicitly verified. All three claims Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F82T1AmplitudeDampingPi2InheritanceRegistration.RegisterF82T1AmplitudeDampingPi2Inheritance"/>
/// (which in turn requires the F81 chain).</para></summary>
public static class F84ThermalAmplitudeDampingPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF84ThermalAmplitudeDampingPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F84ThermalAmplitudeDampingPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f82 = b.Get<F82T1AmplitudeDampingPi2Inheritance>();
            return new F84ThermalAmplitudeDampingPi2Inheritance(ladder, f82);
        });
}
