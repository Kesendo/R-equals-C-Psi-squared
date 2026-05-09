using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F67BondingBellPairPi2Inheritance"/>:
/// F67 bonding-mode-encoded Bell pair as the slowest-mode endpoint of F65's
/// single-excitation spectrum. Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_{−1} = 4</c>
///         (numerator coefficient in α_1 = (4γ₀/(N+1))·sin²(π/(N+1));
///         F25/F65/F73/F76 sibling).</item>
///   <item><see cref="F65XxChainSpectrumPi2Inheritance"/>: F65 source-claim
///         edge — F67 delegates the α_1 closed form (single-excitation rate
///         at k=1) and the bonding-mode amplitude |ψ_1(j)|² to F65, rather
///         than duplicating the formulas.</item>
/// </list>
///
/// <para>Tier consistency: F67 is Tier 1 proven (Absorption Theorem applied
/// to k=1 bonding mode); verified N=3, N=5 in
/// <c>simulations/bell_pair_chain_protection.py</c> with Variant B
/// α_fit/α_1 = 0.9989 (N=3), 0.9963 (N=5). Both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F65XxChainSpectrumPi2InheritanceRegistration.RegisterF65XxChainSpectrumPi2Inheritance"/>
/// (which transitively requires F66 and QubitDimensionalAnchor).</para></summary>
public static class F67BondingBellPairPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF67BondingBellPairPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F67BondingBellPairPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f65 = b.Get<F65XxChainSpectrumPi2Inheritance>();
            return new F67BondingBellPairPi2Inheritance(ladder, f65);
        });
}
