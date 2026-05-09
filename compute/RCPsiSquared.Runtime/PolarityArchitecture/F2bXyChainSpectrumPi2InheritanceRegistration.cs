using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F2bXyChainSpectrumPi2Inheritance"/>:
/// F2b XY chain SE Hamiltonian E_k = 2J·cos(πk/(N+1)). Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c>
///         (HoppingFactor in 2J).</item>
///   <item><see cref="F65XxChainSpectrumPi2Inheritance"/>: F65 sibling edge.
///         F2b's |ψ_k(i)|² eigenvector amplitude squared IS F65's
///         BondingModePopulation; F2b gives the eigenvalues, F65 gives the
///         dissipative rates from the same OBC sine basis. EigenvectorMatchesF65
///         drift check verifies.</item>
/// </list>
///
/// <para>Tier consistency: F2b is Tier 1 proven (standard tight-binding theory
/// for OBC chains); both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F65XxChainSpectrumPi2InheritanceRegistration.RegisterF65XxChainSpectrumPi2Inheritance"/>.</para></summary>
public static class F2bXyChainSpectrumPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF2bXyChainSpectrumPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F2bXyChainSpectrumPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            var f65 = b.Get<F65XxChainSpectrumPi2Inheritance>();
            return new F2bXyChainSpectrumPi2Inheritance(ladder, f65);
        });
}
