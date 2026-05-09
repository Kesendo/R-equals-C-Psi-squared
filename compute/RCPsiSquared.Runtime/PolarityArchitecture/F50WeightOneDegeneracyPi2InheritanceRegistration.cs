using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F50WeightOneDegeneracyPi2Inheritance"/>:
/// F50's d_real(Re = −2γ) = 2N counting result. One typed parent edge:
///
/// <list type="bullet">
///   <item><see cref="Pi2DyadicLadderClaim"/>: provides <c>a_0 = 2</c> twice —
///         as the count multiplier (2 active Pauli letters X, Y) and as the
///         decay rate (Re = −2γ). The same a_0 anchor gives both pieces of
///         F50's formula.</item>
/// </list>
///
/// <para>Tier consistency: F50 is Tier 1 proven (lower + upper bound via SWAP
/// invariance, triangle inequality, and adjacent-transposition generation of
/// S_N); numerically verified N=2..7 in
/// <c>docs/proofs/PROOF_WEIGHT1_DEGENERACY.md</c>. Topology-universal for
/// isotropic Heisenberg (chain, star, ring, complete, tree).</para>
///
/// <para>Requires: <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/>.</para></summary>
public static class F50WeightOneDegeneracyPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF50WeightOneDegeneracyPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F50WeightOneDegeneracyPi2Inheritance>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            return new F50WeightOneDegeneracyPi2Inheritance(ladder);
        });
}
