using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89Path2CardanoClaim"/>: the path-2
/// (SE, DE) S_2-symmetric closed-form characteristic polynomial factorisation.
/// One typed parent edge:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: orbit-closure framework
///         parent (path-2 is one of the pure-path-k topology classes).</item>
/// </list>
///
/// <para>Tier consistency: F89Path2Cardano is Tier 1 derived (sympy-verified
/// symbolic factorisation char(λ) = −(λ+2γ)(λ+6γ)·[cubic], with cubic solvable
/// in radicals via Cardano; numerical roots at q=1.5 verified bit-exactly
/// against populated path-2 fractional rates); parent F89 also Tier 1 derived.</para>
///
/// <para>Requires: <see cref="F89TopologyOrbitClosureRegistration.RegisterF89TopologyOrbitClosure"/>.</para></summary>
public static class F89Path2CardanoClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89Path2CardanoClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89Path2CardanoClaim>(b =>
        {
            var f89 = b.Get<F89TopologyOrbitClosure>();
            return new F89Path2CardanoClaim(f89);
        });
}
