using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89PathKVacSeParsevalClaim"/>: the
/// (vac, SE) self-contribution closed form via Parseval orthogonality. One typed
/// parent edge:
///
/// <list type="bullet">
///   <item><see cref="F89TopologyOrbitClosure"/>: orbit-closure framework
///         parent (defines S_T(t) for ρ_cc + uniform-J multi-bond XY).</item>
/// </list>
///
/// <para>Tier consistency: F89PathKVacSeParseval is Tier 1 derived (Parseval
/// orthogonality of H_B^SE Bloch eigenstates + machine-precision verification
/// 4·10⁻¹⁷ to 6·10⁻¹⁶ across 15 (k, N) pairs); parent F89 also Tier 1 derived.</para>
///
/// <para>Requires: <see cref="F89TopologyOrbitClosureRegistration.RegisterF89TopologyOrbitClosure"/>.</para></summary>
public static class F89PathKVacSeParsevalClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89PathKVacSeParsevalClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89PathKVacSeParsevalClaim>(b =>
        {
            var f89 = b.Get<F89TopologyOrbitClosure>();
            return new F89PathKVacSeParsevalClaim(f89);
        });
}
