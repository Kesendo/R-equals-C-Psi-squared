using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89TopologyOrbitClosure"/>:
/// F89's S_N-orbit closure for ρ_cc spatial-sum coherence under uniform-J
/// multi-bond XY. Two cited parents:
///
/// <list type="bullet">
///   <item><see cref="F73SpatialSumPurityClosurePi2Inheritance"/>: closely
///         related closure for the (vac, SE) coherence block. Same orbit-style
///         argument; F73 yields a closed exponential thanks to a uniform per-
///         element 2γ₀ rate, while F89's (S_1, S_2) block has non-uniform per-
///         element decay (rate 2γ₀ on overlap, 6γ₀ off overlap), so F89 fixes
///         only the orbit-position dependence, not the time dependence.</item>
///   <item><see cref="F71MirrorSymmetryPi2Inheritance"/>: F71's spatial-mirror
///         Z₂ b ↔ N−2−b is a subgroup of the full S_N used in F89.</item>
/// </list>
///
/// <para>Tier consistency: F89, F73, and F71 are all Tier1Derived; no tier downgrade at injection.</para>
///
/// <para>Requires:
/// <see cref="F73SpatialSumPurityClosurePi2InheritanceRegistration.RegisterF73SpatialSumPurityClosurePi2Inheritance"/> +
/// <see cref="F71MirrorSymmetryPi2InheritanceRegistration.RegisterF71MirrorSymmetryPi2Inheritance"/> must
/// have been called first so the parent claims are present in the registry.</para></summary>
public static class F89TopologyOrbitClosureRegistration
{
    public static ClaimRegistryBuilder RegisterF89TopologyOrbitClosure(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89TopologyOrbitClosure>(b =>
        {
            var f73 = b.Get<F73SpatialSumPurityClosurePi2Inheritance>();
            var f71 = b.Get<F71MirrorSymmetryPi2Inheritance>();
            return new F89TopologyOrbitClosure(f73, f71);
        });
}
