using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89TopologyOrbitClosure"/>.
///
/// <para>Cited parent edges (resolved at build time so the runtime ancestor walk
/// can reach them; see <see cref="F89TopologyOrbitClosure"/> for the full F73/F71
/// contrast):
/// <see cref="F73SpatialSumPurityClosurePi2Inheritance"/> (closely related closure
/// for the (vac, SE) coherence block, same orbit-style argument) and
/// <see cref="F71MirrorSymmetryPi2Inheritance"/> (spatial-Z₂ subgroup of S_N).</para>
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
            _ = b.Get<F73SpatialSumPurityClosurePi2Inheritance>();
            _ = b.Get<F71MirrorSymmetryPi2Inheritance>();
            return new F89TopologyOrbitClosure();
        });
}
