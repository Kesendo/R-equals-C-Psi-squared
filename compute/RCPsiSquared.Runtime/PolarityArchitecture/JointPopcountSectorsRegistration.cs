using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="JointPopcountSectors"/>: the U(1)×U(1)
/// foundational block-decomposition. No typed parents — derived directly from
/// XY+Z-dephasing symmetry.</summary>
public static class JointPopcountSectorsRegistration
{
    public static ClaimRegistryBuilder RegisterJointPopcountSectors(
        this ClaimRegistryBuilder builder) =>
        builder.Register<JointPopcountSectors>(_ => new JointPopcountSectors());
}
