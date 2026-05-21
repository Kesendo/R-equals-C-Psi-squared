using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="IbmBlockCpsiHardwareTable"/>: the F86
/// Tier2Verified 32-row table of IBM 2026-04-26 framework-snapshot C_block readouts
/// (Aer / Marrakesh / Kingston / Fez) measured against Theorem 2's universal 1/4 ceiling.
/// Single parent edge:
///
/// <list type="bullet">
///   <item><see cref="QuarterAsBilinearMaxvalClaim"/>: the 1/4 = (1/2)² bilinear-apex
///         maxval, the Theorem 2 ceiling each table row is asserted against.</item>
/// </list>
///
/// <para>Tier consistency: <see cref="IbmBlockCpsiHardwareTable"/> is Tier2Verified under
/// the Tier1Derived <see cref="QuarterAsBilinearMaxvalClaim"/> parent.</para>
///
/// <para>Requires upstream registration: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers QuarterAsBilinearMaxvalClaim).</para></summary>
public static class IbmBlockCpsiHardwareTableRegistration
{
    public static ClaimRegistryBuilder RegisterIbmBlockCpsiHardwareTable(
        this ClaimRegistryBuilder builder) =>
        builder.Register<IbmBlockCpsiHardwareTable>(b =>
            new IbmBlockCpsiHardwareTable(b.Get<QuarterAsBilinearMaxvalClaim>()));
}
