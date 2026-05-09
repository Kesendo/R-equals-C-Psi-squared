using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F79TwoBodyPi2BlockPi2Inheritance"/>:
/// F79 2-body Π²-block decomposition. Three typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="KleinFourCellClaim"/>: F88 4-cell (Π²_Z, Π²_X) refines
///         F79's binary p-classification. F79's Π²-odd splits into Mp/Mm
///         under Π²_X axis.</item>
///   <item><see cref="Pi2OperatorSpaceMirrorClaim"/>: provides 4^N operator
///         space dimension; V_+/V_- have equal dim 4^N/2.</item>
///   <item><see cref="F1Pi2Inheritance"/>: M IS F1's residual operator;
///         F79 reads its 2-body Π²-block structure.</item>
/// </list>
///
/// <para>Tier consistency: F79 is Tier 1 proven (joint analysis with F78);
/// verified N=3..5 across chain, star, disjoint topologies.</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/>
/// (registers KleinFourCellClaim) +
/// <see cref="Pi2OperatorSpaceMirrorRegistration.RegisterPi2OperatorSpaceMirror"/> +
/// <see cref="F1Pi2InheritanceRegistration.RegisterF1Pi2Inheritance"/>.</para></summary>
public static class F79TwoBodyPi2BlockPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF79TwoBodyPi2BlockPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F79TwoBodyPi2BlockPi2Inheritance>(b =>
        {
            var klein = b.Get<KleinFourCellClaim>();
            var mirror = b.Get<Pi2OperatorSpaceMirrorClaim>();
            var f1 = b.Get<F1Pi2Inheritance>();
            return new F79TwoBodyPi2BlockPi2Inheritance(klein, mirror, f1);
        });
}
