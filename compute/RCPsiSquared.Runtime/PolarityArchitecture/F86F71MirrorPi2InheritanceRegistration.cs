using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F86F71MirrorPi2Inheritance"/>:
/// F86 Statement 3 sub-claim (F71 spatial-mirror invariance of Q_peak).
/// Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F71MirrorSymmetryPi2Inheritance"/>: F71 Pi2-Foundation
///         anchor (HalfIntegerMirror N-parity; bond-pair-count ⌊N/2⌋).</item>
///   <item><see cref="F86MirrorGeneralisationLink"/>: F71→F86 lift typed
///         claim (in F71 namespace, registered by RegisterF71Family).</item>
/// </list>
///
/// <para>Tier consistency: F86 Statement 3 is Tier 1 derived, bit-exactly
/// verified at c=2 N=5..7 and c=3 N=5..6. Both parent claims Tier1Derived.</para>
///
/// <para>Requires: <see cref="F71FamilyRegistration.RegisterF71Family"/>
/// (registers F86MirrorGeneralisationLink) +
/// <see cref="F71MirrorSymmetryPi2InheritanceRegistration.RegisterF71MirrorSymmetryPi2Inheritance"/>.</para></summary>
public static class F86F71MirrorPi2InheritanceRegistration
{
    public static ClaimRegistryBuilder RegisterF86F71MirrorPi2Inheritance(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F86F71MirrorPi2Inheritance>(b =>
        {
            var f71 = b.Get<F71MirrorSymmetryPi2Inheritance>();
            var f86Link = b.Get<F86MirrorGeneralisationLink>();
            return new F86F71MirrorPi2Inheritance(f71, f86Link);
        });
}
