using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89TopologyOrbitClosure"/>.
///
/// <para>Parent edges:
/// <see cref="Pi2DyadicLadderClaim"/> (algebraic parent, supplies a_{−1}=4 for
/// both the γ-axis decay coefficient and the J-axis oscillation coefficient in
/// the all-isolated subclass closed form);
/// <see cref="F73SpatialSumPurityClosurePi2Inheritance"/> (closely related closure
/// for the (vac, SE) coherence block, same orbit-style argument and same a_{−1}
/// decay anchor); and
/// <see cref="F71MirrorSymmetryPi2Inheritance"/> (spatial-Z₂ subgroup of S_N).</para>
///
/// <para>Tier consistency: F89, Pi2DyadicLadder, F73, and F71 are all Tier1Derived;
/// no tier downgrade at injection.</para>
///
/// <para>Requires:
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="F73SpatialSumPurityClosurePi2InheritanceRegistration.RegisterF73SpatialSumPurityClosurePi2Inheritance"/> +
/// <see cref="F71MirrorSymmetryPi2InheritanceRegistration.RegisterF71MirrorSymmetryPi2Inheritance"/> must
/// have been called first so the parent claims are present in the registry.</para></summary>
public static class F89TopologyOrbitClosureRegistration
{
    public static ClaimRegistryBuilder RegisterF89TopologyOrbitClosure(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89TopologyOrbitClosure>(b =>
        {
            var ladder = b.Get<Pi2DyadicLadderClaim>();
            _ = b.Get<F73SpatialSumPurityClosurePi2Inheritance>();
            _ = b.Get<F71MirrorSymmetryPi2Inheritance>();
            return new F89TopologyOrbitClosure(ladder);
        });
}
