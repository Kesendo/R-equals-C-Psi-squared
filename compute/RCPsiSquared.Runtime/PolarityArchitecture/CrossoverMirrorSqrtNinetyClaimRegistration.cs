using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="CrossoverMirrorSqrtNinetyClaim"/>:
/// the local crossover mirror (XZ+YZ / ZX+ZY) is the canonical Π turned by
/// √(<see cref="NinetyDegreeMirrorMemoryClaim"/>). One typed parent, the 90° angle-anchor,
/// of which this claim is the continuous square root (the per-site-conjugation face).
///
/// <para>Requires <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> (which registers
/// NinetyDegreeMirrorMemoryClaim).</para></summary>
public static class CrossoverMirrorSqrtNinetyClaimRegistration
{
    public static ClaimRegistryBuilder RegisterCrossoverMirrorSqrtNinetyClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<CrossoverMirrorSqrtNinetyClaim>(b =>
        {
            var ninetyDegree = b.Get<NinetyDegreeMirrorMemoryClaim>();
            return new CrossoverMirrorSqrtNinetyClaim(ninetyDegree);
        });
}
