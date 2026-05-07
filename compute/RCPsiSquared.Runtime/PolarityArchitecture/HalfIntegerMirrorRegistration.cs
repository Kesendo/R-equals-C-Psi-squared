using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Registers <see cref="HalfIntegerMirrorClaim"/> at a specific chain length N.
/// The claim classifies w_XY = N/2 as half-integer (odd N) or integer (even N), inheriting
/// from <see cref="QubitDimensionalAnchorClaim"/> as one more concrete face of the 1/d
/// fixed point.
///
/// <para>Separate from <c>RegisterPi2Family</c> because it is the only Pi2 claim with an
/// N parameter; bundling it would force the family registration into a parameterised API.
/// Keep this opt-in so consumers who do not have a chain length stay parameter-free, and
/// the F1 / F86 families that already carry N can chain it on top.</para></summary>
public static class HalfIntegerMirrorRegistration
{
    public static ClaimRegistryBuilder RegisterHalfIntegerMirror(
        this ClaimRegistryBuilder builder,
        int N) =>
        builder.Register<HalfIntegerMirrorClaim>(b =>
        {
            _ = b.Get<QubitDimensionalAnchorClaim>();
            return new HalfIntegerMirrorClaim(N);
        });
}
