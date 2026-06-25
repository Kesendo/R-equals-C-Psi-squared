using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89MonodromyMirrorClaim"/> (the mirror splits at the Galois
/// boundary). Two typed parent edges: <see cref="F89OcticMonodromyClaim"/> (the S_8 monodromy this mirror
/// analysis acts on) and <see cref="F89BranchLocusPalindromeClaim"/> (the seams' position palindrome this
/// lifts to the braids, and whose Re = −4 fold cannot pass into the Galois action). Both parents Tier 1
/// derived. The builder topo-resolves, so registration order is free, but both parents must be
/// registered.</summary>
public static class F89MonodromyMirrorClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89MonodromyMirrorClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89MonodromyMirrorClaim>(b =>
            new F89MonodromyMirrorClaim(
                b.Get<F89OcticMonodromyClaim>(),
                b.Get<F89BranchLocusPalindromeClaim>()));
}
