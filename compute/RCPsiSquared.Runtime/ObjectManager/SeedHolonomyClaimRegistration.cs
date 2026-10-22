using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Runtime.ObjectManager;

/// <summary>Wiring of <see cref="SeedHolonomyClaim"/> (2026-07-07): the eigenvector holonomy around the
/// (1,2)-block defective seed is the mod-4 memory loop i⁴=1 (M₁ eigenvalues ±i, M₂=−I, M₄=I). One typed
/// parent, <see cref="TwoByTwoExceptionalPointClaim"/>, for one part of its content: the seed toy [[s, 1], [1, −s]]
/// at s = i is a setting of the one closed 2 × 2 (the related structures, the eigenVALUE-swap Monodromy, the algebraic
/// Pi2I4MemoryLoop/NinetyDegreeMirrorMemory, and F86_EP_THROUGH_THE_CLOCK, are named in the claim's anchor
/// and prose children). Requires RegisterTwoByTwoExceptionalPointClaim. Live witness: SeedHolonomyWitness (inspect --root holonomy).</summary>
public static class SeedHolonomyClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSeedHolonomy(this ClaimRegistryBuilder builder) =>
        builder.Register<SeedHolonomyClaim>(b =>
        {
            _ = b.Get<TwoByTwoExceptionalPointClaim>();   // partial scope: the seed toy [[s, 1], [1, −s]] at s = i is a setting of the one closed block
            return new SeedHolonomyClaim();
        });
}
