using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Schicht-1 wiring of <see cref="SeedRungGramClaim"/> (F143).
///
/// <para>Two typed parents, both resolved from the builder: <see cref="ChiralKClaim"/> supplies the
/// involution R the kernel is read in, and <see cref="FrozenDivisorClaim"/> supplies the ⌊N/2⌋ count and
/// the modulus M on the site side, where PROOF_R90_FROZEN_DIVISOR Lemma 5 had the same matrix five days
/// before F143 was minted. The builder resolves topologically, so the order of the Register calls in
/// <c>KnowledgeRegistryFactory</c> is free.</para></summary>
public static class SeedRungGramClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSeedRungGramClaim(this ClaimRegistryBuilder builder) =>
        builder.Register<SeedRungGramClaim>(b =>
        {
            var chiralK = b.Get<ChiralKClaim>();
            var frozenDivisor = b.Get<FrozenDivisorClaim>();
            return new SeedRungGramClaim(chiralK, frozenDivisor);
        });
}
