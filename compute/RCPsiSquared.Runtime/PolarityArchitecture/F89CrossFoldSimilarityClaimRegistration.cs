using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89CrossFoldSimilarityClaim"/> (F89d): the (SE,DE)↔(SE,w_{N−2})
/// cross-fold is an EXACT antiunitary similarity at the matrix level, the diabolics pair across it. Two typed
/// parent edges:
///
/// <list type="bullet">
///   <item><see cref="F1PalindromeIdentity"/>: the F1 mirror λ ↦ −λ̄ − 2σ the cross-fold realises across the
///         two blocks (the antiunitary T = P·K carried on the (SE,DE) block).</item>
///   <item><see cref="F89BranchLocusPalindromeClaim"/>: the spectrum-level cross-fold this claim upgrades to a
///         Jordan-structure-preserving matrix similarity (the diabolic character + gap thus pair).</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (an exact matrix identity + the standard antiunitary-similarity
/// argument that semisimple coalescences map to semisimple). Both parents Tier 1 derived.</para>
///
/// <para>Requires <see cref="F1PalindromeIdentity"/> and <see cref="F89BranchLocusPalindromeClaim"/> to be
/// registered (the builder topo-resolves, so the order of registration is free).</para></summary>
public static class F89CrossFoldSimilarityClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89CrossFoldSimilarityClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89CrossFoldSimilarityClaim>(b =>
            new F89CrossFoldSimilarityClaim(
                b.Get<F1PalindromeIdentity>(),
                b.Get<F89BranchLocusPalindromeClaim>()));
}
