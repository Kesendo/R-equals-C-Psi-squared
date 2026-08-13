using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="FoldCheckerboardParityClaim"/> (F150): the Λ-side coefficient-support
/// parity j + k ≡ d (mod 2) on fold-fixed blocks, and its corollary that disc_Λ(F_res) is even in q and therefore
/// real. Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F89CrossFoldSimilarityClaim"/>: the fold lattice, which supplies the leg this law takes
///         determinants of, and consequence (b) naming BOTH legs as the hypothesis.</item>
///   <item><see cref="BipartiteGaugeCriterionClaim"/>: the Δ = 0 gauge. NOT because it carries the q-side
///         antiparity, which it does not: that is Y1, the pencil reality of L = A + iqK with real ingredients,
///         and both entries call it separate and unconditional. The edge is real for two other reasons: the
///         checkerboard inherits the gauge's Δ = 0 scope, and the gauge's character is what decides whether the
///         odd-j cells are occupied at all, which is exactly the difference between this claim's (1,2)@N=4 row
///         and its (1,3)@N=6 row.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived. On the full characteristic polynomial the determinant identity has no
/// premise at all. The descent to F_res = χ/AT is the one step that is measured rather than proved (the AT factor
/// must carry the same parity), and the claim scopes it as such rather than letting the tier label absorb it. Both
/// parents Tier 1 derived.</para>
///
/// <para>Requires <see cref="F89CrossFoldSimilarityClaim"/> and <see cref="BipartiteGaugeCriterionClaim"/> to be
/// registered (the builder topo-resolves, so the order of registration is free).</para></summary>
public static class FoldCheckerboardParityClaimRegistration
{
    public static ClaimRegistryBuilder RegisterFoldCheckerboardParityClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<FoldCheckerboardParityClaim>(b =>
            new FoldCheckerboardParityClaim(
                b.Get<F89CrossFoldSimilarityClaim>(),
                b.Get<BipartiteGaugeCriterionClaim>()));
}
