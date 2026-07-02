using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="SpectatorIntertwinerClaim"/> (Theorem B of
/// <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c>): the site-summed spectator W(ρ) = Σ_l c_l†ρc_l is an
/// EXACT part-by-part intertwiner of the XY + Z-dephasing Liouvillian, block-shifting (p,q̃) → (p+1,q̃+1) and
/// transporting Jordan chains whenever the eigenvector avoids ker W. Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F89CrossFoldSimilarityClaim"/>: F89d, the exact antiunitary leg that combines with the
///         climbing W-step (plus transpose and Klein full flip) into the containment orbit corollary, the
///         diamond of sectors sharing the defective λ.</item>
///   <item><see cref="AbsorptionTheoremClaim"/>: the rate law behind the block pencil's real part
///         A = −2·diag(n_diff), the diagonal object Lemma 2 intertwines site by site.</item>
/// </list>
///
/// <para>Tier consistency: Tier 1 derived (two exact operator identities + the Jordan-transport lemma,
/// gate-verified at machine zero, SpectatorIntertwinerGateTests commit de4f90a). Both parents Tier 1
/// derived.</para>
///
/// <para>Requires <see cref="F89CrossFoldSimilarityClaim"/> and <see cref="AbsorptionTheoremClaim"/> to be
/// registered (the builder topo-resolves, so the order of registration is free).</para></summary>
public static class SpectatorIntertwinerClaimRegistration
{
    public static ClaimRegistryBuilder RegisterSpectatorIntertwinerClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<SpectatorIntertwinerClaim>(b =>
            new SpectatorIntertwinerClaim(
                b.Get<F89CrossFoldSimilarityClaim>(),
                b.Get<AbsorptionTheoremClaim>()));
}
