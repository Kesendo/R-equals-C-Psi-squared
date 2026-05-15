using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="F89AmplitudeLayerClaim"/>: the typed
/// amplitude-layer bridge from F89c (eigenvalue-layer, AT-locked) to D_k
/// (amplitude-layer, open Tier-1-Derived gap). Two typed parent edges:
///
/// <list type="bullet">
///   <item><see cref="F89UnifiedFaClosedFormClaim"/>: carries the closed-form
///         (P_k, D_k, σ_n, σ-orbit-sum) the Angle A identity decomposes.</item>
///   <item><see cref="F89PathKAtLockMechanismClaim"/>: AT-lock eigenvalue
///         structure λ = −2γ₀ + i·y_n that anchors the F_a witnesses.</item>
/// </list>
///
/// <para>Tier consistency: Tier 2 verified (numerically locked at k=3..6 via
/// the theory probe; path-3 algebraic at (33+14√5)/9). Parents both Tier 1
/// derived; Tier-2 child appropriately downgraded by the open generic-k
/// symbolic derivation gap documented in
/// <c>docs/proofs/PROOF_F89_PATH_D_CLOSED_FORM.md</c> Open Questions.</para>
///
/// <para>Requires: <see cref="F89UnifiedFaClosedFormClaimRegistration.RegisterF89UnifiedFaClosedFormClaim"/>
/// and <see cref="F89PathKAtLockMechanismClaimRegistration.RegisterF89PathKAtLockMechanismClaim"/>.</para></summary>
public static class F89AmplitudeLayerClaimRegistration
{
    public static ClaimRegistryBuilder RegisterF89AmplitudeLayer(
        this ClaimRegistryBuilder builder) =>
        builder.Register<F89AmplitudeLayerClaim>(b =>
        {
            var unified = b.Get<F89UnifiedFaClosedFormClaim>();
            var atLock = b.Get<F89PathKAtLockMechanismClaim>();
            return new F89AmplitudeLayerClaim(unified, atLock);
        });
}
