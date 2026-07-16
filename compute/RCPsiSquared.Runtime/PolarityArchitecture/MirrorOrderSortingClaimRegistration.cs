using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Core.SymmetryFamily;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="MirrorOrderSortingClaim"/> (F131, 2026-07-16):
/// the mirror's order-sorting law. Mirror conjugation reflects a parameter scan,
/// M·G(x₀ + s·δ)·M⁻¹ = G(x₀ + σ_eff·s·δ) with σ_eff = σ_op·χ_M, and a readout of definite
/// mirror parity q reads only the orders q·σ_eff allows: generic / EVEN response / ODD
/// response / IDENTICALLY ZERO. Theorem A = the unitary F71 site-reversal column
/// (unconditional, ⟨O⟩(t) = q·⟨O⟩(−t) for operator-R-even preparation); Theorem B = the ζ²
/// anti-protection law (antiunitary Floquet Θ = T·K, tracking hypotheses).
///
/// <para>Tier1Derived (assembly of proven Tier-1 results; self-check battery at N = 3 in
/// the ctor, moment-level, no eigensolver). Typed parents: <see cref="ChiralKClaim"/>
/// (Θ = T·K, Theorem B's mirror), <see cref="AntilinearTriangleClaim"/> (χ_M, the
/// linear/antilinear character), and the F91 family
/// (<see cref="F71AntiPalindromicGammaSpectralInvariance"/> +
/// <see cref="F92BondAntiPalindromicJSpectralInvariance"/> +
/// <see cref="F93DetuningAntiPalindromicSpectralInvariance"/>, the owned pair-sum
/// invariance on the three scanned axes). Anchor:
/// <c>docs/proofs/PROOF_MIRROR_ORDER_SORTING.md</c> +
/// <c>docs/proofs/PROOF_ZETA2_ANTI_PROTECTION.md</c> +
/// <c>simulations/mirror_order_sorting.py</c>; MirrorWorld adoption
/// <c>compute/MirrorWorld/OrderSorting.cs</c> (run mode sorting N).</para></summary>
public static class MirrorOrderSortingClaimRegistration
{
    public static ClaimRegistryBuilder RegisterMirrorOrderSortingClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<MirrorOrderSortingClaim>(b =>
            new MirrorOrderSortingClaim(
                b.Get<ChiralKClaim>(),
                b.Get<AntilinearTriangleClaim>(),
                b.Get<F71AntiPalindromicGammaSpectralInvariance>(),
                b.Get<F92BondAntiPalindromicJSpectralInvariance>(),
                b.Get<F93DetuningAntiPalindromicSpectralInvariance>()));
}
