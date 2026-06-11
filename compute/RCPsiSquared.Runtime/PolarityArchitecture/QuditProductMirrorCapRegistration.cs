using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="QuditProductMirrorCap"/> (2026-06-11): the operator
/// side of F121. Any per-site mirror W = ⊗q_l intertwining the dissipator palindrome
/// W·L_D = (−L_D − 2Nγ)·W pairs at most (2d)^N of the d^{2N} coherences, full ⟺ d² − 2d = 0
/// ⟺ d = 2; the operator Π_d(ρ) = ρᵀ·Shift^{⊗N} attains the cap with exactly zero residual on
/// the shift-aligned subspace; ord(Π_d) = 2d and ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂ of order 2d² (D₄ at
/// d = 2). Typed parents <see cref="QuditPartialPalindromeCeiling"/> (F121, the combinatorial
/// ceiling; the gap to it is the non-product part) and <see cref="QubitNecessityPi2Inheritance"/>
/// (the d² − 2d = 0 trunk, third appearance). Anchor:
/// <c>docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md</c> §6 +
/// <c>simulations/qudit_product_mirror_cap.py</c>.
///
/// <para>Requires <see cref="QuditPartialPalindromeCeilingRegistration.RegisterQuditPartialPalindromeCeiling"/>
/// and <see cref="QubitNecessityPi2InheritanceRegistration.RegisterQubitNecessityPi2Inheritance"/>
/// earlier in the builder pipeline.</para></summary>
public static class QuditProductMirrorCapRegistration
{
    public static ClaimRegistryBuilder RegisterQuditProductMirrorCap(
        this ClaimRegistryBuilder builder) =>
        builder.Register<QuditProductMirrorCap>(b =>
            new QuditProductMirrorCap(
                b.Get<QuditPartialPalindromeCeiling>(),
                b.Get<QubitNecessityPi2Inheritance>()));
}
