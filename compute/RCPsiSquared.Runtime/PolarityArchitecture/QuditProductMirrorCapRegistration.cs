using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="QuditProductMirrorCap"/> (2026-06-11): the operator
/// side of F121. The former universal product cap is retracted: at d=6,N=2,
/// P_dark⊗P_lit has exact rank 180 &gt; 144=(2d)^N. The operator
/// the restricted map Π_d P_aligned has rank (2d)^N with exactly zero residual;
/// the full permutation Π_d has rank d^(2N), ord(Π_d) = 2d, and ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂ of order 2d² (D₄ at
/// d = 2). Typed parents <see cref="QuditPartialPalindromeCeiling"/> (F121, the combinatorial
/// ceiling) and <see cref="QubitNecessityPi2Inheritance"/>. Anchor:
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
