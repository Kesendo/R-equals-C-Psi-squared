using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="QuditProductMirrorCap"/> (2026-06-11): the operator
/// side of F121. Any per-site mirror pairs at most P(d, N) = max_m (2d)^(N−2m)·(d³ − d²)^m
/// coherences, which is (2d)^N for d ≤ 5 and larger from d = 6 (P_dark⊗P_lit at d = 6, N = 2:
/// rank 180 &gt; 144); full ⟺ d = 2. The operator Π_d(ρ) = ρᵀ·Shift^{⊗N} restricted to the
/// shift-aligned subspace has rank (2d)^N with exactly zero residual; the full permutation Π_d has
/// rank d^(2N), ord(Π_d) = 2d, and ⟨Π_d, D⟩ ≅ Z_d ≀ Z₂ of order 2d² (D₄ at d = 2). Typed parents <see cref="QuditPartialPalindromeCeiling"/> (F121, the combinatorial
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
