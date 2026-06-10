using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="MirrorGroupD4Claim"/> (2026-06-10): the mirror
/// group D₄. Π_Z = R·D (transpose first, ket reflection by X^⊗N second); ⟨R, D⟩ ≅ D₄ of
/// order 8 closes the mirror inventory (the April palindromizers Π_Z/Π_Y, the
/// windowed-converse spine V₄ {I, 𝓕, R, 𝓕R} as Klein subgroup, F114's D, and the fourth
/// diagonal mirror 𝓕D = diag((−1)^{n_Z})); the palindrome splits along the generators
/// (D flips L_H, R carries −2Σγ); the polarity cube axes are the characters of
/// Ad_{Z^⊗N}, Ad_{X^⊗N}, and the transpose θ.
///
/// <para>Tier1Derived (signed-permutation identities, exact; self-check battery at N = 2
/// in the ctor). Typed parents: <see cref="KleinEightCellClaim"/> (the cube whose axes
/// become characters), <see cref="CommutatorDConjugationSign"/> (F114, D's sign law),
/// and <see cref="Pi2KleinV4DephaseSwapGroup"/> (owner of D; its D·Π_Z·D = Π_Y is the
/// dihedral inversion relation). Anchor:
/// <c>docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md</c> +
/// <c>simulations/mirror_inventory_d4.py</c>.</para></summary>
public static class MirrorGroupD4ClaimRegistration
{
    public static ClaimRegistryBuilder RegisterMirrorGroupD4Claim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<MirrorGroupD4Claim>(b =>
            new MirrorGroupD4Claim(
                b.Get<KleinEightCellClaim>(),
                b.Get<CommutatorDConjugationSign>(),
                b.Get<Pi2KleinV4DephaseSwapGroup>()));
}
