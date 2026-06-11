using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="AntilinearTriangleClaim"/> (2026-06-11): the
/// antilinear triangle. The three involutions of operator space, the transpose θ, the
/// entrywise conjugation conj, and the adjoint †, close with id into one Klein four-group
/// († = θ∘conj) graded by linearity ℓ and multiplicativity m; the transport law
/// μ∘L_H∘μ = ℓ(μ)·m(μ)·L_{μ(H)} is the one engine behind five existing proofs (F114 θ-leg,
/// the girth-ladder reversal kill, F112 Lemmas A+B †-leg, F113/F117 Hermitian conjugacy,
/// the K_b mode mirror conj-leg); in the Pauli basis θ = D and † = 𝒦, and
/// ⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ (order 16) is the antilinear double of the mirror group.
///
/// <para>Tier1Derived (one-line algebraic identities, exact; self-check battery at N = 2
/// in the ctor). Typed parents: <see cref="MirrorGroupD4Claim"/> (the D₄ the triangle docks
/// onto and doubles), <see cref="CommutatorDConjugationSign"/> (F114, the θ-leg), and
/// <see cref="LindbladBitBPiBalance"/> (F112, whose Lemmas A+B are the †-leg). The conj-leg
/// claim (ChiralMirrorTrajectoryClaim) lives in RCPsiSquared.Diagnostics and is carried in
/// prose, not as a typed parent. Anchor:
/// <c>docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md</c> +
/// <c>simulations/antilinear_triangle.py</c>.</para></summary>
public static class AntilinearTriangleClaimRegistration
{
    public static ClaimRegistryBuilder RegisterAntilinearTriangleClaim(
        this ClaimRegistryBuilder builder) =>
        builder.Register<AntilinearTriangleClaim>(b =>
            new AntilinearTriangleClaim(
                b.Get<MirrorGroupD4Claim>(),
                b.Get<CommutatorDConjugationSign>(),
                b.Get<LindbladBitBPiBalance>()));
}
