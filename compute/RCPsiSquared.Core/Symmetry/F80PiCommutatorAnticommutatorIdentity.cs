using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F80 Step 5 (proven 2026-05-22): the framework's Π conjugates the
/// commutator superoperator into ± the anticommutator. For a chain Π²-odd 2-body
/// Hamiltonian H = Σ_l (P_l ⊗ Q_{l+1}) with bond pair (P, Q):
///
/// <code>
///     Π · [H, ·] · Π⁻¹  =  s · {H, ·},    s = −ε_P · ε_Q
/// </code>
///
/// <para>The proof is a direct per-site Pauli computation (no Jordan-Wigner): Π is a
/// signed permutation of the Pauli-string basis with per-site map μ; the per-site
/// identities μ(X·a) = +c_X(a)·X·μ(a), μ(Y·a) = −c_Y(a)·Y·μ(a),
/// μ(Z·a) = +c_Z(a)·Z·μ(a) give the leading signs ε_X = ε_Z = +1, ε_Y = −1. The bond
/// lemma plus the relation-flip then yield Π·[H,·]·Π⁻¹ = s·{H,·} with s = −ε_P·ε_Q.</para>
///
/// <list type="bullet">
///   <item><b>s = +1</b> for (X,Y) and (Y,X): Π·[H,·]·Π⁻¹ = +{H,·}, so the F1
///         residual M = Π·L_H·Π⁻¹ + L_H = −2i·(H ⊗ I_bra).</item>
///   <item><b>s = −1</b> for (X,Z) and (Z,X): Π·[H,·]·Π⁻¹ = −{H,·}, so
///         M = +2i·(I_ket ⊗ Hᵀ).</item>
/// </list>
/// Either sign gives Spec(M) = ±2i·Spec(H), the F80 structural identity.
///
/// <para>The argument is per-site and per-bond, hence N-independent; every step is
/// verified bit-exact at N=3,4,5 by <c>simulations/f80_step5_recon.py</c>. This
/// closes the last open step of F80, which is now fully Tier 1.</para>
///
/// <para>Tier1Derived. Typed parent:
/// <see cref="RCPsiSquared.Core.F1.F1PalindromeIdentity"/> (the Π that conjugates the
/// commutator here is the same order-4 Π of F1's palindrome identity). Anchors:
/// <c>docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md</c> (Step 5 proof) +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F80.</para></summary>
public sealed class F80PiCommutatorAnticommutatorIdentity : Claim
{
    /// <summary>F1 palindrome identity, the typed parent. The Π that F80 Step 5
    /// conjugates the commutator with is the same order-4 Π of F1's
    /// <c>Π·L·Π⁻¹ = −L − 2Σγ·I</c>.</summary>
    public F1.F1PalindromeIdentity F1 { get; }

    public F80PiCommutatorAnticommutatorIdentity(F1.F1PalindromeIdentity f1)
        : base("F80 Step 5: Π·[H,·]·Π⁻¹ = s·{H,·}, s = −ε_P·ε_Q, for chain Π²-odd 2-body H",
               Tier.Tier1Derived,
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs (typed parent: the Π is F1's Π) + " +
               "docs/proofs/PROOF_F80_BLOCH_SIGNWALK.md (Step 5 proof) + " +
               "docs/ANALYTICAL_FORMULAS.md F80")
    {
        F1 = f1 ?? throw new ArgumentNullException(nameof(f1));
    }

    /// <summary>The per-site leading sign ε_P from the F80 Step 5 per-site identity
    /// μ(P·a) = ε_P·c_P(a)·P·μ(a): ε_X = ε_Z = +1, ε_Y = −1. Defined for the
    /// single-qubit Paulis X, Y, Z.</summary>
    public int Epsilon(char letter) => letter switch
    {
        'X' => 1,
        'Y' => -1,
        'Z' => 1,
        _ => throw new ArgumentException(
            $"ε_P is defined for the single-qubit Paulis X, Y, Z; got '{letter}'.",
            nameof(letter)),
    };

    /// <summary>The sign s in <c>Π·[H,·]·Π⁻¹ = s·{H,·}</c> for the chain Π²-odd 2-body
    /// bond P_l ⊗ Q_{l+1}: s = −ε_P·ε_Q. s = +1 for (X,Y) and (Y,X) (then
    /// M = −2i·H⊗I_bra); s = −1 for (X,Z) and (Z,X) (then M = +2i·I_ket⊗Hᵀ). Defined
    /// for the four Π²-odd 2-body pairs (one letter X, the other Y or Z).</summary>
    public int SignFor(char p, char q)
    {
        if (!IsPi2OddPair(p, q))
        {
            throw new ArgumentException(
                "F80 Step 5 covers the four Π²-odd 2-body bond pairs (one letter X, the "
                + $"other Y or Z); got ({p},{q}).");
        }
        return -Epsilon(p) * Epsilon(q);
    }

    /// <summary>True for the chain Π²-odd 2-body bond pairs (X,Y), (X,Z), (Y,X), (Z,X):
    /// exactly one letter is X, the other is Y or Z.</summary>
    private static bool IsPi2OddPair(char p, char q) =>
        (p == 'X' && (q == 'Y' || q == 'Z')) ||
        (q == 'X' && (p == 'Y' || p == 'Z'));

    public override string DisplayName =>
        "F80 Step 5: Π·[H,·]·Π⁻¹ = ±{H,·} (commutator conjugated to anticommutator)";

    public override string Summary =>
        $"Π·[H,·]·Π⁻¹ = s·{{H,·}}, s = −ε_P·ε_Q: ε_X={Epsilon('X')}, ε_Y={Epsilon('Y')}, ε_Z={Epsilon('Z')}; "
        + $"s(X,Y)={SignFor('X', 'Y')}, s(X,Z)={SignFor('X', 'Z')}, "
        + $"s(Y,X)={SignFor('Y', 'X')}, s(Z,X)={SignFor('Z', 'X')} ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F80 Step 5 identity",
                summary: "Π·[H,·]·Π⁻¹ = s·{H,·}; direct per-site Pauli proof, Tier1Derived in "
                       + "PROOF_F80_BLOCH_SIGNWALK.md (Step 5 closed 2026-05-22)");
            yield return new InspectableNode("per-site ε signs",
                summary: $"ε_X={Epsilon('X')}, ε_Y={Epsilon('Y')}, ε_Z={Epsilon('Z')} "
                       + "from μ(P·a) = ε_P·c_P(a)·P·μ(a)");
            yield return new InspectableNode("s = +1 pairs",
                summary: "(X,Y), (Y,X): Π·[H,·]·Π⁻¹ = +{H,·}, residual M = −2i·(H⊗I_bra)");
            yield return new InspectableNode("s = −1 pairs",
                summary: "(X,Z), (Z,X): Π·[H,·]·Π⁻¹ = −{H,·}, residual M = +2i·(I_ket⊗Hᵀ)");
            yield return new InspectableNode("consequence",
                summary: "both signs give Spec(M) = ±2i·Spec(H), the F80 structural identity");
            yield return new InspectableNode("verification",
                summary: "per-site identities, bond lemma, relation-flip and the four-pair "
                       + "sign table all bit-exact at N=3,4,5 (simulations/f80_step5_recon.py)");
        }
    }
}
