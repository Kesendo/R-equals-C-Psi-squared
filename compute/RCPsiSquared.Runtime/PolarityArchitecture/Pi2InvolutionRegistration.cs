using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="Pi2InvolutionClaim"/> — the typed
/// statement <c>Π²·L·Π⁻² = L</c> (equivalently <c>[L, Π²] = 0</c>): squaring the
/// F1 palindrome closes the algebra back to identity, and L block-diagonalises
/// in the Π²-eigenbasis.
///
/// <para>This claim is intentionally NOT in <see cref="Pi2FamilyRegistration"/>
/// because it is a <i>consequence</i> of F1 rather than a Pi2 foundation. The
/// Pi2 Family is the foundation; the involution claim is its squared image
/// downstream of F1.</para>
///
/// <para>The Pi2InvolutionClaim is the algebraic backbone of
/// <c>hypotheses/HEISENBERG_RELOADED.md</c> and the
/// <c>docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md</c> proof: every Heisenberg/XXZ
/// Hamiltonian under uniform Z-dephasing satisfies <c>[L, Π²] = 0</c>, and L
/// therefore decomposes into Π²-blocks. Tom 2026-05-09 surfaced the gap: the
/// claim was in Core but not queryable through the Object Manager.</para>
///
/// <para>Edge declared: <see cref="Pi2InvolutionClaim"/> ←
/// <see cref="F1PalindromeIdentity"/>. Squaring F1 directly produces Π²·L·Π⁻² = L
/// (the involution claim's content), so F1 is the natural parent.</para>
///
/// <para>Tier consistency: both Tier1Derived (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="F1Family.F1FamilyRegistration.RegisterF1Family"/>
/// for <see cref="F1PalindromeIdentity"/>.</para></summary>
public static class Pi2InvolutionRegistration
{
    public static ClaimRegistryBuilder RegisterPi2Involution(
        this ClaimRegistryBuilder builder) =>
        builder.Register<Pi2InvolutionClaim>(b =>
        {
            var f1 = b.Get<F1PalindromeIdentity>();   // F1 squared gives Π²·L·Π⁻² = L
            return new Pi2InvolutionClaim(f1);
        });
}
