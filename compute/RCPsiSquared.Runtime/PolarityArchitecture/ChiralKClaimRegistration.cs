using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Wiring of <see cref="ChiralKClaim"/>: the chiral / sublattice symmetry
/// K = diag((−1)^ℓ) = ⊗_{odd i} Z_i, KHK = −H, AZ class BDI, spectrum inversion
/// E_{N+1−k} = −E_k. A typed SIBLING ROOT (no ctor parents) , like
/// <see cref="RCPsiSquared.Core.F1.F1PalindromeIdentity"/> and
/// <see cref="PolynomialFoundationClaim"/>, it derives from no other claim
/// (docs/PI2KB_INHERITANCE_MAP.md). K acts at the Hamiltonian (single-particle) level; F1 acts
/// at the Liouvillian (4^N) level; neither flows from the other.
///
/// <para>Wired 2026-05-30 to bring the second mirror into the Object Manager: the chiral K was
/// prose-and-helper only (deferred on the RegistryWiringAuditTests allowlist), so
/// <c>inspect --claim ChiralKClaim</c> could not render it. Registering it as a root singleton
/// makes the K-mirror a first-class typed-knowledge node alongside the Π palindrome.</para></summary>
public static class ChiralKClaimRegistration
{
    public static ClaimRegistryBuilder RegisterChiralK(this ClaimRegistryBuilder builder) =>
        builder.Register<ChiralKClaim>(_ => new ChiralKClaim());
}
