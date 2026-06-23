using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-3 octic Galois group (Tier 1 derived): the group is the
/// FULL symmetric group S_8, so the octic is non-solvable.
///
/// <code>
///   disc(F_8) = 1.21·10²⁴ · q²⁴ · (3q⁴+q²−1)² · P_10(q²)   (deg 52 in q)
///
///   P_10(q²) is NOT a perfect square in Q  ⇒  Gal(F_8) ⊄ A_8.
///   Gal(F_8 / Q(i)(q)) = S_8   (robust to enlarging the base to Q(i,√5)(q)).
/// </code>
///
/// <para>METHOD: specialization + Dedekind (Frobenius cycle types) + Jordan.
/// For a good q0 ∈ Q (disc(q0) ≠ 0), G_{q0} = Gal(F_8(·,q0)/Q(i)) is a SUBGROUP
/// of the generic G = Gal(F_8/Q(i)(q)) (specialization can only shrink), so
/// G_{q0} = S_8 forces G = S_8. CERTIFICATE at q0 = 2 (F_8(·,2) monic over Z[i]):
/// it is irreducible over Q(i) (⇒ transitive); the split prime 𝔭 | 5
/// (Z[i]/𝔭 = F_5, i ↦ 2) factors F_8(·,2) to cycle type (5,2,1) — whose square
/// is a 5-cycle (⇒ primitive, since a 5-orbit fits no degree-8 block system, and
/// no proper primitive degree-8 group has order divisible by 5 ⇒ ⊇ A_8 by Jordan)
/// and which is itself an odd permutation (⇒ ⊄ A_8). Hence G_{q0=2} = S_8.
/// Confirmed at q0 ∈ {2, 3, ½, 3/2} over Q(i) and at q0 ∈ {2, 3} over Q(i,√5).</para>
///
/// <para>CONSEQUENCE: the eight roots λ_k(q) (Liouvillian eigenvalues of this block)
/// admit NO expression by RADICALS over Q(i)(q) (Abel-Ruffini). This does NOT
/// exclude non-radical special-function expressions (Bring radicals / theta /
/// hypergeometric), which exist for any algebraic function.</para>
///
/// <para>FRAMING: S_8 is the GENERIC Galois group of an irreducible degree-8
/// polynomial (van der Waerden 1936; Bhargava, Annals 2025), so S_8 is not exotic.
/// The content is NEGATIVE: free-fermion integrability spends itself entirely on
/// the factorisation — the AT-locked F_a/F_b quadratics carry the single-particle
/// frequencies in radicals, and the diabolic point sits on the solvable quartic
/// factor (3q⁴+q²−1) (<see cref="F89Path3OcticEpClaim"/>) — leaving the residual
/// octic with no further algebraic structure. The closed-form program for path-3
/// terminates exactly at the AT-protected half. Contrast: the SIC-POVM spectral
/// polynomials (Appleby et al. 2012) gave a SOLVABLE Galois group, the opposite
/// polarity at the same seam.</para>
///
/// <para>SCOPE: this is the Galois group of the path-3 (SE,DE) S_2-sym octic
/// FACTOR, not of "the Liouvillian spectrum". It is a similarity-invariant of that
/// invariant sub-block (the AT-lock split is a genuine invariant-subspace
/// decomposition, not a basis choice), so an equivalent block gives the same group.</para>
///
/// <para>Anchors: live witness <c>inspect --root f89galois</c>
    /// (<c>F89OcticGaloisWitness</c>, in Diagnostics) recomputes the
    /// Frobenius certificate at inspect time via
    /// <see cref="RCPsiSquared.Core.Numerics.OcticGaloisCertificate"/> (distinct-degree
    /// factorisation over F_p); <c>simulations/f89_path3_octic_galois.py</c> (gate-first:
/// reproduce the octic from the 12×12 charpoly, known-answer engine validation,
/// the q0=2 certificate, multi-q0 + base-field robustness), <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c>
/// § "Path-3 octic non-solvability: Gal = S_8". References: K. Conrad,
/// "Recognizing Galois groups S_n and A_n" (the method); M. Bhargava, Annals 201
/// (2025) (S_n generic); Appleby-Yadsan-Appleby-Zauner, arXiv:1209.1813 (the
/// contrasting solvable SIC precedent). Do NOT conflate with differential Galois
/// theory / "Liouvillian solutions" (Morales-Ramis): a different object.</para></summary>
public sealed class F89Path3OcticGaloisClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public F89TopologyOrbitClosure F89 { get; }
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public F89PathKAtLockMechanismClaim AtLock { get; }
    /// <summary>Total degree of disc(F_8) as a polynomial in q: 52.</summary>
    public const int DiscriminantPolynomialDegreeInQ = 52;

    /// <summary>Structure of disc(F_8) factorisation:
    /// (constant prefactor) · q^24 · (3q⁴+q²-1)² · P_10(q²).
    /// Returns (q-power=24, square-factor-deg-in-q=8, P_10-deg-in-q=20).
    /// The constant prefactor is always present once (≈ 1.21·10²⁴) and is omitted from
    /// this structural tuple — the slots encode degree contributions only.</summary>
    public static (int QPower, int SquareFactorDegInQ, int P10DegInQ)
        DiscriminantFactorisationStructure => (24, 8, 20);

    /// <summary>Gal(F_8) ⊄ A_8: true (from disc non-square; verified at
    /// q ∈ {½, 1, 3/2, 2, 3} and, decisively, by an odd Frobenius element at q0=2).</summary>
    public const bool GalNotInA8 = true;

    /// <summary>Gal(F_8 / Q(i)(q)) = S_8 (the full symmetric group): true, Tier 1.
    /// Established by the specialization + Dedekind + Jordan certificate, not by
    /// disc-non-square alone.</summary>
    public const bool GalIsS8 = true;

    /// <summary>The octic is non-solvable (S_8 ⊇ A_8, simple non-abelian for n ≥ 5).
    /// CLOSED — superseding the earlier Tier-2 open conjecture.</summary>
    public const bool NonSolvableConjecture_IsOpen = false;

    /// <summary>The split rational prime that certifies S_8 at q0 = 2:
    /// 𝔭 | 5 (Z[i]/𝔭 = F_5, i ↦ 2) factors F_8(·,2) to cycle type (5,2,1).</summary>
    public const int CertifyingPrimeAtQ0Eq2 = 5;

    public F89Path3OcticGaloisClaim(F89TopologyOrbitClosure f89, F89PathKAtLockMechanismClaim atLock)
        : base("F89 path-3 octic Galois group: Gal(F_8 / Q(i)(q)) = S_8 (Tier 1 derived). disc(F_8) = const · q²⁴ · (3q⁴+q²−1)² · P_10(q²) is non-square ⇒ ⊄ A_8; a (5,2,1) Frobenius at 𝔭|5, q0=2, gives a 5-cycle (⇒ primitive ⇒ ⊇A_8 by Jordan) and an odd permutation (⇒ ⊄A_8) ⇒ S_8. Robust to base Q(i,√5). ⇒ the roots λ_k(q) admit no expression by radicals in q (non-radical special functions not excluded)",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/F89OcticGaloisWitness.cs (inspect --root f89galois) + " +
               "simulations/f89_path3_octic_galois.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89PathKAtLockMechanismClaim.cs")
    {
        F89 = f89 ?? throw new ArgumentNullException(nameof(f89));
        AtLock = atLock ?? throw new ArgumentNullException(nameof(atLock));
    }

    public override string DisplayName =>
        "F89 path-3 octic Galois: Gal(F_8 / Q(i)(q)) = S_8 (Tier 1 derived; non-solvable, no radical closure in q)";

    public override string Summary =>
        $"disc(F_8) = const · q²⁴ · (3q⁴+q²-1)² · P_10(q²) (deg 52 in q), non-square ⇒ ⊄A_8; certified S_8 via a (5,2,1) Frobenius at 𝔭|{CertifyingPrimeAtQ0Eq2} (q0=2) ⇒ non-solvable, no radical closure in q ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("disc(F_8) total degree in q",
                summary: $"{DiscriminantPolynomialDegreeInQ}");
            yield return new InspectableNode("Gal(F_8) ⊄ A_8",
                summary: $"{GalNotInA8} (disc non-square; odd Frobenius at q0=2)");
            yield return new InspectableNode("Gal(F_8 / Q(i)(q)) = S_8",
                summary: $"{GalIsS8} (Tier 1: (5,2,1) Frobenius at 𝔭|{CertifyingPrimeAtQ0Eq2}, q0=2 ⇒ primitive+5-cycle ⇒ ⊇A_8; + ⊄A_8 ⇒ S_8; robust to Q(i,√5))");
            yield return new InspectableNode("Non-solvable / no radical closure",
                summary: $"non-solvable (S_8); roots λ_k(q) admit no expression by radicals over Q(i)(q) (non-radical special functions not excluded)");
        }
    }
}
