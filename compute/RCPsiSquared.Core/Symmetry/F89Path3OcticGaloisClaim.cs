using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-3 octic Galois group property (Tier 1 derived for
/// Gal ⊄ A_8; Tier 2 conjectural for non-solvability):
///
/// <code>
///   disc(F_8) = 1.21·10²⁴ · q²⁴ · (3q⁴+q²−1)² · P_10(q²)
///
///   where P_10(q²) is degree-10 in q² (= 20 in q, even powers only) and
///   is NOT a perfect square in Q (verified at q ∈ {½, 1, 3/2, 2, 3}, all
///   give irrational √disc).
///
///   ⇒ Gal(F_8) ⊄ A_8   (Tier 1, follows from disc-non-square).
/// </code>
///
/// <para>Disc-non-square alone does NOT prove non-solvability (e.g. S_4 is
/// solvable but ⊄ A_4); pinning the exact Galois group requires further
/// resolvent analysis (open). Combined with the verified irreducibility of
/// F_8 over Q[i, √5] (separate Claim) AND the absence of polynomial fits
/// (≤ degree 5 in q) for the per-mode amplitudes, we CONJECTURE (Tier 2)
/// that Gal(F_8) is non-solvable (likely the full S_8).</para>
///
/// <para>Anchors: <c>simulations/_f89_path3_octic_galois.py</c> (disc
/// factorisation + numerical square-root tests), <c>simulations/_f89_path3_octic_amplitude_q_scan.py</c>
/// (no polynomial-q fit), <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c>
/// § "Path-3 octic non-solvability".</para></summary>
public sealed class F89Path3OcticGaloisClaim : Claim
{
    private readonly F89TopologyOrbitClosure _f89;
    private readonly F89PathKAtLockMechanismClaim _atLock;

    /// <summary>Total degree of disc(F_8) as a polynomial in q: 52.</summary>
    public const int DiscriminantPolynomialDegreeInQ = 52;

    /// <summary>Structure of disc(F_8) factorisation:
    /// (constant) · q^24 · (3q⁴+q²-1)² · P_10(q²).
    /// Returns (constant-factor-power=1, q-power=24, square-factor-deg-in-q=8, P_10-deg-in-q=20).</summary>
    public static (int ConstantFactorPower, int QPower, int SquareFactorDegInQ, int P10DegInQ)
        DiscriminantFactorisationStructure => (1, 24, 8, 20);

    /// <summary>Gal(F_8) ⊄ A_8: true (Tier 1, from disc non-square verified
    /// at q ∈ {½, 1, 3/2, 2, 3}).</summary>
    public const bool GalNotInA8 = true;

    /// <summary>The conjecture that Gal(F_8) is non-solvable (likely S_8) is
    /// OPEN (Tier 2). Disc non-square alone does not prove non-solvability;
    /// a complete resolvent analysis is required and not yet performed.</summary>
    public const bool NonSolvableConjecture_IsOpen = true;

    public F89Path3OcticGaloisClaim(F89TopologyOrbitClosure f89, F89PathKAtLockMechanismClaim atLock)
        : base("F89 path-3 octic Galois group: disc(F_8) = const · q²⁴ · (3q⁴+q²−1)² · P_10(q²) is NOT a perfect square in Q[q]; therefore Gal(F_8) ⊄ A_8 (Tier 1). Conjecture (Tier 2 open): Gal is non-solvable (likely S_8)",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_f89_path3_octic_galois.py + " +
               "simulations/_f89_path3_octic_amplitude_q_scan.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89PathKAtLockMechanismClaim.cs")
    {
        _f89 = f89 ?? throw new ArgumentNullException(nameof(f89));
        _atLock = atLock ?? throw new ArgumentNullException(nameof(atLock));
    }

    public override string DisplayName =>
        "F89 path-3 octic Galois: Gal ⊄ A_8 (Tier 1, disc-non-square); non-solvability conjecture Tier 2 open";

    public override string Summary =>
        $"disc(F_8) = const · q²⁴ · (3q⁴+q²-1)² · P_10(q²) (deg 52 in q); P_10 NOT square ⇒ Gal ⊄ A_8 (Tier 1); non-solvability open ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("disc(F_8) total degree in q",
                summary: $"{DiscriminantPolynomialDegreeInQ}");
            yield return new InspectableNode("Gal(F_8) ⊄ A_8 (Tier 1)",
                summary: $"{GalNotInA8} (from disc non-square verified at q ∈ {{½, 1, 3/2, 2, 3}})");
            yield return new InspectableNode("Non-solvability conjecture",
                summary: $"Open (Tier 2): formal Galois-group identification beyond Gal ⊄ A_8 not yet performed");
        }
    }
}
