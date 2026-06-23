using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-3 octic diabolic degeneracy (semisimple, transversal level crossing; Tier 1 derived; analytical
/// from disc factor + numerically verified bit-exact):
///
/// <code>
///   q_EP² = (−1 + √13) / 6 ≈ 0.4343
///   q_EP  ≈ 0.658983
///   λ_EP  = −4γ + 2iJ        (merged eigenvalue at q = q_EP)
/// </code>
///
/// <para>The discriminant of F_8 in λ has a perfect-square factor
/// (3q⁴+q²−1)² which locates the degeneracy. Numerically verified via
/// <c>simulations/f89_path3_ep_locator.py</c>: at q = q_EP, two octic
/// eigenvalues coalesce to machine precision (pair distance ≈ 6.79·10⁻¹⁵)
/// at λ_EP = −4γ + 2iJ.</para>
///
/// <para>The merged-eigenvalue rate Re(λ_EP) = −4γ sits at the AT-spectral
/// midpoint between rate 2γ (overlap) and rate 6γ (no-overlap) of the
/// (SE, DE) sector. The number g_eff = 2/q_EP ≈ 3.034 is the EP-location
/// relation Q_EP = 2/g_eff of the SEPARATE F86a 2-level rate-channel
/// reduction (it fixes the eigenVALUE q_EP/λ_EP); it is NOT a genuine
/// coupling within the octic — the octic's own 2×2 restriction onto the
/// coalescing span is scalar λ·I (off-diagonal ~1e-16), so the eigenVECTORS
/// stay independent (diabolic), not hybridized (defective).</para>
///
/// <para>Re(λ_EP) = −4γ is the J-independent AT-locked rate (the crossing
/// sits at the AT-midpoint for every J); Im(λ_EP) = 2J is the oscillation
/// frequency, which vanishes only in the formal J→0 reading (q_EP itself
/// requires J = q_EP·γ ≠ 0, so there is no degeneracy at J = 0).</para>
///
/// <para><b>Correction (2026-06-21, second-terminal review):</b> the
/// "exceptional point / Jordan-block / defective" character was retracted.
/// Grid-free proof: disc(F_8) carries the EP-condition (3q⁴+q²−1) to even
/// multiplicity 2 (a double zero), so the two eigenvalues cross linearly /
/// analytically through q_EP ⟹ SEMISIMPLE, not defective (a defective
/// √-branch EP forces a simple zero); the octic is irreducible, so this is
/// intrinsic to the double-zero, not a cross-factor crossing. Corroborated
/// artifact-free (g1 = 2 = alg, rank(L−λ_EP·I) = n−2, dep = 0, no generalized
/// eigenvector). Character: a non-defective (semisimple) degeneracy that is
/// obliquely embedded / non-normal (‖P‖ ≈ 3.88 > 1 — ill-conditioned but NOT
/// a Jordan block); NOT "normal". The eigenVALUE double root and
/// λ_EP = −4γ + 2iJ stand; only the eigenVECTOR-coalescence (defective)
/// reading failed. See <c>inspect --root epcharacter</c> (EpCharacter) and
/// <c>simulations/f89_jordan_definitive.py</c>.</para>
/// <para>Anchors: <c>simulations/f89_path3_octic_galois.py</c> (disc
/// factorisation), <c>simulations/f89_path3_ep_locator.py</c> (numerical
/// degeneracy location), <c>simulations/f89_jordan_definitive.py</c> +
/// <c>f89_jordan_corroborate.py</c> (the diabolic-character verdict),
/// <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c>
/// § "Path-3 octic diabolic-degeneracy location".</para></summary>
public sealed class F89Path3OcticEpClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public F89TopologyOrbitClosure F89 { get; }
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    public F89PathKAtLockMechanismClaim AtLock { get; }
    /// <summary>q_EP² = (−1 + √13) / 6: the locus from disc factor (3q⁴+q²−1)² = 0.</summary>
    public static readonly double QEpSquared = (-1.0 + Math.Sqrt(13)) / 6.0;

    /// <summary>q_EP ≈ 0.658983: the diabolic-degeneracy location in q = J/γ.</summary>
    public static readonly double QEp = Math.Sqrt(QEpSquared);

    /// <summary>The merged eigenvalue at q = q_EP: λ_EP = −4γ + 2iJ.
    /// Verified numerically to machine precision (pair distance ≈ 6.79e-15).
    /// The formula's J→0 reading is the real rate λ = −4γ (no oscillation);
    /// it is NOT a degeneracy at J=0 — q_EP fixes J = q_EP·γ ≠ 0.</summary>
    public static Complex MergedEigenvalue(double gamma, double j)
    {
        if (gamma < 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0.");
        if (j < 0) throw new ArgumentOutOfRangeException(nameof(j), j, "J must be ≥ 0.");
        return new Complex(-4 * gamma, 2 * j);
    }

    public F89Path3OcticEpClaim(F89TopologyOrbitClosure f89, F89PathKAtLockMechanismClaim atLock)
        : base("F89 path-3 octic diabolic degeneracy (semisimple) at q² = (−1+√13)/6, q ≈ 0.658983, with merged eigenvalue λ_EP = −4γ + 2iJ; Re(λ_EP) = −4γ sits at the AT-spectral midpoint of rate 2γ (overlap) and rate 6γ (no-overlap)",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/f89_path3_octic_galois.py + " +
               "simulations/f89_path3_ep_locator.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89PathKAtLockMechanismClaim.cs")
    {
        F89 = f89 ?? throw new ArgumentNullException(nameof(f89));
        AtLock = atLock ?? throw new ArgumentNullException(nameof(atLock));
    }

    public override string DisplayName =>
        "F89 path-3 octic diabolic degeneracy: q ≈ 0.659, λ_EP = −4γ + 2iJ at AT-rate-midpoint";

    public override string Summary =>
        $"q² = (−1+√13)/6 ≈ {QEpSquared:F4}, q ≈ {QEp:F6}, λ_EP = −4γ + 2iJ (rate at AT-midpoint between 2γ and 6γ) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("q_EP² = (−1+√13)/6", QEpSquared);
            yield return InspectableNode.RealScalar("q_EP", QEp);
            var lam = MergedEigenvalue(0.05, 0.05 * QEp);
            yield return new InspectableNode("Sample λ_EP at γ=0.05",
                summary: $"({lam.Real}, {lam.Imaginary}i) = (−4γ, 2iJ) at q_EP");
        }
    }
}
