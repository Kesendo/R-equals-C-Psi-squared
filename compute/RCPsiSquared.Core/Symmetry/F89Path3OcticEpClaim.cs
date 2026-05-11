using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-3 octic exceptional point (Tier 1 derived; analytical
/// from disc factor + numerically verified bit-exact):
///
/// <code>
///   q_EP² = (−1 + √13) / 6 ≈ 0.4343
///   q_EP  ≈ 0.658983
///   λ_EP  = −4γ + 2iJ        (merged eigenvalue at q = q_EP)
/// </code>
///
/// <para>The discriminant of F_8 in λ has a perfect-square factor
/// (3q⁴+q²−1)² which locates the EP. Numerically verified via
/// <c>simulations/_f89_path3_ep_locator.py</c>: at q = q_EP, two octic
/// eigenvalues coalesce to machine precision (pair distance ≈ 6.79·10⁻¹⁵)
/// at λ_EP = −4γ + 2iJ.</para>
///
/// <para>The merged-eigenvalue rate Re(λ_EP) = −4γ sits at the AT-spectral
/// midpoint between rate 2γ (overlap) and rate 6γ (no-overlap) of the
/// (SE, DE) sector. The frequency 2J corresponds to a 2-level effective
/// bilinear coupling g_eff = 2/q_EP ≈ 3.034 within the (SE, DE) octic
/// sub-block.</para>
///
/// <para>At J=0 the merged eigenvalue degenerates to λ_EP = −4γ (no
/// oscillation; the EP becomes a real degenerate double root).</para>
///
/// <para>Anchors: <c>simulations/_f89_path3_octic_galois.py</c> (disc
/// factorisation), <c>simulations/_f89_path3_ep_locator.py</c> (numerical
/// EP location), <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c>
/// § "Path-3 octic EP location".</para></summary>
public sealed class F89Path3OcticEpClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    private readonly F89TopologyOrbitClosure _f89;
    // Parent-edge marker for Schicht-1 wiring (consumed by ClaimRegistryBuilder; not used in this class body).
    private readonly F89PathKAtLockMechanismClaim _atLock;

    /// <summary>q_EP² = (−1 + √13) / 6: the locus from disc factor (3q⁴+q²−1)² = 0.</summary>
    public static readonly double QEpSquared = (-1.0 + Math.Sqrt(13)) / 6.0;

    /// <summary>q_EP ≈ 0.658983: the EP location in q = J/γ.</summary>
    public static readonly double QEp = Math.Sqrt(QEpSquared);

    /// <summary>The merged eigenvalue at q = q_EP: λ_EP = −4γ + 2iJ.
    /// Verified numerically to machine precision (pair distance ≈ 6.79e-15).
    /// At J=0 degenerates to a real double root λ = −4γ.</summary>
    public static Complex MergedEigenvalue(double gamma, double j)
    {
        if (gamma < 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0.");
        if (j < 0) throw new ArgumentOutOfRangeException(nameof(j), j, "J must be ≥ 0.");
        return new Complex(-4 * gamma, 2 * j);
    }

    public F89Path3OcticEpClaim(F89TopologyOrbitClosure f89, F89PathKAtLockMechanismClaim atLock)
        : base("F89 path-3 octic exceptional point at q² = (−1+√13)/6, q ≈ 0.658983, with merged eigenvalue λ_EP = −4γ + 2iJ; Re(λ_EP) = −4γ sits at the AT-spectral midpoint of rate 2γ (overlap) and rate 6γ (no-overlap)",
               Tier.Tier1Derived,
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "simulations/_f89_path3_octic_galois.py + " +
               "simulations/_f89_path3_ep_locator.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F89PathKAtLockMechanismClaim.cs")
    {
        _f89 = f89 ?? throw new ArgumentNullException(nameof(f89));
        _atLock = atLock ?? throw new ArgumentNullException(nameof(atLock));
    }

    public override string DisplayName =>
        "F89 path-3 octic EP: q ≈ 0.659, λ_EP = −4γ + 2iJ at AT-rate-midpoint";

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
