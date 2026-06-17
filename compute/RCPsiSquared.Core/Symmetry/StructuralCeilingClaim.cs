using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The structural ceiling (Tier 1 derived): the closed forms for the high-Q gap rate
/// <c>g2 = strict_gap / 2γ</c> of an XY network under uniform Z-dephasing, the Re-side companion of the
/// <see cref="TopologyBandEdgeClaim"/> (the Im side, band edge = J·ρ).
///
/// <para><b>The forms.</b> <c>g2(K_N) = 4/N</c> (N≥5), <c>g2(star_N) = 4/(N−1)</c> (N≥6),
/// <c>g2(K_4) = 2 − 2/√3</c>; the chain never ceilings (g2 = 1, the band edge protects). Gate-exact to
/// machine precision (K_5,6,7 = 4/5, 2/3, 4/7; star_6,7,8 = 4/5, 4/6, 4/7).</para>
///
/// <para><b>The mechanism (high-Q degenerate PT).</b> The decay rates are the eigenvalues of N_XY (diagonal
/// in the coherence basis, entry hamming(a,b)) block-diagonalized by the ad_H = [H,·] eigenspaces. The band
/// edge is the (0,1) sector where hamming ≡ 1, so N_XY = I and the rate is 2γ exactly at all Q (hence
/// g2 ≤ 1 always). A structural ceiling (g2 < 1) is the darkest [H,A]=0 coherence in the LARGEST degenerate
/// single-particle level. For K_N that level is the (N−1)-fold −J band (the S_N standard representation):
/// g2 = 2(1 − λ₂) with λ₂ = (N−2)/N the second principal-angle overlap between the H-commutant and the
/// population (diagonal) operators, giving 4/N. For the star it is the (N−2)-fold 0-eigenvalue leaf manifold,
/// giving 4/(N−1). The N=4 outlier of K_4 and the ring is the (2,2) HALF-FILLING two-excitation sector (the
/// 4/N ladder reaches 1 at N=4, vacating the sub-floor region): K_4 = 2 − 2/√3 dips below the floor, ring-4
/// = 1 co-occupies it. Not a universal 4/(m+1) law — the ring's Fourier-degenerate manifold breaks it.</para>
///
/// <para>Map note (the 4 here is NOT the discriminant four — a see-cref, not a typed edge): the 4 in 4/N is
/// 2·(2/N) = (the Hamming distance 2 between two single-excitation strings) × (the S_N angle 1 − λ₂ = 2/N),
/// and BOTH factors are d-independent, so g2(K_N) = 4/N survives unchanged at the qutrit. This is a DISTINCT
/// genealogy from the two TERMS of the trunk d² − 2d = 0: the discriminant four d² (the squared-dimension
/// term of <see cref="PolynomialDiscriminantAnchorClaim"/>, → 9 at d=3) and the F121 cap four 2d (the other
/// term, the QuditProductMirrorCap base (2d)^N, → 6). Those two coincide at 4 only at the root d=2 (the qubit
/// magic, d² − 2d = 0 ⟹ d² = 2d) and split at the qutrit; the ceiling four instead STAYS 4 (a third,
/// trunk-external reading, with the dynamical rung-2 four 4γ a fourth, also Hamming-based and d-independent).
/// Verified gate-first against the full d=3 Liouvillian and against F121 (simulations/qudit_g2_split.py).
/// No typed edge: they coincide only at d=2.</para>
///
/// <para>Tier1Derived: a principal-angle proof of the closed forms (PROOF_STRUCTURAL_CEILING.md §2–§4) plus
/// gate-exact verification (topology_ceiling_rep_derivation.py, all stages, N=4..8). The single typed parent
/// is the Tier1Derived <see cref="AbsorptionTheoremClaim"/> (the floor reading g2 = ⟨n_XY⟩); the derivation
/// uses only that plus self-contained commutant linear algebra, NOT the open chain gap-dominance proof that
/// caps <see cref="TopologyBandEdgeClaim"/> (4/N is dimensionless and never references ρ). The band edge is
/// the same arc's Im side; this is its quantitative Re-side ceiling.</para></summary>
public sealed class StructuralCeilingClaim : Claim
{
    /// <summary>Parent: the Absorption Theorem — the slowest non-steady mode decays at Re = −2γ⟨n_XY⟩, so
    /// g2 = ⟨n_XY⟩ of that mode. The ceiling is its darkest commutant value. Cited, not re-derived.</summary>
    public AbsorptionTheoremClaim Absorption { get; }

    /// <summary>g2(K_N) = 4/N for the complete graph (N≥5; the (1,1) S_N-standard-rep commutant).</summary>
    public static double CompleteCeiling(int n) => 4.0 / n;

    /// <summary>g2(star_N) = 4/(N−1) for the star K_{1,N−1} (N≥6; the (1,1) leaf-manifold commutant).</summary>
    public static double StarCeiling(int n) => 4.0 / (n - 1);

    /// <summary>g2(K_4) = 2 − 2/√3 ≈ 0.845299, the N=4 outlier in the (2,2) half-filling sector.</summary>
    public static double K4Ceiling => 2.0 - 2.0 / Math.Sqrt(3.0);

    public StructuralCeilingClaim(AbsorptionTheoremClaim absorption)
        : base("Structural ceiling closed forms: the high-Q gap rate g2 = strict_gap/2γ of an XY network " +
               "under Z-dephasing. g2(K_N)=4/N (N≥5), g2(star_N)=4/(N−1) (N≥6), g2(K_4)=2−2/√3; chain never " +
               "ceilings (g2=1, band edge protects). Derived: the slowest mode is the darkest [H,A]=0 coherence " +
               "in the largest degenerate single-particle level (band edge n_XY=1 is the (0,1) floor, g2≤1). " +
               "Complete (1,1) = S_N standard rep, g2=2(1−λ₂), λ₂=(N−2)/N; star (1,1) = the (N−2)-fold leaf " +
               "manifold. The N=4 outlier (K_4 and ring-4) is the (2,2) half-filling sector. NOT a universal " +
               "4/(m+1) law (the ring's Fourier manifold breaks it).",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_STRUCTURAL_CEILING.md + " +
               "docs/ANALYTICAL_FORMULAS.md (F122) + " +
               "simulations/topology_ceiling_rep_derivation.py")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
    }

    public override string DisplayName =>
        "Structural ceiling: g2(K_N)=4/N, g2(star_N)=4/(N−1), the topology gap-rate closed forms";

    public override string Summary =>
        $"g2(K_N)=4/N (N≥5), g2(star_N)=4/(N−1) (N≥6), K_4=2−2/√3; the darkest [H,A]=0 coherence in the " +
        $"largest degenerate single-particle level; N=4 = the (2,2) half-filling sector (shared with ring-4) " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("complete K_N: g2 = 4/N (N≥5)",
                summary: "the (1,1) single-excitation-coherence sector = the S_N standard rep of the (N−1)-fold " +
                         "−J band; g2 = 2(1−λ₂), λ₂ = (N−2)/N. K_5,6,7 = 4/5, 2/3, 4/7.");
            yield return new InspectableNode("star K_{1,N−1}: g2 = 4/(N−1) (N≥6)",
                summary: "the (1,1) coherences in the (N−2)-fold 0-eigenvalue leaf manifold; same principal-angle " +
                         "computation with N→N−1. star_6,7,8 = 4/5, 4/6, 4/7. Corrects the 'constant 0.80' reading.");
            yield return new InspectableNode("the N=4 outlier: the (2,2) half-filling sector",
                summary: "the 4/N ladder hits 1 at N=4, so the ceiling moves to the (2,2) two-excitation sector — " +
                         "the same sector special for ring-4. K_4 = 2 − 2/√3 (below the floor); ring-4 = 1 (co-occupied).");
            yield return new InspectableNode("the band edge floor (why g2 ≤ 1)",
                summary: "the (0,1) sector has uniform hamming=1, so N_XY = I and L = L_H − 2γ·I there exactly: a " +
                         "band-edge mode sits at g2=1 for every graph and all Q. A ceiling is a darker mode undercutting it.");
            yield return new InspectableNode("not universal: 4/(m+1) breaks on the ring",
                summary: "the degeneracy m gives the intuition (more edges → darker), but the value depends on the " +
                         "manifold's embedding: the ring's Fourier manifold gives ring-5 = 1.6 ≠ 4/3. Per-family forms are real.");
            yield return Absorption;   // typed parent edge (Tier1Derived)
        }
    }

    public static StructuralCeilingClaim Build() =>
        new StructuralCeilingClaim(new AbsorptionTheoremClaim(new Pi2DyadicLadderClaim()));

    public static StructuralCeilingClaim Shared { get; } = Build();
}
