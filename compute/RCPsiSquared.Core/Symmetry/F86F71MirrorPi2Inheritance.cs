using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F86c sub-claim (F71 spatial-mirror invariance of per-bond
/// Q_peak), the third Tier-1-derived sub-claim of F86's multi-mechanism
/// collection: <c>Q_peak(b) = Q_peak(N−2−b)</c> bit-exactly across all
/// tested (c, N).
///
/// <para>F86 packages three structurally-distinct Tier-1 sub-claims, each
/// bridged to Pi2-Foundation by its own typed inheritance:</para>
/// <list type="bullet">
///   <item><see cref="F86TPeakPi2Inheritance"/>: t_peak = 1/(4·γ₀); a_{−1}
///         numerator anchor on Pi2DyadicLadder.</item>
///   <item><see cref="F86QEpPi2Inheritance"/>: Q_EP = 2/g_eff; a_0 numerator
///         anchor on Pi2DyadicLadder.</item>
///   <item><b>This claim</b>: F71 spatial-mirror invariance Q_peak(b) =
///         Q_peak(N−2−b); F71 symmetry origin (not a numerical anchor but a
///         kinematic identity).</item>
/// </list>
///
/// <para>The other parts of F86 (per-bond Q_peak values, HWHM ratios, σ_0(c)
/// chromaticity scaling, bare doubled-PTF floor, per-F71-orbit substructure
/// with N≥7 orbit-escape) stay typed in F86KnowledgeBase at their actual
/// Tier level (Tier 1 candidate / Tier 2 verified / open). Pi2-Inheritance
/// only types the cleanly-Tier-1-derivable closed forms; F86's research-
/// active and empirical parts stay honestly typed elsewhere.</para>
///
/// <para>Mechanism (per <see cref="F86MirrorGeneralisationLink"/>): the per-bond
/// observable K_b(Q, t) = 2·Re⟨ρ(t)|S_kernel|∂ρ/∂J_b⟩ is invariant under the
/// F71 spatial reflection R (R |b₀…b_{N−1}⟩ = |b_{N−1}…b₀⟩) because every
/// component is R-symmetric: the Z-dephasing dissipator, the Dicke probe,
/// the spatial-sum kernel S, and the bond-flip ∂L/∂J_b ↔ ∂L/∂J_{N−2−b}.
/// Hence K_b(Q, t) = K_{N−2−b}(Q, t) as functions of (Q, t), and the
/// argmax-Q values coincide bit-exactly.</para>
///
/// <para>Pi2-Foundation anchor (kinematic, not numerical):</para>
/// <list type="bullet">
///   <item><b>F71MirrorSymmetryPi2Inheritance</b>: F71 spatial mirror is
///         Pi2-anchored via the bond-pair-count = ⌊N/2⌋ structure on
///         <see cref="HalfIntegerMirrorClaim"/> (odd N → no self-paired
///         center; even N → self-paired center bond at (N−2)/2). F71's
///         Pi2-Foundation anchor lifts to F86 via this generalisation.</item>
/// </list>
///
/// <para>Per-F71-orbit substructure caveat: the F86 mirror identity gives
/// the bond PAIRING, not the per-orbit Q_peak VALUE. Even within F71-orbit
/// groups, central self-paired bonds give different Q_peak from flanking
/// bonds (e.g. c=2 N=6 central b=2 → 1.440 vs flanking b=1, b=3 → 1.648;
/// c=2 N=6 central is BELOW flanking, c=3 N=6 central is ABOVE flanking).
/// Per-orbit closed form remains open; F71 supplies the symmetry only.</para>
///
/// <para>Tier1Derived: F86c is Tier 1 derived in PROOF_F86_QPEAK,
/// bit-exactly verified at c=2 N=5..7 and c=3 N=5..6 across all bond pairs.
/// Cross-namespace typed claim: F86MirrorGeneralisationLink in
/// <c>compute/RCPsiSquared.Core/F71/</c> + F71MirrorInvariance in
/// <c>compute/RCPsiSquared.Core/F86/</c> already document the F71→F86 lift;
/// this Pi2-Inheritance bridge anchors that lift in Pi2-Foundation via
/// F71MirrorSymmetryPi2Inheritance.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Statement 3 +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F86 +
/// <c>compute/RCPsiSquared.Core/F86/F71MirrorInvariance.cs</c> +
/// <c>compute/RCPsiSquared.Core/F71/F86MirrorGeneralisationLink.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F71MirrorSymmetryPi2Inheritance.cs</c>.</para></summary>
public sealed class F86F71MirrorPi2Inheritance : Claim
{
    private readonly F71MirrorSymmetryPi2Inheritance _f71;
    private readonly F86MirrorGeneralisationLink _f86Link;

    /// <summary>The mirror partner of bond b under F71 spatial reflection: N−2−b.
    /// Delegates to F71MirrorSymmetryPi2Inheritance.MirrorPartner-style logic.</summary>
    public int MirrorPartnerBond(int N, int b) => _f71.MirrorPair(N, b);

    /// <summary>True iff bond b is self-paired under F71 mirror: b = N−2−b ⟺ b = (N−2)/2.
    /// Only integer for even N (no self-paired center for odd N). Delegates to
    /// <see cref="F71MirrorSymmetryPi2Inheritance.IsSelfPaired"/>.</summary>
    public bool IsSelfPairedBond(int N, int b) => _f71.IsSelfPaired(N, b);

    /// <summary>Drift check: the F86 mirror partner sum is N−2 for any valid b.</summary>
    public bool MirrorPartnerSumHolds(int N, int b) => b + _f71.MirrorPair(N, b) == N - 2;

    public F86F71MirrorPi2Inheritance(
        F71MirrorSymmetryPi2Inheritance f71,
        F86MirrorGeneralisationLink f86Link)
        : base("F86c sub-claim: F71 spatial-mirror invariance Q_peak(b) = Q_peak(N−2−b); F71 symmetry origin, kinematic identity",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86_QPEAK.md (Statement 3) + " +
               "docs/ANALYTICAL_FORMULAS.md F86 + " +
               "compute/RCPsiSquared.Core/F86/F71MirrorInvariance.cs + " +
               "compute/RCPsiSquared.Core/F71/F86MirrorGeneralisationLink.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F71MirrorSymmetryPi2Inheritance.cs")
    {
        _f71 = f71 ?? throw new ArgumentNullException(nameof(f71));
        _f86Link = f86Link ?? throw new ArgumentNullException(nameof(f86Link));
    }

    public override string DisplayName =>
        "F86 F71-mirror invariance as Pi2-Foundation F71 symmetry inheritance";

    public override string Summary =>
        $"Q_peak(b) = Q_peak(N−2−b) bit-exactly via F71 spatial mirror; verified c=2 N=5..7, c=3 N=5..6; per-F71-orbit substructure refines bond-class dichotomy ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F86c",
                summary: "Q_peak(b) = Q_peak(N−2−b) bit-exactly under F71 spatial mirror; mechanism: K_b is invariant because L_D, S_kernel, Dicke probe are R-symmetric and ∂L/∂J_b ↔ ∂L/∂J_{N−2−b}");
            yield return new InspectableNode("F71 → F86 lift",
                summary: "F71's c₁(N, b, ρ₀) = c₁(N, N−2−b, ρ₀) generalises to F86's K_CC_pr per-bond Q_peak; same kinematic argument, different observable");
            yield return new InspectableNode("verified envelope",
                summary: "c=2 N=5..7 + c=3 N=5..6 + all bond pairs, max deviation < 10⁻¹⁰");
            yield return new InspectableNode("per-F71-orbit substructure caveat",
                summary: "F71 supplies symmetry, not value; central self-paired bonds give different Q_peak than flanking (c=2 N=6: central 1.440 < flanking 1.648; c=3 N=6: central 1.71 > flanking 1.66, opposite direction). Full per-orbit closed form remains open in F86KnowledgeBase as PerF71OrbitObservation.");
            yield return new InspectableNode("F86 multi-mechanism context",
                summary: "this is the third of three F86 Tier-1 sub-claim Pi2-Inheritances. Companions: F86TPeakPi2Inheritance (a_{-1} anchor on t_peak), F86QEpPi2Inheritance (a_0 anchor on Q_EP). Tier-2 / open / empirical parts stay in F86KnowledgeBase typed family.");
            yield return new InspectableNode("verified examples",
                summary: $"N=5 b=0 mirror = {MirrorPartnerBond(5, 0)} (= b_3); N=6 b=2 self-paired? {IsSelfPairedBond(6, 2)} (b=2 = N/2−1=2, yes); N=7 b=1 mirror = {MirrorPartnerBond(7, 1)} (= b_4)");
        }
    }
}
