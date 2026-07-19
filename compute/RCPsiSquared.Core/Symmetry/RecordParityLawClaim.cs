using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F135, typed (Tier 1 derived): <b>the record parity law</b>. Under
/// H = Σ Δ_ab·Z_aZ_b (Pauli), local Z-dephasing at per-site rates γ_l, initial |+⟩^⊗N, every
/// two-site reduced state is closed-form for all t (Proposition 1, the pair reduction: the
/// Absorption-Theorem substrate |A⟩⟨B|(t) read on a pair page):
///
/// <code>
///   ρ_ab[(z_a,z_b),(z'_a,z'_b)](t) = ¼ · e^{−it·Δ_ab·(z_az_b−z'_az'_b)}
///       · e^{−2t·(γ_a·𝟙[z_a≠z'_a] + γ_b·𝟙[z_b≠z'_b])}
///       · Π_{k∉{a,b}} cos(t·[Δ_ak·(z_a−z'_a) + Δ_bk·(z_b−z'_b)])
/// </code>
///
/// Traced-site γ drops out identically. At the symmetric readout t* = π/(4Δ_S) (all S-bonds
/// at Δ_S, deg(S) ≥ 2, triangle-free at S) whether a witness j ∈ N(S) records the pointer
/// Z_S is the ARITHMETIC of its watcher ratios r_k = Δ_jk/Δ_S:
///
/// <code>
///   I(S:E_j)(t*) = 1 − h₂((1+β_j)/2),   β_j = e^{−2γ_j·t*} · Π_{k∈N(j)∖{S}} |cos((π/2)·r_k)|
/// </code>
///
/// Law A (parity trichotomy): every watcher even → I = 1 exactly (forgiven, watched or not);
/// any watcher odd → I = 0 exactly (blind); non-integer → the generic in-between. Law B:
/// R_perfect = #{even-aligned neighbors} ≤ deg(S), achievable on any triangle-free graph
/// (uniform coupling makes every ratio odd: the leaf law is a parity law, not geometry).
/// Law B′: the uniform 4-cycle's opposite corner holds a perfect distance-2 bit of the
/// ANTI-pointer X_S (invisible to pointer redundancy). Law C (the γ race):
/// D_j(t) = e^{−2γ_j·t}·|sin(2Δ_Sj·t)|·Π|cos(2Δ_jk·t)| exactly, peaked at
/// t_opt = arctan(Δ/γ_j)/(2Δ) &lt; t*; γ_S and traced-site γ never enter.
///
/// <para><b>Typed parent.</b> <see cref="AbsorptionTheoremClaim"/>: Proposition 1 is the
/// absorption substrate (the closed dephasing of every |A⟩⟨B| coherence at rate 2γ·n_diff)
/// evaluated on a pair page; the record radius β_j is the surviving coherence magnitude.
/// Gates: <c>simulations/qd_pointer_opt.py</c> (22/22 vs RK4 at N = 8) and
/// <c>simulations/qd_scout.py</c> (the door). Live: <c>inspect --root record</c>
/// (<c>RecordLawWitness</c>, the pair battery recomputed closed-form-vs-full-state at
/// inspect time).</para></summary>
public sealed class RecordParityLawClaim : Claim
{
    // Parent-edge marker for Schicht-1 wiring: Proposition 1 is the absorption substrate on a pair page.
    public AbsorptionTheoremClaim Absorption { get; }

    public RecordParityLawClaim(AbsorptionTheoremClaim absorption)
        : base("The record parity law (F135): under H = sum Delta*ZZ with local Z-dephasing from " +
               "|+>^N, every two-site reduced state is closed-form (the pair reduction, the " +
               "Absorption substrate on a pair page); at the symmetric readout whether a witness " +
               "records the pointer is the parity arithmetic of its watcher ratios: every watcher " +
               "even -> a perfect record (watched or not), any watcher odd -> exactly blind, " +
               "non-integer -> the explicit in-between 1 - h2((1+beta)/2); R_perfect <= deg(S) " +
               "with equality achievable on any triangle-free graph; the uniform 4-cycle records " +
               "the anti-pointer at distance 2; the gamma race is exact and traced-site gamma " +
               "never enters",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_RECORD_PARITY_LAW.md + " +
               "experiments/QUANTUM_DARWINISM_POINTER_DOOR.md")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
    }

    /// <summary>The readout point t* = π/(4Δ_S) at Δ_S = 1.</summary>
    public static double TStar => Math.PI / 4.0;

    /// <summary>Law A's mutual information from the record radius: I = 1 − h₂((1+β)/2), in bits.</summary>
    public static double LawAInformation(double beta) => 1.0 - Entropy.Binary((1.0 + beta) / 2.0);

    /// <summary>The record radius β_j at t* for watcher ratios rs and witness rate γ_j:
    /// β = e^{−2γ_j·t*}·Π|cos((π/2)·r_k)|.</summary>
    public static double RecordRadius(IEnumerable<double> watcherRatios, double gammaJ = 0.0)
    {
        double beta = Math.Exp(-2.0 * gammaJ * TStar);
        foreach (double r in watcherRatios)
            beta *= Math.Abs(Math.Cos(Math.PI / 2.0 * r));
        return beta;
    }

    public override string DisplayName =>
        "Record parity law (F135): who records is the arithmetic of the angle ratios";

    public override string Summary =>
        "the pair reduction is exact for all t on any graph (Proposition 1, the Absorption substrate " +
        "on a pair page; traced-site γ drops out); at t* = π/(4Δ_S) the watcher-ratio parities decide: " +
        "even → perfect record, odd → exactly blind, non-integer → 1 − h₂((1+β)/2); R_perfect ≤ deg(S); " +
        "the plaquette records the anti-pointer; the γ race D_j(t) is exact " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Proposition 1: the pair reduction",
                summary: "every entry of every two-site reduced state factors into the direct-bond " +
                         "phase, the pair's own dephasing, and one cosine per third site; only γ_a, γ_b " +
                         "survive (watching the writers is free). The Absorption Theorem's |A⟩⟨B| " +
                         "substrate read on a pair page.");
            yield return new InspectableNode("Law A: the parity trichotomy",
                summary: "I(S:E_j)(t*) = 1 − h₂((1+β_j)/2), β_j = e^{−2γ_j t*}·Π|cos((π/2)·r_k)|: " +
                         "every watcher ratio even → I = 1 exactly; any odd → I = 0 exactly; " +
                         "non-integer → generic (r = 3/2 → 0.399124). The leaf law is the " +
                         "empty-watcher column; uniform coupling makes every ratio odd.");
            yield return new InspectableNode("Law B and B′: alignment beats topology; the plaquette",
                summary: "R_perfect = #{even-aligned neighbors} ≤ deg(S), equality on any triangle-free " +
                         "graph; on trees non-neighbors are exactly 0 at t*. The uniform 4-cycle's " +
                         "opposite corner holds a perfect distance-2 bit of the ANTI-pointer X_S " +
                         "(⟨X_SX_j⟩ = 1, zero pointer content): cycles add anti-pointer bits, not " +
                         "pointer records.");
            yield return new InspectableNode("Law C: the γ race",
                summary: "D_j(t) = e^{−2γ_j t}·|sin(2Δ_Sj t)|·Π|cos(2Δ_jk t)| exactly; the " +
                         "distinguishability peak sits at t_opt = arctan(Δ/γ_j)/(2Δ), earlier than t*; " +
                         "at the canonical γ/Δ = 0.05 the record read at t* degrades to I = 0.768040. " +
                         "γ_S never touches D_j.");
            yield return new InspectableNode("typed parent",
                summary: $"AbsorptionTheoremClaim ({Absorption.Tier.Label()}): the closed dephasing of " +
                         "every |A⟩⟨B| coherence, whose pair-page face IS Proposition 1; the record " +
                         "radius is the surviving coherence magnitude.");
            yield return new InspectableNode("live witness (inspect --root record)",
                summary: "RecordLawWitness recomputes the pair battery at inspect time: the F136 " +
                         "classifier's closed-form verdicts against an independent full-state path " +
                         "(2^N closed-form density matrix, partial trace, measured I/correlators), " +
                         "including this law's Law-A in-between and price faces.");
        }
    }
}
