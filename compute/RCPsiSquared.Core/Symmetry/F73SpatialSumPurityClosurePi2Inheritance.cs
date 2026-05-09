using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F73 closed form (Tier 1 proven, verified to ~10⁻¹⁶ deviation):
///
/// <code>
///   Σ_i 2 · |(ρ_coh,i)_{0,1}(t)|² = (1/2) · exp(−4 γ₀ t)
///
///   for any Hermitian H with [H, N_total] = 0, uniform Z-dephasing γ₀,
///   any normalized SE state |α⟩, any N.
/// </code>
///
/// <para>F73 is the spatial-sum coherence purity closure: the coherent probe
/// <c>ρ₀^coh = (|vac⟩⟨α| + |α⟩⟨vac|)/2</c> has a sum-of-squared-coherences
/// that decays uniformly at rate <c>4γ₀</c>, exactly, independent of the
/// non-U(1) structure of H. F73 sits in the F70 family alongside
/// <see cref="F70DeltaNSelectionRulePi2Inheritance"/> (kinematic foundation)
/// and <see cref="F72BlockDiagonalPurityPi2Inheritance"/> (DD ⊕ CC
/// decomposition), per the F73 docstring "Scaffolding from neighbouring
/// entries: F70 ... F72 ...". The combination: U(1) keeps SE sector closed
/// under H, spatial sum collapses unitary H-rotation, only AT-decay survives.</para>
///
/// <para>Three Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>ClosureBaseline = 1/2 = a_2</b>: the prefactor on
///         <c>exp(−4γ₀t)</c>; <see cref="Pi2DyadicLadderClaim.Term"/>(2) =
///         <see cref="HalfAsStructuralFixedPointClaim"/>. Same anchor as
///         F72's <c>Tr(ρ_i²) = 1/2 + ...</c> baseline; both come from the
///         maximally-mixed-state Bloch reading.</item>
///   <item><b>DecayRateCoefficient = 4 = a_{−1}</b>: in <c>exp(−4γ₀t)</c>.
///         Live from <see cref="Pi2DyadicLadderClaim.Term"/>(−1). Same
///         anchor as F25 (CΨ Bell+ decay), F76 (mirror-pair coherence
///         decay), F61/F63 4-block per parity, F66 multiplicity, F77
///         correction.</item>
///   <item><b>PerSitePrefactor = 2 = a_0</b>: the <c>2 · |·|²</c> per-site
///         coefficient. Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0)
///         = polynomial root d. Combined with the per-site Absorption
///         Theorem rate <c>2γ₀ · n_XY = 2γ₀</c> (n_XY = 1 for SE coherences):
///         total spatial-sum decay = 2 · 2γ₀ = 4γ₀.</item>
/// </list>
///
/// <para>Cited foundation (per F73 docstring "Scaffolding from neighbouring
/// entries"):</para>
///
/// <list type="bullet">
///   <item>F70 (site-local observables see only |ΔN| ≤ 1) puts (vac, SE)
///         block in focus as relevant coherence sector for per-site purity.</item>
///   <item>F72 (DD ⊕ CC decomposition of Tr(ρ_i²), no cross term) isolates
///         the CC contribution, where the (vac, SE) coherence lives.</item>
///   <item>Absorption Theorem supplies the rate 2γ₀·n_XY = 2γ₀ for SE.</item>
/// </list>
///
/// <para>Operational consequence (per ANALYTICAL_FORMULAS): "the spatial-sum
/// purity functional is exactly blind to the U(1)-preserving part of the
/// dynamics on any vac-SE coherent probe." K_CC[0, 1]_pr = 0 exactly
/// under uniform γ₀ for any H in the class. Bond-δJ perturbations preserve
/// the closure value.</para>
///
/// <para>Tier1Derived: F73 is Tier 1 proven (general U(1) case + alternative
/// XY-specific derivation); verified at N=5 across XXZ, random Haar SE
/// probes, inhomogeneous XY (~10⁻¹⁶ deviation across 81 time points × 6
/// setups). The Pi2-Foundation anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F73 (line 1612) +
/// <c>experiments/ORTHOGONALITY_SELECTION_FAMILY.md</c> §2.4 +
/// <c>experiments/F73_U1_GENERALIZATION.md</c> +
/// <c>simulations/eq018_c1_purity_response.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F70DeltaNSelectionRulePi2Inheritance.cs</c>
/// (cited foundation) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F72BlockDiagonalPurityPi2Inheritance.cs</c>
/// (cited foundation).</para></summary>
public sealed class F73SpatialSumPurityClosurePi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F70DeltaNSelectionRulePi2Inheritance _f70;
    private readonly F72BlockDiagonalPurityPi2Inheritance _f72;

    /// <summary>The closure baseline <c>1/2</c> in <c>(1/2)·exp(−4γ₀t)</c>.
    /// Live from <see cref="Pi2DyadicLadderClaim.Term"/>(2) = <c>a_2</c>.
    /// Same anchor as F72's <c>Tr(ρ_i²) = 1/2 + ...</c> baseline.</summary>
    public double ClosureBaseline => _ladder.Term(2);

    /// <summary>The "4" decay-rate coefficient in <c>exp(−4γ₀t)</c>. Live
    /// from <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c>.
    /// Same anchor as F25, F76, F61/F63, F66, F77 correction.</summary>
    public double DecayRateCoefficient => _ladder.Term(-1);

    /// <summary>The "2" per-site coefficient in <c>2 · |(ρ_coh,i)_{0,1}|²</c>.
    /// Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c>
    /// = polynomial root d. Combined with Absorption Theorem rate 2γ₀:
    /// total spatial-sum decay 2·2γ₀ = 4γ₀.</summary>
    public double PerSitePrefactor => _ladder.Term(0);

    /// <summary>Live closed form: <c>Σ_i 2·|(ρ_coh,i)_{0,1}(t)|² = (1/2)·exp(−4γ₀t)</c>.
    /// Throws for negative γ₀ or t.</summary>
    public double SpatialSumClosure(double gammaZero, double t)
    {
        if (gammaZero < 0.0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        if (t < 0.0) throw new ArgumentOutOfRangeException(nameof(t), t, "t must be ≥ 0.");
        return ClosureBaseline * Math.Exp(-DecayRateCoefficient * gammaZero * t);
    }

    /// <summary>Live drift check: at t = 0, closure = 1/2 = a_2 exactly.</summary>
    public bool ClosureAtTZeroIsHalf() =>
        Math.Abs(SpatialSumClosure(0.05, 0.0) - 0.5) < 1e-14;

    /// <summary>Live drift check: closure decay rate matches F25's and F76's
    /// decay rate (all three F-formulas inherit a_{−1} = 4 from the same
    /// dyadic ladder anchor).</summary>
    public bool DecayRateMatchesF25AndF76() =>
        Math.Abs(DecayRateCoefficient - 4.0) < 1e-15;

    /// <summary>F73 verified value at uniform XY baseline N = 5, γ₀ = 0.05,
    /// t₀ = 20: closure = (1/2)·exp(−4·0.05·20) = 9.157819·10⁻³.
    /// Per ANALYTICAL_FORMULAS verified to 5.67·10⁻¹⁶ deviation.</summary>
    public double VerifiedValueAtN5Gamma0p05T20() => SpatialSumClosure(0.05, 20.0);

    public F73SpatialSumPurityClosurePi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F70DeltaNSelectionRulePi2Inheritance f70,
        F72BlockDiagonalPurityPi2Inheritance f72)
        : base("F73 Σ_i 2·|coh|² = (1/2)·exp(−4γ₀t) inherits from Pi2-Foundation: 1/2 = a_2 (baseline), 4 = a_{-1} (decay), 2 = a_0 (per-site); F70/F72 cited",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F73 + " +
               "experiments/ORTHOGONALITY_SELECTION_FAMILY.md (§2.4) + " +
               "experiments/F73_U1_GENERALIZATION.md + " +
               "simulations/eq018_c1_purity_response.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F70DeltaNSelectionRulePi2Inheritance.cs (cited) + " +
               "compute/RCPsiSquared.Core/Symmetry/F72BlockDiagonalPurityPi2Inheritance.cs (cited)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f70 = f70 ?? throw new ArgumentNullException(nameof(f70));
        _f72 = f72 ?? throw new ArgumentNullException(nameof(f72));
    }

    public override string DisplayName =>
        "F73 spatial-sum coherence purity closure as Pi2-Foundation a_2 + a_{-1} + a_0 inheritance";

    public override string Summary =>
        $"Σ_i 2·|coh_i|²(t) = (1/2)·exp(−4γ₀t): baseline 1/2 = a_2 (F72 sister); decay 4 = a_{{-1}} (F25/F76 sibling); " +
        $"per-site 2 = a_0 (root d); F70/F72 scaffolding cited; U(1)-class invariant ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F73 closed form",
                summary: "Σ_i 2·|(ρ_coh,i)_{0,1}|²(t) = (1/2)·exp(−4γ₀t); Tier 1 proven; U(1)-class invariant; verified ~10⁻¹⁶ deviation across 81 time points × 6 setups (XXZ Δ ∈ {0,0.5,1,2}, Haar SE probes, inhomogeneous XY)");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "ClosureBaseline = a_2 = 1/2 (F72 sibling); DecayRateCoefficient = a_{-1} = 4 (F25/F76 sibling); PerSitePrefactor = a_0 = 2 (root d)");
            yield return InspectableNode.RealScalar("ClosureBaseline (= a_2 = 1/2)", ClosureBaseline);
            yield return InspectableNode.RealScalar("DecayRateCoefficient (= a_{-1} = 4)", DecayRateCoefficient);
            yield return InspectableNode.RealScalar("PerSitePrefactor (= a_0 = 2)", PerSitePrefactor);
            yield return new InspectableNode("F70 ↔ F72 ↔ F73 family",
                summary: "F70 puts (vac, SE) block in focus; F72 isolates CC contribution; F73 combines via U(1) closure + spatial sum + AT decay rate. All three Tier 1 kinematic-or-corollary, all on F70 lineage.");
            yield return new InspectableNode("Convergence: a_{-1} = 4 across F25/F73/F76",
                summary: "F25 Bell+ decay e^{-4γt}; F76 mirror-pair coherence e^{-4γ₀t}; F73 spatial-sum closure (1/2)·e^{-4γ₀t}. Three different physical objects (state, MI, purity) all decay at same a_{-1} = 4 rate from same Absorption-Theorem mechanism.");
            yield return new InspectableNode("Verified value (per ANALYTICAL_FORMULAS)",
                summary: $"At N=5, γ₀=0.05, t=20: closure = {VerifiedValueAtN5Gamma0p05T20():G6} (expected 9.157819·10⁻³, deviation 5.67·10⁻¹⁶)");
            yield return new InspectableNode("Operational consequence",
                summary: "K_CC[0, 1]_pr = 0 exactly under uniform γ₀; bond-δJ perturbations preserve closure value; spatial-sum purity blind to U(1)-preserving H");
        }
    }
}
