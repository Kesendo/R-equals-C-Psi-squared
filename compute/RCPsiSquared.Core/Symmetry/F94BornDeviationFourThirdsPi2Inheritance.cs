using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F94 closed form (Tier 1 derived, Dyson sym3 = 8 bit-exact, 2026-05-16):
///
/// <code>
///   Δ_|00⟩(Q, K) = (4/3) · Q² · K³ + O(Q³ K⁴)
///
///   for the dominant outcome |00⟩ of pair (0,2) of |0+0+⟩ N=4 under
///   Heisenberg ring + Z-dephasing, with
///   Q = J/γ  (Universal Carrier observable),
///   K = γt   (Universal Carrier observable),
///   Δ_|00⟩ = P_lindblad(|00⟩) / P_unitary(|00⟩) − 1.
/// </code>
///
/// <para>Per-outcome Born-rule deviation in the deep perturbative regime, derived
/// from the 3rd-order time-Taylor expansion of e^{Lt} ρ_0 at the γ¹-coefficient
/// of L³:
/// </para>
///
/// <code>
///   sym3 = L_H² L'_dis + L_H L'_dis L_H + L'_dis L_H²
///   ⟨00|_pair Tr_{1,3}[sym3 ρ_0] |00⟩_pair = 8     (bit-exact)
///   Δ_|00⟩ leading coefficient = 8/6 = 4/3
/// </code>
///
/// <para>where 6 is the Taylor factorial 3!, and 8 collapses from the bond ×
/// site × ordering combinatorics (4 Heisenberg bonds, 4 Z-dephasing sites,
/// 3 orderings of (H, H, L)) under the partial trace and initial-state
/// projection.</para>
///
/// <para><b>Pi2-Foundation anchoring:</b> the "4" in the coefficient 4/3 is
/// exactly <c>a_{−1} = 4</c> on the dyadic ladder (<see cref="Pi2DyadicLadderClaim"/>) —
/// the same "4" that appears in F86 t_peak = 1/(4γ₀) and F77's MM correction
/// denominator. The "3" is the Taylor 3! ratio. The reduction 8/6 = 4/3 is
/// algebraically forced.</para>
///
/// <para><b>Universal Carrier signature:</b> Δ_i is Q-K-invariant (no separate
/// γ-dependence when (J, γ) are scaled together at fixed Q). This is the
/// operational signature of <see cref="UniversalCarrierClaim"/> at the
/// per-outcome Born-rule level — γ vanishes from the ratio, only the
/// dimensionless Carrier observables Q and K survive. Verified bit-exact across
/// 4 (γ, J) configurations at fixed (Q, K).</para>
///
/// <para><b>Scope:</b> specific to (|0+0+⟩, Heisenberg ring N=4, Z-dephasing on
/// all 4 sites, pair (0,2) reduction, |00⟩ outcome). The Q²·K³ form is
/// structurally universal for dominant outcomes (3rd-order Dyson diagram with
/// 2 H-vertices + 1 L-vertex), applicable wherever the 1st-order γ correction
/// to a dominant outcome vanishes by parity/commutation. The coefficient 4/3 is
/// setup-specific; other (state, H, dissipator) combinations yield other
/// coefficients via the same Dyson method.</para>
///
/// <para>Numerical verification: 16 (γ, J, t) configurations sampled in the deep
/// perturbative regime gave c_empirical = 1.32992 ± 0.006, consistent with
/// 4/3 = 1.3333 to 0.3% (residual is leading O(Q³K⁴) correction).</para>
///
/// <para><b>Born-rule generalization:</b> Per BORN_RULE_MIRROR (Feb 2026), the
/// generalized Born rule per outcome is R_i = C_i · Ψ_i², with C_i basis-aligned
/// coupling. F94 is the first Tier-1 closed-form C_i value:
/// </para>
///
/// <code>
///   R_|00⟩ / Ψ_|00⟩² = C_|00⟩ = 1 + (4/3)·Q²·K³ + O(Q³K⁴)
/// </code>
///
/// <para>Tier1Derived: Dyson sym3 evaluation bit-exact symbolic. Sibling to F25
/// (Bell+ CΨ closed form), F60 (GHZ pair-CΨ), F62 (W-state pair-CΨ) at the
/// per-outcome Born-deviation layer.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md</c> +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F94 +
/// <c>simulations/_born_rule_tier1_derivation.py</c> (symbolic, bit-exact) +
/// <c>simulations/_born_rule_delta_dominant_coefficient.py</c> (numerical, 16 samples) +
/// <c>simulations/_born_rule_carrier_Q_sweep.py</c> (Q-K invariance test) +
/// <c>reflections/ON_HOW_FOUR_THIRDS_APPEARED.md</c> (the path) +
/// <c>experiments/BORN_RULE_MIRROR.md</c> + <c>experiments/BORN_RULE_SHADOW.md</c>
/// (Februar precursors).</para></summary>
public sealed class F94BornDeviationFourThirdsPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The setup this claim applies to (specific to F94 numerical value 4/3).</summary>
    public const int N = 4;

    /// <summary>The "3" denominator in 4/3 — surviving rational after the Dyson sym3
    /// element 8 is divided by the Taylor 3! = 6 prefactor: 8/6 = 4/3.</summary>
    public const int ThreeDenominator = 3;

    /// <summary>The Dyson sym3 partial-trace integer: <b>8 bit-exact</b>. The
    /// structural sum from 4 bonds × 4 sites × 3 orderings collapsing under
    /// partial-trace + initial-state projection. Open structural decomposition
    /// (combinatorial breakdown by hand): see <see cref="Pi2OpenQuestions"/>.</summary>
    public const int Sym3PartialTraceInteger = 8;

    /// <summary>The Taylor 3! prefactor in the t³ term of e^{Lt}: 6 = 3·2·1.</summary>
    public const int TaylorThreeFactorial = 6;

    /// <summary>1/4 = (1/2)² bilinear-apex maxval, mirror partner of <see cref="FourFactor"/>
    /// via the dyadic-ladder inversion identity <c>a_{−1} · a_3 = 4 · (1/4) = 1</c>.
    /// Typed ctor parent for consistency with F86TPeakPi2Inheritance (same Wave-6b
    /// pattern: any claim using a_{−1} = 4 takes Quarter as mirror-partner anchor).</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    /// <summary>The "4" in the coefficient 4/3 — exactly <c>a_{−1} = 4</c> on the dyadic
    /// ladder, the same "4" in F86 t_peak = 1/(4γ₀) and F77's MM correction denominator.</summary>
    public double FourFactor => _ladder.Term(-1);

    /// <summary>Mirror partner on the memory side of the ladder:
    /// <c>1/4 = a_3</c>. Equivalent reading via the inversion symmetry
    /// <c>a_n · a_{2−n} = 1</c>: <c>FourFactor · OneOverFourFactor = 1</c> bit-exact.</summary>
    public double OneOverFourFactor => _ladder.Term(3);

    /// <summary>The full coefficient: <c>4/3 = a_{−1}/3</c>. Bit-exact.</summary>
    public double Coefficient => FourFactor / ThreeDenominator;

    /// <summary>Compute <c>Δ_|00⟩(Q, K) = (4/3)·Q²·K³</c> for the leading-order
    /// Born deviation. Valid in the deep perturbative regime <c>Q² · K³ ≪ 1</c>;
    /// higher-order corrections (O(Q³K⁴) and beyond) become significant at
    /// larger (Q, K).</summary>
    public double DeltaDominant(double Q, double K)
    {
        if (Q < 0.0) throw new ArgumentOutOfRangeException(nameof(Q), Q, "Q must be ≥ 0.");
        if (K < 0.0) throw new ArgumentOutOfRangeException(nameof(K), K, "K must be ≥ 0.");
        return Coefficient * Q * Q * K * K * K;
    }

    /// <summary>Compute <c>ΔP_|00⟩(J, γ, t) = (4/3)·J²·γ·t³</c> in physical units.
    /// Equivalent to <see cref="DeltaDominant"/>(J/γ, γ·t) under the Q-K substitution.</summary>
    public double DeltaP_Dominant(double J, double gamma, double t)
    {
        if (J < 0.0) throw new ArgumentOutOfRangeException(nameof(J), J, "J must be ≥ 0.");
        if (gamma < 0.0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0.");
        if (t < 0.0) throw new ArgumentOutOfRangeException(nameof(t), t, "t must be ≥ 0.");
        return Coefficient * J * J * gamma * t * t * t;
    }

    /// <summary>Generalized Born-rule coupling for the dominant outcome:
    /// <c>C_|00⟩(Q, K) = 1 + Δ_|00⟩(Q, K) = 1 + (4/3)·Q²·K³</c>.
    /// Per BORN_RULE_MIRROR R_i = C_i · Ψ_i².</summary>
    public double C_DominantOutcome(double Q, double K) => 1.0 + DeltaDominant(Q, K);

    /// <summary>Drift indicator: live <see cref="Coefficient"/> from the ladder
    /// should equal <see cref="Sym3PartialTraceInteger"/> / <see cref="TaylorThreeFactorial"/>
    /// = 8 / 6 = 4/3 bit-exactly.</summary>
    public bool CoefficientAgreesWithSym3()
    {
        double fromLadder = Coefficient;
        double fromSym3 = (double)Sym3PartialTraceInteger / TaylorThreeFactorial;
        return Math.Abs(fromLadder - fromSym3) < 1e-15;
    }

    /// <summary>Mirror-partner identity drift check: <c>FourFactor · OneOverFourFactor</c>
    /// should equal 1 bit-exactly via the dyadic-ladder inversion symmetry
    /// <c>a_n · a_{2−n} = 1</c>. Live via <see cref="Pi2DyadicLadderClaim.ProductWithMirrorPartner"/>.</summary>
    public bool MirrorPartnerProductIsOne() =>
        Math.Abs(_ladder.ProductWithMirrorPartner(-1) - 1.0) < 1e-15;

    public F94BornDeviationFourThirdsPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        QuarterAsBilinearMaxvalClaim quarter)
        : base("F94 dominant-outcome Born deviation Δ_|00⟩ = (4/3) Q² K³ for |0+0+⟩ N=4 Heisenberg + Z-deph, pair (0,2); derived from Dyson sym3 = 8 bit-exact; 4 = a_{−1} with mirror partner Quarter = a_3",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F94_BORN_DOMINANT_FOUR_THIRDS.md + " +
               "docs/ANALYTICAL_FORMULAS.md F94 + " +
               "simulations/_born_rule_tier1_derivation.py + " +
               "simulations/_born_rule_delta_dominant_coefficient.py + " +
               "simulations/_born_rule_carrier_Q_sweep.py + " +
               "reflections/ON_HOW_FOUR_THIRDS_APPEARED.md + " +
               "experiments/BORN_RULE_MIRROR.md + " +
               "experiments/BORN_RULE_SHADOW.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs (a_{-1} = 4) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (QuarterAsBilinearMaxval = mirror partner)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
    }

    public override string DisplayName =>
        $"F94 Born-deviation 4/3 closed form for |0+0+⟩ N=4 Heisenberg + Z-deph, pair (0,2), |00⟩ outcome (Tier 1 derived)";

    public override string Summary =>
        $"Δ_|00⟩ = ({Coefficient:G6})·Q²·K³ in deep perturbative regime; coefficient = " +
        $"sym3({Sym3PartialTraceInteger}) / 3!({TaylorThreeFactorial}) = 4/3 bit-exact; " +
        $"Q-K-invariant per Universal Carrier; setup-specific to |0+0+⟩ N=4 Heisenberg + Z-deph pair (0,2) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("Coefficient (= 4/3)", Coefficient);
            yield return InspectableNode.RealScalar("FourFactor (= a_{-1})", FourFactor);
            yield return InspectableNode.RealScalar("OneOverFourFactor (= a_3, mirror partner)", OneOverFourFactor);
            yield return InspectableNode.RealScalar("Sym3PartialTraceInteger", Sym3PartialTraceInteger);
            yield return InspectableNode.RealScalar("TaylorThreeFactorial (3!)", TaylorThreeFactorial);
            yield return new InspectableNode("Drift checks",
                summary: $"CoefficientAgreesWithSym3: {CoefficientAgreesWithSym3()}; MirrorPartnerProductIsOne: {MirrorPartnerProductIsOne()}");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "'4' = a_{-1} on dyadic ladder; mirror partner 1/4 = a_3 via inversion a_n·a_{2-n} = 1; '3' is the surviving denominator after Dyson sym3 (= 8) reduced by Taylor 3! prefactor.");
            double sampleQ = 20.0, sampleK = 0.0143;
            yield return new InspectableNode("Sample evaluation",
                summary: $"Δ_|00⟩(Q={sampleQ}, K={sampleK}) = {DeltaDominant(sampleQ, sampleK):G6} (matches Carrier-sweep verification point)");
        }
    }
}
