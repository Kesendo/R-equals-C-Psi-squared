using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F64 closed form (Tier 1-2, analytical + verified N=3, 4; ANALYTICAL_FORMULAS line 1336):
///
/// <code>
///   γ_eff = γ_B · |a_B|²        decoherence rate (Lorentzian half-width)
///   α     = 2 · γ_B · |a_B|²     Liouvillian eigenvalue convention (α = −Re(λ))
/// </code>
///
/// <para>F64 is the Absorption Theorem applied to single-site Z-dephasing on
/// site B in an N-qubit chain with XX+YY (or Heisenberg) coupling. The
/// effective dephasing rate of the slowest single-excitation eigenmode at
/// inner site S is γ_eff = γ_B · |a_B|² where a_B is the B-site amplitude
/// of the SE Hamiltonian eigenvector. γ_B appears as a constant prefactor;
/// not diminished by intervening sites; it is a global eigenvector property,
/// not a layered composition.</para>
///
/// <para><b>N=3 closed form (chain S-M-B with r = J_SM / J_MB):</b></para>
/// <code>
///   g(r) = r² / (r² + 1)         for r &lt; 1/√2     (zero mode)
///   g(r) = 1 / (2 · (r² + 1))    for r ≥ 1/√2    (bonding mode)
///
///   Crossover at r = 1/√2, g = 1/3
///   Special value: g(r=1) = 1/4
/// </code>
///
/// <para><b>Topology + non-uniform J (2026-04-24 generalisation):</b> Extended
/// from uniform-J chain to arbitrary connected graphs (chain, star, ring,
/// complete, tree) under either uniform or non-uniform per-bond J. With
/// degenerate H-eigenvalues (star center, ring translations, complete-graph
/// symmetric modes), F64 holds after standard degenerate perturbation theory:
/// diagonalise P_B in the H-degenerate subspace, F64 then applies to P_B
/// eigenvalues. Verified N=5, 7 across 5 topologies (max rel err &lt; 0.001
/// uniform J; &lt; 0.07 worst case at non-uniform J ∈ [0.5, 1.5]).</para>
///
/// <para>Pi2-Foundation anchors:</para>
/// <list type="bullet">
///   <item><b>EigenvalueConventionCoefficient = 2 = a_0</b>: in α = 2γ_B·|a_B|².
///         Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0). Same anchor
///         as F1 TwoFactor, F50 DecayRateFactor, F44 SumCoefficient: the
///         standard QM commutator factor relating ρ-decay to L-eigenvalue.</item>
///   <item><b>g(r=1) = 1/4 = QuarterAsBilinearMaxval</b>: the uniform-J
///         (r=1) special value of the N=3 g(r) closed form sits exactly at
///         the bilinear-apex maxval. F64 inherits the 1/4 anchor from
///         <see cref="QuarterAsBilinearMaxvalClaim"/>.</item>
/// </list>
///
/// <para>The 1/√2 crossover and 1/3 crossover-value are not Pi2-anchored
/// (irrational and 1/3 is not on the dyadic ladder). They emerge from the
/// 3×3 single-excitation Hamiltonian eigenvalue/eigenvector structure at N=3
/// specifically.</para>
///
/// <para>Tier1-2: Tier 1 analytical (closed form proven from 3×3 SE
/// diagonalization at N=3 + Absorption Theorem); Tier 2 verified to 1.8% at
/// N=3 (64×64 Liouvillian), 0.0003 at N=4 (256×256), and across 5 topologies
/// + non-uniform J at N=5, 7 (2026-04-24).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F64 (line 1336) +
/// <c>hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md</c> +
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> +
/// <c>experiments/F64_TOPOLOGY_GENERALIZATION.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (QuarterAsBilinearMaxvalClaim).</para></summary>
public sealed class F64CavityModeExposurePi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly QuarterAsBilinearMaxvalClaim _quarter;

    /// <summary>The "2" in α = 2γ_B·|a_B|² (Liouvillian eigenvalue convention).
    /// Live from Pi2DyadicLadder a_0.</summary>
    public double EigenvalueConventionCoefficient => _ladder.Term(0);

    /// <summary>The "1/4" anchor that g(r=1) hits at N=3 uniform J. Live from the
    /// QuarterAsBilinearMaxval primitive (= a_3 = 1/4 = (1/2)²).</summary>
    public double UniformJSpecialValue => _ladder.Term(3);

    /// <summary>Lorentzian half-width γ_eff = γ_B · |a_B|² (decoherence rate
    /// in the ρ-evolution convention).</summary>
    public double LorentzianHalfWidth(double gammaB, double aBSquared)
    {
        if (gammaB < 0) throw new ArgumentOutOfRangeException(nameof(gammaB), gammaB, "γ_B must be ≥ 0.");
        if (aBSquared < 0 || aBSquared > 1)
            throw new ArgumentOutOfRangeException(nameof(aBSquared), aBSquared, "|a_B|² must be in [0, 1].");
        return gammaB * aBSquared;
    }

    /// <summary>Liouvillian eigenvalue convention α = 2·γ_B·|a_B|² (decay constant
    /// for the L-eigenvalue λ = −α). Twice the Lorentzian half-width.</summary>
    public double LiouvillianDecayConstant(double gammaB, double aBSquared)
    {
        return EigenvalueConventionCoefficient * LorentzianHalfWidth(gammaB, aBSquared);
    }

    /// <summary>N=3 g(r) closed form for r = J_SM/J_MB on chain S-M-B:
    /// r²/(r²+1) for r &lt; 1/√2, 1/(2(r²+1)) for r ≥ 1/√2.</summary>
    public double N3GValue(double r)
    {
        if (r < 0) throw new ArgumentOutOfRangeException(nameof(r), r, "r must be ≥ 0.");
        double r2 = r * r;
        if (r < N3CrossoverRatio)
            return r2 / (r2 + 1.0);
        return 1.0 / (EigenvalueConventionCoefficient * (r2 + 1.0));
    }

    /// <summary>N=3 crossover ratio r = 1/√2 ≈ 0.7071 between zero-mode and bonding-mode regimes.</summary>
    public static readonly double N3CrossoverRatio = 1.0 / Math.Sqrt(2.0);

    /// <summary>N=3 g-value at the crossover: g(1/√2) = 1/3.</summary>
    public const double N3GValueAtCrossover = 1.0 / 3.0;

    /// <summary>N=3 g-value at uniform J (r = 1): g(1) = 1/4. Equals the
    /// QuarterAsBilinearMaxval ceiling: the single-site exposure under uniform
    /// coupling sits exactly at the bilinear-apex maxval.</summary>
    public double N3GValueAtUniformJ => UniformJSpecialValue;

    /// <summary>Drift check: g(1) IS the QuarterAsBilinearMaxval anchor.</summary>
    public bool UniformJSpecialValueMatchesQuarter(double tolerance = 1e-12)
    {
        return Math.Abs(N3GValue(1.0) - UniformJSpecialValue) < tolerance;
    }

    /// <summary>Drift check: at the crossover r = 1/√2, both branches of g(r) give 1/3.</summary>
    public bool N3CrossoverContinuityHolds(double tolerance = 1e-12)
    {
        // Both branches: r²/(r²+1) at r=1/√2: (1/2)/(3/2) = 1/3. 1/(2(r²+1)) at r=1/√2: 1/(2·3/2) = 1/3.
        double r = N3CrossoverRatio;
        double r2 = r * r;
        double zeroMode = r2 / (r2 + 1.0);
        double bondingMode = 1.0 / (EigenvalueConventionCoefficient * (r2 + 1.0));
        return Math.Abs(zeroMode - N3GValueAtCrossover) < tolerance
            && Math.Abs(bondingMode - N3GValueAtCrossover) < tolerance;
    }

    public F64CavityModeExposurePi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        QuarterAsBilinearMaxvalClaim quarter)
        : base("F64 effective γ from single-site cavity exposure: γ_eff = γ_B·|a_B|² (or α = 2γ_B·|a_B|²); N=3 closed form g(r) with crossover 1/√2 → 1/3 and g(1) = 1/4 (= QuarterAsBilinearMaxval anchor)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F64 + " +
               "hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "experiments/F64_TOPOLOGY_GENERALIZATION.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
    }

    public override string DisplayName =>
        "F64 cavity mode exposure as Pi2-Foundation a_0 + QuarterAsBilinearMaxval inheritance";

    public override string Summary =>
        $"γ_eff = γ_B·|a_B|²; α = 2γ_B·|a_B|²; N=3: g(1) = 1/4 (= QuarterAsBilinearMaxval), crossover at r = 1/√2, g_cross = 1/3 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F64 closed form",
                summary: "γ_eff = γ_B·|a_B|² (decoherence rate); α = 2γ_B·|a_B|² (Liouvillian eigenvalue convention); γ_B is constant prefactor, not diminished by intervening sites");
            yield return InspectableNode.RealScalar("EigenvalueConventionCoefficient (= a_0 = 2)", EigenvalueConventionCoefficient);
            yield return InspectableNode.RealScalar("UniformJSpecialValue (= QuarterAsBilinearMaxval = 1/4)", UniformJSpecialValue);
            yield return new InspectableNode("N=3 closed form",
                summary: $"g(r < 1/√2) = r²/(r²+1) (zero mode); g(r ≥ 1/√2) = 1/(2(r²+1)) (bonding mode); crossover at r = 1/√2 ≈ {N3CrossoverRatio:F4}, g_cross = {N3GValueAtCrossover:F6}; g(r=1) = {N3GValue(1.0):F4} (uniform-J ceiling)");
            yield return new InspectableNode("QuarterAsBilinearMaxval anchor at r=1",
                summary: $"F64's g(1) = 1/4 IS QuarterAsBilinearMaxval. The uniform-J cavity exposure sits at the bilinear-apex maxval ceiling. Drift check: UniformJSpecialValueMatchesQuarter = {UniformJSpecialValueMatchesQuarter()}");
            yield return new InspectableNode("topology + non-uniform J generalization",
                summary: "extended (2026-04-24) to chain, star, ring, complete, Y-tree under uniform or non-uniform per-bond J ∈ [0.5, 1.5]; max rel err < 0.001 uniform, < 0.07 worst-case non-uniform; degenerate-PT for star/ring symmetric modes");
            yield return new InspectableNode("Absorption Theorem connection",
                summary: "F64 IS the Absorption Theorem α = 2γ·⟨n_XY⟩ applied to single-site exposure (⟨n_XY⟩_B = |a_B|²); same '2' as F33's WeightOneRateCoefficient and F50's DecayRateFactor");
            yield return new InspectableNode("verified examples",
                summary: $"r=0.5: g = {N3GValue(0.5):F4} (zero mode); r=1: g = {N3GValue(1.0):F4} (bonding mode at uniform J); r=2: g = {N3GValue(2.0):F4}");
        }
    }
}
