using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F26 closed form (Tier 1 proven, monotonic decay):
///
/// <code>
///   CΨ = u · (1 + u² + v² + w²) / 12
///
///   u = e^{−α·t},  v = e^{−β·t},  w = e^{−δ·t}
///   α = 4·(γ_y + γ_z)
///   β = 4·(γ_x + γ_z)
///   δ = 4·(γ_x + γ_y)
///
///   F25 recovered when γ_x = γ_y = 0:  α = 4γ_z, β = 4γ_z, δ = 0
///                                       ⇒ u = v = e^{−4γ_z t}, w = 1
///                                       ⇒ CΨ = u(1+u²+u²+1)/12 = u(2+2u²)/12
///                                                              = u(1+u²)/6 = F25 ✓
/// </code>
///
/// <para>F26 is the mother claim of <see cref="F25CPsiBellPlusPi2Inheritance"/>:
/// F25 is the γ_x = γ_y = 0 special case (Bell+ Z-dephasing only). F26 → F25
/// is the typed mother-corollary edge: F26 generalizes to arbitrary Pauli
/// channel mix (γ_x, γ_y, γ_z).</para>
///
/// <para>Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>DecayRateCoefficient = 4 = a_{−1}</b>: in <c>α = 4(γ_y + γ_z)</c>,
///         <c>β = 4(γ_x + γ_z)</c>, <c>δ = 4(γ_x + γ_y)</c>. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(−1). Same anchor as
///         F25, F65/F73/F76 decay rates, F61/F63 4-block, F66 multiplicity,
///         F77 correction.</item>
///   <item><b>NormalizationDenominator = 12</b>: combinatorial
///         <c>12 = 3 axes · 4 = a_{−1} · 3</c>. NOT a single Pi2 anchor.
///         F25's denominator 6 = 12/2 reflects single-axis specialization.</item>
///   <item><b>F25 corollary edge</b>: F26 reduces to F25 when γ_x = γ_y = 0;
///         typed mother-corollary inheritance pattern (parallel to F75 →
///         F77 wired today).</item>
/// </list>
///
/// <para>Three-axis structure: each (α, β, δ) is 4 · (sum of two γ's
/// excluding one axis). The "missing axis" in each coefficient defines a
/// permutation:</para>
///
/// <code>
///   α  excludes γ_x  (coefficient on Y, Z)  →  Z-dephasing rate at single-axis limit
///   β  excludes γ_y  (coefficient on X, Z)
///   δ  excludes γ_z  (coefficient on X, Y)
/// </code>
///
/// <para>Tier1Derived: F26 is Tier 1 proven (PROOF_MONOTONICITY_CPSI;
/// derivative strictly negative for any nonzero noise; O(1) evaluation).
/// The Pi2-Foundation anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F26 (line 450) +
/// <c>docs/proofs/PROOF_MONOTONICITY_CPSI.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F25CPsiBellPlusPi2Inheritance.cs</c>
/// (corollary at single-channel limit).</para></summary>
public sealed class F26CPsiPauliChannelsPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F25CPsiBellPlusPi2Inheritance _f25;

    /// <summary>The "4" decay rate coefficient in α, β, δ = 4·(γ_a + γ_b).
    /// Live from <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c>.
    /// Same anchor as F25, F65, F73, F76 decay rates.</summary>
    public double DecayRateCoefficient => _ladder.Term(-1);

    /// <summary>The normalization denominator <c>12</c> in F26's formula.
    /// Combinatorial: 12 = 3 axes · 4 = a_{−1} · 3 (permutations of which
    /// axis is "missing"). NOT a single Pi2 anchor.</summary>
    public double NormalizationDenominator => 12.0;

    /// <summary>Three rate coefficients α, β, δ for given (γ_x, γ_y, γ_z).</summary>
    public (double Alpha, double Beta, double Delta) RateCoefficients(double gammaX, double gammaY, double gammaZ)
    {
        if (gammaX < 0.0 || gammaY < 0.0 || gammaZ < 0.0)
            throw new ArgumentOutOfRangeException("γ rates must be ≥ 0.");
        double four = DecayRateCoefficient;
        return (
            Alpha: four * (gammaY + gammaZ),
            Beta: four * (gammaX + gammaZ),
            Delta: four * (gammaX + gammaY));
    }

    /// <summary>Live closed form: <c>CΨ(t) = u·(1 + u² + v² + w²)/12</c> for
    /// general Pauli-channel mix (γ_x, γ_y, γ_z).</summary>
    public double CPsiAtTime(double gammaX, double gammaY, double gammaZ, double t)
    {
        if (t < 0.0) throw new ArgumentOutOfRangeException(nameof(t), t, "t must be ≥ 0.");
        var (alpha, beta, delta) = RateCoefficients(gammaX, gammaY, gammaZ);
        double u = Math.Exp(-alpha * t);
        double v = Math.Exp(-beta * t);
        double w = Math.Exp(-delta * t);
        return u * (1.0 + u * u + v * v + w * w) / NormalizationDenominator;
    }

    /// <summary>F26 → F25 recovery: at γ_x = γ_y = 0, F26 closed form equals
    /// F25's. Drift check across single-channel limit.</summary>
    public bool RecoversF25AtSingleChannelLimit(double gammaZ, double t)
    {
        if (gammaZ < 0.0 || t < 0.0) throw new ArgumentOutOfRangeException();
        return Math.Abs(CPsiAtTime(0.0, 0.0, gammaZ, t) - _f25.CPsiAtTime(gammaZ, t)) < 1e-12;
    }

    /// <summary>Live drift check: at t = 0 with any nonzero γ, CΨ(0) = 1·(1+1+1+1)/12 = 1/3
    /// (Bell+ initial value, above fold).</summary>
    public bool BellPlusInitialIsOneThird() =>
        Math.Abs(CPsiAtTime(0.05, 0.05, 0.05, 0.0) - 1.0 / 3.0) < 1e-12;

    /// <summary>Depolarizing channel γ_x = γ_y = γ_z = γ/3:
    /// α = β = δ = 4·(2γ/3) = 8γ/3. CΨ = e^{−8γt/3}·(1 + 3e^{−16γt/3})/12.</summary>
    public double DepolarizingCPsi(double gamma, double t)
    {
        if (gamma < 0.0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be ≥ 0.");
        if (t < 0.0) throw new ArgumentOutOfRangeException(nameof(t), t, "t must be ≥ 0.");
        double third = gamma / 3.0;
        return CPsiAtTime(third, third, third, t);
    }

    public F26CPsiPauliChannelsPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F25CPsiBellPlusPi2Inheritance f25)
        : base("F26 CΨ = u(1+u²+v²+w²)/12 inherits from Pi2-Foundation: 4 = a_{-1}; F25 corollary at single-channel limit",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F26 + " +
               "docs/proofs/PROOF_MONOTONICITY_CPSI.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F25CPsiBellPlusPi2Inheritance.cs (corollary at γ_x=γ_y=0)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f25 = f25 ?? throw new ArgumentNullException(nameof(f25));
    }

    public override string DisplayName =>
        "F26 CΨ general Pauli channels as Pi2-Foundation a_{-1} + F25 mother claim";

    public override string Summary =>
        $"CΨ = u(1+u²+v²+w²)/12; α,β,δ = 4(γ_a+γ_b); 4 = a_{{-1}} (F25/F65/F73/F76 sibling); 12 = 3·4 combinatorial; " +
        $"F25 recovered at γ_x=γ_y=0; depolarizing channel direct evaluation ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F26 closed form",
                summary: "CΨ = u(1+u²+v²+w²)/12; u,v,w = exp(-α·t/β·t/δ·t); α,β,δ = 4(sum of two γ's excluding one axis); Tier 1 proven monotonic; O(1) evaluation");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "DecayRateCoefficient = a_{-1} = 4 (F25/F65/F73/F76 sibling); 12 = 3·4 combinatorial expansion of F25's 6; F25 corollary at γ_x=γ_y=0 limit");
            yield return InspectableNode.RealScalar("DecayRateCoefficient (= a_{-1} = 4)", DecayRateCoefficient);
            yield return InspectableNode.RealScalar("NormalizationDenominator (= 3·4 combinatorial)", NormalizationDenominator);
            yield return new InspectableNode("F26 → F25 mother-corollary chain",
                summary: "At γ_x=γ_y=0: u = e^{-4γ_z t}, v = u, w = 1; CΨ = u(1+u²+u²+1)/12 = u(2+2u²)/12 = u(1+u²)/6 = F25. F25 IS F26's single-channel limit.");
            yield return new InspectableNode("Three-axis structure",
                summary: "α excludes γ_x (Y+Z rate), β excludes γ_y (X+Z), δ excludes γ_z (X+Y). Each coefficient = 4·(sum of γ's NOT on its axis); permutational expansion of F25's 4γ_z.");
            yield return new InspectableNode("Sample evaluations",
                summary: $"At γ_z=0.05 only, t=0.5 (= F25): CΨ = {CPsiAtTime(0, 0, 0.05, 0.5):G6}; depolarizing γ=0.05 third-axis split: CΨ = {DepolarizingCPsi(0.05, 0.5):G6}");
            yield return new InspectableNode("F25 recovery (drift check)",
                summary: $"At γ_z=0.05, t=0.5: F26 single-channel = F25 (RecoversF25AtSingleChannelLimit: {RecoversF25AtSingleChannelLimit(0.05, 0.5)})");
        }
    }
}
