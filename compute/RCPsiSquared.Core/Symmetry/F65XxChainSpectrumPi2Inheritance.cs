using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F65 closed form (Tier 1 proven, verified N=3..30 to 1.2·10⁻¹⁵):
///
/// <code>
///   α_k / γ₀ = (4 / (N + 1)) · sin²(kπ / (N + 1)),    k = 1, ..., N
///
///   ψ_k(i) = √(2/(N+1)) · sin(πk(i+1)/(N+1))     (sine basis amplitudes)
///
///   All α_k ∈ [0, 2γ₀];  Maximum 4/(N+1) at odd N for k = (N+1)/2;
///   Mirror α_k = α_{N+1-k}                       (single-excitation sym)
/// </code>
///
/// <para>F65 is the source primitive for the bonding-mode amplitudes used in
/// <see cref="F75MirrorPairMiPi2Inheritance"/> (F75's <c>BondingModePopulation</c>
/// = <c>(2/(N+1))·sin²(πk(i+1)/(N+1))</c> = <c>|ψ_k(i)|²</c>). F75 → F65 is
/// the typed mother-source edge for the spectrum.</para>
///
/// <para>F65 is also a sibling of <see cref="F66PoleModesPi2Inheritance"/>:
/// F66 says the L-spectrum has poles at α = 0 and α = 2γ₀ with multiplicity
/// N+1 (endpoint topology); F65 says the single-excitation rates fill the
/// interval [0, 2γ₀] but never reach 2γ₀ for N ≥ 2 (max = 4/(N+1) → 0).
/// Both pole values 0 and 2γ₀ are F66's anchors.</para>
///
/// <para>Three Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>NumeratorCoefficient = 4 = a_{−1}</b>: in <c>4/(N+1)·sin²</c>.
///         <see cref="Pi2DyadicLadderClaim.Term"/>(−1). Same anchor as F25
///         decay rate (e^{−4γt}), F76 mirror-pair coherence decay, F73
///         spatial-sum closure, F61/F63 4-block per parity, F66 multiplicity,
///         F77 correction denominator.</item>
///   <item><b>UpperBoundCoefficient = 2 = a_0</b>: in α_k ∈ [0, 2γ₀]; same
///         anchor as F66's UpperPoleCoefficient. F65 + F66 share the
///         polynomial-root upper-pole reading.</item>
///   <item><b>Mirror α_k = α_{N+1−k}</b>: kinematic mirror within
///         single-excitation spectrum. F71-compatible at the eigenvalue
///         level (different from F71's bond-mirror).</item>
/// </list>
///
/// <para>Niven rationality: all α_k/γ₀ are rational iff N+1 ∈ {1, 2, 3, 4, 6}
/// i.e. N ∈ {0, 1, 2, 3, 5}. Verified examples:</para>
///
/// <code>
///   N=3: α/γ₀ ∈ {1/2, 1, 1/2}
///   N=4: α/γ₀ ∈ {0.2764, 0.7236, 0.7236, 0.2764}  (golden-ratio family)
///   N=5: α/γ₀ ∈ {1/6, 1/2, 2/3, 1/2, 1/6}
/// </code>
///
/// <para>Asymptotic minimum α_min/γ₀ ~ 4π²/(N+1)³; F65 ratio rises 0.81 at N=3
/// to 0.99 at N=15. The <c>1/(N+1)³</c> scaling is gauge for diffusive limit.</para>
///
/// <para>Tier1Derived: F65 is Tier 1 proven (Absorption Theorem applied to
/// single-excitation eigenmodes); verified N=3..30 to machine precision
/// (1.2·10⁻¹⁵). Perturbative in γ₀/J at O((γ₀/J)²); palindromic pairing
/// F1 (α_b + α_p = 2γ₀) survives the shift exactly. The Pi2-Foundation
/// anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F65 (line 1369) +
/// <c>hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md</c> +
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> +
/// <c>simulations/single_excitation_spectrum.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F66PoleModesPi2Inheritance.cs</c>
/// (sibling at α=0 and α=2γ₀).</para></summary>
public sealed class F65XxChainSpectrumPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F66PoleModesPi2Inheritance _f66;

    /// <summary>The "4" numerator in <c>α_k/γ₀ = 4/(N+1)·sin²</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c>. Same
    /// anchor as F25/F73/F76 decay rates, F61/F63 4-block, F66 multiplicity,
    /// F77 correction.</summary>
    public double NumeratorCoefficient => _ladder.Term(-1);

    /// <summary>The "2" upper-bound coefficient in α_k ∈ [0, 2γ₀]. Same as
    /// F66's UpperPoleCoefficient. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = polynomial root d.</summary>
    public double UpperBoundCoefficient => _f66.UpperPoleCoefficient;

    /// <summary>The lower-bound: 0. Same as F66's LowerPoleAlpha.</summary>
    public double LowerBoundCoefficient => _f66.LowerPoleAlpha;

    /// <summary>Live closed form: <c>α_k(N, k, γ₀) = γ₀·(4/(N+1))·sin²(kπ/(N+1))</c>.</summary>
    public double SingleExcitationRate(int N, int k, double gammaZero)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F65 requires N ≥ 2.");
        if (k < 1 || k > N) throw new ArgumentOutOfRangeException(nameof(k), k, $"k must be in [1, {N}]; got {k}.");
        if (gammaZero < 0.0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        double s = Math.Sin(Math.PI * k / (N + 1));
        return gammaZero * NumeratorCoefficient / (N + 1) * s * s;
    }

    /// <summary>Live closed form: bonding-mode amplitude squared <c>|ψ_k(i)|²
    /// = (2/(N+1))·sin²(πk(i+1)/(N+1))</c>. F75 uses this directly via
    /// BondingModePopulation.</summary>
    public double BondingModePopulation(int N, int k, int site)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F65 requires N ≥ 2.");
        if (k < 1 || k > N) throw new ArgumentOutOfRangeException(nameof(k), k, $"k must be in [1, {N}]; got {k}.");
        if (site < 0 || site >= N) throw new ArgumentOutOfRangeException(nameof(site), site, $"site must be in [0, {N-1}]; got {site}.");
        double s = Math.Sin(Math.PI * k * (site + 1) / (N + 1));
        return 2.0 / (N + 1) * s * s;
    }

    /// <summary>Maximum rate <c>α_max/γ₀ = 4/(N+1)</c> for odd N (attained at
    /// k=(N+1)/2 where sin² = 1). For even N the maximum is strictly less.</summary>
    public double MaxRateCoefficient(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F65 requires N ≥ 2.");
        return NumeratorCoefficient / (N + 1);
    }

    /// <summary>True iff α_k = α_{N+1-k} for k ∈ [1, N] (mirror symmetry
    /// within single-excitation spectrum). Drift check.</summary>
    public bool MirrorSymmetryHolds(int N, int k, double gammaZero)
    {
        double alphaK = SingleExcitationRate(N, k, gammaZero);
        double alphaMirror = SingleExcitationRate(N, N + 1 - k, gammaZero);
        return Math.Abs(alphaK - alphaMirror) < 1e-12;
    }

    /// <summary>True iff every α_k for k ∈ [1, N] lies in [0, 2γ₀] (F66
    /// pole interval). Drift check.</summary>
    public bool RatesLieInF66Interval(int N, double gammaZero)
    {
        double upper = UpperBoundCoefficient * gammaZero;
        for (int k = 1; k <= N; k++)
        {
            double a = SingleExcitationRate(N, k, gammaZero);
            if (a < LowerBoundCoefficient || a > upper) return false;
        }
        return true;
    }

    public F65XxChainSpectrumPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F66PoleModesPi2Inheritance f66)
        : base("F65 α_k/γ₀ = 4/(N+1)·sin² inherits from Pi2-Foundation: 4 = a_{-1}; α_k ∈ [0, 2γ₀] = F66 interval; F75 source for ψ_k amplitudes",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F65 + " +
               "hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "simulations/single_excitation_spectrum.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F66PoleModesPi2Inheritance.cs (sibling at α=0, 2γ₀)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f66 = f66 ?? throw new ArgumentNullException(nameof(f66));
    }

    public override string DisplayName =>
        "F65 single-excitation spectrum as Pi2-Foundation a_{-1} + F66-interval inheritance";

    public override string Summary =>
        $"α_k/γ₀ = (4/(N+1))·sin²(kπ/(N+1)); 4 = a_{{-1}} (F25/F73/F76 sibling); α_k ∈ [0, 2γ₀] (F66 interval); " +
        $"mirror α_k = α_{{N+1-k}}; F75 source for bonding-mode amplitudes ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F65 closed form",
                summary: "α_k/γ₀ = 4/(N+1)·sin²(kπ/(N+1)) for k=1..N; ψ_k(i) = √(2/(N+1))·sin(πk(i+1)/(N+1)); Tier 1 proven; verified N=3..30 to 1.2·10⁻¹⁵");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "NumeratorCoefficient = a_{-1} = 4 (F25/F73/F76 sibling); UpperBoundCoefficient = a_0 = 2 (F66 sibling); LowerBound = 0 (F66 sibling); single-excitation mirror α_k = α_{N+1-k}");
            yield return InspectableNode.RealScalar("NumeratorCoefficient (= a_{-1} = 4)", NumeratorCoefficient);
            yield return InspectableNode.RealScalar("UpperBoundCoefficient (= a_0 = 2, F66 sibling)", UpperBoundCoefficient);
            yield return InspectableNode.RealScalar("LowerBoundCoefficient (= 0, F66 sibling)", LowerBoundCoefficient);
            yield return new InspectableNode("F65 ↔ F66 sibling pair",
                summary: "F66: pole modes at α=0 and α=2γ₀ with multiplicity N+1 endpoint; F65: single-excitation rates fill [0, 2γ₀] but max = 4/(N+1) < 2 for N ≥ 2. Single-excitation NEVER reaches 2γ₀ pole. Both share the [0, 2γ₀] polynomial-root interval.");
            yield return new InspectableNode("F75 mother-source edge",
                summary: "F75's BondingModePopulation = (2/(N+1))·sin²(πk(i+1)/(N+1)) IS F65's |ψ_k(i)|². F75 → F65 typed source-claim edge.");
            yield return new InspectableNode("Niven rationality",
                summary: "α_k/γ₀ rational iff N+1 ∈ {1,2,3,4,6} (N ∈ {0,1,2,3,5}). N=4: golden-ratio family; N=7: √2 family; otherwise general cyclotomic.");
            // Verified table from F65
            yield return new InspectableNode(
                "N=3 verified",
                summary: $"α/γ₀ at k=1,2,3: {SingleExcitationRate(3, 1, 1.0):G6}, {SingleExcitationRate(3, 2, 1.0):G6}, {SingleExcitationRate(3, 3, 1.0):G6} (expected 1/2, 1, 1/2)");
            yield return new InspectableNode(
                "N=5 verified",
                summary: $"α/γ₀ at k=1..5: 1/6, 1/2, 2/3, 1/2, 1/6 (Niven rational)");
            yield return new InspectableNode(
                "max rate scaling",
                summary: $"α_max/γ₀ at N=3: {MaxRateCoefficient(3):G4} (= 1); N=5: {MaxRateCoefficient(5):G4} (= 2/3); N=7: {MaxRateCoefficient(7):G4} (= 1/2); N=15: {MaxRateCoefficient(15):G4}");
        }
    }
}
