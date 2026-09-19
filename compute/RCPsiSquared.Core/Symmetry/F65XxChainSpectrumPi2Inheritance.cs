using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F65 first-order endpoint-rate coefficient comb (Tier 1 proven,
/// verified N=3..30 to 1.2·10⁻¹⁵). For the uniform open XX chain with one
/// dephased endpoint:
/// The model is the uniform open XX chain with one dephased endpoint.
///
/// <code>
///   a_k = (4/(N+1))·sin²(kπ/(N+1)),                    k = 1, ..., N
///   α_k^full = γ₀·a_k + O(γ₀³/J²)
///   α_k^full/γ₀ = a_k + O((γ₀/J)²)
///
///   ψ_k(i) = √(2/(N+1)) · sin(πk(i+1)/(N+1))     (sine basis amplitudes)
///
///   All a_k ∈ [0, 2];  Maximum 4/(N+1) at odd N for k = (N+1)/2;
///   Mirror a_k = a_{N+1-k}                       (single-excitation sym)
/// </code>
///
/// <para>F65 is the source primitive for the bonding-mode amplitudes used in
/// <see cref="F75MirrorPairMiPi2Inheritance"/> (F75's <c>BondingModePopulation</c>
/// = <c>(2/(N+1))·sin²(πk(i+1)/(N+1))</c> = <c>|ψ_k(i)|²</c>). F75 → F65 is
/// the typed mother-source edge for the spectrum.</para>
///
/// <para>F65 is also a sibling of <see cref="F66PoleModesPi2Inheritance"/>:
/// F66 says the L-spectrum has poles at α = 0 and α = 2γ₀ with multiplicity
/// N+1 (endpoint topology); F65 says the first-order single-excitation
/// coefficients lie in [0, 2] but never reach 2 for N ≥ 2
/// (a_max ≤ 4/(N+1) → 0; equality holds iff N is odd, and a_max is strictly
/// smaller for even N).
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
///   <item><b>UpperBoundCoefficient = 2 = a_0</b>: in a_k ∈ [0, 2]; same
///         anchor as F66's UpperPoleCoefficient. F65 + F66 share the
///         polynomial-root upper-pole reading.</item>
///   <item><b>Mirror a_k = a_{N+1−k}</b>: kinematic mirror within
///         single-excitation spectrum. F71-compatible at the eigenvalue
///         level (different from F71's bond-mirror).</item>
/// </list>
///
/// <para>Niven rationality: all first-order coefficients a_k are rational iff
/// N+1 ∈ {1, 2, 3, 4, 6}. In the physical domain N ≥ 2 this is
/// N ∈ {2, 3, 5}; the formal extension is N ∈ {0, 1, 2, 3, 5}. Exact Niven rationality
/// belongs to the first-order coefficient comb, not to finite-γ₀/J full-L
/// rates. Verified examples:</para>
///
/// <code>
///   N=3: a ∈ {1/2, 1, 1/2}
///   N=4: a ∈ {0.2764, 0.7236, 0.7236, 0.2764}  (golden-ratio family)
///   N=5: a ∈ {1/6, 1/2, 2/3, 1/2, 1/6}
/// </code>
///
/// <para>Asymptotic minimum a_min ~ 4π²/(N+1)³; F65 ratio rises 0.81 at N=3
/// to 0.99 at N=15. The <c>1/(N+1)³</c> scaling is gauge for diffusive limit.</para>
///
/// <para>Tier1Derived: F65 is Tier 1 proven (Absorption Theorem applied to
/// single-excitation eigenmodes); verified N=3..30 to machine precision
/// (1.2·10⁻¹⁵). At finite γ₀/J the relative full-L rate shift is O((γ₀/J)²),
/// equivalently the absolute shift δα_k = O(γ₀³/J²), so no exact finite-γ₀/J
/// full-L rationality is claimed. Palindromic pairing F1
/// (α_b + α_p = 2γ₀) survives the shift exactly. The Pi2-Foundation
/// anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F65 +
/// <c>hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md</c> +
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> +
/// <c>simulations/single_excitation_spectrum.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F66PoleModesPi2Inheritance.cs</c>
/// (sibling at α=0 and α=2γ₀).</para></summary>
public sealed class F65XxChainSpectrumPi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// currently not typed for this Claim.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed bit_a-twin sibling, if one exists. Currently null
    /// (no bit_a twin is typed for this Claim; this is an open slot in the
    /// cubic-architecture coverage).</summary>
    public Claim? BitATwin => null;

    // Absorption-Theorem descendant; bit_a twin holds by the Hadamard X↔Z duality (PROOF_BIT_A_TWIN_VIA_HADAMARD.md).
    public BitATwinClassification BitATwinStatus => BitATwinClassification.CoveredByHadamardDuality;
    public Pi2DyadicLadderClaim Ladder { get; }
    public F66PoleModesPi2Inheritance F66 { get; }
    /// <summary>The "4" numerator in <c>a_k = 4/(N+1)·sin²</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c>. Same
    /// anchor as F25/F73/F76 decay rates, F61/F63 4-block, F66 multiplicity,
    /// F77 correction.</summary>
    public double NumeratorCoefficient => Ladder.Term(-1);

    /// <summary>The "2" upper bound for the first-order coefficient a_k. Same as
    /// F66's UpperPoleCoefficient. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = polynomial root d.</summary>
    public double UpperBoundCoefficient => F66.UpperPoleCoefficient;

    /// <summary>The lower-bound: 0. Same as F66's LowerPoleAlpha.</summary>
    public double LowerBoundCoefficient => F66.LowerPoleAlpha;

    /// <summary>Exact first-order endpoint coefficient
    /// <c>a_k = (4/(N+1))·sin²(kπ/(N+1))</c>.</summary>
    public double FirstOrderRateCoefficient(int N, int k)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F65 requires N ≥ 2.");
        if (k < 1 || k > N) throw new ArgumentOutOfRangeException(nameof(k), k, $"k must be in [1, {N}]; got {k}.");
        double s = Math.Sin(Math.PI * k / (N + 1));
        return NumeratorCoefficient / (N + 1) * s * s;
    }

    /// <summary>Returns the first-order term γ₀·a_k. It is not the finite-γ₀/J
    /// full-L rate <c>α_k^full = γ₀·a_k + O(γ₀³/J²)</c>.</summary>
    public double FirstOrderRateTerm(int N, int k, double gammaZero)
    {
        if (gammaZero < 0.0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return gammaZero * FirstOrderRateCoefficient(N, k);
    }

    /// <summary>Compatibility wrapper. Returns only the first-order term
    /// γ₀·a_k, not a finite-γ₀/J full-L rate. New code should use
    /// <see cref="FirstOrderRateTerm"/> or <see cref="FirstOrderRateCoefficient"/>.</summary>
    public double SingleExcitationRate(int N, int k, double gammaZero) =>
        FirstOrderRateTerm(N, k, gammaZero);

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

    /// <summary>Maximum first-order coefficient <c>a_max = 4/(N+1)</c> for odd N (attained at
    /// k=(N+1)/2 where sin² = 1). For even N the maximum is strictly less.</summary>
    public double MaxFirstOrderCoefficient(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F65 requires N ≥ 2.");
        double oddNMaximum = NumeratorCoefficient / (N + 1);
        if ((N & 1) == 1) return oddNMaximum;

        double centralOffset = Math.Cos(Math.PI / (2.0 * (N + 1)));
        return oddNMaximum * centralOffset * centralOffset;
    }

    /// <summary>Compatibility wrapper for <see cref="MaxFirstOrderCoefficient"/>.</summary>
    public double MaxRateCoefficient(int N) => MaxFirstOrderCoefficient(N);

    /// <summary>True iff the returned first-order terms obey
    /// γ₀a_k = γ₀a_{N+1-k} for k ∈ [1, N]. Drift check.</summary>
    public bool MirrorSymmetryHolds(int N, int k, double gammaZero)
    {
        double alphaK = FirstOrderRateTerm(N, k, gammaZero);
        double alphaMirror = FirstOrderRateTerm(N, N + 1 - k, gammaZero);
        return Math.Abs(alphaK - alphaMirror) < 1e-12;
    }

    /// <summary>True iff every returned first-order term γ₀a_k lies in
    /// [0, 2γ₀], the interval delimited by F66's exact poles. Drift check.</summary>
    public bool RatesLieInF66Interval(int N, double gammaZero)
    {
        double upper = UpperBoundCoefficient * gammaZero;
        for (int k = 1; k <= N; k++)
        {
            double a = FirstOrderRateTerm(N, k, gammaZero);
            if (a < LowerBoundCoefficient || a > upper) return false;
        }
        return true;
    }

    public F65XxChainSpectrumPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F66PoleModesPi2Inheritance f66)
        : base("F65 first-order coefficient comb a_k = (4/(N+1))·sin²(kπ/(N+1)) inherits from Pi2-Foundation: 4 = a_{-1}; γ₀a_k lies in the F66 interval [0, 2γ₀]; F75 source for ψ_k amplitudes",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F65 + " +
               "hypotheses/PRIMORDIAL_GAMMA_CONSTANT.md + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "simulations/single_excitation_spectrum.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F66PoleModesPi2Inheritance.cs (sibling at α=0, 2γ₀)")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        F66 = f66 ?? throw new ArgumentNullException(nameof(f66));
    }

    public override string DisplayName =>
        "F65 first-order endpoint-rate coefficient comb as Pi2-Foundation a_{-1} + F66-interval inheritance";

    public override string Summary =>
        $"a_k = (4/(N+1))·sin²(kπ/(N+1)); α_k^full = γ₀·a_k + O(γ₀³/J²), equivalently α_k^full/γ₀ = a_k + O((γ₀/J)²); " +
        $"exact Niven rationality belongs to the first-order coefficient comb; mirror a_k = a_{{N+1-k}}; F75 source ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F65 first-order coefficient comb",
                summary: "a_k = (4/(N+1))·sin²(kπ/(N+1)) for k=1..N; α_k^full = γ₀·a_k + O(γ₀³/J²), equivalently α_k^full/γ₀ = a_k + O((γ₀/J)²); ψ_k(i) = √(2/(N+1))·sin(πk(i+1)/(N+1)); Tier 1 proven; verified N=3..30 to 1.2·10⁻¹⁵");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "NumeratorCoefficient = a_{-1} = 4 (F25/F73/F76 sibling); UpperBoundCoefficient = a_0 = 2 (F66 sibling); LowerBound = 0 (F66 sibling); first-order mirror a_k = a_{N+1-k}");
            yield return InspectableNode.RealScalar("NumeratorCoefficient (= a_{-1} = 4)", NumeratorCoefficient);
            yield return InspectableNode.RealScalar("UpperBoundCoefficient (= a_0 = 2, F66 sibling)", UpperBoundCoefficient);
            yield return InspectableNode.RealScalar("LowerBoundCoefficient (= 0, F66 sibling)", LowerBoundCoefficient);
            yield return new InspectableNode("F65 ↔ F66 sibling pair",
                summary: "F66: exact pole modes at α=0 and α=2γ₀ with multiplicity N+1 endpoint; F65: first-order terms γ₀a_k lie in [0, 2γ₀] but a_max ≤ 4/(N+1) < 2 for N ≥ 2; equality holds iff N is odd, and a_max is strictly smaller for even N. The first-order comb never reaches the upper pole coefficient.");
            yield return new InspectableNode("F75 mother-source edge",
                summary: "F75's BondingModePopulation = (2/(N+1))·sin²(πk(i+1)/(N+1)) IS F65's |ψ_k(i)|². F75 → F65 typed source-claim edge.");
            yield return new InspectableNode("Niven rationality",
                summary: "Exact Niven rationality belongs to the first-order coefficient comb: a_k rational for every k iff N+1 ∈ {1,2,3,4,6} (physical N ∈ {2,3,5}; formal extension {0,1,2,3,5}); no exact finite-γ₀/J full-L rationality is claimed. N=4: golden-ratio family; N=7: √2 family; otherwise general cyclotomic.");
            // Verified table from F65
            yield return new InspectableNode(
                "N=3 verified",
                summary: $"a_k at k=1,2,3: {FirstOrderRateCoefficient(3, 1):G6}, {FirstOrderRateCoefficient(3, 2):G6}, {FirstOrderRateCoefficient(3, 3):G6} (expected 1/2, 1, 1/2)");
            yield return new InspectableNode(
                "N=5 verified",
                summary: $"a_k at k=1..5: 1/6, 1/2, 2/3, 1/2, 1/6 (Niven rational)");
            yield return new InspectableNode(
                "max rate scaling",
                summary: $"a_max at N=3: {MaxFirstOrderCoefficient(3):G4} (= 1); N=5: {MaxFirstOrderCoefficient(5):G4} (= 2/3); N=7: {MaxFirstOrderCoefficient(7):G4} (= 1/2); N=15: {MaxFirstOrderCoefficient(15):G4}");
        }
    }
}
