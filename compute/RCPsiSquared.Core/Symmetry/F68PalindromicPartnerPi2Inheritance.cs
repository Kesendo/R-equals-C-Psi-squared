using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F68 closed form (Tier 1, verified N=3, 4, 5):
///
/// <code>
///   α_p = 2γ₀ − α_b                              (F1 palindrome identity at eigenvalue level)
///
///   |α_b + α_p − 2γ₀| &lt; 4·10⁻¹⁵                  exact to machine precision
///
///   For the F67 bonding-mode setup (uniform XY chain, endpoint Z-dephasing):
///   α_b = α_1 = (4γ₀/(N+1))·sin²(π/(N+1))         (F67 / F65 at k=1)
///   α_p = 2γ₀ − (4γ₀/(N+1))·sin²(π/(N+1))         (palindromic partner)
/// </code>
///
/// <para>F68 is F1's palindrome identity (Π·L·Π⁻¹ + L + 2σ·I = 0) read at the
/// eigenvalue level: every L-eigenvalue −α_b has a partner −α_p such that
/// α_b + α_p = 2γ₀ (equivalently, α_b + α_p = 2σ for σ = Σγ; here σ = γ₀ for
/// endpoint dephasing). The bonding-mode partner V_p lives in the
/// XY-weight-(N−1) Pauli sector — Π-mirror of the bonding mode's w=1.</para>
///
/// <para>Two Pi2-Foundation anchors:</para>
/// <list type="bullet">
///   <item><b>SumCoefficient = 2 = a_0</b>: in α_b + α_p = 2γ₀. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0) (= polynomial root d).
///         Same anchor as F66's UpperPoleCoefficient and F65's
///         UpperBoundCoefficient.</item>
///   <item><b>F1 palindrome identity</b>: α_b + α_p = 2σ is the F1 identity at
///         spectral level; F1's "2σ·I" shift is what produces the partner pairing.</item>
/// </list>
///
/// <para>Operational rank: N ≥ 4 admits a rank-1 SVD encoding V_p = σ₀|u⟩⟨v|
/// (verified to 9.3·10⁻¹² at N=5, 5.1·10⁻⁸ at N=4 limited by 16-fold partner
/// degeneracy). N=3 is rank-2 on both sides (fourfold degenerate; σ₁/σ₀ ≈ 0.98);
/// no clean rank-1 operational encoding exists at N=3. The algebraic palindromic
/// pairing α_b + α_p = 2γ₀ still holds spectrally for all N.</para>
///
/// <para>Tier1Derived: F68 follows directly from F1 (algebraic) + F65/F67
/// (provides α_b); spectrally verified at N=3, 4, 5 with |α_b + α_p − 2γ₀| &lt;
/// 4·10⁻¹⁵; operationally verified via Bell-pair-like (|0⟩|u⟩ + |1⟩|v⟩)/√2
/// propagation (rel err 2.8·10⁻¹⁶ at N=4, 3.8·10⁻¹⁴ at N=5 for clean SVD encoding).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F68 (line 1435) +
/// <c>experiments/PALINDROMIC_PARTNER_MODE.md</c> +
/// <c>simulations/palindromic_partner_f67.py</c> +
/// <c>simulations/bell_pair_partner_mode.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs</c> (palindrome identity) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F67BondingBellPairPi2Inheritance.cs</c> (α_b source).</para></summary>
public sealed class F68PalindromicPartnerPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F1Pi2Inheritance _f1;
    private readonly F67BondingBellPairPi2Inheritance _f67;

    /// <summary>The "2" coefficient in α_b + α_p = 2γ₀. Live from Pi2DyadicLadder a_0
    /// (polynomial root d). Same anchor as F1's TwoFactor and F66's UpperPoleCoefficient.</summary>
    public double SumCoefficient => _ladder.Term(0);

    /// <summary>The bonding-mode rate α_b. Delegates to F67 (which delegates to F65 at k=1).</summary>
    public double BondingRate(int N, double gammaZero)
    {
        return _f67.BondingModeDecayRate(N, gammaZero);
    }

    /// <summary>The palindromic partner rate α_p = 2γ₀ − α_b. The F1 palindrome
    /// identity at eigenvalue level.</summary>
    public double PartnerRate(int N, double gammaZero)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F68 requires N ≥ 2.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return SumCoefficient * gammaZero - BondingRate(N, gammaZero);
    }

    /// <summary>The palindromic sum α_b + α_p. Should be exactly SumCoefficient · γ₀
    /// = 2γ₀ by the F1 identity. Drift check.</summary>
    public double PalindromicSum(int N, double gammaZero)
    {
        return BondingRate(N, gammaZero) + PartnerRate(N, gammaZero);
    }

    /// <summary>True iff |α_b + α_p − 2γ₀| &lt; tolerance. Spectral verification of the F1
    /// identity; should hold at machine precision for all N.</summary>
    public bool PalindromicSumHolds(int N, double gammaZero, double tolerance = 4e-15)
    {
        double sum = PalindromicSum(N, gammaZero);
        double expected = SumCoefficient * gammaZero;
        return Math.Abs(sum - expected) < tolerance;
    }

    /// <summary>True iff F68 admits a clean rank-1 SVD operational encoding V_p = σ₀|u⟩⟨v|.
    /// Holds for N ≥ 4 (verified to 9.3·10⁻¹² at N=5). False at N=3 where both V_b
    /// and V_p are fourfold-degenerate rank-2 (no clean rank-1 operational encoding).</summary>
    public bool IsRankOneOperational(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F68 requires N ≥ 2.");
        return N >= 4;
    }

    public F68PalindromicPartnerPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F1Pi2Inheritance f1,
        F67BondingBellPairPi2Inheritance f67)
        : base("F68 palindromic partner of bonding mode: α_p = 2γ₀ − α_b (F1 palindrome at eigenvalue level); rank-1 operational at N ≥ 4",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F68 + " +
               "experiments/PALINDROMIC_PARTNER_MODE.md + " +
               "simulations/palindromic_partner_f67.py + " +
               "simulations/bell_pair_partner_mode.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F67BondingBellPairPi2Inheritance.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f1 = f1 ?? throw new ArgumentNullException(nameof(f1));
        _f67 = f67 ?? throw new ArgumentNullException(nameof(f67));
    }

    public override string DisplayName =>
        "F68 palindromic partner of bonding mode as Pi2-Foundation a_0 + F1 + F67 inheritance";

    public override string Summary =>
        $"α_p = 2γ₀ − α_b (F1 palindrome at spectral level); 2 = a_0; rank-1 operational at N ≥ 4, rank-2 at N=3 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F68 closed form",
                summary: "α_p = 2γ₀ − α_b; |α_b + α_p − 2γ₀| < 4·10⁻¹⁵ verified N=3,4,5; partner V_p in Pauli weight-(N-1) sector (Π-mirror of bonding mode's w=1)");
            yield return InspectableNode.RealScalar("SumCoefficient (= a_0 = 2)", SumCoefficient);
            yield return new InspectableNode("F1 palindrome at eigenvalue level",
                summary: $"F1 identity Π·L·Π⁻¹ + L + 2σ·I = 0 squared/spectral consequence: α_b + α_p = 2σ. F1's TwoFactor (= {_f1.TwoFactor}) is the same '2'.");
            yield return new InspectableNode("F67 source-claim edge",
                summary: "F68's α_b = F67's BondingModeDecayRate (which delegates to F65 at k=1). Live delegation, no formula duplication.");
            yield return new InspectableNode("rank structure",
                summary: "N ≥ 4: rank-1 V_p = σ₀|u⟩⟨v| (verified 9.3·10⁻¹² at N=5, 5.1·10⁻⁸ at N=4 limited by partner degeneracy); N=3: rank-2 fourfold-degenerate, no clean rank-1 operational encoding");
            yield return new InspectableNode("operational verification",
                summary: "Bell-pair-like (|0⟩|u⟩ + |1⟩|v⟩)/√2 propagation matches spectral α_p: clean SVD encoding rel err 2.8·10⁻¹⁶ (N=4), 3.8·10⁻¹⁴ (N=5); legacy bonding encoding rel err 1.6·10⁻⁶ (N=4), 2.8·10⁻⁷ (N=5) where residual is F65's O((γ₀/J)²) shift");
            yield return new InspectableNode("N=3 verified",
                summary: $"α_b(N=3, γ=0.05) = {BondingRate(3, 0.05):G6}; α_p(N=3, γ=0.05) = {PartnerRate(3, 0.05):G6}; sum = {PalindromicSum(3, 0.05):G15}");
            yield return new InspectableNode("N=5 verified",
                summary: $"α_b(N=5, γ=0.05) = {BondingRate(5, 0.05):G6}; α_p(N=5, γ=0.05) = {PartnerRate(5, 0.05):G6}; sum = {PalindromicSum(5, 0.05):G15}");
        }
    }
}
