using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F41 closed form (Tier 1, corollary of D10):
///
/// <code>
///   t_Pi = 2π / ω_min = π / (2·J · sin²(π/(2N)))
///
///   ω_min = 4·J · sin²(π/(2N))     slowest palindromic SFF modulation frequency
/// </code>
///
/// <para>F41 is the period of the slowest palindromic modulation in the spectral
/// form factor (SFF) of a Heisenberg chain under Z-dephasing. The palindromic
/// modulation is F1's structural consequence: the Liouvillian eigenvalues come
/// in pairs λ ↔ −λ − 2σ (F1's Π·L·Π⁻¹ = −L − 2σ·I), so the SFF has a periodic
/// modulation whose period equals the slowest pair difference.</para>
///
/// <para><b>Asymptotic scaling (large N):</b></para>
/// <code>
///   sin(π/(2N)) → π/(2N)  for N → ∞
///   sin²(π/(2N)) → π²/(4N²)
///   t_Pi → π / (2·J · π²/(4N²)) = 2N² / (π·J)
/// </code>
/// <para>So t_Pi grows as 2N²/(π·J), or equivalently ~N²/π² in J-units. The
/// palindromic modulation is a short-time effect at any finite N (cf. F42's
/// timescale-separation: t_Pi/t_H ~ N²/4^N → 0 for N → ∞).</para>
///
/// <para>Pi2-Foundation anchor:</para>
/// <list type="bullet">
///   <item><b>HoppingCoefficient = 2 = a_0</b>: in 2·J denominator. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0). Same anchor as F1's
///         TwoFactor (the "2σ" shift), F66's UpperPoleCoefficient, F50's
///         DegeneracyFactor and DecayRateFactor.</item>
/// </list>
///
/// <para>The "2" appears twice in F41: once as 2J in the denominator (commutator-
/// hopping factor), and structurally as 2π in the period definition (full
/// oscillation = 2π). Both = a_0 on the dyadic ladder. The π/(2N) argument
/// scaling and the 4·J in ω_min are derivable from these.</para>
///
/// <para>F1 connection: F41's existence (a finite palindromic period) follows from
/// F1's palindrome identity at the spectral level. Without F1, eigenvalue pairs
/// wouldn't pair as λ ↔ −λ − 2σ, and there would be no palindromic modulation
/// in the SFF. F41 reads F1 in the time domain.</para>
///
/// <para>Tier1Derived: F41 is Tier 1 corollary of D10's w=1 dispersion derivation;
/// confirmed by FFT peak matching &lt;1% for N=2..4, 6 in
/// <c>experiments/SPECTRAL_FORM_FACTOR.md</c>. Pi2-Foundation anchoring is
/// algebraic-trivial composition through the 2·J coefficient.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F41 (line 838) +
/// <c>experiments/SPECTRAL_FORM_FACTOR.md</c> +
/// <c>docs/proofs/derivations/D10_W1_DISPERSION.md</c> (parent derivation) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs</c> (palindrome
/// identity at spectral level) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F41PalindromicTimePi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F1Pi2Inheritance _f1;

    /// <summary>The "2" hopping coefficient in 2·J denominator. Live from
    /// Pi2DyadicLadder a_0. Same anchor as F1's TwoFactor.</summary>
    public double HoppingCoefficient => _ladder.Term(0);

    /// <summary>Slowest palindromic SFF modulation frequency:
    /// <c>ω_min = 4·J · sin²(π/(2N))</c>. Equals 2·HoppingCoefficient·J·sin²(π/(2N)).</summary>
    public double MinFrequency(int N, double J)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F41 requires N ≥ 2.");
        if (J <= 0) throw new ArgumentOutOfRangeException(nameof(J), J, "J must be > 0.");
        double s = Math.Sin(Math.PI / (HoppingCoefficient * N));
        return HoppingCoefficient * HoppingCoefficient * J * s * s;
    }

    /// <summary>Palindromic period: <c>t_Pi = π / (2·J · sin²(π/(2N)))</c>.
    /// Period of the slowest palindromic modulation in the SFF.</summary>
    public double PalindromicTime(int N, double J)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F41 requires N ≥ 2.");
        if (J <= 0) throw new ArgumentOutOfRangeException(nameof(J), J, "J must be > 0.");
        double s = Math.Sin(Math.PI / (HoppingCoefficient * N));
        return Math.PI / (HoppingCoefficient * J * s * s);
    }

    /// <summary>Asymptotic palindromic time at large N: <c>t_Pi → 2N² / (π·J)</c>.
    /// The exact PalindromicTime(N, J) approaches this value as N → ∞.</summary>
    public double AsymptoticPalindromicTime(int N, double J)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F41 requires N ≥ 2.");
        if (J <= 0) throw new ArgumentOutOfRangeException(nameof(J), J, "J must be > 0.");
        return HoppingCoefficient * N * N / (Math.PI * J);
    }

    /// <summary>Drift check: <c>t_Pi · ω_min = 2π</c> exactly (full-period definition).</summary>
    public bool PeriodFrequencyProductIsTwoPi(int N, double J, double tolerance = 1e-12)
    {
        double product = PalindromicTime(N, J) * MinFrequency(N, J);
        return Math.Abs(product - 2.0 * Math.PI) < tolerance;
    }

    public F41PalindromicTimePi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F1Pi2Inheritance f1)
        : base("F41 palindromic time t_Pi = π/(2J·sin²(π/(2N))) as Pi2-Foundation a_0 + F1 inheritance (palindrome identity at time-domain level)",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F41 + " +
               "experiments/SPECTRAL_FORM_FACTOR.md + " +
               "docs/proofs/derivations/D10_W1_DISPERSION.md + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f1 = f1 ?? throw new ArgumentNullException(nameof(f1));
    }

    public override string DisplayName =>
        "F41 palindromic time as Pi2-Foundation a_0 + F1 inheritance";

    public override string Summary =>
        $"t_Pi = π/(2J·sin²(π/(2N))); 2 = a_0; ω_min = 4J·sin²(π/(2N)); F1 palindrome at time-domain level; asymptotic t_Pi → 2N²/(π·J) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F41 closed form",
                summary: "t_Pi = π/(2·J·sin²(π/(2N))); ω_min = 4·J·sin²(π/(2N)); period of slowest palindromic SFF modulation; FFT-confirmed <1% at N=2..4, 6");
            yield return InspectableNode.RealScalar("HoppingCoefficient (= a_0 = 2)", HoppingCoefficient);
            yield return new InspectableNode("F1 palindrome at time domain",
                summary: $"the palindromic modulation exists because F1 pairs eigenvalues λ ↔ −λ − 2σ; F41 reads the period of the slowest such pair-difference. F1's TwoFactor (= {_f1.TwoFactor}) is the same '2' as F41's HoppingCoefficient.");
            yield return new InspectableNode("asymptotic scaling",
                summary: "sin(π/(2N)) → π/(2N) for N → ∞; t_Pi → 2N²/(π·J); the palindromic modulation is a short-time effect at finite N; t_Pi/t_H ~ N²/4^N → 0 (F42 timescale separation)");
            yield return new InspectableNode("N=3, J=1 verified",
                summary: $"sin²(π/6) = 1/4; ω_min = 4·1·(1/4) = 1; t_Pi = π/(2·1·(1/4)) = 2π ≈ {PalindromicTime(3, 1.0):G6}");
            yield return new InspectableNode("N=5, J=1 verified",
                summary: $"sin²(π/10) ≈ 0.0955; ω_min ≈ 0.382; t_Pi ≈ {PalindromicTime(5, 1.0):G6}; asymptotic = {AsymptoticPalindromicTime(5, 1.0):G6}");
            yield return new InspectableNode("N=20 large-N",
                summary: $"t_Pi ≈ {PalindromicTime(20, 1.0):G6}; asymptotic 2·400/π = {AsymptoticPalindromicTime(20, 1.0):G6}; convergence ~1/N² to leading order");
        }
    }
}
