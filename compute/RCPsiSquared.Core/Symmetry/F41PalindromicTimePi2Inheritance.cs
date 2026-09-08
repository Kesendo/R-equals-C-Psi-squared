using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F41 closed form (Tier 1, corollary of D10):
///
/// <code>
///   t_Pi = 2π / ω_min = π / (4·J · sin²(π/(2N)))
///
///   ω_min = 4·J · (1 − cos(π/N)) = 8·J · sin²(π/(2N))
///                                  k=1 frequency of the D10 (0,1) block
/// </code>
///
/// <para>F41 is the full period of the k=1 frequency in D10's (0,1)
/// coherence block. F1 pairs that eigenvalue with a frequency-negated partner,
/// so the pair contributes <c>2cos(ω_min t)</c> to the frequency trace
/// amplitude. The SFF is the squared modulus of the full trace amplitude and
/// can therefore contain doubled and cross frequencies; F41 is not by itself
/// the period of the complete SFF.</para>
///
/// <para><b>Asymptotic scaling (large N):</b></para>
/// <code>
///   sin(π/(2N)) → π/(2N)  for N → ∞
///   sin²(π/(2N)) → π²/(4N²)
///   t_Pi → π / (4·J · π²/(4N²)) = N² / (π·J)
/// </code>
/// <para>So t_Pi grows as N²/(π·J), or equivalently N²/π at J = 1.
/// Comparison with the producer's multiplicity-dependent raw density scale is
/// a finite-N descriptive ratio, not a physical short/long-time division.</para>
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
/// <para>F1 connection: D10 supplies ω_min; F1 supplies the partner at −ω_min.
/// Together they give the cosine term in the trace amplitude. Neither identity
/// promotes the producer's raw multiset density scale to a physical timescale.</para>
///
/// <para>Tier1Derived: F41 is a Tier 1 corollary of D10's (0,1) coherence-block
/// dispersion. The separate finite-N producer associates an FFT candidate with
/// ω_min at N=2..4 and N=6, but not at N=5 or N=7. Pi2-Foundation anchoring is
/// algebraic-trivial composition through the 2·J coefficient.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F41 (line 838) +
/// <c>experiments/SPECTRAL_FORM_FACTOR.md</c> +
/// <c>docs/proofs/derivations/D10_W1_DISPERSION.md</c> (parent derivation) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs</c> (palindrome
/// identity at spectral level) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F41PalindromicTimePi2Inheritance : Claim, IZ2AxisClaim
{

    /// <summary>The F1² / Π²_Z axis (bit_b parity, n_Y + n_Z mod 2). The
    /// canonical Pi²-Inheritance axis. The bit_a-twin (Π²_X / F61 axis) is
    /// currently not typed for this Claim.</summary>
    public Z2Axis Z2Axis => Z2Axis.BitB;

    /// <summary>The typed bit_a-twin sibling, if one exists. Currently null
    /// (no bit_a twin is typed for this Claim; this is an open slot in the
    /// cubic-architecture coverage).</summary>
    public Claim? BitATwin => null;
    public Pi2DyadicLadderClaim Ladder { get; }
    public F1Pi2Inheritance F1 { get; }

    /// <summary>The "2" hopping coefficient, the amplitude 2J that XX + YY gives
    /// each bond. Live from Pi2DyadicLadder a_0; same anchor as F1's TwoFactor.
    /// F41 prints its powers rather than a_0 itself: ω_min carries
    /// a_0³ = 8 = a_{−2}, and t_Pi = 2π/ω_min carries a_0² = 4 = a_{−1}.
    /// Three rungs, one anchor.</summary>
    public double HoppingCoefficient => Ladder.Term(0);

    /// <summary>The k=1 frequency of D10's (0,1) coherence block:
    /// <c>ω_min = 4·J · (1 − cos(π/N)) = 8·J · sin²(π/(2N))</c>, the k=1 mode of
    /// the D10 (0,1) coherence block dispersion. Equals a_0²·J·(1 − cos(π/N)).</summary>
    public double MinFrequency(int N, double J)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F41 requires N ≥ 2.");
        if (J <= 0) throw new ArgumentOutOfRangeException(nameof(J), J, "J must be > 0.");
        // omega_min is the k=1 mode of the D10 (0,1) block dispersion omega_k = 4J(1 - cos(pi k / N)),
        // and 1 - cos(pi/N) = 2 sin^2(pi/(2N)), so omega_min = 8 J sin^2(pi/(2N)).
        return HoppingCoefficient * HoppingCoefficient * J * (1.0 - Math.Cos(Math.PI / N));
    }

    /// <summary>Palindromic-pair trace-amplitude period:
    /// <c>t_Pi = π / (4·J · sin²(π/(2N)))</c>.</summary>
    public double PalindromicTime(int N, double J)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F41 requires N ≥ 2.");
        if (J <= 0) throw new ArgumentOutOfRangeException(nameof(J), J, "J must be > 0.");
        double s = Math.Sin(Math.PI / (HoppingCoefficient * N));
        return Math.PI / (HoppingCoefficient * HoppingCoefficient * J * s * s);
    }

    /// <summary>Asymptotic palindromic time at large N: <c>t_Pi → N² / (π·J)</c>.
    /// The exact PalindromicTime(N, J) approaches this value as N → ∞.</summary>
    public double AsymptoticPalindromicTime(int N, double J)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F41 requires N ≥ 2.");
        if (J <= 0) throw new ArgumentOutOfRangeException(nameof(J), J, "J must be > 0.");
        return (double)N * N / (Math.PI * J);
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
        : base("F41 palindromic-pair trace-amplitude period t_Pi = π/(4J·sin²(π/(2N))); D10 supplies ω_min and F1 supplies its frequency-negated partner",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F41 + " +
               "experiments/SPECTRAL_FORM_FACTOR.md + " +
               "docs/proofs/derivations/D10_W1_DISPERSION.md + " +
               "compute/RCPsiSquared.Core/Symmetry/F1Pi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        F1 = f1 ?? throw new ArgumentNullException(nameof(f1));
    }

    public override string DisplayName =>
        "F41 palindromic-pair trace-amplitude period as D10 + F1 inheritance";

    public override string Summary =>
        $"t_Pi = π/(4J·sin²(π/(2N))) in the Pauli normalisation H = J·Σ(XX+YY+ZZ); ω_min = 8J·sin²(π/(2N)); D10 supplies ω_min, F1 supplies the −ω_min partner, and the pair contributes 2cos(ω_min·t) to the trace amplitude; asymptotic t_Pi → N²/(π·J). The ring and star Im-max claims use the spin normalisation J·Σ S_i·S_j, which is this J divided by 4 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F41 closed form",
                summary: "t_Pi = π/(4·J·sin²(π/(2N))); ω_min = 8·J·sin²(π/(2N)) = 4·J·(1 − cos(π/N)); full period of the k=1 D10 frequency and its F1 cosine pair in the trace amplitude");
            yield return InspectableNode.RealScalar("HoppingCoefficient (= a_0 = 2)", HoppingCoefficient);
            yield return new InspectableNode("F1 palindrome at time domain",
                summary: $"D10 supplies ω_min and F1 pairs it with −ω_min; together they contribute 2cos(ω_min·t) to the trace amplitude. Squaring the full sum can create doubled and cross frequencies. F1's TwoFactor (= {F1.TwoFactor}) is the same '2' as F41's HoppingCoefficient.");
            yield return new InspectableNode("asymptotic scaling",
                summary: "sin(π/(2N)) → π/(2N) for N → ∞; t_Pi → N²/(π·J). Any ratio to the producer's raw multiset density scale is descriptive and finite-N, not a physical time-regime boundary");
            yield return new InspectableNode("N=3, J=1 verified",
                summary: $"sin²(π/6) = 1/4; ω_min = 8·1·(1/4) = 2; t_Pi = π/(4·1·(1/4)) = π ≈ {PalindromicTime(3, 1.0):G6}");
            yield return new InspectableNode("N=5, J=1 verified",
                summary: $"sin²(π/10) ≈ 0.0955; ω_min = 8·0.0955 ≈ 0.764; t_Pi ≈ {PalindromicTime(5, 1.0):G6}; asymptotic = {AsymptoticPalindromicTime(5, 1.0):G6}");
            yield return new InspectableNode("N=20 large-N",
                summary: $"t_Pi ≈ {PalindromicTime(20, 1.0):G6}; asymptotic 400/π = {AsymptoticPalindromicTime(20, 1.0):G6}; convergence ~1/N² to leading order");
        }
    }
}
