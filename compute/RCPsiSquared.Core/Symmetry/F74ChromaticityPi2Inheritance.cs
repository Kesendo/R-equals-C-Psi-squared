using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F74 closed form (Tier 1 combinatorial; ANALYTICAL_FORMULAS line 1652):
///
/// <code>
///   c(n, N) = min(n, N − 1 − n) + 1     chromaticity of (n, n+1) coherence block
///   HD ∈ {1, 3, 5, ..., 2c − 1}          odd Hamming-distance ladder
///   pure rates = 2γ₀ · HD                Absorption Theorem at J = 0
/// </code>
///
/// <para>F74 is the combinatorial counting result for the (n, n+1) popcount
/// coherence block in an N-qubit system under uniform Z-dephasing. At J = 0
/// the block contains exactly c distinct pure dephasing rates, one per allowed
/// odd Hamming distance between popcount-n and popcount-(n+1) basis states.</para>
///
/// <para><b>Mirror structure of the chromaticity:</b></para>
/// <list type="bullet">
///   <item><b>Mono-chromatic (c = 1):</b> at n = 0 (vacuum-SE block) and n = N−1
///         (XOR-side); single pure rate 2γ₀. F73's spatial-sum coherence closure
///         applies to the n=0 monochromatic case.</item>
///   <item><b>Maximum chromaticity:</b> for odd N, c_max = (N+1)/2 at the unique
///         center block n = (N−1)/2. For even N, c_max = N/2 at the two adjacent
///         center blocks n = N/2 − 1 and n = N/2.</item>
///   <item><b>Mirror axis at n = (N−1)/2:</b> c(n, N) = c(N−1−n, N); the
///         chromaticity is palindromic in n. Same axis as F43's sector-pairing
///         (w → N−w with axis at N/2 for the full sector); F74 lives at the
///         (n, n+1)-coherence-block level with axis at (N−1)/2.</item>
/// </list>
///
/// <para>Pi2-Foundation anchor:</para>
/// <list type="bullet">
///   <item><b>RateCoefficient = 2 = a_0</b>: in 2γ₀·HD pure rates. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(0). Same anchor as F1
///         TwoFactor, F50 DecayRateFactor, F3 RateCoefficient, F43 XorRateCoefficient.</item>
/// </list>
///
/// <para>The min(n, N−1−n) structure and the +1 offset are combinatorial
/// (no Pi2 anchor). The HD = 2k+1 ladder gives only odd values because
/// popcount-n ↔ popcount-(n+1) coherences require an odd parity of differing
/// bit positions.</para>
///
/// <para>Tier1Derived: F74 is Tier 1 combinatorial; proof in
/// <c>docs/proofs/derivations/D6</c>-related material; verified N=3..7.
/// Pi2-Foundation anchoring is composition through Pi2DyadicLadder.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F74 (line 1652) +
/// <c>experiments/Q_SCALE_THREE_BANDS.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F74ChromaticityPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The "2" in pure rate 2γ₀·HD. Live from Pi2DyadicLadder a_0.</summary>
    public double RateCoefficient => _ladder.Term(0);

    /// <summary>Chromaticity c(n, N) = min(n, N−1−n) + 1.</summary>
    public int Chromaticity(int n, int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F74 requires N ≥ 2.");
        if (n < 0 || n >= N) throw new ArgumentOutOfRangeException(nameof(n), n, $"n must be in [0, {N - 1}]; got {n}.");
        return Math.Min(n, N - 1 - n) + 1;
    }

    /// <summary>Odd Hamming-distance ladder for the (n, n+1) coherence block:
    /// {1, 3, 5, ..., 2c − 1}.</summary>
    public IReadOnlyList<int> HammingDistanceValues(int n, int N)
    {
        int c = Chromaticity(n, N);
        var hds = new int[c];
        for (int k = 0; k < c; k++) hds[k] = 2 * k + 1;
        return hds;
    }

    /// <summary>Pure dephasing rate at HD: 2γ₀ · HD (Absorption Theorem at J = 0).</summary>
    public double PureRate(double gammaZero, int hammingDistance)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        if (hammingDistance < 1) throw new ArgumentOutOfRangeException(nameof(hammingDistance), hammingDistance, "HD must be ≥ 1.");
        return RateCoefficient * gammaZero * hammingDistance;
    }

    /// <summary>True iff the (n, n+1) block is monochromatic (c = 1, single pure rate).
    /// Holds at n = 0 (vac-SE) and n = N−1 (XOR-side).</summary>
    public bool IsMonochromatic(int n, int N) => Chromaticity(n, N) == 1;

    /// <summary>The maximum chromaticity for chain length N: (N+1)/2 for odd N,
    /// N/2 for even N. Attained at n = (N−1)/2 (odd N, unique) or
    /// n = N/2 − 1 and n = N/2 (even N, pair).</summary>
    public int MaxChromaticity(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F74 requires N ≥ 2.");
        return (N + 1) / 2;
    }

    /// <summary>True iff (n, n+1) is a center block achieving MaxChromaticity.</summary>
    public bool IsCenterBlock(int n, int N) => Chromaticity(n, N) == MaxChromaticity(N);

    /// <summary>Drift check: chromaticity is palindromic in n: c(n, N) = c(N − 1 − n, N).</summary>
    public bool MirrorPalindromicityHolds(int N)
    {
        for (int n = 0; n < N; n++)
            if (Chromaticity(n, N) != Chromaticity(N - 1 - n, N)) return false;
        return true;
    }

    public F74ChromaticityPi2Inheritance(Pi2DyadicLadderClaim ladder)
        : base("F74 chromaticity c(n, N) = min(n, N−1−n) + 1; pure rates 2γ₀·HD with HD ∈ {1, 3, ..., 2c−1}; mono at n=0/N−1, c_max at center; mirror palindromic in n",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F74 + " +
               "experiments/Q_SCALE_THREE_BANDS.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F74 chromaticity as Pi2-Foundation a_0 (in 2γ₀·HD pure rates) inheritance";

    public override string Summary =>
        $"c(n, N) = min(n, N−1−n) + 1; pure rates 2γ₀·HD; 2 = a_0; mirror axis n = (N−1)/2; mono-chromatic at n=0 and n=N−1 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F74 closed form",
                summary: "c(n, N) = min(n, N−1−n) + 1; HD ∈ {1, 3, ..., 2c−1}; pure rates 2γ₀·HD via Absorption Theorem at J = 0; combinatorial proof");
            yield return InspectableNode.RealScalar("RateCoefficient (= a_0 = 2)", RateCoefficient);
            yield return new InspectableNode("verified table N=5",
                summary: $"c(0,5)={Chromaticity(0,5)} (mono); c(1,5)={Chromaticity(1,5)}; c(2,5)={Chromaticity(2,5)} (center, max); c(3,5)={Chromaticity(3,5)}; c(4,5)={Chromaticity(4,5)} (mono)");
            yield return new InspectableNode("verified table N=6",
                summary: $"c(0,6)={Chromaticity(0,6)}; c(1,6)={Chromaticity(1,6)}; c(2,6)={Chromaticity(2,6)} (center pair); c(3,6)={Chromaticity(3,6)} (center pair); c(4,6)={Chromaticity(4,6)}; c(5,6)={Chromaticity(5,6)}");
            yield return new InspectableNode("max chromaticity",
                summary: $"MaxChromaticity(3) = {MaxChromaticity(3)}; MaxChromaticity(5) = {MaxChromaticity(5)}; MaxChromaticity(7) = {MaxChromaticity(7)}");
            yield return new InspectableNode("monochromatic blocks",
                summary: "n = 0 (vacuum-SE block) and n = N−1 are always c = 1; F73's spatial-sum coherence closure applies to the n=0 monochromatic case");
            yield return new InspectableNode("F43 sibling",
                summary: "F43's sector-pairing K_freq(w) = K_freq(N−w) lives at the full sector level (axis at N/2). F74's chromaticity palindrome lives at the (n, n+1)-coherence-block level (axis at (N−1)/2). Both are F1's Π palindrome at different decomposition layers.");
        }
    }
}
