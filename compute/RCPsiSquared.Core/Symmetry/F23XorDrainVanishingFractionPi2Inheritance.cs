using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F23 closed form (Tier 1, combinatorial proof):
///
/// <code>
///   fraction(XOR) = (N + 1) / 4^N
/// </code>
///
/// <para>F23 is the structural counting result for the XOR drain sector
/// fraction in the 4^N Pauli operator space on N qubits. The XOR drain is the
/// sector that hosts GHZ-fragility for small N; F23 says this sector has
/// vanishing measure at macroscopic N. Numerical examples:</para>
///
/// <list type="bullet">
///   <item>N = 3: 6.25%</item>
///   <item>N = 5: 0.59%</item>
///   <item>N = 8: 0.014%</item>
///   <item>N = 20: ≈ 10⁻¹¹</item>
/// </list>
///
/// <para>The (N+1) numerator counts the XOR-drain eigenmodes, one per popcount
/// sector (the Pauli strings of XY-weight N that span the drain number 2^N);
/// the 4^N denominator is the full Pauli-operator-space dimension. F23 typifies
/// the contrast: linear growth in N (numerator) against exponential decay
/// (denominator), making the XOR drain's measure go to zero exponentially fast.
/// The drain's share of operator space vanishes with N while GHZ_N's coherence sits
/// entirely in it at every N (F22): what shrinks is the drain's measure, not GHZ's exposure.</para>
///
/// <para>Pi2-Foundation anchors:</para>
/// <list type="bullet">
///   <item><b>BaseFactor = 4 = a_{−1}</b>: in 4^N denominator. Live from
///         <see cref="Pi2DyadicLadderClaim.Term"/>(−1). Same anchor as F25
///         decay rate (e^{−4γt}), F65/F76 numerators, F73 spatial-sum closure,
///         F61/F63 4-block per parity.</item>
///   <item><b>OperatorSpaceDim = 4^N</b>: the full Pauli operator-space
///         dimension d² where d = 2 (qubit). Sibling of
///         <see cref="Pi2OperatorSpaceMirrorClaim"/> which typifies the
///         4^N ↔ 1/4^N mirror pair (a_n · a_{2−n} = 1 at n = −(2N−1)).</item>
/// </list>
///
/// <para>F23 is the "operator-space-dimension counting result" sibling of
/// other 4^N-based claims: F61/F63 use 4-block per parity (4^N total split
/// into 4-blocks × N parities); F73 uses 4 · γ_T1 decay coefficient; F23
/// uses 4^N as the denominator of fraction-counting.</para>
///
/// <para>Tier1Derived: the (N+1) count is the eigenvectors X^⊗N·P_k, one per popcount
/// sector: F4's kernel dimension for a connected chain, Π_c(|c|+1) over components
/// (PROOF_F4_KERNEL_DIMENSION_BY_COMPONENTS.md), carried onto the −2Σγ eigenspace by
/// F1 where F1 holds, and the same vectors directly under a uniform field; the table is in
/// <c>experiments/N_INFINITY_PALINDROME.md</c>. Combined with the 4^N total it gives
/// the closed form directly. Valid for any N under Z-dephasing on the
/// XY/Heisenberg chain with every bond non-zero and every site dephased (γ_l > 0),
/// whose count the N+1 is (a zero bond raises it to Π_c(|c|+1); a zero γ on one seat
/// leaves N+1 at N = 3, 4, while γ on the N = 3 middle seat alone gives 12 on XY and 10
/// on Heisenberg; a longitudinal field keeps it only if uniform); number conservation alone does not
/// give it (XY plus a non-uniform longitudinal field keeps [H, ΣZ] = 0 and leaves 2
/// modes at N = 3, 4), a pure Ising ZZ chain gives 2^N and a generic Hermitian H
/// none (GLOSSARY, the −2Σγ row).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F23 +
/// <c>experiments/N_INFINITY_PALINDROME.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>
/// (sibling).</para></summary>
public sealed class F23XorDrainVanishingFractionPi2Inheritance : Claim, IZ2AxisClaim
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
    public Pi2OperatorSpaceMirrorClaim Mirror { get; }

    /// <summary>The "4" base in 4^N. Live from Pi2DyadicLadder a_{−1}. Same anchor
    /// as F25/F65/F73/F76 decay rate numerators.</summary>
    public double BaseFactor => Ladder.Term(-1);

    /// <summary>Live computation of the operator-space dimension 4^N. For N=1..6
    /// matches Pi2OperatorSpaceMirrorClaim.PairAt(N).OperatorSpace exactly; for
    /// N ≥ 7 extends the same pattern via BaseFactor^N.</summary>
    public double OperatorSpaceDim(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F23 requires N ≥ 1.");
        return Math.Pow(BaseFactor, N);
    }

    /// <summary>The XOR-drain eigenmode count: (N + 1), one per popcount sector.</summary>
    public int XorDrainCount(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F23 requires N ≥ 1.");
        return N + 1;
    }

    /// <summary>The XOR-drain fraction <c>(N + 1) / 4^N</c> of all 4^N Pauli strings.
    /// Vanishes exponentially: 6.25% at N=3, 0.59% at N=5, 0.014% at N=8, ≈10⁻¹¹ at N=20.</summary>
    public double XorDrainFraction(int N)
    {
        return XorDrainCount(N) / OperatorSpaceDim(N);
    }

    /// <summary>Whether the XOR drain is "macroscopically negligible" at N.
    /// Threshold = 10⁻⁶ (one part per million); first reached at N = 12
    /// ((N+1)/4^N = 2.9·10⁻⁶ at N = 11, 7.8·10⁻⁷ at N = 12). The threshold is a
    /// convention; what it reads is the drain's share of operator space, not GHZ's
    /// exposure to it, which is total at every N.</summary>
    public bool IsMacroscopicallyNegligible(int N, double threshold = 1e-6)
    {
        if (threshold <= 0) throw new ArgumentOutOfRangeException(nameof(threshold), threshold, "threshold must be > 0.");
        return XorDrainFraction(N) < threshold;
    }

    /// <summary>Drift check: F23's denominator agrees with Pi2OperatorSpaceMirror's
    /// pinned table for N ∈ [1, 6]. Returns true iff <c>OperatorSpaceDim(N)</c>
    /// matches <c>Mirror.PairAt(N).OperatorSpace</c>.</summary>
    public bool MatchesMirrorTable(int N, double tolerance = 1e-12)
    {
        if (N < 1 || N > 6) throw new ArgumentOutOfRangeException(nameof(N), N, "Mirror table covers N=1..6 only.");
        var pair = Mirror.PairAt(N);
        if (pair == null) return false;
        return Math.Abs(OperatorSpaceDim(N) - pair.OperatorSpace) < tolerance;
    }

    public F23XorDrainVanishingFractionPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        Pi2OperatorSpaceMirrorClaim mirror)
        : base("F23 XOR drain vanishing fraction: fraction(XOR) = (N+1)/4^N; 4 = a_{-1}; 4^N is Pi2OperatorSpaceMirror sibling; vanishes exponentially at macroscopic N",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F23 + " +
               "experiments/N_INFINITY_PALINDROME.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs")
    {
        Ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        Mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
    }

    public override string DisplayName =>
        "F23 XOR drain (N+1)/4^N as Pi2-Foundation a_{-1} + OperatorSpaceMirror inheritance";

    public override string Summary =>
        $"fraction(XOR) = (N+1)/4^N; 4 = a_{{-1}} (= {BaseFactor}); vanishes exponentially: 6.25% at N=3, 0.59% at N=5, ≈10⁻¹¹ at N=20; the drain's share of operator space vanishes with N, GHZ_N's exposure to it does not ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F23 closed form",
                summary: "fraction(XOR) = (N+1)/4^N; the (N+1) numerator counts XOR-drain eigenmodes, one per popcount sector (the eigenvectors X^N P_k; F4's kernel dimension carried onto -2*Sigma_gamma by F1 where F1 holds; the table is in N_INFINITY_PALINDROME); 4^N denominator is full Pauli operator-space");
            yield return InspectableNode.RealScalar("BaseFactor (= a_{-1} = 4)", BaseFactor);
            yield return new InspectableNode("Pi2OperatorSpaceMirror sibling",
                summary: "F23's 4^N denominator IS Pi2OperatorSpaceMirror's OperatorSpace; both share the same a_{-1}^N = (4)^N = d²·d²·...·d² = (2^N)² ladder anchor at index n = -(2N-1)");
            yield return new InspectableNode("vanishing-measure reading",
                summary: "GHZ_N's coherence lives in the XOR drain at every N; what vanishes with N (below one ppm from N = 12) is the drain's share of operator space, what a random operator puts there. F23 quantifies the drain's measure, not GHZ's exposure.");
            yield return new InspectableNode("N=3 verified",
                summary: $"XorDrainCount={XorDrainCount(3)}, OperatorSpaceDim={OperatorSpaceDim(3)}, fraction={XorDrainFraction(3):P4} (6.25%)");
            yield return new InspectableNode("N=5 verified",
                summary: $"XorDrainCount={XorDrainCount(5)}, OperatorSpaceDim={OperatorSpaceDim(5)}, fraction={XorDrainFraction(5):P4} (0.59%)");
            yield return new InspectableNode("N=8 verified",
                summary: $"XorDrainCount={XorDrainCount(8)}, OperatorSpaceDim={OperatorSpaceDim(8)}, fraction={XorDrainFraction(8):P6} (0.014%)");
            yield return new InspectableNode("N=20 verified",
                summary: $"XorDrainCount={XorDrainCount(20)}, OperatorSpaceDim={OperatorSpaceDim(20):G6}, fraction={XorDrainFraction(20):G4}");
        }
    }
}
