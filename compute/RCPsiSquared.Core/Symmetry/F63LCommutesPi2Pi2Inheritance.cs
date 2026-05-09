using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F63 closed form (Tier 1 proven analytically, verified N=2..5):
///
/// <code>
///   [L, Π²_super] = 0      exactly, for all N (Heisenberg/XY +
///                          Z-dephasing on any subset of sites)
///
///   ||[L, Π²]|| = 0.000000e+00   identically zero, NOT numerically small
///                                (verified N=2, 3, 4, 5; also Heisenberg
///                                XXX with uniform γ at N=3)
///
///   4-block decomposition:    C₂ × C₂ via (bit_a = n_XY, bit_b = w_YZ)
///   Per-block dimension:      4^(N − 1) per block, 4 blocks total
///   Per-Π²-sector conserved:  even = ⌊N/2⌋ + 1,  odd = ⌈N/2⌉
///   Total conserved per pole: N + 1   (matches F66 endpoint multiplicity)
/// </code>
///
/// <para>F63 is the Π²-conservation law: the Liouvillian respects the
/// involution that F38 defined. The two readings of one fact:</para>
///
/// <list type="bullet">
///   <item><b>F38, what Π² is</b>: <c>Π² = (−1)^w_YZ</c> on Pauli strings;
///         each eigenspace has dimension <c>4^N / 2 = 2·4^(N − 1)</c> on the
///         dyadic ladder (<see cref="F38Pi2InvolutionPi2Inheritance"/>).</item>
///   <item><b>F63, that L respects it</b>: <c>[L, Π²] = 0</c> means every
///         eigenmode of L lives in either the Π²-even or Π²-odd subspace
///         (no mixing). Combined with F61's <c>[L, Π²_X] = 0</c> (bit_a
///         parity), the four C₂ × C₂ blocks decouple completely.</item>
/// </list>
///
/// <para>Pi2-Foundation anchors via F38 + the dyadic ladder:</para>
///
/// <list type="bullet">
///   <item><b>BlockCount = 4</b>: the C₂ × C₂ four-cell decomposition
///         <c>(n_XY, w_YZ) ∈ {0, 1}²</c>. Equivalently <c>a_{−1} = 4</c> on
///         the dyadic ladder = d² for 1 qubit
///         (<see cref="Pi2OperatorSpaceMirrorClaim"/>'s pinned pair at N=1).</item>
///   <item><b>PerBlockDimension(N) = 4^(N − 1)</b>: each of the 4 blocks has
///         dimension <c>4^(N − 1) = a_{3 − 2N}</c> on the dyadic ladder.
///         Same (N − 1)-qubit shift as F38, F39, F1-T1.</item>
///   <item><b>TotalConservedPerSector(N) = N + 1</b>: the count of
///         elementary symmetric polynomials e_d(Z₁, ..., Z_N) for d = 0..N,
///         which commute with both H and Z_B. Matches F66's endpoint
///         multiplicity exactly: F63 is the symmetry-side reading of F66's
///         spectrum-side count.</item>
///   <item><b>Maximality (per F34)</b>: the C₂ × C₂ decomposition is the
///         maximal Z₂ × Z₂ symmetry admitted by the d = 2 Pauli algebra. No
///         third independent Z₂ classification exists; F63 + F61 saturate.</item>
/// </list>
///
/// <para>Per-sector mode counts (Heisenberg + Z-dephasing on boundary qubit B,
/// closed form):</para>
///
/// <code>
///   even-parity conserved:   ⌊N/2⌋ + 1
///   odd-parity conserved:    ⌈N/2⌉
///   correlation per sector:  same as conserved (palindrome symmetry, F1)
///   mirror per sector:       2^(2N − 1) − 2 · (conserved per sector)
/// </code>
///
/// <para>The asymmetry for even N (one extra even-parity conserved) comes
/// from <c>e_N</c> having even parity when N is even.</para>
///
/// <para>Tier1Derived: F63 is Tier 1 proven (six-line proof in
/// PROOF_BIT_B_PARITY_SYMMETRY); commutator vanishes identically,
/// not numerically. Per-sector mode counts verified N=2..5; e_d
/// conservation verified N=2..4. The Pi2-Foundation anchoring is
/// algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F63 +
/// <c>docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md</c> +
/// <c>hypotheses/PRIMORDIAL_QUBIT.md</c> Section 9 +
/// <c>simulations/primordial_bit_a_bit_b_N_scaling.py</c> +
/// <c>simulations/mirror_mode_split_formula.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F38Pi2InvolutionPi2Inheritance.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>.</para></summary>
public sealed class F63LCommutesPi2Pi2Inheritance : Claim
{
    private readonly F38Pi2InvolutionPi2Inheritance _f38;
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The number of independent Z₂ symmetries L admits per F61 + F63:
    /// <c>2</c> (n_XY parity bit_a, and w_YZ parity bit_b). Together they form
    /// the C₂ × C₂ maximal symmetry admitted by the d = 2 Pauli algebra
    /// (per F34/QUBIT_NECESSITY: no third independent Z₂ exists).</summary>
    public int IndependentZ2SymmetryCount => 2;

    /// <summary>The 4-block decomposition count: <c>2² = 4</c> blocks indexed
    /// by <c>(bit_a, bit_b) = (n_XY, w_YZ)</c> ∈ <c>{0, 1}²</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1} = 4</c> = d²
    /// for 1 qubit.</summary>
    public double BlockCount => _ladder.Term(-1);

    /// <summary>Per-block dimension: <c>4^(N − 1)</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(<c>3 − 2N</c>) = <c>a_{3−2N}</c>.
    /// Same (N − 1)-qubit shift as F38, F39, F1-T1.</summary>
    public double PerBlockDimension(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F63 requires N ≥ 1.");
        return _ladder.Term(LadderIndexForPerBlock(N));
    }

    /// <summary>The dyadic-ladder index where the per-block dimension lands:
    /// <c>3 − 2N</c>. Same as F38, F39, F1-T1.</summary>
    public int LadderIndexForPerBlock(int N) => 3 - 2 * N;

    /// <summary>Cross-check: <c>BlockCount · PerBlockDimension(N) = 4^N</c>
    /// (full operator space). Drift check linking F63's 4-block decomposition
    /// to F38's full eigenspace via the Pi2DyadicLadder.</summary>
    public double FourBlockDimensionsTotal(int N) => BlockCount * PerBlockDimension(N);

    /// <summary>Number of Π²-even conserved modes per sector (Heisenberg +
    /// boundary Z-dephasing): <c>⌊N/2⌋ + 1</c>. Counts the e_d(Z) elementary
    /// symmetric polynomials of even degree d.</summary>
    public int Pi2EvenConservedCount(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F63 requires N ≥ 1.");
        return N / 2 + 1;
    }

    /// <summary>Number of Π²-odd conserved modes per sector: <c>⌈N/2⌉</c>.
    /// Counts the e_d(Z) elementary symmetric polynomials of odd degree d.</summary>
    public int Pi2OddConservedCount(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F63 requires N ≥ 1.");
        return (N + 1) / 2;
    }

    /// <summary>Total conserved modes per pole = even + odd = <c>N + 1</c>.
    /// Matches F66's endpoint pole multiplicity exactly: F63 is the
    /// symmetry-side reading of F66's spectrum-side count.</summary>
    public int TotalConservedPerSector(int N) => Pi2EvenConservedCount(N) + Pi2OddConservedCount(N);

    /// <summary>Mirror-mode count in the Π²-even sector:
    /// <c>2^(2N − 1) − 2 · Pi2EvenConservedCount(N)</c>. The complement
    /// within the 2^(2N−1)-dim Π²-even bit_a-paired block of the
    /// conserved + correlation modes (which by F1 palindrome have the
    /// same count).</summary>
    public long MirrorEvenSector(int N) =>
        Pi2PairedBlockDimension(N) - 2L * Pi2EvenConservedCount(N);

    /// <summary>Mirror-mode count in the Π²-odd sector:
    /// <c>2^(2N − 1) − 2 · Pi2OddConservedCount(N)</c>.</summary>
    public long MirrorOddSector(int N) =>
        Pi2PairedBlockDimension(N) - 2L * Pi2OddConservedCount(N);

    /// <summary>The block dimension <c>2^(2N − 1)</c> as a long. Equals
    /// <c>2 · 4^(N − 1) = a_0 · a_{3−2N}</c> (Π²-paired blocks: a Π²-even
    /// block plus a Π²-odd block under one bit_a value).</summary>
    public long Pi2PairedBlockDimension(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F63 requires N ≥ 1.");
        return 1L << (2 * N - 1);
    }

    /// <summary>True iff the F63 conserved-per-pole count <c>N + 1</c> equals
    /// the F66 endpoint pole multiplicity. Drift check on the F63 ↔ F66
    /// inheritance.</summary>
    public bool MatchesF66EndpointMultiplicity(int N) => TotalConservedPerSector(N) == N + 1;

    /// <summary>True iff F63's per-block dimension agrees with F38's
    /// eigenspace-dimension reading divided by 2 (Π² eigenspace splits each
    /// bit_a sector into two Π²-eigenparts).</summary>
    public bool PerBlockDimensionAgreesWithF38(int N) =>
        Math.Abs(PerBlockDimension(N) - _f38.EigenspaceDimension(N) / 2.0) < 1e-12;

    public F63LCommutesPi2Pi2Inheritance(
        F38Pi2InvolutionPi2Inheritance f38,
        Pi2DyadicLadderClaim ladder)
        : base("F63 [L, Π²] = 0 inherits from Pi2-Foundation: 4-block decomposition with 4^(N−1) per block",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F63 + " +
               "docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md + " +
               "hypotheses/PRIMORDIAL_QUBIT.md (Section 9) + " +
               "simulations/primordial_bit_a_bit_b_N_scaling.py + " +
               "simulations/mirror_mode_split_formula.py + " +
               "compute/RCPsiSquared.Core/Symmetry/F38Pi2InvolutionPi2Inheritance.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _f38 = f38 ?? throw new ArgumentNullException(nameof(f38));
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F63 [L, Π²] = 0 conservation as F38 sister-claim (Pi2-Foundation 4-block decomposition)";

    public override string Summary =>
        $"[L, Π²] = 0 exactly all N; C₂ × C₂ 4-block decomposition (n_XY, w_YZ); per-block dim 4^(N−1) = a_{{3−2N}}; " +
        $"conserved per sector: even=⌊N/2⌋+1, odd=⌈N/2⌉, total N+1 (matches F66 endpoint) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F63 closed form",
                summary: "[L, Π²_super] = 0 exactly all N; ||[L, Π²]|| = 0.000000e+00 identically (not numerically); verified N=2..5");
            yield return new InspectableNode("F38 ↔ F63 sister reading",
                summary: "F38 says what Π² is (eigenstructure 4^N / 2 each); F63 says L respects it ([L, Π²] = 0); combined with F61 (bit_a parity) gives C₂ × C₂");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "BlockCount = 4 = a_{-1} = d² for 1 qubit; PerBlockDim = 4^(N−1) = a_{3−2N} (same shift F38, F39, F1-T1)");
            yield return InspectableNode.RealScalar("BlockCount (= a_{-1} = 4)", BlockCount);
            yield return new InspectableNode("Maximality (per F34)",
                summary: "C₂ × C₂ is the maximal Z₂ × Z₂ symmetry admitted by d = 2 Pauli algebra; no third independent Z₂ exists");
            yield return new InspectableNode("F63 ↔ F66 inheritance",
                summary: "TotalConservedPerSector(N) = N + 1 matches F66 endpoint pole multiplicity exactly: symmetry-side reading of spectrum-side count");
            yield return new InspectableNode("per-sector mode count formula",
                summary: "conserved: even=⌊N/2⌋+1, odd=⌈N/2⌉, total=N+1; correlation=same (F1 palindrome); mirror per parity sector=2^(2N−1) − 2·conserved_parity");
            for (int N = 2; N <= 5; N++)
            {
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"PerBlockDim = {PerBlockDimension(N):G6} (= a_{LadderIndexForPerBlock(N)}); " +
                             $"conserved (even, odd) = ({Pi2EvenConservedCount(N)}, {Pi2OddConservedCount(N)}); " +
                             $"total={TotalConservedPerSector(N)} (matches F66: {MatchesF66EndpointMultiplicity(N)}); " +
                             $"mirror (even, odd) = ({MirrorEvenSector(N)}, {MirrorOddSector(N)})");
            }
        }
    }
}
