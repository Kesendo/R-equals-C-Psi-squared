using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F62 closed form (Tier 1, analytical, verified N=2..10):
///
/// <code>
///   CΨ(0) = 2(N² − 4N + 8) / (3 N³)
///
///   ρ_ab = diag((N−2)/N, 1/N, 1/N, 0)
///        + (1/N)|01⟩⟨10| + (1/N)|10⟩⟨01|
///
///   Tr(ρ_ab²) = (N² − 4N + 8)/N²,  L1 = 2/N,  Ψ = 2/(3N)
///
///   For N ≥ 3: CΨ(0) &lt; 1/4 (W_N born below the fold).
/// </code>
///
/// <para>F62 is the W-state companion to <see cref="F60GhzBornBelowFoldPi2Inheritance"/>:
/// both are pair-CΨ closed forms for canonical N-qubit states that fail to cross
/// the bilinear-apex fold for N ≥ 3. The two states sit at opposite ends of the
/// popcount spectrum:</para>
///
/// <code>
///   F60: GHZ_N = (|0...0⟩ + |1...1⟩)/√2,  popcount {0, N}
///        CΨ(0) = 1/(2^N − 1)              [Hilbert-space dimension scaling]
///
///   F62: W_N   = (1/√N) Σ_j |1_j⟩,        popcount 1
///        CΨ(0) = 2(N² − 4N + 8)/(3 N³)    [polynomial in N]
/// </code>
///
/// <para>Both fall below 1/4 for N ≥ 3, but with different mechanisms:</para>
///
/// <list type="bullet">
///   <item><b>F60 mechanism</b>: GHZ has only TWO basis states, so the
///         off-diagonal weight 1/2 is split among 2^N − 1 mass-positions of
///         the d²-dim Hilbert space; CΨ scales as 1/(2^N − 1).</item>
///   <item><b>F62 mechanism</b>: W_N has N basis states, each with weight
///         1/N; pair-reduced ρ has off-diagonal 1/N; CΨ scales as 1/N.
///         Cubic decay 1/N³ comes from L1 (2/N), Ψ (2/(3N)), and Tr(ρ²)
///         all carrying 1/N factors.</item>
/// </list>
///
/// <para>Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>FoldPosition = 1/4 = a_3</b>: the bilinear-apex maxval; same
///         anchor as F57, Dicke, F60. F62's "below the fold" assertion
///         compares CΨ(0) against <see cref="QuarterAsBilinearMaxvalClaim"/>.</item>
///   <item><b>"2" coefficient = a_0</b>: the numerator factor "2" in
///         <c>2(N² − 4N + 8)</c>. Live from <see cref="Pi2DyadicLadderClaim.Term"/>(0)
///         = polynomial root d. Same anchor as F1, F66, F86 Q_EP, F60
///         (mirror partner of F60's 1/2 off-diagonal: a_0 ↔ a_2).</item>
///   <item><b>F61 cited in F62 proof</b>: per ANALYTICAL_FORMULAS F62
///         "Combined with the Parity Selection Rule (F61), this proves that
///         single-excitation states on Heisenberg chains under Z-dephasing
///         never cross CΨ = 1/4." F61 is the bit_a parity claim: W-states
///         have purely even n_XY content, so SE optimisers cannot reach
///         odd-n_XY modes (per F61's accessibility ceiling). F62 + F61
///         together close the "single-excitation regime is structurally
///         outside the framework" reading.</item>
///   <item><b>"3 N³" denominator</b>: combinatorial (3 from L1 normalization
///         2/3, N³ from the three 1/N factors). NOT Pi2-anchored.</item>
/// </list>
///
/// <para>Tier1Derived: F62 is Tier 1 analytical; verified numerically
/// N=2..10 within machine precision. The Pi2-Foundation anchoring is
/// algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F62 (line 1255) +
/// <c>simulations/cpsi_wn_analytical.py</c> +
/// <c>experiments/CUSP_LENS_CONNECTION.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (QuarterAsBilinearMaxvalClaim) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F60GhzBornBelowFoldPi2Inheritance.cs</c>
/// (sibling claim) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F61BitAParityPi2Inheritance.cs</c>
/// (cited parity selection rule).</para></summary>
public sealed class F62WStateBornBelowFoldPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The fold position <c>1/4 = a_3</c> on the dyadic ladder.
    /// Same anchor as F57, Dicke, F60; F62's "below the fold" assertion
    /// compares CΨ(0) against <see cref="QuarterAsBilinearMaxvalClaim"/>.</summary>
    public double FoldPosition => _ladder.Term(3);

    /// <summary>The "2" numerator coefficient in F62's closed form. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = <c>a_0</c> = polynomial
    /// root d. Same anchor as F1, F66, F86 Q_EP.</summary>
    public double NumeratorTwoCoefficient => _ladder.Term(0);

    /// <summary>The off-diagonal element of W_N's reduced pair density matrix:
    /// <c>1/N</c>. NOT Pi2-anchored as constant (N-dependent), but structurally
    /// the per-site weight from the equal-superposition normalization.</summary>
    public double PairOffDiagonalElement(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F62 requires N ≥ 2.");
        return 1.0 / N;
    }

    /// <summary>Live closed form: <c>CΨ(0) = 2(N² − 4N + 8) / (3 N³)</c>.</summary>
    public double CPsiAtZeroForWState(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F62 requires N ≥ 2.");
        return NumeratorTwoCoefficient * (N * N - 4 * N + 8) / (3.0 * N * N * N);
    }

    /// <summary>True iff W_N is born below the fold:
    /// <c>CΨ(0) &lt; 1/4 ⇔ N ≥ 3</c>. Proof per F62 cubic monotonicity
    /// 3N³ − 8N² + 32N − 64 > 0 for N ≥ 3 (= 41 at N=3).</summary>
    public bool IsBornBelowFold(int N) => CPsiAtZeroForWState(N) < FoldPosition;

    /// <summary>The smallest N at which W_N is born below the fold:
    /// <c>N = 3</c>. At N = 2, W_2 = Bell+ (above fold, crosses under
    /// dynamics). Same threshold as F60 (GHZ_N).</summary>
    public int SmallestNBelowFold => 3;

    /// <summary>Cross-check: at N = 2 (W_2 as Bell+), CΨ(0) = 1/3 ≈ 0.333,
    /// above the fold. Drift indicator on the "Bell+ as edge case" reading.</summary>
    public bool BellPlusAboveFold() => CPsiAtZeroForWState(2) > FoldPosition;

    public F62WStateBornBelowFoldPi2Inheritance(Pi2DyadicLadderClaim ladder)
        : base("F62 CΨ(0)_W_N = 2(N²−4N+8)/(3N³) inherits from Pi2-Foundation: 2 = a_0; fold = a_3 = 1/4; F61 parity selection rule cited",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F62 + " +
               "simulations/cpsi_wn_analytical.py + " +
               "experiments/CUSP_LENS_CONNECTION.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (QuarterAsBilinearMaxval) + " +
               "compute/RCPsiSquared.Core/Symmetry/F60GhzBornBelowFoldPi2Inheritance.cs (sibling) + " +
               "compute/RCPsiSquared.Core/Symmetry/F61BitAParityPi2Inheritance.cs (cited)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "F62 CΨ(0)_W_N as Pi2-Foundation inheritance (W-state companion to F60 GHZ)";

    public override string Summary =>
        $"CΨ(0)_W_N = 2(N²−4N+8)/(3N³): fold = a_3 = 1/4; numerator 2 = a_0; off-diagonal 1/N (combinatorial); " +
        $"born below fold for N ≥ 3 (cubic monotonicity); F61 parity rule completes the proof ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F62 closed form",
                summary: "CΨ(0)_W_N = 2(N²−4N+8)/(3N³); Tier 1 analytical; verified numerically N=2..10 within machine precision");
            yield return new InspectableNode("F60 ↔ F62 sibling pair",
                summary: "F60 (GHZ, popcount {0,N}): CΨ(0) = 1/(2^N − 1) [exponential in N]; F62 (W, popcount 1): CΨ(0) = 2(N²−4N+8)/(3N³) [polynomial in N]; both below fold for N ≥ 3");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "FoldPosition = a_3 = 1/4 (same as F57, Dicke, F60); NumeratorTwoCoefficient = a_0 = 2 (same as F1, F66, F86 Q_EP); 1/N off-diagonal combinatorial");
            yield return InspectableNode.RealScalar("FoldPosition (= a_3 = 1/4)", FoldPosition);
            yield return InspectableNode.RealScalar("NumeratorTwoCoefficient (= a_0 = 2)", NumeratorTwoCoefficient);
            yield return new InspectableNode("F61 parity-rule citation",
                summary: "F62 + F61 together: SE states have purely even n_XY (F61); W_N is SE; W_N born below fold for N ≥ 3 (F62); SE optimisers structurally cannot cross 1/4 on Heisenberg + Z-dephasing");
            yield return new InspectableNode("operational consequence",
                summary: "single-excitation regime is outside framework's quantum band for N ≥ 3, regardless of γ. Multi-excitation encodings required (Dicke, mirror-pair, etc.)");
            // Verified table from ANALYTICAL_FORMULAS F62 at N ∈ {2, 3, 5, 10}
            foreach (int N in new[] { 2, 3, 5, 10 })
            {
                double cpsi = CPsiAtZeroForWState(N);
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"CΨ(0) = 2({N}²−4·{N}+8)/(3·{N}³) = {cpsi:G6}; " +
                             $"above fold: {cpsi > FoldPosition} ({(cpsi > FoldPosition ? "Bell+ regime, W_2 = Bell+" : "born classical")})");
            }
        }
    }
}
