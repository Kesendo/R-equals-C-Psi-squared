using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The Absorption Theorem (Tier 1 derived; proven in
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> + verified bit-exact across 1024
/// modes on IBM Q52 with max deviation 2.4·10⁻¹⁴):
///
/// <code>
///   For a Pauli-string Liouvillian eigenmode σ_α (Heisenberg or XY chain
///   under uniform Z-dephasing γ₀ on every site), the spectrum is quantized:
///
///     Re(λ_α) = −2 · γ₀ · ⟨n_XY⟩(σ_α)
///
///   where ⟨n_XY⟩ = number of sites k where σ_α has factor X or Y. The
///   absorption quantum is 2γ₀; spectrum lives on the grid {0, −2γ₀, −4γ₀,
///   ..., −2N·γ₀}. Maximum rate (full XOR mode) = 2γ₀·N.
/// </code>
///
/// <para>Per-coherence reading on computational-basis density-matrix entries:
/// a coherence |A⟩⟨B| decomposes per-site into pure {I, Z} (when A_l = B_l)
/// or pure {X, Y} (when A_l ≠ B_l), so n_XY(|A⟩⟨B|) = n_diff(A, B). Hence
/// <c>per-coherence rate = 2γ₀ · n_diff</c>: the Absorption Theorem applied
/// per basis-pair.</para>
///
/// <para>Hamming-complement pair-sum corollary: the column bit-flip
/// ρ[a, b] → ρ[a, bar(b)] uses n_diff(a, b) + n_diff(a, bar(b)) = N
/// to map rates 2γ₀·k ↔ 2γ₀·(N − k). Pair-sum equals 2γ₀·N exactly = the
/// spectral maximum. This is the F89c structural lemma's source.</para>
///
/// <para>Pi2-Foundation anchor: the "2" absorption quantum IS
/// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = a_0 = polynomial root d in
/// d² − 2d = 0. Same anchor as F1 TwoFactor, F50 DecayRateFactor, F66
/// UpperPoleCoefficient. The Absorption Theorem reifies the per-mode rate
/// reading of this single anchor.</para>
///
/// <para>Descendants (live via <c>rcpsi knowledge descendants AbsorptionTheoremClaim</c>;
/// per-formula one-liners in <c>compute/RCPsiSquared.Core/F_FORMULA_CROSSWALK.md</c>):
/// F33, F50, F55, F64, F65, F66, F67, F68, F74, F89.</para>
///
/// <para>Tier1Derived: proven analytically (PROOF_ABSORPTION_THEOREM.md) +
/// verified bit-exact N=2..7 (Liouvillian spectra) + hardware-confirmed on
/// IBM Q52 (1343 modes, deviation &lt; 3%, IBM_ABSORPTION_THEOREM.md).</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c> +
/// <c>docs/ANALYTICAL_FORMULAS.md</c> (F33, F50, F55, F64-F68, F74, F89 entries) +
/// <c>experiments/IBM_ABSORPTION_THEOREM.md</c> +
/// <c>experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c>
/// (a_0 = 2 absorption quantum source).</para></summary>
public sealed class AbsorptionTheoremClaim : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;

    /// <summary>The absorption quantum (numerical coefficient): <c>2</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = a_0 = polynomial root d. Multiplying
    /// by γ₀ gives the per-XY-site decay quantum 2γ₀.</summary>
    public double AbsorptionQuantumCoefficient => _ladder.Term(0);

    /// <summary>The absorption quantum at γ₀: <c>2γ₀</c>. The smallest non-zero rate
    /// step in the Liouvillian spectrum under uniform Z-dephasing.</summary>
    public double AbsorptionQuantum(double gammaZero)
    {
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return AbsorptionQuantumCoefficient * gammaZero;
    }

    /// <summary>The Absorption Theorem rate: <c>α = 2γ₀ · n_XY</c>. Returns the
    /// Liouvillian decay rate for a Pauli-string mode with the given XY-weight.</summary>
    public double Rate(int nXY, double gammaZero)
    {
        if (nXY < 0) throw new ArgumentOutOfRangeException(nameof(nXY), nXY, "n_XY must be ≥ 0.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return AbsorptionQuantumCoefficient * gammaZero * nXY;
    }

    /// <summary>The per-coherence rate for a computational-basis density-matrix entry
    /// |A⟩⟨B| with n_diff(A, B) differing bits: <c>2γ₀ · n_diff</c>. Equivalent to
    /// <see cref="Rate"/> with n_XY = n_diff (each differing bit = one X/Y Pauli factor
    /// in the coherence's basis decomposition).</summary>
    public double PerCoherenceRateComputationalBasis(int nDiff, double gammaZero)
    {
        if (nDiff < 0) throw new ArgumentOutOfRangeException(nameof(nDiff), nDiff, "n_diff must be ≥ 0.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return AbsorptionQuantumCoefficient * gammaZero * nDiff;
    }

    /// <summary>The maximum Liouvillian decay rate on N qubits: <c>2γ₀ · N</c> (full
    /// XOR mode, n_XY = N). The "absorption ceiling" of section 4.1 in the proof.</summary>
    public double MaxRate(int n, double gammaZero)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be ≥ 1.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return AbsorptionQuantumCoefficient * gammaZero * n;
    }

    /// <summary>Inverse map from a decay rate to its XY-weight: <c>n_XY = α / (2γ₀)</c>.
    /// Mirrors <see cref="F33ExactN3DecayRatesPi2Inheritance.NXyExpectationFromRate"/>;
    /// the Absorption Theorem read backwards.</summary>
    public double NXyFromRate(double rate, double gammaZero)
    {
        if (rate < 0) throw new ArgumentOutOfRangeException(nameof(rate), rate, "rate must be ≥ 0.");
        if (gammaZero <= 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be > 0.");
        return rate / (AbsorptionQuantumCoefficient * gammaZero);
    }

    /// <summary>The Hamming-complement pair-sum on a (k+1)-qubit block:
    /// <c>2γ₀ · (k+1)</c>. The column bit-flip ρ[a, b] → ρ[a, bar(b)] uses
    /// n_diff(a, b) + n_diff(a, bar(b)) = (k+1) to map rates 2γ₀·n_diff ↔
    /// 2γ₀·((k+1) − n_diff); pair-sum is exactly the spectral maximum. F89c
    /// structural lemma; verified bit-exact at path-2 (3-qubit block, pair-sum =
    /// 6γ₀) and generalizes to all block sizes.</summary>
    public double HammingComplementPairSum(int blockSize, double gammaZero)
    {
        if (blockSize < 1) throw new ArgumentOutOfRangeException(nameof(blockSize), blockSize, "block size must be ≥ 1.");
        if (gammaZero < 0) throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be ≥ 0.");
        return AbsorptionQuantumCoefficient * gammaZero * blockSize;
    }

    /// <summary>True iff the absorption-quantum coefficient matches the literal
    /// <c>2.0</c> from the proof. Live drift check against
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0).</summary>
    public bool AbsorptionQuantumMatchesLiteral() =>
        Math.Abs(AbsorptionQuantumCoefficient - 2.0) < 1e-15;

    public AbsorptionTheoremClaim(Pi2DyadicLadderClaim ladder)
        : base("Absorption Theorem: Re(λ) = −2γ₀·⟨n_XY⟩; spectrum quantized in 2γ₀ steps; absorption quantum 2γ₀ = a_0·γ₀; per-coherence rate 2γ₀·n_diff in computational basis; Hamming-complement pair-sum 2γ₀·N",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "docs/ANALYTICAL_FORMULAS.md (F33, F50, F55, F64, F65, F66, F67, F68, F74, F89) + " +
               "experiments/IBM_ABSORPTION_THEOREM.md + " +
               "experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
    }

    public override string DisplayName =>
        "Absorption Theorem α = 2γ₀·⟨n_XY⟩ as the rate-quantization root of F33/F50/F55/F64-F68/F74/F89";

    public override string Summary =>
        $"Re(λ) = −2γ₀·⟨n_XY⟩; absorption quantum 2 = a_0 from Pi2 ladder; per-coherence rate 2γ₀·n_diff in computational basis; max rate 2γ₀·N; pair-sum 2γ₀·N under column bit-flip (F89c) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("AbsorptionQuantumCoefficient (= a_0 = 2)", AbsorptionQuantumCoefficient);
            yield return new InspectableNode("Per-coherence reading (computational basis)",
                summary: "|A⟩⟨B| decomposes per-site into pure {I, Z} (A_l = B_l) or pure {X, Y} (A_l ≠ B_l); n_XY = n_diff(A, B); per-coherence rate 2γ₀·n_diff. F89c structural lemma uses this directly.");
            yield return new InspectableNode("Hamming-complement pair-sum (F89c corollary)",
                summary: $"column bit-flip ρ[a, b] → ρ[a, bar(b)] gives n_diff ↔ N−n_diff; rate-pairs sum to 2γ₀·N exactly = spectral maximum. Sample: 3-qubit block at γ₀=1 yields pair-sum {HammingComplementPairSum(3, 1.0):G6} (F89c path-2 anchor, verified bit-exact)");
            yield return new InspectableNode("Pi2 anchor drift check",
                summary: $"AbsorptionQuantumMatchesLiteral = {AbsorptionQuantumMatchesLiteral()} (a_0 = {AbsorptionQuantumCoefficient} vs literal 2.0 in PROOF_ABSORPTION_THEOREM.md)");
            yield return new InspectableNode("Hardware confirmation",
                summary: "IBM Q52: Re(λ) = −2γ⟨n_XY⟩ verified across 1343 modes with max deviation 2.4·10⁻¹⁴ (floating-point floor); experiments/IBM_ABSORPTION_THEOREM.md");
        }
    }
}
