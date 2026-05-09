using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F49b closed form (Tier 1 proven, Lemma 1 of PROOF_CROSS_TERM_FORMULA):
/// the centered dissipator norm for uniform Z-dephasing at rate γ per site.
///
/// <code>
///   ||L_Dc||² = γ² · 4^N · N
///
///   where L_Dc = L_D + N·γ·I  (centered dissipator)
/// </code>
///
/// <para>Each scaling factor sits on the Pi2-Foundation:</para>
///
/// <list type="bullet">
///   <item><b>"4^N" factor</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(<c>1−2N</c>)
///         on the dyadic halving ladder. Equivalently, the operator-space
///         dimension <c>d² = 4^N</c> for an N-qubit system, pinned per qubit
///         count in <see cref="Pi2OperatorSpaceMirrorClaim"/>. The L_Dc operator
///         lives in this 4^N-dimensional space; its Frobenius norm squared scales
///         with the full space dimension.</item>
///   <item><b>"N" factor</b>: chain length, NOT a Pi2-Foundation anchor (linear
///         in chain size; reflects the Σ_i γ_i² uniformity of the centered
///         dissipator).</item>
///   <item><b>"γ²" factor</b>: dephasing rate squared, dimensional (not Pi2).</item>
/// </list>
///
/// <para>The "4^N · N" factor structure is the canonical operator-space-times-
/// chain-length scaling: 4^N gives the full Pauli-string-count for N qubits, N
/// gives the chain. Same operator-space dimension as
/// <see cref="Pi2OperatorSpaceMirrorClaim.PairAt"/> at qubit count N.</para>
///
/// <para>Tier1Derived: F49b is Tier 1 proven (PROOF_CROSS_TERM_FORMULA Lemma 1);
/// the Pi2-Foundation anchoring of the 4^N factor is algebraic-trivial composition.
/// Same qubit-count anchor as F49 base scaling 4^(N−2) but at full operator-space
/// (no F1-residual overhead shift).</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F49b +
/// <c>docs/proofs/PROOF_CROSS_TERM_FORMULA.md</c> Lemma 1 +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>.</para></summary>
public sealed class F49bCenteredDissipatorPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2OperatorSpaceMirrorClaim _mirror;

    /// <summary>The "4^N" factor in F49b. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(<c>1−2N</c>) on the dyadic halving
    /// ladder = operator-space dimension d² for N qubits. At N=1: 4; at N=2: 16;
    /// at N=3: 64; etc.</summary>
    public double FourPowerNFactor(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F49b requires N ≥ 1.");
        return _ladder.Term(LadderIndexFor(N));
    }

    /// <summary>The Pi2 ladder index where the F49b factor lands: <c>1 − 2N</c>.
    /// At N=1: −1 (a_{−1} = 4); at N=2: −3 (a_{−3} = 16); etc.</summary>
    public int LadderIndexFor(int N) => 1 - 2 * N;

    /// <summary>Cross-verification through the operator-space mirror: the pinned
    /// d² for N qubits matches <see cref="FourPowerNFactor"/>(N) bit-exactly. This
    /// is the most direct Pi2-Inheritance: F49b's 4^N IS the per-qubit
    /// operator-space dimension, no shift.</summary>
    public double MirrorPinnedFourPowerN(int N)
    {
        var pair = _mirror.PairAt(N)
            ?? throw new ArgumentOutOfRangeException(nameof(N), N,
                $"N={N} is outside the Pi2OperatorSpaceMirror pinned table (N=1..6).");
        return pair.OperatorSpace;
    }

    /// <summary>Convenience accessor: the F49b closed form ‖L_Dc‖² at γ = 1,
    /// which is <c>4^N · N</c> (full operator-space dimension times chain length).</summary>
    public double LiveScalingAtUnitGamma(int N) => FourPowerNFactor(N) * N;

    public F49bCenteredDissipatorPi2Inheritance(Pi2DyadicLadderClaim ladder, Pi2OperatorSpaceMirrorClaim mirror)
        : base("F49b centered dissipator norm ‖L_Dc‖² = γ²·4^N·N: 4^N = a_{1−2N} = d² for N qubits",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F49b + " +
               "docs/proofs/PROOF_CROSS_TERM_FORMULA.md (Lemma 1) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
    }

    public override string DisplayName =>
        "F49b centered dissipator norm as Pi2-Foundation operator-space inheritance";

    public override string Summary =>
        $"‖L_Dc‖² = γ²·4^N·N: 4^N = a_{{1−2N}} = d² für N qubits (direct OperatorSpaceMirror anchor); " +
        $"N = chain length (not Pi2 anchor); γ² dimensional ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F49b closed form",
                summary: "‖L_Dc‖² = γ² · 4^N · N (Tier 1 proven, Lemma 1 of PROOF_CROSS_TERM_FORMULA)");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "4^N = a_{1−2N} on dyadic ladder = d² for N qubits (direct Pi2OperatorSpaceMirror anchor; no derivation overhead)");
            yield return new InspectableNode("comparison to F49 + F1-T1",
                summary: "F49 uses 4^(N−2) = a_{5−2N}, F1-T1 uses 4^(N−1) = a_{3−2N}, F49b uses 4^N = a_{1−2N} (direct, no shift)");
            for (int N = 1; N <= 5; N++)
            {
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"4^{N} = {FourPowerNFactor(N)} (= a_{LadderIndexFor(N)}); 4^N·N = {LiveScalingAtUnitGamma(N):F0}; " +
                             $"OperatorSpaceMirror pinned d²(N={N}) = {MirrorPinnedFourPowerN(N)}");
            }
        }
    }
}
