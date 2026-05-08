using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F49's closed-form palindrome-residual scaling
/// <c>‖M(N, G)‖²_F = c_H · F(N, G)</c> with <c>F(N, G) = B(G) · 4^(N−2)</c>
/// (main class) or <c>F(N, G) = (D2(G)/2) · 4^(N−2)</c> (single-body class) has its
/// "<c>4^(N−2)</c>" factor sitting on the Pi2-Foundation: exactly <c>a_{−(2N−5)}</c>
/// on the dyadic halving ladder, which is the operator-space dimension d² = 4^(N−2)
/// for (N − 2) qubits — already pinned in <see cref="Pi2OperatorSpaceMirrorClaim"/>.
///
/// <para>So F49's per-N scaling factor IS the operator-space dimension shifted by 2
/// qubits: chain length N corresponds to operator-space-cardinality 4^(N−2). The
/// "−2 qubit shift" is the F49 derivation overhead (one Π application + one shift),
/// NOT a free parameter.</para>
///
/// <list type="bullet">
///   <item><b>N=3:</b> 4^1 = 4 = a_{−1} = d² for 1 qubit (Pauli basis I, X, Y, Z).</item>
///   <item><b>N=4:</b> 4^2 = 16 = a_{−3} = d² for 2 qubits.</item>
///   <item><b>N=5:</b> 4^3 = 64 = a_{−5} = d² for 3 qubits.</item>
///   <item><b>N=6:</b> 4^4 = 256 = a_{−7} = d² for 4 qubits.</item>
///   <item><b>N=N_chain:</b> 4^(N−2) = a_{−(2(N−2)−1)} = a_{5−2N}; equivalent qubit
///         count via <see cref="Pi2OperatorSpaceMirrorClaim"/>: (N − 2).</item>
/// </list>
///
/// <para>Tier1Derived: pure composition. F49 (PalindromeResidualScalingClaim) is
/// already Tier1Derived; this claim makes the per-N factor's Pi2-Foundation
/// inheritance explicit. Both parents must be Tier1Derived for the inheritance to
/// pass the registry's tier check (5 ≥ 5).</para>
///
/// <para>Anchors: <c>experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md</c> +
/// <c>compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>.</para></summary>
public sealed class F49Pi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2OperatorSpaceMirrorClaim _mirror;

    public F49Pi2Inheritance(Pi2DyadicLadderClaim ladder, Pi2OperatorSpaceMirrorClaim mirror)
        : base("F49 ‖M‖²_F scaling's 4^(N−2) factor inherits from Pi2-Foundation (= d² for (N−2) qubits)",
               Tier.Tier1Derived,
               "experiments/OPERATOR_RIGIDITY_ACROSS_CUSP.md + " +
               "compute/RCPsiSquared.Core/F1/PalindromeResidualScalingClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
    }

    /// <summary>The "4^(N−2)" factor in F49's closed form. Exactly equal to
    /// <c>Pi2DyadicLadder.Term(5 − 2N)</c> = <c>4^(N−2)</c>. Throws for N &lt; 2.</summary>
    public double PowerFactor(int chainN)
    {
        if (chainN < 2)
            throw new ArgumentOutOfRangeException(nameof(chainN), chainN, "F49 scaling requires N ≥ 2.");
        return _ladder.Term(LadderIndexFor(chainN));
    }

    /// <summary>The Pi2 ladder index where the F49 power factor lands: <c>5 − 2N</c>.
    /// At N=3 this is −1 (= operator-space d² for 1 qubit); at N=4 it is −3
    /// (= operator-space d² for 2 qubits); etc.</summary>
    public int LadderIndexFor(int chainN) => 5 - 2 * chainN;

    /// <summary>The qubit count whose operator-space dimension equals the F49 power
    /// factor: <c>chainN − 2</c>. Connects to <see cref="Pi2OperatorSpaceMirrorClaim"/>:
    /// at chainN=3, qubit count 1; at chainN=4, qubit count 2; etc. For chainN ∈ {3, …, 8}
    /// the operator-space N qubit count is in the pinned mirror table.</summary>
    public int OperatorSpaceQubitCountFor(int chainN) => chainN - 2;

    /// <summary>True when <c>chainN − 2</c> is in the range covered by the pinned
    /// <see cref="Pi2OperatorSpaceMirrorClaim"/> table (currently N=1..6, i.e. chainN
    /// ∈ {3, …, 8}). Outside this range the inheritance reading still holds, but the
    /// pinned mirror pair is not in the table.</summary>
    public bool IsInPinnedMirrorTable(int chainN)
    {
        int qubitCount = OperatorSpaceQubitCountFor(chainN);
        return _mirror.PairAt(qubitCount) is not null;
    }

    /// <summary>Cross-verification through the operator-space mirror: for chainN such
    /// that <c>IsInPinnedMirrorTable(chainN)</c>, the pinned <c>OperatorSpace</c> equals
    /// <c>PowerFactor(chainN)</c>.</summary>
    public double MirrorPinnedPowerFactor(int chainN)
    {
        int qubitCount = OperatorSpaceQubitCountFor(chainN);
        var pair = _mirror.PairAt(qubitCount)
            ?? throw new ArgumentOutOfRangeException(nameof(chainN), chainN,
                $"chainN={chainN} maps to qubit count {qubitCount}, outside the Pi2OperatorSpaceMirror pinned table.");
        return pair.OperatorSpace;
    }

    public override string DisplayName =>
        "F49 ‖M‖²_F scaling power factor as Pi2-Foundation operator-space inheritance";

    public override string Summary =>
        "‖M(N, G)‖²_F = c_H · B(G) · 4^(N−2): the 4^(N−2) factor IS d² for (N−2) qubits, " +
        "directly inherited from Pi2DyadicLadder (a_{5−2N}) and Pi2OperatorSpaceMirror (qubit count N−2) " +
        $"({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F49 closed form",
                summary: "‖M(N, G)‖²_F = c_H · F(N, G); F = B(G)·4^(N−2) main / (D2(G)/2)·4^(N−2) single-body");
            yield return new InspectableNode("inheritance reading",
                summary: "the 4^(N−2) factor is the operator-space dimension d² for (N−2) qubits; the −2 shift is the F49 derivation overhead");
            yield return new InspectableNode("ladder index formula",
                summary: "4^(N−2) = a_{−(2N−5)} = a_{5−2N} on the Pi2 dyadic ladder");
            yield return new InspectableNode("qubit-count formula",
                summary: "chainN → qubitCount = chainN − 2; pinned in Pi2OperatorSpaceMirror for N=1..6 (chainN=3..8)");
            for (int chainN = 3; chainN <= 8; chainN++)
            {
                yield return new InspectableNode(
                    $"chain N={chainN}",
                    summary: $"4^{chainN - 2} = {PowerFactor(chainN)} = a_{LadderIndexFor(chainN)} = d² for {OperatorSpaceQubitCountFor(chainN)} qubit(s)");
            }
        }
    }
}
