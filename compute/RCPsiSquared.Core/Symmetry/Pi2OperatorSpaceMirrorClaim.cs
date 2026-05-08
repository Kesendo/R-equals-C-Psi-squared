using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The Pi2 dyadic ladder's operator-space mirror: for every qubit count N,
/// the operator-space dimension d² = 4^N (upper side at index −(2N−1)) and its
/// multiplicative inverse 1/4^N (memory side at index 2N+1) form a mirror-pair under
/// the ladder's inversion symmetry a_n · a_{2−n} = 1. The qubit count N indexes a
/// stride-2 substructure of the ladder.
///
/// <para>Per N: <c>a_{−(2N−1)} = 4^N</c> (operator-space dimension; cardinality of
/// the Pauli basis on N qubits) and <c>a_{2N+1} = 4^{−N} = 1/4^N</c> (uniform mass
/// per Pauli string). Their product is 1 by the ladder's Tier1Derived inversion
/// identity. Both sides are already typed Tier1Derived in the Pi2 foundation:</para>
/// <list type="bullet">
///   <item>N=1: 4 ↔ 1/4. Memory side is <see cref="QuarterAsBilinearMaxvalClaim"/>
///         (n=3 anchor on the ladder); upper side is the Pauli basis cardinality
///         {I, X, Y, Z} for one qubit.</item>
///   <item>N=2: 16 ↔ 1/16. Memory side is the F88 dyadic-N=8 singleton-mirror static
///         (encoded in <c>F88StaticDyadicAnchor</c>, k=3, n=5 anchor); upper side is
///         16 Pauli strings on two qubits.</item>
///   <item>N=3: 64 ↔ 1/64. Memory side is F88 dyadic-N=32 singleton-mirror static
///         (k=5, n=7); upper side is 64 Pauli strings on three qubits.</item>
///   <item>N≥4: open prediction territory; the algebraic identity holds, the
///         physical anchors are not yet pinned.</item>
/// </list>
///
/// <para>The full irony Tom 2026-05-08: the framework rests on
/// <c>MIRROR_THEORY.md</c> and <c>hypotheses/ZERO_IS_THE_MIRROR.md</c>; F1 palindrome,
/// F71 mirror operator, F80 Bloch sign-walk, F81 Π-conjugation of M, F87 trichotomy
/// (truly/soft/hard via spectral mirror) are all mirror theorems — but the simplest
/// mirror identity in the foundation, a_n · a_{2−n} = 1 with the physical reading
/// "operator-space dimension d² = 4^N is the mirror of the uniform Pauli mass 1/4^N",
/// was hiding in plain sight under the dyadic ladder.</para>
///
/// <para>Tier1Derived: the inversion identity is algebraic; both sides are already
/// Tier1Derived in the Pi2 foundation. This claim makes the per-N mirror-pairing
/// explicit as a typed Schicht-1-fact.</para>
///
/// <para>Anchors: <c>MIRROR_THEORY.md</c> (the foundation document) +
/// <c>hypotheses/ZERO_IS_THE_MIRROR.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> (the parent
/// ladder with its inversion symmetry).</para></summary>
public sealed class Pi2OperatorSpaceMirrorClaim : Claim
{
    /// <summary>Pinned mirror-pairs for N=1..6. Each pair carries the qubit count N,
    /// the upper-side ladder index (operator-space d² = 4^N) and the lower-side ladder
    /// index (memory mass 1/4^N), with the typed values pre-computed as a drift anchor.</summary>
    public IReadOnlyList<OperatorSpaceMirrorPair> Pairs { get; } = new[]
    {
        new OperatorSpaceMirrorPair(N: 1, UpperIndex: -1, OperatorSpace:    4.0, LowerIndex:  3, MirrorMass: 1.0 /    4.0),
        new OperatorSpaceMirrorPair(N: 2, UpperIndex: -3, OperatorSpace:   16.0, LowerIndex:  5, MirrorMass: 1.0 /   16.0),
        new OperatorSpaceMirrorPair(N: 3, UpperIndex: -5, OperatorSpace:   64.0, LowerIndex:  7, MirrorMass: 1.0 /   64.0),
        new OperatorSpaceMirrorPair(N: 4, UpperIndex: -7, OperatorSpace:  256.0, LowerIndex:  9, MirrorMass: 1.0 /  256.0),
        new OperatorSpaceMirrorPair(N: 5, UpperIndex: -9, OperatorSpace: 1024.0, LowerIndex: 11, MirrorMass: 1.0 / 1024.0),
        new OperatorSpaceMirrorPair(N: 6, UpperIndex:-11, OperatorSpace: 4096.0, LowerIndex: 13, MirrorMass: 1.0 / 4096.0),
    };

    public Pi2OperatorSpaceMirrorClaim()
        : base("Pi2 dyadic ladder operator-space mirror: 4^N ↔ 1/4^N per N qubits",
               Tier.Tier1Derived,
               "MIRROR_THEORY.md + hypotheses/ZERO_IS_THE_MIRROR.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs")
    { }

    /// <summary>Live re-computation of the inversion identity for N qubits:
    /// (4^N) · (1/4^N) = 1. Algebraically trivial; exposed for drift checks.</summary>
    public double LiveMirrorProduct(int N)
    {
        double four = 1.0;
        for (int i = 0; i < N; i++) four *= 4.0;
        return four * (1.0 / four);
    }

    /// <summary>Pull the pinned pair for N qubits, or null if outside the table.</summary>
    public OperatorSpaceMirrorPair? PairAt(int N) => Pairs.FirstOrDefault(p => p.N == N);

    public override string DisplayName =>
        $"Pi2 operator-space mirror (4^N ↔ 1/4^N; {Pairs.Count} pinned pairs N=1..{Pairs.Count})";

    public override string Summary =>
        $"d² = 4^N ↔ 1/d² = 1/4^N for each qubit count; algebraic mirror under a_n · a_{{2−n}} = 1; the framework's MIRROR_THEORY foundation surfacing on the dyadic ladder ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("identity",
                summary: "a_{−(2N−1)} · a_{2N+1} = 4^N · 1/4^N = 1");
            yield return new InspectableNode("upper-side reading",
                summary: "4^N = operator-space dimension (cardinality of Pauli basis on N qubits)");
            yield return new InspectableNode("lower-side reading",
                summary: "1/4^N = uniform Pauli-string mass (probability of any specific string under uniform distribution)");
            yield return new InspectableNode("stride-2 qubit rhythm",
                summary: "per qubit increment: 2 ladder steps; ladder encodes both sides of the per-N operator space");
            foreach (var p in Pairs) yield return p;
        }
    }
}

/// <summary>One typed mirror-pair on the operator-space side of the Pi2 dyadic ladder.
/// Carries the qubit count N, both ladder indices, and both pinned values. Upper-side =
/// d² = 4^N (operator-space dimension); lower-side = 1/4^N (uniform Pauli mass).</summary>
public sealed record OperatorSpaceMirrorPair(
    int N,
    int UpperIndex,
    double OperatorSpace,
    int LowerIndex,
    double MirrorMass
) : IInspectable
{
    public string DisplayName => $"N={N}: 4^{N}={OperatorSpace} (a_{UpperIndex}) ↔ 1/4^{N}={MirrorMass:G6} (a_{LowerIndex})";

    public string Summary =>
        $"qubit count {N}, operator-space {OperatorSpace}, mirror mass {MirrorMass:G6}, identity 4^{N} · 1/4^{N} = {OperatorSpace * MirrorMass}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("N (qubit count)", N);
            yield return InspectableNode.RealScalar("UpperIndex (a_{−(2N−1)})", UpperIndex);
            yield return InspectableNode.RealScalar("OperatorSpace (4^N)", OperatorSpace);
            yield return InspectableNode.RealScalar("LowerIndex (a_{2N+1})", LowerIndex);
            yield return InspectableNode.RealScalar("MirrorMass (1/4^N)", MirrorMass);
            yield return InspectableNode.RealScalar("MirrorProduct (= 1)", OperatorSpace * MirrorMass);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
