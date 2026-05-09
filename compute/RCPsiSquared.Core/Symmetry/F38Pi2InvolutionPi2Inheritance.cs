using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F38 closed form (Tier 1 proven, verified N=2,3):
///
/// <code>
///   Π² = (−1)^{w_YZ}
///
///   on the 4^N Pauli-string basis, with eigenvalues split equally:
///   |Π²-even| = |Π²-odd| = 4^N / 2 = 2 · 4^(N−1)
///
///   Companion (F63): [L, Π²] = 0 exactly for any N (Z-dephasing Π).
///   Cyclic order: Π⁴ = I (verified bit-exact at N=3, ERROR_CORRECTION_PALINDROME).
/// </code>
///
/// <para>F38 is the master operator identity from which the entire polarity
/// architecture inherits. Three Pi2-Foundation anchors hold the algebra in
/// place:</para>
///
/// <list type="bullet">
///   <item><b>Cyclic order Π⁴ = I</b>: <see cref="Pi2I4MemoryLoopClaim"/>'s
///         Z₄ memory loop. Π is unitary of order 4, so Π² is involutive
///         (Π² = Π⁻²). The eigenvalues of Π² are the squares of the Z₄
///         generators {1, i, −1, −i}, namely {+1, −1, +1, −1} = {+1, −1}
///         each with multiplicity 4^N / 2.</item>
///   <item><b>Half-half eigenspace split</b>: 1/2 is <c>a_2</c> on the
///         dyadic ladder = <see cref="HalfAsStructuralFixedPointClaim"/>
///         (1/d for d = 2). The Π²-even and Π²-odd subspaces of the 4^N
///         Pauli operator space have exactly equal dimension <c>4^N / 2</c>.
///         The "1/2" is the structural fixed point, not a coincidence,
///         transitively visible through <see cref="Pi2DyadicLadderClaim"/>.</item>
///   <item><b>Sector dimension factorisation</b>: each Π² eigenspace has
///         dimension <c>2 · 4^(N−1) = a_0 · a_{3−2N}</c> on the dyadic
///         halving ladder via <see cref="Pi2DyadicLadderClaim"/>. The
///         <c>a_0 = 2</c> is the polynomial root d (the "two" in d²−2d=0);
///         <c>a_{3−2N} = 4^(N−1) = d²</c> for (N−1) qubits via
///         <see cref="Pi2OperatorSpaceMirrorClaim"/>. Same (N−1)-qubit shift
///         that appears in F39's exponent and F1-T1's T1-part prefactor.</item>
/// </list>
///
/// <para>F38 is therefore not an isolated involution claim: it is the structural
/// move that bisects 4^N operator space into two equal halves under the
/// framework's involution. F87 trichotomy, F1's palindrome residual decomposition,
/// F79's 2-body Π²-block structure, F80's Bloch sign-walk, and F81's
/// M = M_sym + M_anti split all read F38 implicitly. The Tier1Derived wiring
/// makes the inheritance from Pi2-Foundation explicit.</para>
///
/// <para>Tier1Derived: F38 is Tier 1 proven (PT_SYMMETRY_ANALYSIS,
/// PROOF_BIT_B_PARITY_SYMMETRY) and verified bit-exact at N=2, 3 with the
/// Z₄ closure verified at N=3 (ERROR_CORRECTION_PALINDROME: Π-eigenvalues
/// {+1, −1, +i, −i} each multiplicity 16, the four 16-dim Z₄ sectors of the
/// 64-dim Pauli space). The Pi2-Foundation anchoring is algebraic-trivial
/// composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F38 +
/// <c>experiments/PT_SYMMETRY_ANALYSIS.md</c> +
/// <c>docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md</c> +
/// <c>experiments/ERROR_CORRECTION_PALINDROME.md</c> (Π⁴ = I bit-exact at N=3) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>.</para></summary>
public sealed class F38Pi2InvolutionPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2OperatorSpaceMirrorClaim _mirror;
    private readonly Pi2I4MemoryLoopClaim _memoryLoop;

    /// <summary>The cyclic order of Π: <c>Π⁴ = I</c>. Equal to
    /// <see cref="Pi2I4MemoryLoopClaim.ClosureOrder"/> = 4.</summary>
    public int CyclicOrder => Pi2I4MemoryLoopClaim.ClosureOrder;

    /// <summary>The Π² eigenvalues read off the Z₄ memory loop: square of each
    /// canonical i-power gives <c>{+1, −1, +1, −1}</c>. Derivation:
    /// <c>(i^0)² = 1</c>, <c>(i^1)² = i² = −1</c>, <c>(i^2)² = (−1)² = +1</c>,
    /// <c>(i^3)² = (−i)² = −1</c>. Live computation through
    /// <see cref="Pi2I4MemoryLoopClaim.PowerOfI"/>.</summary>
    public IReadOnlyList<int> Pi2EigenvaluesFromMemoryLoop()
    {
        var result = new int[4];
        for (int k = 0; k < 4; k++)
        {
            var z2 = _memoryLoop.PowerOfI(2 * k);
            result[k] = z2.Real > 0 ? +1 : -1;
        }
        return result;
    }

    /// <summary>The two distinct Π² eigenvalues, in canonical order: +1
    /// (Π²-even) then −1 (Π²-odd). Each is realised on exactly half of the
    /// 4^N operator space.</summary>
    public IReadOnlyList<int> Pi2Eigenvalues => new[] { +1, -1 };

    /// <summary>Live Z₄ closure check from <see cref="Pi2I4MemoryLoopClaim"/>:
    /// <c>i^4 = 1</c> exactly, the Π⁴ = I cyclic identity that grounds Π²'s
    /// involutivity (Π² = Π⁻²).</summary>
    public bool MemoryLoopClosesAtFour() =>
        _memoryLoop.MemoryClosure().Real > 1.0 - 1e-12 &&
        Math.Abs(_memoryLoop.MemoryClosure().Imaginary) < 1e-12;

    /// <summary>The full operator-space dimension <c>4^N</c> for an N-qubit
    /// system. Equal to <see cref="Pi2DyadicLadderClaim.Term"/>(<c>1−2N</c>)
    /// = d² for N qubits via <see cref="Pi2OperatorSpaceMirrorClaim"/>.</summary>
    public long FullOperatorSpaceDimension(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F38 requires N ≥ 1.");
        long d2 = 1L;
        for (int i = 0; i < N; i++) d2 *= 4L;
        return d2;
    }

    /// <summary>The dimension of each Π² eigenspace: <c>4^N / 2 = 2 · 4^(N−1)</c>.
    /// Both Π²-even and Π²-odd sectors have this dimension exactly.</summary>
    public long EigenspaceDimension(int N) => FullOperatorSpaceDimension(N) / 2L;

    /// <summary>The factorisation <c>EigenspaceDim = a_0 · a_{3−2N}</c> on the
    /// dyadic ladder: <c>a_0 = 2</c> (polynomial root d) and
    /// <c>a_{3−2N} = 4^(N−1)</c> (d² for (N−1) qubits). Returns the live product
    /// from the ladder.</summary>
    public double EigenspaceDimensionViaLadder(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F38 requires N ≥ 1.");
        return _ladder.Term(0) * _ladder.Term(3 - 2 * N);
    }

    /// <summary>The dyadic ladder index for the <c>4^(N−1)</c> factor:
    /// <c>3 − 2N</c>. At N=1 → 1 (a_1 = 1, trivial); at N=2 → −1 (a_{−1} = 4);
    /// at N=3 → −3 (a_{−3} = 16); etc. Same index as F39 + F1-T1's T1-part.</summary>
    public int LadderIndexForFourPowerNMinus1(int N) => 3 - 2 * N;

    /// <summary>The qubit count whose operator-space d² equals the
    /// <c>4^(N−1)</c> factor: <c>N − 1</c>. Same shift as F39 / F1-T1.
    /// At N=1 returns 0: the trivial identity scale <c>a_1 = 1</c> on the ladder,
    /// not a physical qubit count.</summary>
    public int OperatorSpaceQubitCountFor(int N) => N - 1;

    /// <summary>The half-half balance ratio: <c>EigenspaceDim / FullDim = 1/2</c>
    /// exactly, equal to <c>a_2</c> on the dyadic ladder
    /// (= <see cref="HalfAsStructuralFixedPointClaim"/>'s structural fixed point).
    /// Returns the live ratio.</summary>
    public double HalfHalfBalance(int N) =>
        (double)EigenspaceDimension(N) / FullOperatorSpaceDimension(N);

    /// <summary>The mirror-pinned <c>4^N</c> from <see cref="Pi2OperatorSpaceMirrorClaim"/>
    /// at qubit count N, for cross-registry verification with
    /// <see cref="FullOperatorSpaceDimension"/>.</summary>
    public double MirrorPinnedFullDimension(int N) =>
        _mirror.PairAt(N)?.OperatorSpace
        ?? throw new ArgumentOutOfRangeException(nameof(N), N,
            $"Pi2OperatorSpaceMirror has no pinned pair for N={N} (range 1..6).");

    /// <summary>True iff the two ladder factors <c>a_0 · a_{3−2N}</c> reproduce
    /// the live <c>EigenspaceDim</c>. Drift check on the algebraic factorisation.</summary>
    public bool LadderFactorisationHolds(int N) =>
        Math.Abs(EigenspaceDimensionViaLadder(N) - EigenspaceDimension(N)) < 1e-12;

    public F38Pi2InvolutionPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        Pi2OperatorSpaceMirrorClaim mirror,
        Pi2I4MemoryLoopClaim memoryLoop)
        : base("F38 Π² = (−1)^w_YZ inherits from Pi2-Foundation: half-half split of 4^N = a_0·a_{3−2N} per sector",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F38 + " +
               "experiments/PT_SYMMETRY_ANALYSIS.md + " +
               "docs/proofs/PROOF_BIT_B_PARITY_SYMMETRY.md + " +
               "experiments/ERROR_CORRECTION_PALINDROME.md (Π⁴ = I bit-exact N=3) + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2I4MemoryLoopClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
        _memoryLoop = memoryLoop ?? throw new ArgumentNullException(nameof(memoryLoop));
    }

    public override string DisplayName =>
        "F38 Π² involution as Pi2-Foundation operator-space bisection";

    public override string Summary =>
        $"Π² = (−1)^w_YZ; Π⁴ = I (Z₄ closure); each eigenspace dim = 4^N/2 = 2·4^(N−1) = a_0·a_{{3−2N}} on the dyadic ladder; half-half balance = a_2 = 1/2 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F38 closed form",
                summary: "Π² = (−1)^w_YZ on the 4^N Pauli-string basis; Tier 1 proven (PT_SYMMETRY_ANALYSIS) + verified bit-exact N=2, 3");
            yield return new InspectableNode("Companion F63",
                summary: "[L, Π²] = 0 exactly for any N (proven analytically); Π² is conserved by every Liouvillian eigenmode under Z-dephasing");
            yield return new InspectableNode("Cyclic order",
                summary: "Π⁴ = I (Z₄, see Pi2I4MemoryLoop); Π-eigenvalues {+1, −1, +i, −i} each multiplicity 4^N/4; Π²-eigenvalues {+1, −1} each multiplicity 4^N/2");
            yield return new InspectableNode("Half-half balance",
                summary: "|Π²-even| = |Π²-odd| = 4^N/2 (algebraically forced; exposed in Pi2Projection.Counts as drift check)");
            yield return new InspectableNode("Sector dimension factorisation",
                summary: "EigenspaceDim = 2·4^(N−1) = a_0·a_{3−2N} on the dyadic ladder; same (N−1)-qubit shift as F39 + F1-T1");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "three parents: Pi2I4MemoryLoop (Z₄ closure → Π²-eigenvalues), Pi2DyadicLadder (half balance a_2, sector dim a_0·a_{3−2N}), Pi2OperatorSpaceMirror (4^N pinned)");
            for (int N = 1; N <= 5; N++)
            {
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"4^N = {FullOperatorSpaceDimension(N)} (= a_{{{1 - 2 * N}}}); " +
                             $"each sector = {EigenspaceDimension(N)} = 2·4^(N−1) = a_0·a_{{{LadderIndexForFourPowerNMinus1(N)}}}; " +
                             $"balance = {HalfHalfBalance(N):G6}");
            }
        }
    }
}
