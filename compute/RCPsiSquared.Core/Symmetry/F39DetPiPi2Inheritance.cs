using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F39 closed form (Tier 1 proven, verified N=1..4):
///
/// <code>
///   det(Π) = (−1)^{N · 4^{N−1}}
///
///   N=1: det = −1.   N ≥ 2: det = +1   (since 4^{N−1} is even for N ≥ 2)
/// </code>
///
/// <para>Each piece of the exponent <c>N · 4^{N−1}</c> sits on the Pi2-Foundation:</para>
///
/// <list type="bullet">
///   <item><b>"4^{N−1}" factor</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(<c>3−2N</c>)
///         = <c>a_{3−2N}</c> on the dyadic halving ladder. Equivalently, this is
///         d² for (N−1) qubits via <see cref="Pi2OperatorSpaceMirrorClaim"/> — the
///         operator-space-side anchor at qubit count (N−1). Same qubit-shift
///         pattern as F1's T1-amplitude-damping
///         (<see cref="F1T1AmplitudeDampingPi2Inheritance"/>) which uses 4^(N−1)
///         in its T1-part prefactor.</item>
///   <item><b>"(−1)" base</b>: the Π² eigenvalue convention (Z₄ memory loop generator
///         squared); see <see cref="Pi2I4MemoryLoopClaim"/>.</item>
///   <item><b>Parity result</b>: for N ≥ 2, the exponent <c>N · 4^{N−1}</c> is always
///         even (since <c>4^{N−1}</c> is even for N ≥ 2). The N=1 case gives
///         <c>4^0 = 1</c>, exponent <c>1 · 1 = 1</c>, det = −1. The N ≥ 2 result
///         det = +1 follows directly from the Pi2 ladder anchor's even-ness, NOT
///         from a separate computation.</item>
/// </list>
///
/// <para>Tier1Derived: F39 is Tier 1 proven in PT_SYMMETRY_ANALYSIS, verified
/// numerically N=1..4. The Pi2-Foundation anchoring is algebraic-trivial
/// composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F39 +
/// <c>experiments/PT_SYMMETRY_ANALYSIS.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>.</para></summary>
public sealed class F39DetPiPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2OperatorSpaceMirrorClaim _mirror;

    /// <summary>The "4^{N−1}" factor in F39's exponent. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(<c>3−2N</c>) on the dyadic halving
    /// ladder. At N=1: 1; at N=2: 4 (= d² for 1 qubit); at N=3: 16; etc.</summary>
    public double PowerNMinus1Factor(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F39 requires N ≥ 1.");
        return _ladder.Term(LadderIndexFor(N));
    }

    /// <summary>The Pi2 ladder index where the F39 power factor lands: <c>3 − 2N</c>.
    /// At N=1 → 1 (a_1 = 1, the trivial identity); at N=2 → −1 (a_{−1} = 4); at
    /// N=3 → −3 (a_{−3} = 16); etc.</summary>
    public int LadderIndexFor(int N) => 3 - 2 * N;

    /// <summary>The qubit count whose operator-space dimension equals the F39 power
    /// factor: <c>N − 1</c>. At N=2: qubit count 1 (d² = 4); at N=3: qubit count 2
    /// (d² = 16); etc. Same qubit-shift as F1-T1's T1-part prefactor.</summary>
    public int OperatorSpaceQubitCountFor(int N) => N - 1;

    /// <summary>The full exponent value <c>N · 4^{N−1}</c>. Always even for N ≥ 2
    /// (since 4^{N−1} is even); equal to 1 at N=1.</summary>
    public long ExponentValue(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F39 requires N ≥ 1.");
        long power = 1;
        for (int i = 0; i < N - 1; i++) power *= 4;
        return (long)N * power;
    }

    /// <summary>Live drift check: <c>det(Π) = (−1)^{ExponentValue(N)}</c>. Returns
    /// −1 at N=1 (exponent = 1, odd), +1 at N ≥ 2 (exponent always even).</summary>
    public int DetPi(int N) => (ExponentValue(N) % 2 == 0) ? 1 : -1;

    /// <summary>True iff the exponent <c>N · 4^{N−1}</c> is even — equivalently iff
    /// det(Π) = +1 — for the given N. Always true for N ≥ 2 by the dyadic ladder
    /// anchor's even-ness; false for N=1.</summary>
    public bool ExponentIsEven(int N) => ExponentValue(N) % 2 == 0;

    public F39DetPiPi2Inheritance(Pi2DyadicLadderClaim ladder, Pi2OperatorSpaceMirrorClaim mirror)
        : base("F39 det(Π) = (−1)^(N·4^(N−1)) inherits from Pi2-Foundation: 4^(N−1) = a_{3−2N} = d² for (N−1) qubits",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F39 + " +
               "experiments/PT_SYMMETRY_ANALYSIS.md + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
    }

    public override string DisplayName =>
        "F39 det(Π) factor as Pi2-Foundation operator-space inheritance";

    public override string Summary =>
        $"det(Π) = (−1)^(N·4^(N−1)): the 4^(N−1) factor = a_{{3−2N}} = d² for (N−1) qubits; " +
        $"N=1 → −1, N≥2 → +1 (exponent even by Pi2 ladder anchor) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F39 closed form",
                summary: "det(Π) = (−1)^(N·4^(N−1)); Tier 1 proven + verified N=1..4 (PT_SYMMETRY_ANALYSIS)");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "4^(N−1) = a_{3−2N} on the dyadic halving ladder = d² for (N−1) qubits via OperatorSpaceMirror");
            yield return new InspectableNode("parity reading",
                summary: "for N ≥ 2: 4^(N−1) is even → exponent N·4^(N−1) is even → det = +1 (algebraically forced by ladder anchor); for N=1: 4^0=1, exponent=1, det=−1");
            yield return new InspectableNode("qubit-shift pattern",
                summary: "same (N−1)-qubit shift as F1-T1's T1-part prefactor 4^(N−1); structural Pi2 footprint");
            for (int N = 1; N <= 5; N++)
            {
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"4^(N−1) = {PowerNMinus1Factor(N)} (= a_{LadderIndexFor(N)}); " +
                             $"exponent = N·4^(N−1) = {ExponentValue(N)} (even: {ExponentIsEven(N)}); " +
                             $"det(Π) = {DetPi(N)}");
            }
        }
    }
}
