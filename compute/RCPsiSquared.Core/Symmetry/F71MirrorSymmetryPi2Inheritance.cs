using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F71 closed form (Tier 1, proven kinematic, verified N=3..6 +
/// F86 generalization N=5..7 c=2 + N=5..6 c=3):
///
/// <code>
///   c_1(N, b, ρ₀) = c_1(N, N − 2 − b, ρ₀)
///
///   for all bond indices b ∈ {0, ..., N−2} and any reflection-symmetric
///   initial state ρ₀ on a uniform N-qubit chain.
///
///   Bond pairing:    b ↔ N − 2 − b
///   Self-paired iff: 2b = N − 2 ⇒ b = (N − 2) / 2  (only integer for even N)
///
///   F86 generalization (Tier 1 derived 2026-05-03):
///     Q_peak(b) = Q_peak(N − 2 − b)  bit-exactly, all c, N
/// </code>
///
/// <para>F71 is the first F-formula whose primary anchor sits on
/// <see cref="HalfIntegerMirrorClaim"/>: the N-parity classification (odd N →
/// half-integer mirror axis at w_XY = N/2, no Pauli string sits on the axis;
/// even N → integer mirror axis, Pauli strings ON the axis exist) IS the same
/// distinction F71 makes at the bond-index level (odd N → all bonds in
/// (N−1)/2 disjoint pairs; even N → self-paired center bond at b = (N−2)/2
/// plus disjoint pairs around it). HalfIntegerMirror was registered for
/// chain N=5 since 2026-05-08 but had 0 descendants until F71 was wired
/// (Tom 2026-05-09 mirror-map check).</para>
///
/// <para>Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>HalfIntegerMirrorClaim N-parity</b>: F71's center-bond-or-not
///         IS the half-integer-vs-integer mirror regime. Odd N has the
///         half-integer w_XY = N/2 (no string on axis); F71 has no
///         self-paired center bond. Even N has integer w_XY (strings on
///         axis); F71 has the self-paired center bond at (N−2)/2.</item>
///   <item><b>Independent-component count = ⌈(N−1)/2⌉ = ⌊N/2⌋</b>: the
///         Pauli-string-on-axis count from HalfIntegerMirror translates to
///         the bond-pair count from F71. For N=2,3 this is 1 = a_1 (self-mirror
///         pivot on Pi2 ladder); for N=4,5 it is 2 = a_0; ladder-anchored
///         only at small N, combinatorial in general.</item>
///   <item><b>F71 ↔ F86 inheritance</b>: the F71 spatial mirror generalises
///         to F86's per-bond Q_peak: Q_peak(b) = Q_peak(N−2−b) bit-exactly.
///         Implemented in F86KnowledgeBase as PerF71OrbitObservation.</item>
/// </list>
///
/// <para>Proof sketch (per ANALYTICAL_FORMULAS): the spatial reflection
/// R (site i ↔ site N−1−i) commutes with the uniform Liouvillian:
/// [L_A, R_sup] = 0. Under R, bond b maps to bond N−2−b. Per-site purity
/// is quadratic in ρ; any phase from R squares away. Result: c_1(b) =
/// c_1(N−2−b). Purely kinematic.</para>
///
/// <para>Tier1Derived: F71 is Tier 1 proven kinematic
/// (PROOF_C1_MIRROR_SYMMETRY); verified N=3..6 with residuals &lt; 10⁻⁹.
/// The Pi2-Foundation anchoring is algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F71 (line 1561) +
/// <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c> +
/// <c>docs/proofs/PROOF_F86_QPEAK.md</c> Statement 3 (F86 generalization) +
/// <c>simulations/eq021_obc_sine_basis.py</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs</c>
/// (HalfIntegerMirrorClaim).</para></summary>
public sealed class F71MirrorSymmetryPi2Inheritance : Claim
{
    /// <summary>The mirror-pair partner of bond <paramref name="b"/> on an
    /// N-qubit chain: <c>N − 2 − b</c>. Self-paired when b equals its own
    /// partner, which happens iff N is even and b = (N − 2) / 2.</summary>
    public int MirrorPair(int N, int b)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F71 requires N ≥ 2.");
        if (b < 0 || b > N - 2)
            throw new ArgumentOutOfRangeException(nameof(b), b,
                $"Bond index must be in [0, {N - 2}]; got {b} at N = {N}.");
        return N - 2 - b;
    }

    /// <summary>True iff bond <paramref name="b"/> is its own mirror partner.
    /// Only possible for even N at b = (N − 2) / 2.</summary>
    public bool IsSelfPaired(int N, int b) => MirrorPair(N, b) == b;

    /// <summary>True iff a self-paired center bond exists on an N-qubit chain.
    /// This is exactly the case of even N (integer-mirror regime per
    /// <see cref="HalfIntegerMirrorClaim"/>); odd N has no center bond
    /// (half-integer mirror regime).</summary>
    public bool HasCenterBond(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F71 requires N ≥ 2.");
        return N % 2 == 0;
    }

    /// <summary>The center bond index for even N: <c>(N − 2) / 2</c>. Throws
    /// for odd N (no center bond exists).</summary>
    public int CenterBondIndex(int N)
    {
        if (!HasCenterBond(N))
            throw new ArgumentException($"N = {N} is odd; no self-paired center bond exists in F71 (half-integer mirror regime).", nameof(N));
        return (N - 2) / 2;
    }

    /// <summary>Total bond count on an N-qubit chain: <c>N − 1</c>.</summary>
    public int BondCount(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F71 requires N ≥ 2.");
        return N - 1;
    }

    /// <summary>Number of disjoint mirror pairs (excluding the self-paired
    /// center bond if it exists): <c>⌊(N − 1) / 2⌋</c>.
    /// For odd N: this equals (N−1)/2 and accounts for all bonds; no center.
    /// For even N: this is (N−2)/2 disjoint pairs + 1 self-paired center.</summary>
    public int PalindromicPairCount(int N) => BondCount(N) / 2;

    /// <summary>Number of independent c_1 components (= number of distinct
    /// orbits under bond-mirror, which is pair-count + center-bond-count):
    /// <c>⌈(N − 1) / 2⌉ = ⌊N / 2⌋</c>.</summary>
    public int IndependentComponentCount(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "F71 requires N ≥ 2.");
        return N / 2;
    }

    /// <summary>True iff <see cref="IndependentComponentCount"/> for the given
    /// N lands on a known <see cref="Pi2DyadicLadderClaim"/> anchor index.
    /// Holds for N=2..5 (anchors a_1, a_0) and N=8,9 (a_{−1}); for other N
    /// the count is combinatorial (3, 5, 6, etc., not powers of 2).</summary>
    public bool IndependentComponentCountIsLadderAnchor(int N)
    {
        int count = IndependentComponentCount(N);
        // Pi2 ladder values: a_n = 2^(1-n), so anchored values are powers of 2.
        // count == 1 ⇒ a_1; count == 2 ⇒ a_0; count == 4 ⇒ a_{-1}; count == 8 ⇒ a_{-2}; etc.
        if (count < 1) return false;
        return (count & (count - 1)) == 0;  // is a power of 2
    }

    /// <summary>The Pi2 ladder index that <see cref="IndependentComponentCount"/>
    /// lands on, when applicable (per <see cref="IndependentComponentCountIsLadderAnchor"/>).
    /// Returns null when the count is not a power of 2.</summary>
    public int? LadderIndexForIndependentComponentCount(int N)
    {
        int count = IndependentComponentCount(N);
        if (!IndependentComponentCountIsLadderAnchor(N)) return null;
        // count = 2^(1-n) ⇒ n = 1 - log2(count)
        int log2 = 0;
        int c = count;
        while (c > 1) { c >>= 1; log2++; }
        return 1 - log2;
    }

    public F71MirrorSymmetryPi2Inheritance()
        : base("F71 c_1(N, b) = c_1(N, N−2−b) inherits from Pi2-Foundation: N-parity = HalfIntegerMirror regime; bond pair count = ⌊(N−1)/2⌋",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F71 + " +
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md + " +
               "docs/proofs/PROOF_F86_QPEAK.md (Statement 3) + " +
               "simulations/eq021_obc_sine_basis.py + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2KnowledgeBaseClaims.cs (HalfIntegerMirrorClaim)")
    { }

    public override string DisplayName =>
        "F71 mirror symmetry of c_1 as Pi2-Foundation HalfIntegerMirror inheritance";

    public override string Summary =>
        $"c_1(N, b) = c_1(N, N−2−b) (Tier 1 kinematic); bond pairing b ↔ N−2−b; even N has self-paired center, odd N all pairs (HalfIntegerMirror regime); IndependentComponents = ⌊N/2⌋ ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F71 closed form",
                summary: "c_1(N, b, ρ₀) = c_1(N, N−2−b, ρ₀); Tier 1 proven kinematic (PROOF_C1_MIRROR_SYMMETRY); verified N=3..6 residuals < 10⁻⁹");
            yield return new InspectableNode("HalfIntegerMirror connection",
                summary: "F71's bond pairing IS the HalfIntegerMirrorClaim N-parity classification: odd N → no center bond (half-integer w_XY); even N → self-paired center (integer w_XY)");
            yield return new InspectableNode("F86 generalization (Tier 1 derived 2026-05-03)",
                summary: "Q_peak(b) = Q_peak(N−2−b) bit-exact; verified N=5..7 c=2 + N=5..6 c=3 < 10⁻¹⁰; implemented as PerF71OrbitObservation in F86KnowledgeBase");
            yield return new InspectableNode("Pi2 ladder partial overlap",
                summary: "IndependentComponentCount(N) = ⌊N/2⌋ lands on Pi2 ladder anchor for N=2,3 (1 = a_1), N=4,5 (2 = a_0), N=8,9 (4 = a_{-1}); combinatorial elsewhere");
            // Sample bond-pair structure
            for (int N = 2; N <= 7; N++)
            {
                string regime = HasCenterBond(N) ? "even, integer-mirror, self-paired center" : "odd, half-integer-mirror, all pairs";
                int? ladderIdx = LadderIndexForIndependentComponentCount(N);
                string ladderInfo = ladderIdx.HasValue ? $", lands on a_{{{ladderIdx.Value}}}" : "";
                yield return new InspectableNode(
                    $"N={N}",
                    summary: $"{regime}; {BondCount(N)} bonds, {PalindromicPairCount(N)} pairs, {IndependentComponentCount(N)} independent components{ladderInfo}");
            }
        }
    }
}
