using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The dyadic halving ladder forced by the Pi2 foundation: a_n = 2^(1−n),
/// bidirectional in n ∈ ℤ.
///
/// <para>The trunk d²−2d=0 forces d=2 (<see cref="QubitDimensionalAnchorClaim"/>); the
/// Bloch decomposition ρ = (I + r·σ)/2 forces the 1/2 baseline shift
/// (<see cref="HalfAsStructuralFixedPointClaim"/>); the bilinear apex max p·(1−p)=1/4 at
/// p=1/2 forces the square (<see cref="QuarterAsBilinearMaxvalClaim"/>). Iterating the
/// halving operation by algebraic continuation gives the geometric sequence
/// 2, 1, 1/2, 1/4, 1/8, 1/16, ... Each term is the previous halved.</para>
///
/// <para>The ladder is bidirectional and inversion-symmetric: a_n · a_{2−n} = 1 for every
/// integer n (immediate from 2^{1−n} · 2^{1−(2−n)} = 2^0). The fixpoint of the inversion
/// is n=1 where a_1 = 1, the trivial identity scale. Tom's reading 2026-05-08:
/// "es muss auch eine Regel geben wo es wieder auf 1/1 springt, wenn der Bruch zu groß wird"
/// — the rule IS the inversion symmetry. There is no escape into ever-smaller fractions:
/// every n ≥ 2 (memory side) has a mirror at 2−n ≤ 0 (upper / operator-space side).</para>
///
/// <para>The upper side (n ≤ −1) carries the operator-space dimensions
/// d² = 4^N for an N-qubit system: a_{−(2N−1)} = 2^{2N} = 4^N. So the ladder unifies
/// "memory" 1/2^k anchors with "operator-space" 4^N dimensions through the inversion. See
/// <c>memory/feedback_d2_operator_space.md</c> for the project's d² operator-space
/// reading.</para>
///
/// <para>This is not a speculative pattern; the first non-trivial indices (n = 0, 2, 3)
/// are already typed Tier1Derived Claims. The ladder makes the inheritance explicit:
/// 1/4 is the square of 1/2 is the inverse of d. There is no alien hiding here; the
/// structure is forced at the base by d=2 and inherits algebraically.</para>
///
/// <para><b>n=0 carries two co-existing structural readings (2026-05-12):</b> the value
/// a₀ = 2 is BOTH the qubit operator-space dimension d=2 (from d²−2d=0; per
/// <see cref="QubitDimensionalAnchorClaim"/>) AND the absorption quantum coefficient 2 in
/// Re(λ) = −2γ·⟨n_XY⟩ (the Liouvillian eigenvalue grid step under uniform Z-dephasing;
/// per <see cref="AbsorptionTheoremClaim"/>). Same number, two roles: the qubit's
/// algebraic dimension and the dynamics' rate-quantum coefficient. The ladder makes the
/// doppelte Rolle visible: a₀ is the linchpin between Hilbert-space algebra (where it
/// counts basis dimensions) and Liouville-space dynamics (where it scales the carrier
/// γ₀ into the absorption quantum 2γ₀ per <see cref="UniversalCarrierClaim"/>). Both
/// readings sit at the same ladder index; neither dominates.</para>
///
/// <para>Tier: Tier1Derived. The closed form a_n = 2^(1−n) is trivial; the lineage from
/// the existing Pi2 foundation Claims is documented per known anchor entry. The open
/// part — does each n≥4 entry have a physical anchor in the framework? — is documented
/// as <see cref="OpenAnchorIndices"/> (Tier1Candidate prediction territory: searchable
/// with the same algebraic forcing logic).</para>
///
/// <para>Anchors: <c>docs/EXCLUSIONS.md</c> (d²−2d=0) +
/// <c>reflections/ON_THE_HALF.md</c> (the half as anchor + ladder) +
/// <c>docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md</c> (1/4 = (1/2)²).</para></summary>
public sealed class Pi2DyadicLadderClaim : Claim
{
    /// <summary>Pinned table of indices where the ladder term has a typed Tier1Derived
    /// Claim in the Pi2 foundation. Multiple entries per index are allowed when the same
    /// ladder term carries multiple co-existing structural roles (currently: n=0 has two,
    /// the qubit-dimension reading and the absorption-quantum-coefficient reading; same
    /// value 2.0, two co-existing meanings, both Tier 1 derived).</summary>
    public IReadOnlyList<DyadicAnchor> KnownAnchors { get; } = new[]
    {
        new DyadicAnchor(N: 0, Value: 2.0,
            ClaimType: typeof(QubitDimensionalAnchorClaim),
            ClaimName: "QubitDimensionalAnchorClaim",
            Role: "root: d=2 from d²−2d=0 (operator-space dimension of one qubit)"),
        new DyadicAnchor(N: 0, Value: 2.0,
            ClaimType: typeof(AbsorptionTheoremClaim),
            ClaimName: "AbsorptionTheoremClaim",
            Role: "absorption quantum coefficient: 2γ in Re(λ) = −2γ·⟨n_XY⟩ (Liouvillian eigenvalue grid step under uniform Z-dephasing)"),
        new DyadicAnchor(N: 2, Value: 0.5,
            ClaimType: typeof(HalfAsStructuralFixedPointClaim),
            ClaimName: "HalfAsStructuralFixedPointClaim",
            Role: "polarity baseline: ρ = (I + r·σ)/2 baseline 1/d, three-faces fixed point"),
        new DyadicAnchor(N: 3, Value: 0.25,
            ClaimType: typeof(QuarterAsBilinearMaxvalClaim),
            ClaimName: "QuarterAsBilinearMaxvalClaim",
            Role: "bilinear maxval: max p·(1−p) = (1/2)² = 1/4"),
    };

    /// <summary>Indices n ∈ {1, 4, 5, ...} where the ladder predicts a value but no typed
    /// Claim has been hinged yet. n=1 is the trivial identity scale (a_1=1); n≥4 are
    /// open predictions of the algebraic continuation.</summary>
    public IReadOnlyList<int> OpenAnchorIndices { get; } = new[] { 1, 4, 5, 6, 7, 8 };

    public Pi2DyadicLadderClaim()
        : base("Dyadic halving ladder a_n = 2^(1−n) (Pi2 foundation continuation)",
               Tier.Tier1Derived,
               "docs/EXCLUSIONS.md:251 (d²−2d=0) + reflections/ON_THE_HALF.md + docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md (1/4 = (1/2)²)")
    { }

    /// <summary>Closed form: <c>a_n = 2^(1−n)</c>, defined for all integer indices.
    /// Returns 4.0 at n=−1, 2.0 at n=0, 1.0 at n=1, 0.5 at n=2, 0.25 at n=3, etc.
    /// The ladder extends bidirectionally: positive n into the halving "memory" side
    /// (1/2, 1/4, 1/8, ...); negative n into the doubling "operator-space" side
    /// (4, 8, 16, ...).</summary>
    public double Term(int n) => Math.Pow(2.0, 1 - n);

    /// <summary>The mirror-partner index: <c>2 − n</c>. Encodes the inversion symmetry
    /// <c>a_n · a_{2−n} = 1</c> (every term times its mirror gives unity). Self-mirror
    /// fixpoint is <c>n = 1</c> where <c>a_1 = 1</c> and the partner index equals n itself.
    /// </summary>
    /// <example>
    /// MirrorPartnerIndex(0) = 2  (a_0=2 ↔ a_2=1/2)<br/>
    /// MirrorPartnerIndex(3) = -1 (a_3=1/4 ↔ a_{-1}=4)<br/>
    /// MirrorPartnerIndex(1) = 1  (self)
    /// </example>
    public int MirrorPartnerIndex(int n) => 2 - n;

    /// <summary>The product of <c>a_n · a_{MirrorPartnerIndex(n)}</c>. Algebraically equals
    /// 1.0 for every integer n (consequence of <c>2^{1−n} · 2^{1−(2−n)} = 2^0 = 1</c>);
    /// computed live as a drift check.</summary>
    public double ProductWithMirrorPartner(int n) => Term(n) * Term(MirrorPartnerIndex(n));

    /// <summary>The unique self-mirror index where <c>a_n = 1</c> and the partner is n
    /// itself. Always 1 by the inversion symmetry; exposed as a constant for documentation
    /// and tests.</summary>
    public int SelfMirrorIndex => 1;

    /// <summary>True iff index <paramref name="n"/> appears in <see cref="KnownAnchors"/>.</summary>
    public bool IsKnownAnchorIndex(int n) => KnownAnchors.Any(a => a.N == n);

    /// <summary>The first typed anchor at index <paramref name="n"/> if any exists, else
    /// null. When the index carries multiple roles (e.g. n=0 has both QubitDimensional
    /// and AbsorptionQuantum readings), this returns the first listed; use
    /// <see cref="AnchorsAt"/> to get all roles.</summary>
    public DyadicAnchor? AnchorAt(int n) => KnownAnchors.FirstOrDefault(a => a.N == n);

    /// <summary>All typed anchors at index <paramref name="n"/> (plural). At n=0 returns
    /// both QubitDimensionalAnchorClaim (d=2 = operator-space dim) and
    /// AbsorptionTheoremClaim (a₀=2 = absorption quantum coefficient); at n=2 and n=3
    /// currently returns one each. Returns empty enumerable when n has no typed anchor.
    /// </summary>
    public IEnumerable<DyadicAnchor> AnchorsAt(int n) => KnownAnchors.Where(a => a.N == n);

    /// <summary>The operator-space dimension <c>d² = 4^N</c> for an N-qubit system lands
    /// on the upper (negative-n) side of the ladder at index <c>−(2N − 1)</c>:
    /// <c>a_{−(2N−1)} = 2^{1−(−(2N−1))} = 2^{2N} = 4^N</c>. Returns the ladder index for
    /// the d² anchor of an N-qubit system.</summary>
    public int OperatorSpaceIndexForN(int N) => -(2 * N - 1);

    public override string DisplayName =>
        $"Pi2 dyadic halving ladder (a_n = 2^(1−n); {KnownAnchors.Count} known anchors)";

    public override string Summary =>
        $"a_n = 2^(1−n) bidirectional, inversion-symmetric a_n · a_{{2−n}} = 1; {KnownAnchors.Count} typed anchors at 3 distinct indices n ∈ {{0, 2, 3}} (n=0 carries TWO co-existing roles: qubit-dimension d=2 AND absorption-quantum coefficient 2γ); n=1 self-mirror pivot, n≥4 memory side, n≤−1 operator-space side d²=4^N at n=−(2N−1) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("closed form",
                summary: "a_n = 2^(1−n) for all n ∈ ℤ (bidirectional)");
            yield return new InspectableNode("forcing logic",
                summary: "d²−2d=0 → d=2 (root); ρ = (I + r·σ)/2 → 1/d=1/2 (shift); max p·(1−p)=1/4 at p=1/2 → square; iterate halving");
            yield return new InspectableNode("inversion symmetry",
                summary: "a_n · a_{2−n} = 1 (algebraic; every term has a multiplicative mirror); self-mirror at n=1, a_1=1");
            yield return new InspectableNode("upper side (n ≤ −1)",
                summary: "a_{−1}=4, a_{−2}=8, a_{−3}=16, ...; doubling powers; operator-space d²=4^N lands at n=−(2N−1)");
            foreach (var a in KnownAnchors)
                yield return new InspectableNode(
                    $"n={a.N} (anchor)",
                    summary: $"a_{a.N} = {a.Value:G6} ← {a.ClaimName} ({a.Role})");
            yield return new InspectableNode("n=1 (self-mirror pivot)",
                summary: "a_1 = 1.0; identity scale, fixpoint of the inversion a_n ↔ a_{2−n}");
            yield return new InspectableNode("n≥4 (memory side, open prediction)",
                summary: "a_4=0.125, a_5=0.0625, ...; algebra-continued; F88 dyadic-N singleton-mirror lands here (commit 6ef2cb9)");
        }
    }
}

/// <summary>One typed entry on the Pi2 dyadic halving ladder: index n, value 2^(1−n),
/// the <see cref="ClaimType"/> in the Pi2 foundation that anchors this term, and a
/// short role string for the ledger view.</summary>
public sealed record DyadicAnchor(
    int N,
    double Value,
    Type ClaimType,
    string ClaimName,
    string Role);
