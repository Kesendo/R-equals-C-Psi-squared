using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The dyadic halving ladder forced by the Pi2 foundation: a_n = 2^(1−n).
///
/// <para>The trunk d²−2d=0 forces d=2 (<see cref="QubitDimensionalAnchorClaim"/>); the
/// Bloch decomposition ρ = (I + r·σ)/2 forces the 1/2 baseline shift
/// (<see cref="HalfAsStructuralFixedPointClaim"/>); the bilinear apex max p·(1−p)=1/4 at
/// p=1/2 forces the square (<see cref="QuarterAsBilinearMaxvalClaim"/>). Iterating the
/// halving operation by algebraic continuation gives the geometric sequence
/// 2, 1, 1/2, 1/4, 1/8, 1/16, ... Each term is the previous halved.</para>
///
/// <para>This is not a speculative pattern; the first three non-trivial entries
/// (n = 0, 2, 3) are already typed Tier1Derived Claims. The ladder makes the inheritance
/// explicit: 1/4 is the square of 1/2 is the inverse of d. There is no alien hiding here;
/// the structure is forced at the base by d=2 and inherits algebraically.</para>
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
    /// Claim in the Pi2 foundation.</summary>
    public IReadOnlyList<DyadicAnchor> KnownAnchors { get; } = new[]
    {
        new DyadicAnchor(N: 0, Value: 2.0,
            ClaimType: typeof(QubitDimensionalAnchorClaim),
            ClaimName: "QubitDimensionalAnchorClaim",
            Role: "root: d=2 from d²−2d=0"),
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
    public IReadOnlyList<int> OpenAnchorIndices => new[] { 1, 4, 5, 6, 7, 8 };

    public Pi2DyadicLadderClaim()
        : base("Dyadic halving ladder a_n = 2^(1−n) (Pi2 foundation continuation)",
               Tier.Tier1Derived,
               "docs/EXCLUSIONS.md:251 (d²−2d=0) + reflections/ON_THE_HALF.md + docs/proofs/PROOF_BLOCK_CPSI_QUARTER.md (1/4 = (1/2)²)")
    { }

    /// <summary>Closed form: <c>a_n = 2^(1−n)</c>. Returns 2.0 for n=0, 1.0 for n=1,
    /// 0.5 for n=2, 0.25 for n=3, 0.125 for n=4, etc.</summary>
    public double Term(int n)
    {
        if (n < 0)
            throw new ArgumentOutOfRangeException(nameof(n), n, "Pi2DyadicLadder is indexed from n=0 (the d=2 root).");
        return Math.Pow(2.0, 1 - n);
    }

    /// <summary>True iff index <paramref name="n"/> appears in <see cref="KnownAnchors"/>.</summary>
    public bool IsKnownAnchorIndex(int n) => KnownAnchors.Any(a => a.N == n);

    /// <summary>The typed anchor at index <paramref name="n"/> if one exists, else null.</summary>
    public DyadicAnchor? AnchorAt(int n) => KnownAnchors.FirstOrDefault(a => a.N == n);

    public override string DisplayName =>
        $"Pi2 dyadic halving ladder (a_n = 2^(1−n); {KnownAnchors.Count} known anchors)";

    public override string Summary =>
        $"a_n = 2^(1−n): 2, 1, 1/2, 1/4, 1/8, ...; {KnownAnchors.Count} typed anchors at n ∈ {{0, 2, 3}}; n=1 trivial, n≥4 open ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("closed form",
                summary: "a_n = 2^(1−n) = 2 · (1/2)^n");
            yield return new InspectableNode("forcing logic",
                summary: "d²−2d=0 → d=2 (root); ρ = (I + r·σ)/2 → 1/d=1/2 (shift); max p·(1−p)=1/4 at p=1/2 → square; iterate halving");
            foreach (var a in KnownAnchors)
                yield return new InspectableNode(
                    $"n={a.N} (anchor)",
                    summary: $"a_{a.N} = {a.Value:G6} ← {a.ClaimName} ({a.Role})");
            yield return new InspectableNode("n=1 (trivial)",
                summary: "a_1 = 1.0; identity scale, no typed Claim hinged");
            yield return new InspectableNode("n≥4 (open prediction)",
                summary: "a_4=0.125, a_5=0.0625, ...; algebra-continued; open: which framework anchors land here?");
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
