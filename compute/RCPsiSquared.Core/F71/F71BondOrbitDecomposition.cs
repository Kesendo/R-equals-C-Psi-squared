using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F71;

/// <summary>F71 spatial-mirror bond pairing: bond b ↔ bond N−2−b. For an N-qubit chain
/// (N−1 bonds), this gives ⌈(N−1)/2⌉ orbits. If N is even, the middle bond b=(N−2)/2 is
/// self-paired (its mirror partner is itself).
///
/// <para>Per F71-mirror invariance, any F71-symmetric chain observable evaluated per-bond
/// is identical for paired bonds. F86 generalisation: per-bond Q_peak is bit-exact mirror-
/// invariant (Statement 3 of <c>PROOF_F86_QPEAK</c>).</para>
/// </summary>
public sealed class F71BondOrbitDecomposition : Claim
{
    public int N { get; }
    public int NumBonds => N - 1;
    public int NumOrbits => (NumBonds + 1) / 2;
    public bool HasSelfPairedCentralOrbit => NumBonds % 2 == 1;
    public IReadOnlyList<F71BondOrbit> Orbits { get; }

    public F71BondOrbitDecomposition(int N)
        : base("F71 bond-mirror orbit decomposition",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md + ANALYTICAL_FORMULAS F71 + PROOF_F86_QPEAK Statement 3")
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"N must be ≥ 2; got {N}");
        this.N = N;
        var orbits = new List<F71BondOrbit>();
        int numBonds = N - 1;
        for (int b = 0; b < (numBonds + 1) / 2; b++)
        {
            int mirror = numBonds - 1 - b;
            if (b == mirror) orbits.Add(new F71BondOrbit(b, null));
            else orbits.Add(new F71BondOrbit(b, mirror));
        }
        Orbits = orbits;
    }

    public override string DisplayName => $"F71 orbits (N={N}, {NumBonds} bonds, {NumOrbits} orbits)";

    public override string Summary =>
        HasSelfPairedCentralOrbit
            ? $"{NumOrbits} orbits including a self-paired central bond at b={(NumBonds - 1) / 2}"
            : $"{NumOrbits} orbit pairs (N−1 = {NumBonds} is even, no self-paired)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return InspectableNode.RealScalar("N", N);
            yield return InspectableNode.RealScalar("number of bonds", NumBonds);
            yield return InspectableNode.RealScalar("number of orbits", NumOrbits);
            foreach (var orbit in Orbits) yield return orbit;
        }
    }
}

/// <summary>One F71 bond orbit: either a (b, N−2−b) pair when distinct, or a self-paired
/// central bond when b = N−2−b. The latter exists iff numBonds = N−1 is odd, i.e. N is
/// even.</summary>
public sealed class F71BondOrbit : IInspectable
{
    public int BondA { get; }
    public int? BondB { get; }
    public bool IsSelfPaired => BondB is null;

    public F71BondOrbit(int bondA, int? bondB)
    {
        BondA = bondA;
        BondB = bondB;
    }

    public string DisplayName => IsSelfPaired
        ? $"orbit{{b={BondA}}} (self-paired)"
        : $"orbit{{b={BondA} ↔ b={BondB}}}";

    public string Summary => IsSelfPaired
        ? $"central bond b={BondA} fixed by mirror"
        : $"bond pair (b={BondA}, b={BondB})";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return InspectableNode.RealScalar("bond A", BondA);
            if (BondB is not null) yield return InspectableNode.RealScalar("bond B (mirror)", BondB.Value);
            yield return new InspectableNode("self-paired", summary: IsSelfPaired ? "yes" : "no");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
