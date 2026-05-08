using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F71-orbit-grouped JW k-mode profile for one c=2 block: bridges the JW track
/// (T1-T3) to the F71 mirror-pair structure {b, N−2−b}.
///
/// <para>The witness exposed here is that the per-bond k-mode profile from
/// <see cref="C2BondKModeProfile"/> is <b>F71-mirror-invariant</b>:</para>
///
/// <list type="bullet">
///   <item>K_90 within each orbit is bit-identical (integer comparison, deviation = 0).</item>
///   <item>K_99 within each orbit is bit-identical.</item>
///   <item>RowL1Profile within each orbit agrees to floating-point precision
///   (<see cref="F71MirrorRowL1Tolerance"/> = 1e-12).</item>
///   <item>TopThreeKIndices within each orbit are identical (deterministic sort on
///   identical row-L1 values + ascending-k tiebreak).</item>
/// </list>
///
/// <para>The mathematical statement underlying these witnesses: the OBC sine-mode parity
/// <c>ψ_k(N−1−j) = (−1)^(k+1) ψ_k(j)</c> implies <c>C_{N−2−b}(k1, k2) = (−1)^(k1+k2) C_b(k1, k2)</c>.
/// The sign cancels in absolute values, so <c>|C_{N−2−b}| = |C_b|</c> matrix-element-wise,
/// and every quantity derived from <c>|C_b|</c> (row-L1 profile, K_90, K_99, top-K indices)
/// is automatically F71-mirror-invariant. This is a textbook XY-JW + F71-algebra theorem;
/// the <see cref="F71MirrorRowL1Tolerance"/> bound only catches floating-point drift.</para>
///
/// <para><b>Class-level Tier: <c>Tier2Verified</c>.</b> The underlying theorem is Tier1
/// algebraic content (textbook XY-JW + F71 sine-mode parity); this primitive's runtime
/// witness is the numerical residual at construction. Setting tier to Tier2 mirrors
/// <see cref="PerF71OrbitKTable"/>'s convention (the analogous K-resonance witness table
/// also keeps Tier2 despite F71 mirror invariance being algebraically exact for Q_peak
/// per Statement 3 of <c>PROOF_F86_QPEAK</c>).</para>
///
/// <para>Why this matters: at N=8 the JW K_90 already organises by F71 orbit
/// ({orbit 0, orbit 3} = 7 vs {orbit 1, orbit 2} = 8 — boundary + center vs intermediate),
/// even though the slow K-resonance Q_peak/HWHM observable does not yet show the split
/// loudly. At N≥9 the same orbital reorganisation becomes loud in Q_peak — the F71-orbit
/// "Ausbruch". The JW lens detects the precursor one N-step earlier than the K-resonance
/// lens; this primitive gives that observation a typed-knowledge home.</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c> for the F71 mirror-pair algebra.</para>
/// </summary>
public sealed class PerF71OrbitKModeTable : Claim
{
    /// <summary>Floating-point tolerance for the row-L1 within-orbit deviation. The
    /// algebra gives bit-exact equality; this bound only catches FP drift.</summary>
    public const double F71MirrorRowL1Tolerance = 1e-12;

    public CoherenceBlock Block { get; }

    /// <summary>Composition: T3 per-bond k-mode profile. Provides the bond-level witnesses
    /// that this table groups by F71 orbit.</summary>
    public C2BondKModeProfile Profile { get; }

    /// <summary>One <see cref="OrbitKModeWitness"/> per F71 orbit, in orbit-decomposition
    /// order. Self-paired orbits hold the single bond's witnesses; mirror-pair orbits
    /// hold the (validated-identical) bond-A witnesses.</summary>
    public IReadOnlyList<OrbitKModeWitness> OrbitWitnesses { get; }

    /// <summary>Maximum |K_90_A − K_90_B| over all 2-bond orbits. Expected = 0 (bit-exact).</summary>
    public int MaxK90WithinOrbitDeviation { get; }

    /// <summary>Maximum |K_99_A − K_99_B| over all 2-bond orbits. Expected = 0 (bit-exact).</summary>
    public int MaxK99WithinOrbitDeviation { get; }

    /// <summary>Maximum element-wise |rowL1_A[k] − rowL1_B[k]| over all 2-bond orbits and
    /// all k = 1..N. Expected &lt; <see cref="F71MirrorRowL1Tolerance"/>.</summary>
    public double MaxRowL1WithinOrbitDeviation { get; }

    /// <summary>Public factory: validates c=2, builds T3 <see cref="C2BondKModeProfile"/>,
    /// then walks the F71 orbit decomposition and joins each orbit with the per-bond
    /// witnesses, throwing if the F71-mirror invariance is violated (which would indicate
    /// a numerical regression, since mirror-paired bonds are bit-identical at the algebraic
    /// level).</summary>
    public static PerF71OrbitKModeTable Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"PerF71OrbitKModeTable applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var profile = C2BondKModeProfile.Build(block);
        var orbits = new F71BondOrbitDecomposition(block.N).Orbits;
        var bondWitnesses = profile.Bonds;

        var orbitWitnesses = new List<OrbitKModeWitness>(orbits.Count);
        int maxK90Dev = 0;
        int maxK99Dev = 0;
        double maxRowL1Dev = 0;
        foreach (var orbit in orbits)
        {
            if (orbit.IsSelfPaired)
            {
                var bw = bondWitnesses[orbit.BondA];
                orbitWitnesses.Add(new OrbitKModeWitness(
                    Orbit: orbit,
                    K90: bw.K90,
                    K99: bw.K99,
                    TopThreeKIndices: bw.TopThreeKIndices,
                    RowL1Profile: bw.RowL1Profile,
                    K90WithinOrbitDeviation: 0,
                    K99WithinOrbitDeviation: 0,
                    RowL1WithinOrbitDeviation: 0));
            }
            else
            {
                var bA = bondWitnesses[orbit.BondA];
                var bB = bondWitnesses[orbit.BondB!.Value];
                int k90Dev = Math.Abs(bA.K90 - bB.K90);
                int k99Dev = Math.Abs(bA.K99 - bB.K99);
                double rowL1Dev = MaxAbsDifference(bA.RowL1Profile, bB.RowL1Profile);

                if (k90Dev != 0 || k99Dev != 0 || rowL1Dev > F71MirrorRowL1Tolerance)
                    throw new InvalidOperationException(
                        $"F71 mirror invariance violated at orbit {{b={orbit.BondA}, b={orbit.BondB}}}: " +
                        $"|ΔK_90|={k90Dev}, |ΔK_99|={k99Dev}, max|ΔrowL1|={rowL1Dev:E2}");

                maxK90Dev = Math.Max(maxK90Dev, k90Dev);
                maxK99Dev = Math.Max(maxK99Dev, k99Dev);
                maxRowL1Dev = Math.Max(maxRowL1Dev, rowL1Dev);

                orbitWitnesses.Add(new OrbitKModeWitness(
                    Orbit: orbit,
                    K90: bA.K90,
                    K99: bA.K99,
                    TopThreeKIndices: bA.TopThreeKIndices,
                    RowL1Profile: bA.RowL1Profile,
                    K90WithinOrbitDeviation: k90Dev,
                    K99WithinOrbitDeviation: k99Dev,
                    RowL1WithinOrbitDeviation: rowL1Dev));
            }
        }

        return new PerF71OrbitKModeTable(block, profile, orbitWitnesses, maxK90Dev, maxK99Dev, maxRowL1Dev);
    }

    private static double MaxAbsDifference(IReadOnlyList<double> a, IReadOnlyList<double> b)
    {
        if (a.Count != b.Count)
            throw new ArgumentException($"profile length mismatch: {a.Count} vs {b.Count}");
        double m = 0;
        for (int i = 0; i < a.Count; i++)
            m = Math.Max(m, Math.Abs(a[i] - b[i]));
        return m;
    }

    private PerF71OrbitKModeTable(
        CoherenceBlock block,
        C2BondKModeProfile profile,
        IReadOnlyList<OrbitKModeWitness> orbitWitnesses,
        int maxK90Dev,
        int maxK99Dev,
        double maxRowL1Dev)
        : base("c=2 F71-orbit JW k-mode witness table",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md")
    {
        Block = block;
        Profile = profile;
        OrbitWitnesses = orbitWitnesses;
        MaxK90WithinOrbitDeviation = maxK90Dev;
        MaxK99WithinOrbitDeviation = maxK99Dev;
        MaxRowL1WithinOrbitDeviation = maxRowL1Dev;
    }

    public override string DisplayName =>
        $"F71-orbit JW k-mode table (N={Block.N}, {OrbitWitnesses.Count} orbits)";

    public override string Summary =>
        $"{OrbitWitnesses.Count} F71 orbits at c=2 N={Block.N}; F71-mirror invariant " +
        $"(max-Δ K_90 = {MaxK90WithinOrbitDeviation}, max-Δ rowL1 = {MaxRowL1WithinOrbitDeviation:G3}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumOrbits", summary: OrbitWitnesses.Count.ToString());
            yield return Profile;
            yield return InspectableNode.RealScalar("MaxK90WithinOrbitDeviation", MaxK90WithinOrbitDeviation);
            yield return InspectableNode.RealScalar("MaxK99WithinOrbitDeviation", MaxK99WithinOrbitDeviation);
            yield return InspectableNode.RealScalar("MaxRowL1WithinOrbitDeviation", MaxRowL1WithinOrbitDeviation, "G3");
            yield return InspectableNode.Group("OrbitWitnesses",
                OrbitWitnesses.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One F71-orbit's JW k-mode witnesses: K_90, K_99, top-three k-indices, row-L1
/// profile (all F71-mirror-invariant per the algebra), plus the within-orbit deviations
/// recorded at construction (zero for self-paired orbits, ≤ tolerance for mirror-pair
/// orbits).</summary>
public sealed record OrbitKModeWitness(
    F71BondOrbit Orbit,
    int K90,
    int K99,
    IReadOnlyList<int> TopThreeKIndices,
    IReadOnlyList<double> RowL1Profile,
    int K90WithinOrbitDeviation,
    int K99WithinOrbitDeviation,
    double RowL1WithinOrbitDeviation
) : IInspectable
{
    public string DisplayName =>
        Orbit.IsSelfPaired
            ? $"OrbitKModeWitness {{b={Orbit.BondA}}} (self-paired)"
            : $"OrbitKModeWitness {{b={Orbit.BondA} ↔ b={Orbit.BondB}}}";

    public string Summary =>
        $"K_90 = {K90}, K_99 = {K99}, top-3 = [{string.Join(",", TopThreeKIndices)}]";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return Orbit;
            yield return InspectableNode.RealScalar("K_90", K90);
            yield return InspectableNode.RealScalar("K_99", K99);
            yield return new InspectableNode("top-3 k indices",
                summary: $"[{string.Join(",", TopThreeKIndices)}]");
            yield return new InspectableNode("row-L1 profile (k = 1..N)",
                summary: string.Join(" ", RowL1Profile.Select(v => v.ToString("F4"))));
            yield return InspectableNode.RealScalar("ΔK_90 within orbit", K90WithinOrbitDeviation);
            yield return InspectableNode.RealScalar("ΔK_99 within orbit", K99WithinOrbitDeviation);
            yield return InspectableNode.RealScalar("ΔrowL1 within orbit", RowL1WithinOrbitDeviation, "G3");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
