using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track, T3 of 4: per-bond k-mode profile distilled
/// from <see cref="C2BlockJwDecomposition"/> by row-summing |C_b[k1, k2]| over k2.
///
/// <para>For each bond b ∈ [0, N−2] this primitive reduces the N×N coefficient matrix
/// <c>C_b[k1, k2]</c> from T2 to a length-N row-L1 profile</para>
///
/// <para><c>rowL1[k] = Σ_{k2 = 0..N-1} |C_b[k, k2]|</c> for k = 0..N−1.</para>
///
/// <para>From this profile we derive three per-bond summary statistics:</para>
/// <list type="bullet">
///   <item><c>K_90(b)</c> — minimum number of k-modes (sorted by row-L1 descending) needed
///   to capture 90% of the bond's total L1 mass.</item>
///   <item><c>K_99(b)</c> — same but at 99%; monotone <c>K_99(b) ≥ K_90(b)</c>.</item>
///   <item><c>TopThreeKIndices(b)</c> — the 3 dominant k-mode 1-indexed indices, sorted by
///   row-L1 descending. For <c>N &lt; 3</c> all available indices are returned in row-L1
///   descending order.</item>
/// </list>
///
/// <para>Aggregate <see cref="EndpointK90Mean"/> / <see cref="InteriorK90Mean"/> average
/// <c>K_90</c> over the two bond classes; both return <see cref="double.NaN"/> when no bond
/// of that class exists. <see cref="MinInteriorK90"/> returns the smallest <c>K_90</c> among
/// Interior bonds (the innermost-Interior reference for the Endpoint comparison).</para>
///
/// <para>The Endpoint-vs-Interior K_90 asymmetry direction flips with N: at N = 5..7 the
/// Endpoint K_90 mean dominates the Interior K_90 mean, but at N = 8 the flanking Interior
/// bonds (b ∈ {1, 2, 4, 5}) develop more k-uniform profiles than the Endpoints, while the
/// innermost-Interior bond stays the most localized. Comparing the Endpoint mean against
/// <see cref="MinInteriorK90"/> preserves the structural claim "boundary bonds are at least
/// as spread as the deepest-Interior bond" across N = 5..8.</para>
///
/// <para><b>Class-level Tier: <c>Tier2Verified</c>.</b> The row-L1 distillation is a
/// numerical post-processing of T2's algebraically derived coefficient matrix; the
/// distillation itself adds no analytic content beyond T2's textbook XY-JW identity. The
/// open Tier 1 promotion path is the closed-form Endpoint HWHM/Q_peak from the JW spectrum,
/// which is the open work tracked alongside the JW track in <c>F86OpenQuestions</c> (Item 1'
/// Direction (b'')). Until the analytical closed form is derived (F90 has numerical Tier-1;
/// analytical closed form remains open), the per-bond k-mode profile here is
/// the structural lens that highlights which sine modes carry the bond's spectral weight.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track).</para>
///
/// <para>F90 status (2026-05-11): the F86 c=2 ↔ F89 path-(N−1) bridge identity
/// achieves numerical Tier-1 for Direction (b'') via per-bond Hellmann-Feynman
/// (bit-exact 20/22 bonds at N=5..8). The JW-track primitives in this file remain
/// active as the alternative analytical route toward the closed-form HWHM_left/Q_peak
/// constants; the per-bond numerical answer itself is no longer the open piece.
/// See <c>docs/proofs/PROOF_F90_F86C2_BRIDGE.md</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F90F86C2BridgeIdentity.cs</c>.</para>
/// </summary>
public sealed class C2BondKModeProfile : Claim
{
    /// <summary>Composition: T2 per-bond bilinear coefficients C_b(k1, k2). T2's sine-mode
    /// basis (T1) is reachable through <c>Decomposition.Modes</c>.</summary>
    public C2BlockJwDecomposition Decomposition { get; }

    /// <summary>One <see cref="BondKModeWitness"/> per bond, in bond-index order, exposing
    /// row-L1 profile + K_90 / K_99 / TopThreeKIndices.</summary>
    public IReadOnlyList<BondKModeWitness> Bonds { get; }

    /// <summary>Mean of <see cref="BondKModeWitness.K90"/> over bonds tagged
    /// <see cref="BondClass.Endpoint"/>; <see cref="double.NaN"/> if no Endpoint bond.</summary>
    public double EndpointK90Mean { get; }

    /// <summary>Mean of <see cref="BondKModeWitness.K90"/> over bonds tagged
    /// <see cref="BondClass.Interior"/>; <see cref="double.NaN"/> if no Interior bond.</summary>
    public double InteriorK90Mean { get; }

    /// <summary>Minimum <see cref="BondKModeWitness.K90"/> among bonds tagged
    /// <see cref="BondClass.Interior"/> (the most localized Interior bond, structurally the
    /// innermost-Interior reference). <see cref="double.NaN"/> if no Interior bond.</summary>
    public double MinInteriorK90 { get; }

    /// <summary>Public factory: validates c=2, builds T2 <see cref="C2BlockJwDecomposition"/>,
    /// then distills each bond's coefficient matrix to a row-L1 profile and computes K_90,
    /// K_99, top-three k-modes. Aggregates K_90 by <see cref="BondClass"/>.</summary>
    public static C2BondKModeProfile Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"C2BondKModeProfile applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var decomp = C2BlockJwDecomposition.Build(block);
        int N = block.N;

        var bondWitnesses = new BondKModeWitness[decomp.Bonds.Count];
        for (int b = 0; b < decomp.Bonds.Count; b++)
        {
            var bondCoeffs = decomp.Bonds[b];
            var rowL1 = new double[N];
            for (int k1 = 0; k1 < N; k1++)
            {
                double sum = 0;
                for (int k2 = 0; k2 < N; k2++)
                    sum += Math.Abs(bondCoeffs.Coefficients[k1, k2]);
                rowL1[k1] = sum;
            }

            // Sort k-mode indices (1-indexed) by row-L1 descending; stable for ties via
            // ascending k as secondary key.
            var sortedIndices = Enumerable.Range(1, N)
                .OrderByDescending(k => rowL1[k - 1])
                .ThenBy(k => k)
                .ToArray();

            double total = rowL1.Sum();
            int k90 = MinModesForFraction(sortedIndices, rowL1, total, 0.90);
            int k99 = MinModesForFraction(sortedIndices, rowL1, total, 0.99);

            int topCount = Math.Min(3, N);
            var topThree = sortedIndices.Take(topCount).ToArray();

            bondWitnesses[b] = new BondKModeWitness(
                Bond: bondCoeffs.Bond,
                BondClass: bondCoeffs.BondClass,
                K90: k90,
                K99: k99,
                TopThreeKIndices: topThree,
                RowL1Profile: rowL1);
        }

        double endpointMean = MeanK90(bondWitnesses, BondClass.Endpoint);
        double interiorMean = MeanK90(bondWitnesses, BondClass.Interior);
        double interiorMin = MinK90(bondWitnesses, BondClass.Interior);

        return new C2BondKModeProfile(decomp, bondWitnesses, endpointMean, interiorMean, interiorMin);
    }

    /// <summary>Minimum number of modes (taken in <paramref name="sortedIndices"/> order)
    /// needed for cumulative row-L1 mass to reach <paramref name="fraction"/> · <paramref name="total"/>.
    /// Returns 1 for the trivial total = 0 case (all-zero profile).</summary>
    private static int MinModesForFraction(
        IReadOnlyList<int> sortedIndices,
        IReadOnlyList<double> rowL1,
        double total,
        double fraction)
    {
        if (total <= 0) return 1;
        double threshold = fraction * total;
        double accum = 0;
        for (int i = 0; i < sortedIndices.Count; i++)
        {
            accum += rowL1[sortedIndices[i] - 1];
            if (accum >= threshold) return i + 1;
        }
        return sortedIndices.Count;
    }

    private static double MeanK90(IReadOnlyList<BondKModeWitness> bonds, BondClass cls)
    {
        var matching = bonds.Where(w => w.BondClass == cls).ToArray();
        if (matching.Length == 0) return double.NaN;
        return matching.Average(w => (double)w.K90);
    }

    private static double MinK90(IReadOnlyList<BondKModeWitness> bonds, BondClass cls)
    {
        var matching = bonds.Where(w => w.BondClass == cls).ToArray();
        if (matching.Length == 0) return double.NaN;
        return matching.Min(w => (double)w.K90);
    }

    private C2BondKModeProfile(
        C2BlockJwDecomposition decomposition,
        IReadOnlyList<BondKModeWitness> bonds,
        double endpointK90Mean,
        double interiorK90Mean,
        double minInteriorK90)
        : base("c=2 bond k-mode profile: row-L1 distillation of T2 coefficients",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track)")
    {
        Decomposition = decomposition;
        Bonds = bonds;
        EndpointK90Mean = endpointK90Mean;
        InteriorK90Mean = interiorK90Mean;
        MinInteriorK90 = minInteriorK90;
    }

    public override string DisplayName =>
        $"c=2 bond k-mode profile (N={Decomposition.Block.N}, bonds={Decomposition.Block.NumBonds})";

    public override string Summary =>
        $"K_90 mean: Endpoint = {FormatMean(EndpointK90Mean)}, Interior = {FormatMean(InteriorK90Mean)} ({Tier.Label()})";

    private static string FormatMean(double value) =>
        double.IsNaN(value) ? "n/a" : value.ToString("F2");

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Decomposition.Block.N.ToString());
            yield return new InspectableNode("NumBonds", summary: Decomposition.Block.NumBonds.ToString());
            yield return Decomposition;
            yield return InspectableNode.RealScalar("EndpointK90Mean", EndpointK90Mean, "F2");
            yield return InspectableNode.RealScalar("InteriorK90Mean", InteriorK90Mean, "F2");
            yield return InspectableNode.RealScalar("MinInteriorK90", MinInteriorK90, "F2");
            yield return InspectableNode.Group("Bonds", Bonds, Bonds.Count);
        }
    }
}

/// <summary>Per-bond k-mode profile witness: row-L1 distillation of the T2 N×N coefficient
/// matrix plus the K_90 / K_99 / top-three summary statistics.
///
/// <para><see cref="RowL1Profile"/> is the length-N row-L1 vector indexed 0..N−1 with
/// <c>RowL1Profile[k − 1] = Σ_{k2} |C_b[k − 1, k2]|</c> for k = 1..N.</para>
///
/// <para><see cref="K90"/> / <see cref="K99"/> are the minimum number of modes (sorted by
/// row-L1 descending) needed for cumulative L1 mass to reach 90% / 99% of the bond's total.</para>
///
/// <para><see cref="TopThreeKIndices"/> contains the 3 dominant 1-indexed k-mode indices,
/// sorted by row-L1 descending; for <c>N &lt; 3</c> all available indices are returned.</para>
/// </summary>
public sealed record BondKModeWitness(
    int Bond,
    BondClass BondClass,
    int K90,
    int K99,
    IReadOnlyList<int> TopThreeKIndices,
    IReadOnlyList<double> RowL1Profile
) : IInspectable
{
    public string DisplayName => $"bond {Bond} k-mode profile ({BondClass})";

    public string Summary =>
        $"K_90 = {K90}, K_99 = {K99}, top-3 = [{string.Join(",", TopThreeKIndices)}]";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("bond", summary: Bond.ToString());
            yield return new InspectableNode("bond class", summary: BondClass.ToString());
            yield return new InspectableNode("K_90", summary: K90.ToString());
            yield return new InspectableNode("K_99", summary: K99.ToString());
            yield return new InspectableNode("top-3 k indices",
                summary: $"[{string.Join(",", TopThreeKIndices)}]");
            yield return new InspectableNode("row-L1 profile (k = 1..N)",
                summary: string.Join(" ", RowL1Profile.Select(v => v.ToString("F4"))));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
