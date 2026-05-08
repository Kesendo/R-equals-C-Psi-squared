using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.F71;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F71-orbit-grouped EP-resonance Q-peak witness table for one c=2 block: bridges
/// the JW track (T7 step B) to the F71 mirror-pair structure {b, N−2−b}.
///
/// <para>The witness exposed here is that the per-bond L(Q) Q-peak scan from
/// <see cref="C2BondLQPeakScan"/> is <b>F71-mirror-invariant</b>:</para>
///
/// <list type="bullet">
///   <item>Q_peak within each orbit agrees to floating-point precision
///   (<see cref="F71MirrorQPeakTolerance"/> = 1e-8).</item>
///   <item>‖xB‖_max within each orbit agrees to floating-point precision
///   (<see cref="F71MirrorKMaxTolerance"/> = 1e-8).</item>
/// </list>
///
/// <para>The mathematical statement underlying these witnesses: per F71-mirror invariance
/// (Statement 3 of <c>PROOF_F86_QPEAK</c>), the per-bond Q-resonance peak position and
/// height are bit-exact equal under bond b ↔ bond N−2−b at the algebraic level. The
/// underlying L(Q)-eigenbasis projection and its Frobenius norm depend only on the
/// F71-symmetric structure of <c>D + j·M_h_total</c>, which carries the same
/// OBC sine-mode parity that fixed K-mode invariants in
/// <see cref="PerF71OrbitKModeTable"/>. The <see cref="F71MirrorQPeakTolerance"/> /
/// <see cref="F71MirrorKMaxTolerance"/> bounds only catch FP drift from EVD ordering
/// and parabolic-refinement arithmetic.</para>
///
/// <para><b>Class-level Tier: <c>Tier2Verified</c>.</b> Numerical orbit aggregation;
/// the underlying F71-mirror algebraic invariance is Tier1 (OBC sine-mode parity per
/// PROOF_C1_MIRROR_SYMMETRY), the runtime witness here is the numerical residual at
/// construction. Mirrors <see cref="PerF71OrbitKModeTable"/>'s Tier convention for the
/// JW-track orbit-aggregation primitives.</para>
///
/// <para>Why this matters: for the EP-resonance lens, the F71-orbit aggregation is the
/// natural reduction of N−1 per-bond witnesses to ⌈(N−1)/2⌉ orbit-level witnesses
/// (Endpoint + ⌊(N−1)/2⌋−1 Interior orbits + optional self-paired central orbit at even
/// N). The Endpoint vs Interior split is the load-bearing structure for the
/// BareDoubledPtfXPeak ≈ 2.197 floor and the per-bond `g_eff_probe(Endpoint) ≈
/// g_eff_probe(Interior)/1.6` candidate from F86 Item 1' Direction (a').</para>
///
/// <para>Anchors: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// <c>docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md</c> for the F71 mirror-pair algebra.</para>
/// </summary>
public sealed class PerF71OrbitLQPeakTable : Claim
{
    /// <summary>Floating-point tolerance for the within-orbit Q_peak deviation. Algebra
    /// gives bit-exact equality; this bound only catches FP drift from EVD ordering and
    /// parabolic refinement.</summary>
    public const double F71MirrorQPeakTolerance = 1e-8;

    /// <summary>Floating-point tolerance for the within-orbit ‖xB‖_max deviation.</summary>
    public const double F71MirrorKMaxTolerance = 1e-8;

    public CoherenceBlock Block { get; }

    /// <summary>Composition: T7 step B per-bond L(Q) Q-peak scan. Provides the bond-level
    /// witnesses that this table groups by F71 orbit.</summary>
    public C2BondLQPeakScan PeakScan { get; }

    /// <summary>One <see cref="OrbitLQPeakWitness"/> per F71 orbit, in orbit-decomposition
    /// order. Self-paired orbits hold the single bond's witnesses; mirror-pair orbits
    /// hold the (validated-identical) bond-A witnesses.</summary>
    public IReadOnlyList<OrbitLQPeakWitness> OrbitWitnesses { get; }

    /// <summary>Maximum |Q_peak_A − Q_peak_B| over all 2-bond orbits. Expected
    /// &lt; <see cref="F71MirrorQPeakTolerance"/>.</summary>
    public double MaxQPeakWithinOrbitDeviation { get; }

    /// <summary>Maximum |‖xB‖_max_A − ‖xB‖_max_B| over all 2-bond orbits. Expected
    /// &lt; <see cref="F71MirrorKMaxTolerance"/>.</summary>
    public double MaxKMaxWithinOrbitDeviation { get; }

    /// <summary>Public factory: validates c=2, builds T7 <see cref="C2BondLQPeakScan"/>,
    /// then walks the F71 orbit decomposition and joins each orbit with the per-bond
    /// witnesses, throwing if F71-mirror invariance is violated (which would indicate
    /// a numerical regression beyond the FP-drift tolerance).</summary>
    public static PerF71OrbitLQPeakTable Build(CoherenceBlock block, IReadOnlyList<double>? qGrid = null)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"PerF71OrbitLQPeakTable applies only to the c=2 stratum; got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        var peakScan = C2BondLQPeakScan.Build(block, qGrid);
        var orbits = new F71BondOrbitDecomposition(block.N).Orbits;
        var bondWitnesses = peakScan.Bonds;

        var orbitWitnesses = new List<OrbitLQPeakWitness>(orbits.Count);
        double maxQPeakDev = 0;
        double maxKMaxDev = 0;
        foreach (var orbit in orbits)
        {
            if (orbit.IsSelfPaired)
            {
                var bw = bondWitnesses[orbit.BondA];
                orbitWitnesses.Add(new OrbitLQPeakWitness(
                    Orbit: orbit,
                    QPeak: bw.QPeak,
                    KMax: bw.XbNormAtPeak,
                    HwhmLeft: bw.HwhmLeft,
                    HwhmRight: bw.HwhmRight,
                    HwhmLeftOverQPeak: bw.HwhmLeftOverQPeak,
                    QPeakWithinOrbitDeviation: 0,
                    KMaxWithinOrbitDeviation: 0));
            }
            else
            {
                var bA = bondWitnesses[orbit.BondA];
                var bB = bondWitnesses[orbit.BondB!.Value];
                double qPeakDev = Math.Abs(bA.QPeak - bB.QPeak);
                double kMaxDev = Math.Abs(bA.XbNormAtPeak - bB.XbNormAtPeak);

                bool qPeakViolated = qPeakDev > F71MirrorQPeakTolerance;
                bool kMaxViolated = kMaxDev > F71MirrorKMaxTolerance;
                if (qPeakViolated || kMaxViolated)
                {
                    string violator = (qPeakViolated, kMaxViolated) switch
                    {
                        (true, true) => "QPeak + KMax tolerances",
                        (true, false) => "QPeak tolerance",
                        (false, true) => "KMax tolerance",
                        _ => "(unreachable)"
                    };
                    throw new InvalidOperationException(
                        $"F71 mirror invariance violated at orbit {{b={orbit.BondA}, b={orbit.BondB}}} ({violator}): " +
                        $"|ΔQ_peak|={qPeakDev:E2} (tol {F71MirrorQPeakTolerance:G3}), " +
                        $"|ΔKMax|={kMaxDev:E2} (tol {F71MirrorKMaxTolerance:G3})");
                }

                maxQPeakDev = Math.Max(maxQPeakDev, qPeakDev);
                maxKMaxDev = Math.Max(maxKMaxDev, kMaxDev);

                orbitWitnesses.Add(new OrbitLQPeakWitness(
                    Orbit: orbit,
                    QPeak: bA.QPeak,
                    KMax: bA.XbNormAtPeak,
                    HwhmLeft: bA.HwhmLeft,
                    HwhmRight: bA.HwhmRight,
                    HwhmLeftOverQPeak: bA.HwhmLeftOverQPeak,
                    QPeakWithinOrbitDeviation: qPeakDev,
                    KMaxWithinOrbitDeviation: kMaxDev));
            }
        }

        return new PerF71OrbitLQPeakTable(block, peakScan, orbitWitnesses, maxQPeakDev, maxKMaxDev);
    }

    private PerF71OrbitLQPeakTable(
        CoherenceBlock block,
        C2BondLQPeakScan peakScan,
        IReadOnlyList<OrbitLQPeakWitness> orbitWitnesses,
        double maxQPeakDev,
        double maxKMaxDev)
        : base("c=2 F71-orbit EP-resonance Q-peak witness table",
               Tier.Tier2Verified,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + " +
               "docs/proofs/PROOF_C1_MIRROR_SYMMETRY.md")
    {
        Block = block;
        PeakScan = peakScan;
        OrbitWitnesses = orbitWitnesses;
        MaxQPeakWithinOrbitDeviation = maxQPeakDev;
        MaxKMaxWithinOrbitDeviation = maxKMaxDev;
    }

    public override string DisplayName =>
        $"F71-orbit L(Q) peak table (N={Block.N}, {OrbitWitnesses.Count} orbits)";

    public override string Summary =>
        $"{OrbitWitnesses.Count} F71 orbits at c=2 N={Block.N}; F71-mirror invariant " +
        $"(max-Δ Q_peak = {MaxQPeakWithinOrbitDeviation:G3}, max-Δ KMax = {MaxKMaxWithinOrbitDeviation:G3}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("NumOrbits", summary: OrbitWitnesses.Count.ToString());
            yield return PeakScan;
            yield return InspectableNode.RealScalar("MaxQPeakWithinOrbitDeviation", MaxQPeakWithinOrbitDeviation, "G3");
            yield return InspectableNode.RealScalar("MaxKMaxWithinOrbitDeviation", MaxKMaxWithinOrbitDeviation, "G3");
            yield return InspectableNode.Group("OrbitWitnesses",
                OrbitWitnesses.Cast<IInspectable>().ToArray());
        }
    }
}

/// <summary>One F71-orbit's L(Q) Q-peak witnesses: parabolic-refined Q_peak, peak height
/// ‖xB‖_max, optional left/right HWHM (linear interpolation between bracketing grid
/// points), plus the within-orbit deviations recorded at construction (zero for
/// self-paired orbits, ≤ tolerance for mirror-pair orbits).</summary>
public sealed record OrbitLQPeakWitness(
    F71BondOrbit Orbit,
    double QPeak,
    double KMax,
    double? HwhmLeft,
    double? HwhmRight,
    double? HwhmLeftOverQPeak,
    double QPeakWithinOrbitDeviation,
    double KMaxWithinOrbitDeviation
) : IInspectable
{
    public string DisplayName =>
        Orbit.IsSelfPaired
            ? $"OrbitLQPeakWitness {{b={Orbit.BondA}}} (self-paired)"
            : $"OrbitLQPeakWitness {{b={Orbit.BondA} ↔ b={Orbit.BondB}}}";

    public string Summary =>
        $"Q_peak = {QPeak:F4}, ‖xB‖_max = {KMax:F4}, " +
        $"HWHM_left = {HwhmLeft?.ToString("F4") ?? "n/a"}, " +
        $"HWHM_left/Q_peak = {HwhmLeftOverQPeak?.ToString("F4") ?? "n/a"}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return Orbit;
            yield return InspectableNode.RealScalar("Q_peak", QPeak, "F4");
            yield return InspectableNode.RealScalar("‖xB‖_max", KMax, "F4");
            if (HwhmLeft is not null) yield return InspectableNode.RealScalar("HWHM_left", HwhmLeft.Value, "F4");
            if (HwhmRight is not null) yield return InspectableNode.RealScalar("HWHM_right", HwhmRight.Value, "F4");
            if (HwhmLeftOverQPeak is not null)
                yield return InspectableNode.RealScalar("HWHM_left/Q_peak", HwhmLeftOverQPeak.Value, "F4");
            yield return InspectableNode.RealScalar("ΔQ_peak within orbit", QPeakWithinOrbitDeviation, "G3");
            yield return InspectableNode.RealScalar("ΔKMax within orbit", KMaxWithinOrbitDeviation, "G3");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
