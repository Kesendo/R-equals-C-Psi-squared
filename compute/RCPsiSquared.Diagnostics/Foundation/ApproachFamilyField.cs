using System.Globalization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's telescope onto the family of approach shapes: how the approach to the
/// cusp ¼ depends on the start. Sweeps the partial-entanglement initial state |ψ(α)⟩ = cosα|00⟩ +
/// sinα|11⟩ over s = sin2α and reads the family through <see cref="OddHarmonicApproach"/>. Five readings:
/// the carrier (the slowest rate 4γ, shared by every member); the family of starts (CΨ(0)=s/3, the start
/// height is the entanglement); the crossing threshold (crosses ¼ iff s &gt; 3/4); the shape parameter
/// (the harmonic fraction s²/2, the fast-mode content growing quadratically); and the carrier collapse
/// (every member → e^(−4γt) late, the shape is the early harmonic, the slowing is ours). Closed-form,
/// N-free; the Bell+ member reproduces F25.</summary>
public sealed class ApproachFamilyField : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private readonly double _gamma;
    private readonly double[] _sLadder;
    private readonly double _tMaxFactor;

    public ApproachFamilyField(double gamma = 0.5, double sLo = 0.3, double sHi = 1.0,
        int sPoints = 8, double tMaxFactor = 6.0)
    {
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        if (sLo <= 0 || sHi > 1.0 || sHi <= sLo)
            throw new ArgumentOutOfRangeException(nameof(sHi), $"need 0 < sLo < sHi ≤ 1; got [{sLo}, {sHi}]");
        if (sPoints < 2) throw new ArgumentOutOfRangeException(nameof(sPoints), $"need at least two s points; got {sPoints}");
        if (tMaxFactor <= 1.0) throw new ArgumentOutOfRangeException(nameof(tMaxFactor), $"tMaxFactor must exceed 1; got {tMaxFactor}");
        _gamma = gamma;
        _tMaxFactor = tMaxFactor;
        _sLadder = new double[sPoints];
        for (int i = 0; i < sPoints; i++) _sLadder[i] = sLo + (sHi - sLo) * i / (sPoints - 1);
    }

    public string DisplayName =>
        $"ApproachFamilyField (the family of approach shapes, γ={_gamma.ToString("0.###", Inv)}, s∈[{_sLadder[0].ToString("0.##", Inv)}, {_sLadder[^1].ToString("0.##", Inv)}], {_sLadder.Length} members)";

    public string Summary
    {
        get
        {
            int nCross = _sLadder.Count(OddHarmonicApproach.Crosses);
            return $"|ψ(α)⟩=cosα|00⟩+sinα|11⟩ swept over s=sin2α: CΨ(0)=s/3, crosses ¼ iff s>3/4 " +
                   $"({nCross}/{_sLadder.Length} members cross), harmonic fraction s²/2 (Bell+ ½); every member shares " +
                   $"the carrier 4γ={OddHarmonicApproach.CarrierRate(_gamma).ToString("0.###", Inv)} and collapses onto it. The shape is the early harmonic; the carrier is universal.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1. The carrier (the shared, slowest mode).
            yield return new InspectableNode(
                displayName: "the carrier (universal, the slowest mode)",
                summary: $"every member ends as w₀·e^(−4γt), the carrier rate 4γ={OddHarmonicApproach.CarrierRate(_gamma).ToString("0.###", Inv)} (the HD=1 mode, the eigenvalue −γ₀). The fast 12γ={OddHarmonicApproach.HarmonicRate(_gamma).ToString("0.###", Inv)} harmonic dies early; the carrier is what survives, shared by the whole family.");

            // 2. The family of starts: CΨ(0) = s/3.
            var starts = _sLadder.Select(OddHarmonicApproach.InitialCpsi).ToArray();
            yield return new InspectableNode(
                displayName: "the family of starts (CΨ(0) = s/3)",
                summary: $"the start height is the entanglement: CΨ(0)=s/3 from {starts[0].ToString("0.###", Inv)} to {starts[^1].ToString("0.###", Inv)} across the sweep.",
                payload: new InspectablePayload.Curve("CΨ(0) vs s", _sLadder, starts, "s = sin2α", "CΨ(0)"));

            // 3. The crossing threshold (s = 3/4): CΨ(0) − ¼ crosses zero there (no NaN, unlike t_cross).
            var margin = _sLadder.Select(s => OddHarmonicApproach.InitialCpsi(s) - OddHarmonicApproach.Cusp).ToArray();
            double bellCross = OddHarmonicApproach.CrossingTime(OddHarmonicApproach.BellPlusS, _gamma);
            yield return new InspectableNode(
                displayName: "the crossing threshold (s > 3/4)",
                summary: $"the approach reaches ¼ only above the entanglement threshold s=3/4 (CΨ(0)>¼; the curve CΨ(0)−¼ crosses zero there). s=3/4 starts exactly on the cusp; below it the state never reaches ¼. The Bell+ member crosses at t={bellCross.ToString("0.###", Inv)}.",
                payload: new InspectablePayload.Curve("CΨ(0) − ¼ vs s", _sLadder, margin, "s = sin2α", "CΨ(0) − ¼ (zero at 3/4)"));

            // 4. The shape parameter (the harmonic fraction = s²/2).
            var harm = _sLadder.Select(OddHarmonicApproach.HarmonicFraction).ToArray();
            yield return new InspectableNode(
                displayName: "the shape parameter (harmonic fraction s²/2)",
                summary: $"the fast 12γ mode carries a fraction s²/2 of the start, growing quadratically: {harm[0].ToString("0.###", Inv)} to {harm[^1].ToString("0.###", Inv)}. Only strong entanglement excites it; Bell+ (s=1) is the 50/50 member.",
                payload: new InspectablePayload.Curve("harmonic fraction vs s", _sLadder, harm, "s = sin2α", "w₁/(w₀+w₁) = s²/2"));

            // 5. The carrier collapse (the Bell+ member's two-exponential), the slowing-is-ours tie.
            var tGrid = TimeGrid();
            var bell = tGrid.Select(t => OddHarmonicApproach.Cpsi(OddHarmonicApproach.BellPlusS, _gamma, t)).ToArray();
            yield return new InspectableNode(
                displayName: "the carrier collapse (the slowing is ours)",
                summary: "every member runs parallel to the carrier 4γ at late time and collapses onto it; the shape is the early 12γ harmonic transient. The Bell+ member is (1/6)e^(−4γt)+(1/6)e^(−12γt) = F25. Same reading as spiral_slowing: the carrier is steady, the slowing is the observable's.",
                payload: new InspectablePayload.Curve("Bell+ member CΨ(t)", tGrid, bell, "t", "CΨ (crosses ¼)"));
        }
    }

    private double[] TimeGrid()
    {
        double tMax = _tMaxFactor * OddHarmonicApproach.CrossingTime(OddHarmonicApproach.BellPlusS, _gamma);
        const int points = 200;
        var grid = new double[points];
        for (int i = 0; i < points; i++) grid[i] = tMax * i / (points - 1);
        return grid;
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
