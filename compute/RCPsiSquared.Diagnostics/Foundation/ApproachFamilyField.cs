using System.Globalization;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Object Manager's view of a named family of free two-qubit Z-dephasing curves: how the
/// scalar-quarter crossing depends on the start. Sweeps the partial-entanglement initial state |ψ(α)⟩ = cosα|00⟩ +
/// sinα|11⟩ over s = sin2α and reads the family through <see cref="OddHarmonicApproach"/>. Five readings:
/// the shared slow exponential (4γ within this family); the family of starts (s is the initial
/// pure-state concurrence while CΨ(0)=s/3 is one third of it); the crossing threshold (a temporal
/// downward crossing occurs iff γ &gt; 0 and s &gt; 3/4); the shape parameter
/// (for s&gt;0, the cubic-term fraction s²/2, continuously extended to 0 at s=0); and the late-time exponent
/// (every nonzero member has a 4γ late-time term, while the 12γ term is transient). Closed-form,
/// N-free; the Bell+ member reproduces F25. Wired into the typed graph (Core) as
/// <c>ApproachFamilyCarrierClaim</c> with the Absorption Theorem and F25 as its two typed
/// parents; C2 and Two Readings remain non-ancestral prose comparisons. Render it with
/// <c>inspect --claim ApproachFamilyCarrierClaim</c>.</summary>
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
        if (sLo < 0 || sHi > 1.0 || sHi <= sLo)
            throw new ArgumentOutOfRangeException(nameof(sHi), $"need 0 ≤ sLo < sHi ≤ 1; got [{sLo}, {sHi}]");
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
            int nCross = _sLadder.Count(s => OddHarmonicApproach.HasDownwardCrossing(s, _gamma));
            return $"|ψ(α)⟩=cosα|00⟩+sinα|11⟩ swept over s=sin2α: s is the initial pure-state concurrence, " +
                   $"and CΨ(0)=s/3 is exactly one third of it. A temporal downward crossing occurs iff γ>0 and s>3/4 " +
                   $"({nCross}/{_sLadder.Length} members cross), cubic-term fraction s²/2 for s>0 (Bell+ ½; continuous value 0 at s=0); every nonzero member " +
                   $"(0<s≤1) contains the late-time 4γ carrier exponential={OddHarmonicApproach.CarrierRate(_gamma).ToString("0.###", Inv)}, " +
                   "while s=0 has w₀=w₁=0 and no late-time term; " +
                   "at γ=0 every curve is constant and none crosses. The weights set the shape and the cubic 12γ contribution fades first.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // 1. The shared late-time exponential in this named family.
            yield return new InspectableNode(
                displayName: "the shared slow exponential (4γ within this family)",
                summary: $"every nonzero member contains w₀·e^(−4γt), with rate 4γ={OddHarmonicApproach.CarrierRate(_gamma).ToString("0.###", Inv)}. " +
                         $"The 12γ={OddHarmonicApproach.HarmonicRate(_gamma).ToString("0.###", Inv)} term dies earlier. " +
                         "At s=0 both weights vanish, so that endpoint has no late-time term. " +
                         "This is a closed-form property of this setup, not an additional spectral-object claim.");

            // 2. The family of starts: CΨ(0) = s/3.
            var starts = _sLadder.Select(OddHarmonicApproach.InitialCpsi).ToArray();
            yield return new InspectableNode(
                displayName: "the family of starts (CΨ(0) = s/3)",
                summary: $"s is the initial pure-state concurrence; the distinct linear CΨ readout is exactly one third of it: CΨ(0)=s/3, from {starts[0].ToString("0.###", Inv)} to {starts[^1].ToString("0.###", Inv)} across the sweep.",
                payload: new InspectablePayload.Curve("CΨ(0) vs initial concurrence s", _sLadder, starts, "s = sin2α (initial pure-state concurrence)", "CΨ(0) = s/3"));

            // 3. The crossing threshold (s = 3/4): CΨ(0) − ¼ crosses zero there (no NaN, unlike t_cross).
            var margin = _sLadder.Select(s => OddHarmonicApproach.InitialCpsi(s) - OddHarmonicApproach.Cusp).ToArray();
            double bellCross = OddHarmonicApproach.CrossingTime(OddHarmonicApproach.BellPlusS, _gamma);
            yield return new InspectableNode(
                displayName: "the temporal downward-crossing threshold (γ > 0 and s > 3/4)",
                summary: $"a genuine temporal downward crossing occurs iff γ>0 and s>3/4. At s=3/4 the curve touches ¼ at t=0 but does not cross; below it the curve starts below ¼. At γ=0 CΨ(t)=s/3 is constant and does not cross. For this positive-γ sweep, the Bell+ member crosses at t={bellCross.ToString("0.###", Inv)}.",
                payload: new InspectablePayload.Curve("CΨ(0) − ¼ vs s", _sLadder, margin, "s = sin2α", "CΨ(0) − ¼ (zero at 3/4)"));

            // 4. The shape parameter (the cubic-term fraction = s²/2).
            var harm = _sLadder.Select(OddHarmonicApproach.HarmonicFraction).ToArray();
            yield return new InspectableNode(
                displayName: "the shape parameter (cubic-term fraction s²/2)",
                summary: $"for every s>0 the cubic 12γ term is present and carries the relative weight s²/2 of the start: {harm[0].ToString("0.###", Inv)} to {harm[^1].ToString("0.###", Inv)}. This share grows quadratically and becomes appreciable as s grows; the plotted shape parameter has continuous value 0 at s=0, where the total start is zero. Bell+ (s=1) is the 50/50 member.",
                payload: new InspectablePayload.Curve("cubic-term fraction vs s", _sLadder, harm, "s = sin2α", "w₁/(w₀+w₁) = s²/2 for s>0"));

            // 5. The exact Bell+ specialization.
            var tGrid = TimeGrid();
            var bell = tGrid.Select(t => OddHarmonicApproach.Cpsi(OddHarmonicApproach.BellPlusS, _gamma, t)).ToArray();
            yield return new InspectableNode(
                displayName: "the Bell+ specialization (exactly F25)",
                summary: "the Bell+ member is (1/6)e^(−4γt)+(1/6)e^(−12γt)=F25. Its 4γ term controls late time and its 12γ term is transient; no broader slowing mechanism is inferred.",
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
