using System;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The Symphony first movement: ONE evolution, many lenses, one timeline, the events axis.
/// The assertions guard the contract — the trajectory is built exactly once and every lens reads it;
/// CΨ falls monotonically under pure dephasing (the proven dCΨ/dt &lt; 0); the events list is sorted;
/// all lens curves share the one t grid; and the N=2 Bell+ trajectory reproduces the F25 closed form
/// CΨ(t) = f(1+f²)/6, f = e^{−4γt}.</summary>
public class SymphonyTests
{
    private static System.Collections.Generic.List<IInspectable> Children(Symphony s) =>
        ((IInspectable)s).Children.ToList();

    [Fact]
    public void Constructor_RejectsBadArgs()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new Symphony(n: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Symphony(n: Symphony.MaxN + 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Symphony(gamma: 0.0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new Symphony(tPoints: 1));
    }

    [Fact]
    public void CarrierPair_DefaultsToZeroOne_AndValidates()
    {
        Assert.Equal((0, 1), new Symphony(n: 3).CarrierPair);
        Assert.Equal((1, 2), new Symphony(n: 3, carrierPair: (1, 2)).CarrierPair);
        // out of range and duplicate sites are rejected
        Assert.Throws<ArgumentOutOfRangeException>(() => new Symphony(n: 3, carrierPair: (0, 3)));
        Assert.Throws<ArgumentException>(() => new Symphony(n: 3, carrierPair: (1, 1)));
    }

    [Fact]
    public void OneEvolution_BuiltOnce_NoMatterHowManyLensesRead()
    {
        var s = new Symphony(n: 3);
        // Touch every lens, the summary, and the events repeatedly.
        _ = s.Summary;
        _ = Children(s).Select(c => c.Summary).ToList();
        _ = s.States;
        _ = s.TimeGrid;
        _ = s.FirstQuarterCrossing();
        _ = s.Summary;
        Assert.Equal(1, s.EvolveCount);   // ONE evolution, shared by all lenses
    }

    [Fact]
    public void AllLensCurves_ShareTheOneTimeGrid()
    {
        var s = new Symphony(n: 3, tPoints: 41);
        Assert.Equal(41, s.TimeGrid.Count);
        Assert.Equal(41, s.States.Count);

        // Pull the curve payloads off the lens nodes and check every one is 41 long on the same grid.
        var curves = Children(s)
            .SelectMany(LensCurves)
            .ToList();
        Assert.NotEmpty(curves);
        foreach (var curve in curves)
        {
            Assert.Equal(41, curve.X.Count);
            Assert.Equal(41, curve.Y.Count);
            for (int i = 0; i < 41; i++)
                Assert.Equal(s.TimeGrid[i], curve.X[i], 12);
        }
    }

    private static System.Collections.Generic.IEnumerable<InspectablePayload.Curve> LensCurves(IInspectable node)
    {
        if (node.Payload is InspectablePayload.Curve c) yield return c;
        foreach (var child in node.Children)
            foreach (var inner in LensCurves(child))
                yield return inner;
    }

    [Fact]
    public void Cpsi_MonotoneDecreasing_UnderPureDephasing_N2Bell()
    {
        // dCΨ/dt < 0 is proven for local Markovian channels (PROOF_MONOTONICITY_CPSI / F25). The clean
        // pure-dephasing case is N=2 Bell+: |00⟩ and |11⟩ both sit in the XY hopping kernel (the hop
        // connects |01⟩↔|10⟩ only), so the Hamiltonian does nothing and only Z-dephasing acts — CΨ
        // falls monotonically along F25. (At N=3 the |110⟩ component CAN hop to |101⟩, so Bell+ is no
        // longer a kernel state and the Hamiltonian briefly raises CΨ; that is real physics, not a
        // dephasing trajectory — so the monotone assertion belongs to the genuine pure-dephasing case.)
        var s = new Symphony(n: 2, gamma: 0.1, initialState: InitialStateKind.BellPair, tMax: 8.0, tPoints: 80);
        var cpsi = s.States.Select(Symphony.Cpsi).ToList();
        for (int i = 1; i < cpsi.Count; i++)
            Assert.True(cpsi[i] <= cpsi[i - 1] + 1e-12,
                $"CΨ rose at step {i}: {cpsi[i - 1]} -> {cpsi[i]}");
        Assert.True(cpsi[^1] < cpsi[0], "CΨ should strictly decay overall");
    }

    [Fact]
    public void Events_SortedByTime()
    {
        var s = new Symphony(n: 3);
        var eventsNode = Children(s).Single(c => c.DisplayName == "events");
        var times = eventsNode.Children
            .Select(c => ParseTime(c.DisplayName))
            .ToList();
        Assert.NotEmpty(times);
        for (int i = 1; i < times.Count; i++)
            Assert.True(times[i] >= times[i - 1], $"events out of order at {i}: {times[i - 1]} then {times[i]}");
    }

    private static double ParseTime(string displayName)
    {
        // displayName looks like "t=0.376, K=0.0376"
        var tTok = displayName.Split(',')[0].Trim();   // "t=0.376"
        return double.Parse(tTok.Substring(2), System.Globalization.CultureInfo.InvariantCulture);
    }

    [Fact]
    public void N2_BellPair_ReproducesF25ClosedForm()
    {
        // F25: CΨ(t) = f·(1+f²)/6, f = e^{−4γt}, for Bell+ under Z-dephasing at N=2. CΨ(0) = 1/3.
        double gamma = 0.1;
        var s = new Symphony(n: 2, gamma: gamma, initialState: InitialStateKind.BellPair, tMax: 6.0, tPoints: 60);
        Assert.Equal(1.0 / 3.0, Symphony.Cpsi(s.States[0]), 9);
        for (int i = 0; i < s.TimeGrid.Count; i++)
        {
            double t = s.TimeGrid[i];
            double f = Math.Exp(-4.0 * gamma * t);
            double expected = f * (1.0 + f * f) / 6.0;
            Assert.Equal(expected, Symphony.Cpsi(s.States[i]), 8);
        }
    }

    [Fact]
    public void N2_BellPair_CrossesQuarter_AtF25Dose()
    {
        // F25 crossing: f* = 0.8612, K = γ·t_cross = 0.0374. Our grid + interpolation should land close.
        var s = new Symphony(n: 2, gamma: 0.1, initialState: InitialStateKind.BellPair, tMax: 6.0, tPoints: 120);
        var cross = s.FirstQuarterCrossing();
        Assert.NotNull(cross);
        Assert.Equal(0.0374, cross!.Value.Dose, 3);   // K = γ·t_cross
    }

    [Fact]
    public void N3_BellPair_StartsBelowQuarter_NoSpuriousCrossing()
    {
        // The global Ψ-normalization /(d−1) puts CΨ(0) = 1/(d−1) = 1/7 ≈ 0.143 < ¼ at N=3 (d=8): a
        // state that begins below ¼ has NO crossing. The witness must say so honestly, not report a
        // t=0 fold.
        var s = new Symphony(n: 3, initialState: InitialStateKind.BellPair);
        Assert.Equal(1.0 / 7.0, Symphony.Cpsi(s.States[0]), 9);
        Assert.True(Symphony.Cpsi(s.States[0]) < 0.25);
        Assert.Null(s.FirstQuarterCrossing());
        Assert.DoesNotContain("¼ crossing at", s.Summary);
    }

    [Fact]
    public void Palindrome_Holds_ForTrulyXY()
    {
        // The palindrome lens reads the same MirrorSystem F1 verdict; XY is truly, so it holds.
        var s = new Symphony(n: 3, hType: HamiltonianType.XY);
        var lens = Children(s).Single(c => c.DisplayName == "lens: palindrome");
        Assert.Contains("holds = True", lens.Summary);
    }

    [Fact]
    public void LightContent_BellPair_StartsAtTwo_DiagonalIsDark()
    {
        // Bell+ on sites (0,1): coherences |00..0⟩⟨11 0..0| have popcount(i⊕j) = 2; the purity-weighted
        // mean light at t=0 is 2·(coherence weight) / (total weight). For pure Bell+ embedded the four
        // nonzero ρ entries are 1/2 each: two on the diagonal (light 0), two off (light 2) -> mean = 1.
        var s = new Symphony(n: 3, initialState: InitialStateKind.BellPair);
        Assert.Equal(1.0, Symphony.LightContent(s.States[0]), 9);
        // Light decays toward the dark {I,Z} diagonal as coherences dephase.
        Assert.True(Symphony.LightContent(s.States[^1]) < Symphony.LightContent(s.States[0]));
    }

    [Fact]
    public void Witness_SurfacesAllLensAndEventChildren()
    {
        var labels = Children(new Symphony(n: 3)).Select(c => c.DisplayName).ToList();
        Assert.Contains("score", labels);
        Assert.Contains("lens: palindrome", labels);
        Assert.Contains("lens: quarter (CΨ)", labels);
        Assert.Contains("lens: dose (K)", labels);
        Assert.Contains("lens: light", labels);
        Assert.Contains("events", labels);
    }

    [Fact]
    public void Witness_RendersToJson()
    {
        var json = InspectionJsonExporter.ToJson(new Symphony(n: 2));
        Assert.Contains("one evolution", json);
        Assert.Contains("events", json);
    }

    [Fact]
    public void LocalCpsi_AtN2_EqualsGlobalCpsi_BitExact()
    {
        // At N=2, tracing onto (0,1) keeps both qubits: the reduced state IS the full state,
        // so the local lens reproduces the global CΨ exactly (and thus F25).
        var s = new Symphony(n: 2, gamma: 0.1, initialState: InitialStateKind.BellPair, tMax: 6.0, tPoints: 40);
        foreach (var rho in s.States)
            Assert.Equal(Symphony.Cpsi(rho), s.LocalCpsi(rho), 12);
    }

    [Fact]
    public void LocalCpsi_N3BellPair_StartsAtOneThird()
    {
        // The reduced pair (0,1) of Bell+ is the Bell+ 4×4 state: local CΨ(0) = 1/3, ABOVE ¼,
        // even though the global CΨ(0) = 1/(d−1) = 1/7 is below ¼.
        var s = new Symphony(n: 3, initialState: InitialStateKind.BellPair);
        Assert.Equal(1.0 / 3.0, s.LocalCpsi(s.States[0]), 9);
    }

    [Fact]
    public void LocalCpsi_SingleExcitation_HasNoPairCoherence()
    {
        // |100⟩ reduced onto (0,1) is |10⟩⟨10|, diagonal: zero coherence, local CΨ ≡ 0.
        var s = new Symphony(n: 3, initialState: InitialStateKind.SingleExcitation);
        Assert.Equal(0.0, s.LocalCpsi(s.States[0]), 12);
    }

    [Fact]
    public void LocalCpsi_ConfigurablePair_ReadsTheChosenSites()
    {
        // Bell+ lives on (0,1): pair (0,1) sees the coherence (CΨ=1/3), pair (1,2) sees a diagonal
        // mixture (CΨ=0). The lens reads the pair it is told to.
        var onCarrier = new Symphony(n: 3, carrierPair: (0, 1), initialState: InitialStateKind.BellPair);
        var offCarrier = new Symphony(n: 3, carrierPair: (1, 2), initialState: InitialStateKind.BellPair);
        Assert.Equal(1.0 / 3.0, onCarrier.LocalCpsi(onCarrier.States[0]), 9);
        Assert.Equal(0.0, offCarrier.LocalCpsi(offCarrier.States[0]), 9);
    }

    [Fact]
    public void QuarterCrossings_TagDirection_DownThenUp()
    {
        // A hand-built curve that dips below ¼ and recovers ABOVE ¼ must yield one down then one up crossing.
        double[] curve = { 0.40, 0.30, 0.20, 0.18, 0.30, 0.40 };  // down (idx 1→2), up (idx 3→4)
        var dirs = Symphony.QuarterCrossingDirections(curve);
        Assert.Equal(new[] { -1, +1 }, dirs);
    }
}
