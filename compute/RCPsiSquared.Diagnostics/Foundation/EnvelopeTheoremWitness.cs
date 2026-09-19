using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live finite rise atlas associated with the historically named
/// <c>CpsiEnvelopeTheoremClaim</c>. It reports what selected N/Q/K windows resolve. A zero sampled
/// rise is not an absence theorem, and a positive sampled rise is numerical evidence only after the
/// stated reporting bar and refinement controls are applied.
///
/// <para>It reuses <see cref="Symphony"/> as the evolve-CΨ engine and <see cref="QuarterEnvelope.Of"/>
/// to read the envelope, exactly as the Symphony tests do — it does not re-implement the propagation.
/// The named J=5, γ=0.01, tMax=25 grid resolves above-bar carrier-pair rises.</para>
///
/// <para>Guard: N in 3..<see cref="Symphony.MaxN"/>. At N=2 the carrier pair is the full state, so
/// local and global readings coincide. That identity is all this guard asserts: the autonomous N=2
/// successive-local-maxima question remains unproved, and the five exact examples do not refute it.
/// The finite global/local comparison needs a bath.</para>
///
/// <para>Children: the finite global reading, the local carrier-pair reading, and the named-state
/// triptych control (SingleExcitation = two named grid readings; BondingMode = initial H-eigenstate
/// plus no above-bar rise on its named grid; Bell+ = above-bar reading on the named grid), plus the local CΨ(t) curve
/// payload.</para></summary>
public sealed class EnvelopeTheoremWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>A reporting convention for the finite atlas. A predecessor rise above this value is
    /// resolved on this named grid; refinement comparisons are reported separately. The bar does not
    /// separate physical signal from noise.</summary>
    public const double RiseReportingBar = 1e-3;

    private const double JStrong = 5.0, GammaStrong = 0.01, TMaxStrong = 25.0;
    private const int FinePoints = 1600, CoarsePoints = 400;

    public int N { get; }

    public EnvelopeTheoremWitness(int n = 3)
    {
        if (n < 3 || n > Symphony.MaxN)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"the witness needs N in 3..{Symphony.MaxN}: at N=2 the carrier pair IS the full state " +
                "(the partial trace is the identity), so local ≡ global. The finite global/local comparison " +
                $"needs a larger system; the autonomous N=2 peak question remains unproved. Got {n}.");
        N = n;
    }

    private bool _built;
    private EnvelopeReading _globalBell, _localBell, _singleCoarse, _singleFine, _localBonding;
    private double[] _payloadT = Array.Empty<double>(), _payloadCurve = Array.Empty<double>();

    private void Ensure()
    {
        if (_built) return;
        _built = true;

        var bell = new Symphony(n: N, j: JStrong, gamma: GammaStrong,
            initialState: InitialStateKind.BellPair, tMax: TMaxStrong, tPoints: FinePoints);
        var t = bell.TimeGrid.ToArray();
        var localCurve = bell.States.Select(bell.LocalCpsi).ToArray();
        _globalBell = GlobalEnvelope(bell);   // the same core the static GlobalReading uses — no drift
        _localBell = QuarterEnvelope.Of(localCurve, t, riseTol: RiseReportingBar);
        _payloadT = t;
        _payloadCurve = localCurve;

        _singleCoarse = LocalReading(InitialStateKind.SingleExcitation, CoarsePoints);
        _singleFine = LocalReading(InitialStateKind.SingleExcitation, FinePoints);
        _localBonding = LocalReading(InitialStateKind.BondingMode, FinePoints);
    }

    private EnvelopeReading LocalReading(InitialStateKind init, int points)
    {
        var s = new Symphony(n: N, j: JStrong, gamma: GammaStrong,
            initialState: init, tMax: TMaxStrong, tPoints: points);
        return QuarterEnvelope.Of(s.States.Select(s.LocalCpsi).ToArray(), s.TimeGrid.ToArray(),
            riseTol: RiseReportingBar);
    }

    /// <summary>The full-state (global) CΨ envelope on an already-evolved engine — the exact detector the
    /// witness reads for its finite global child. A shared core so the live witness and the atlas sweep
    /// cannot drift.</summary>
    private static EnvelopeReading GlobalEnvelope(Symphony bell) =>
        QuarterEnvelope.Of(bell.States.Select(Symphony.Cpsi).ToArray(), bell.TimeGrid.ToArray(),
            riseTol: RiseReportingBar);

    /// <summary>The global CΨ envelope reading for a Bell+ carrier at (<paramref name="n"/>,
    /// <paramref name="j"/>, <paramref name="gamma"/>) over [0, <paramref name="tMax"/>] on
    /// <paramref name="points"/> grid points — the witness's own detector, parameterised so the
    /// envelope_n4_rise finite atlas reuses it verbatim instead of re-deriving it.</summary>
    public static EnvelopeReading GlobalReading(int n, double j, double gamma, double tMax, int points) =>
        GlobalEnvelope(new Symphony(n: n, j: j, gamma: gamma,
            initialState: InitialStateKind.BellPair, tMax: tMax, tPoints: points));

    /// <summary>The number of predecessor rises exceeding <see cref="RiseReportingBar"/> in this
    /// finite global CΨ window. Zero means none was resolved on this grid, not that none exists.</summary>
    public static int GlobalRiseCount(int n, double j, double gamma, double tMax, int points) =>
        GlobalReading(n, j, gamma, tMax, points).RiseCount;

    public EnvelopeReading GlobalBell { get { Ensure(); return _globalBell; } }
    public EnvelopeReading LocalBell { get { Ensure(); return _localBell; } }
    public EnvelopeReading SingleCoarse { get { Ensure(); return _singleCoarse; } }
    public EnvelopeReading SingleFine { get { Ensure(); return _singleFine; } }
    public EnvelopeReading LocalBonding { get { Ensure(); return _localBonding; } }

    public string DisplayName =>
        $"EnvelopeTheoremWitness (finite rise atlas, N={N}, J={JStrong.ToString("0.#", Inv)}, γ={GammaStrong.ToString("0.###", Inv)})";

    public string Summary
    {
        get
        {
            Ensure();
            string globalClause = _globalBell.RiseCount == 0
                ? $"At N={N}, no global predecessor rise above the reporting bar was resolved (RiseCount = 0)"
                : $"At N={N}, this finite grid resolves {_globalBell.RiseCount} global predecessor rises above the reporting bar";
            return $"finite CΨ rise atlas (historical typed home: CpsiEnvelopeTheoremClaim): {globalClause}. " +
                   "That finite reading neither proves absence nor decides the still-unproved autonomous N=2 peak-sequence claim. " +
                   "The reduced carrier pair has no general monotonicity guarantee; " +
                   $"this named grid reports local RiseCount = {_localBell.RiseCount} and raw max positive Δ = " +
                   $"{_localBell.MaxRiseMagnitude.ToString("0.#####", Inv)} against the reporting bar " +
                   $"{RiseReportingBar.ToString("0.###", Inv)}. The named grids invite a wider atlas; they do not draw a theorem boundary.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            Ensure();

            yield return new InspectableNode("the finite global reading",
                summary: $"full-state CΨ at N={N}, J={JStrong.ToString("0.#", Inv)}, γ={GammaStrong.ToString("0.###", Inv)}, " +
                         $"1600 pts: predecessor rises above {RiseReportingBar.ToString("0.###", Inv)} = {_globalBell.RiseCount}; " +
                         $"max Δ = {_globalBell.MaxRiseMagnitude.ToString("0.#####", Inv)}. " +
                         "No-rise rows are finite null samples; rise rows are finite numerical evidence. " +
                         "PROOF_MONOTONICITY_CPSI retracts the old universal package and leaves the autonomous N=2 peak question unproved.");

            yield return new InspectableNode("the local carrier-pair reading",
                summary: $"the reduced carrier-pair CΨ has no theorem; this named grid reports RiseCount = " +
                         $"{_localBell.RiseCount}, max Δ = {_localBell.MaxRiseMagnitude.ToString("0.#####", Inv)} " +
                         $"against the reporting bar {RiseReportingBar.ToString("0.###", Inv)}. Compare a separately evolved refined grid. " +
                         "This is one named N/Q/K window" +
                         (N < Symphony.MaxN ? $"; try --N {N + 1} to add another atlas row." : "."));

            // The grid comparison is N-dependent: at N=3 SingleExcitation's sub-bar rise is absent on
            // the finer grid; at N=4 the named localized-state row also rises above the bar. Report which
            // finite case is live rather than hardcoding N=3 or inferring a mechanism.
            bool singleBelowBarOnBoth = _singleCoarse.RiseCount == 0 && _singleFine.RiseCount == 0;
            string singleClause = singleBelowBarOnBoth
                ? $"SingleExcitation = sub-bar grid readings (max Δ {_singleCoarse.MaxRiseMagnitude.ToString("0.#####", Inv)} at 400 pts " +
                  $"and {_singleFine.MaxRiseMagnitude.ToString("0.#####", Inv)} at 1600 pts, both < bar; RiseCount " +
                  $"{_singleFine.RiseCount} at 1600)"
                : $"SingleExcitation = an above-bar finite row at this N (max Δ {_singleFine.MaxRiseMagnitude.ToString("0.#####", Inv)} at 1600 pts " +
                  $"{(_singleFine.MaxRiseMagnitude >= RiseReportingBar ? "> bar" : "< bar")}, RiseCount {_singleFine.RiseCount} at 1600: " +
                  "this is only the stated grid classification)";
            yield return new InspectableNode("the named-state control (above bar / two grids / H-eigenstate)",
                summary: $"Bell+ = rise above the reporting bar on this named grid (max Δ {_localBell.MaxRiseMagnitude.ToString("0.#####", Inv)} > bar); " +
                         $"{singleClause}; " +
                         $"BondingMode = initial H-eigenstate with no rise above the reporting bar on its named grid " +
                         $"(RiseCount {_localBonding.RiseCount}). " +
                         "These are finite numerical classifications; the two-grid comparison supplies no cause or convergence claim, " +
                         "the reporting bar does not separate physical signal from noise, and no general state-class law is inferred.");

            yield return new InspectableNode("the local CΨ(t) finite curve",
                summary: $"the carrier-pair CΨ over t for Bell+ in the named J={JStrong.ToString("0.#", Inv)}, " +
                         $"γ={GammaStrong.ToString("0.###", Inv)} finite window.",
                payload: new InspectablePayload.Curve("local CΨ(t)", _payloadT, _payloadCurve, "t", "local CΨ"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
