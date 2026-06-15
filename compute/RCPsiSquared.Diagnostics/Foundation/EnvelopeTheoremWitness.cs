using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the CΨ Envelope Theorem (the typed claim is
/// <c>CpsiEnvelopeTheoremClaim</c>). The proof PREDICTS the full-state envelope is non-increasing; this
/// witness CONFIRMS it live and shows the theorem's BOUNDARY — the reduced carrier pair escapes it (its
/// beat envelope rises, the freedom). Two independent computations meeting: the theorem vs the live
/// evolution.
///
/// <para>It reuses <see cref="Symphony"/> as the evolve-CΨ engine and <see cref="QuarterEnvelope.Of"/>
/// to read the envelope, exactly as the Symphony tests do — it does not re-implement the propagation.
/// Strong-coupling regime J=5, γ=0.01, tMax=25, where the freedom is loud.</para>
///
/// <para>Guard: N in 3..<see cref="Symphony.MaxN"/>. NOT N=2: there the carrier pair IS the full state
/// (the partial trace is the identity), so local ≡ global and the freedom vanishes by construction, not
/// physics. The N=2-proven case is the claim's domain; the freedom needs a bath.</para>
///
/// <para>Children: the theorem (global, RiseCount 0), the freedom (local pair, rises above the
/// genuineness bar; grows with N), and the state-class triptych control (SingleExcitation = sub-bar
/// artifacts that vanish under refinement; BondingMode = silent H-eigenstate; vs Bell+ = genuine), plus
/// the local CΨ(t) curve payload.</para></summary>
public sealed class EnvelopeTheoremWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>A real beating rise exceeds this; grid-clip artifacts sit far below (genuine ≥ 8.9e-3,
    /// SingleExcitation artifact ≈ 5.5e-4). The convention-robust separator between signal and noise.</summary>
    public const double GenuinenessBar = 1e-3;

    private const double JStrong = 5.0, GammaStrong = 0.01, TMaxStrong = 25.0;
    private const int FinePoints = 1600, CoarsePoints = 400;

    public int N { get; }

    public EnvelopeTheoremWitness(int n = 3)
    {
        if (n < 3 || n > Symphony.MaxN)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"the witness needs N in 3..{Symphony.MaxN}: at N=2 the carrier pair IS the full state " +
                "(the partial trace is the identity), so local ≡ global and the freedom vanishes by " +
                "construction, not physics — the identity case. The N=2-proven theorem is the claim's " +
                $"domain; the freedom needs a bath (N≥3). Got {n}.");
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
        _localBell = QuarterEnvelope.Of(localCurve, t);
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
        return QuarterEnvelope.Of(s.States.Select(s.LocalCpsi).ToArray(), s.TimeGrid.ToArray());
    }

    /// <summary>The full-state (global) CΨ envelope on an already-evolved engine — the exact detector the
    /// witness reads for its theorem child. A shared core so the live witness and the boundary sweep
    /// cannot drift.</summary>
    private static EnvelopeReading GlobalEnvelope(Symphony bell) =>
        QuarterEnvelope.Of(bell.States.Select(Symphony.Cpsi).ToArray(), bell.TimeGrid.ToArray());

    /// <summary>The global CΨ envelope reading for a Bell+ carrier at (<paramref name="n"/>,
    /// <paramref name="j"/>, <paramref name="gamma"/>) over [0, <paramref name="tMax"/>] on
    /// <paramref name="points"/> grid points — the witness's own detector, parameterised so the
    /// envelope_n4_rise boundary sweep reuses it verbatim instead of re-deriving it.</summary>
    public static EnvelopeReading GlobalReading(int n, double j, double gamma, double tMax, int points) =>
        GlobalEnvelope(new Symphony(n: n, j: j, gamma: gamma,
            initialState: InitialStateKind.BellPair, tMax: tMax, tPoints: points));

    /// <summary>The number of predecessor-rises in the global CΨ envelope at (N, J, γ): 0 means the
    /// full-state envelope is non-increasing (the N=2 theorem's behaviour), &gt;0 means it RISES.</summary>
    public static int GlobalRiseCount(int n, double j, double gamma, double tMax, int points) =>
        GlobalReading(n, j, gamma, tMax, points).RiseCount;

    public EnvelopeReading GlobalBell { get { Ensure(); return _globalBell; } }
    public EnvelopeReading LocalBell { get { Ensure(); return _localBell; } }
    public EnvelopeReading SingleCoarse { get { Ensure(); return _singleCoarse; } }
    public EnvelopeReading SingleFine { get { Ensure(); return _singleFine; } }
    public EnvelopeReading LocalBonding { get { Ensure(); return _localBonding; } }

    public string DisplayName =>
        $"EnvelopeTheoremWitness (the theorem live, N={N}, J={JStrong.ToString("0.#", Inv)}, γ={GammaStrong.ToString("0.###", Inv)})";

    public string Summary
    {
        get
        {
            Ensure();
            string globalClause = _globalBell.IsNonIncreasing
                ? $"At N={N} the live evolution CONFIRMS it (global RiseCount = 0, non-increasing)"
                : $"At N={N} the live evolution REFUTES it (global RiseCount = {_globalBell.RiseCount} rises): a candidate " +
                  "falsification of the over-broad 'verified N=3-5' paraphrase (the proof is Tier-1 for the 2-qubit case only)";
            return $"the CΨ Envelope Theorem live (typed home: CpsiEnvelopeTheoremClaim): the proof predicts the " +
                   $"full-state envelope is non-increasing. {globalClause}. The reduced carrier pair has no theorem and " +
                   $"its beat envelope rises (local RiseCount = {_localBell.RiseCount}, max Δ = " +
                   $"{_localBell.MaxRiseMagnitude.ToString("0.#####", Inv)} > {GenuinenessBar.ToString("0.###", Inv)}, the " +
                   "freedom, beating). Two independent computations meeting (theorem vs live evolution).";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            Ensure();

            yield return new InspectableNode("the theorem (global, live)",
                summary: $"full-state CΨ at N={N}, J={JStrong.ToString("0.#", Inv)}, γ={GammaStrong.ToString("0.###", Inv)}, " +
                         $"1600 pts: envelope RiseCount = {_globalBell.RiseCount} " +
                         $"({(_globalBell.IsNonIncreasing ? "non-increasing ✓ — consistent with the N=2 theorem (N≥3 is not guaranteed)" : "RISES — the N≥3 full-state envelope is OPEN; at N≥4 strong coupling the internal J-coupling is the Part-6 coherence injector, NOT a falsification of the N=2 theorem")}). " +
                         "Proven Tier-1 for N=2; N≥3 OPEN (PROOF_MONOTONICITY_CPSI Part 5 / CpsiEnvelopeTheoremClaim). " +
                         "The rise boundary is charted (EnvelopeBoundaryTests): an N≥4 floor, Q_c(4)≈27, Q_c(5)≈45.");

            yield return new InspectableNode("the freedom (local carrier pair)",
                summary: $"the reduced carrier-pair CΨ has no theorem and its beat envelope RISES: RiseCount = " +
                         $"{_localBell.RiseCount}, max Δ = {_localBell.MaxRiseMagnitude.ToString("0.#####", Inv)} " +
                         $"(> the genuineness bar {GenuinenessBar.ToString("0.###", Inv)}), grid-sensitive: verify with ≥4× t-points. " +
                         "The freedom grows with N (a richer bath)" +
                         (N < Symphony.MaxN ? $"; try --N {N + 1} for a larger rise." : "."));

            // The artifact control is N-dependent: at N=3 SingleExcitation's apparent rises are sub-bar grid
            // noise that vanish under refinement; at N=4 the richer bath lets even the localized state beat,
            // so the bar no longer cleanly separates it. Report which case is live rather than hardcoding N=3.
            bool singleIsArtifact = _singleCoarse.MaxRiseMagnitude < GenuinenessBar
                                    && _singleFine.MaxRiseMagnitude < GenuinenessBar
                                    && _singleFine.RiseCount == 0;
            string singleClause = singleIsArtifact
                ? $"SingleExcitation = artifacts only (max Δ {_singleCoarse.MaxRiseMagnitude.ToString("0.#####", Inv)} at 400 pts " +
                  $"and {_singleFine.MaxRiseMagnitude.ToString("0.#####", Inv)} at 1600 pts, both < bar; RiseCount " +
                  $"{_singleFine.RiseCount} at 1600: they vanish under refinement)"
                : $"SingleExcitation = ALSO beating at this N (max Δ {_singleFine.MaxRiseMagnitude.ToString("0.#####", Inv)} at 1600 pts " +
                  $"{(_singleFine.MaxRiseMagnitude >= GenuinenessBar ? "> bar" : "< bar")}, RiseCount {_singleFine.RiseCount} at 1600: " +
                  "the richer bath lets even the localized state beat, so the bar no longer cleanly separates it here)";
            yield return new InspectableNode("the state-class control (beats / artifacts / silent)",
                summary: $"Bell+ = genuine beating (max Δ {_localBell.MaxRiseMagnitude.ToString("0.#####", Inv)} > bar); " +
                         $"{singleClause}; " +
                         $"BondingMode = silent (RiseCount {_localBonding.RiseCount}: an H-eigenstate's pair CΨ decays without beating). " +
                         "Genuine rise needs a state that is not an H-mode and survives refinement.");

            yield return new InspectableNode("the local CΨ(t) (the beating)",
                summary: $"the carrier-pair CΨ over t for Bell+ at the strong-coupling regime; the late-time beat " +
                         "envelope is where the freedom's rises sit.",
                payload: new InspectablePayload.Curve("local CΨ(t)", _payloadT, _payloadCurve, "t", "local CΨ"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
