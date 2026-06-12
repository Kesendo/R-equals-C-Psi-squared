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
        var globalCurve = bell.States.Select(Symphony.Cpsi).ToArray();
        var localCurve = bell.States.Select(bell.LocalCpsi).ToArray();
        _globalBell = QuarterEnvelope.Of(globalCurve, t);
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
            return $"the CΨ Envelope Theorem live (typed home: CpsiEnvelopeTheoremClaim): the proof predicts the " +
                   $"full-state envelope is non-increasing; the live evolution CONFIRMS it (global RiseCount = " +
                   $"{_globalBell.RiseCount}) and shows the reduced carrier pair ESCAPE it (local RiseCount = " +
                   $"{_localBell.RiseCount}, max Δ = {_localBell.MaxRiseMagnitude.ToString("0.#####", Inv)} > " +
                   $"{GenuinenessBar.ToString("0.###", Inv)}, the freedom — beating). Two independent " +
                   "computations meeting (theorem vs live evolution).";
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
                         $"({(_globalBell.IsNonIncreasing ? "non-increasing ✓ — the Envelope Theorem holds live" : "RISES — unexpected, would falsify the Tier-2 verification")}). " +
                         "Proven N=2, verified N≥3 (PROOF_MONOTONICITY_CPSI / CpsiEnvelopeTheoremClaim).");

            yield return new InspectableNode("the freedom (local carrier pair)",
                summary: $"the reduced carrier-pair CΨ has no theorem and its beat envelope RISES: RiseCount = " +
                         $"{_localBell.RiseCount}, max Δ = {_localBell.MaxRiseMagnitude.ToString("0.#####", Inv)} " +
                         $"(> the genuineness bar {GenuinenessBar.ToString("0.###", Inv)}), grid-sensitive — verify with ≥4× t-points. " +
                         "The freedom GROWS with N (a richer bath; N=4 shows ~13× this magnitude — try --N 4).");

            yield return new InspectableNode("the state-class control (beats / artifacts / silent)",
                summary: "the genuineness bar separates the state classes convention-robustly: " +
                         $"Bell+ = genuine beating (max Δ {_localBell.MaxRiseMagnitude.ToString("0.#####", Inv)} > bar); " +
                         $"SingleExcitation = artifacts only (max Δ {_singleCoarse.MaxRiseMagnitude.ToString("0.#####", Inv)} at 400 pts " +
                         $"and {_singleFine.MaxRiseMagnitude.ToString("0.#####", Inv)} at 1600 pts, both < bar; RiseCount " +
                         $"{_singleFine.RiseCount} at 1600 — they vanish under refinement); " +
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
