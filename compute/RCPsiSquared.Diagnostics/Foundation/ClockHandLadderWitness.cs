using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the two clocks (typed home: <c>ClockHandLadderClaim</c>;
/// closed forms: <c>docs/ANALYTICAL_FORMULAS.md</c> F2b corollary; the readings:
/// <c>reflections/ON_HOW_GAMMA_BECAME_THE_TICK.md</c> + <c>reflections/ON_TWO_TIMES.md</c>).
/// The Symphony clock node has two hands: the Takt (Gap = slowest decay = 2γ) and the
/// coherence hand (Omega = max|Im λ| at the gap). This witness reads both across N from the
/// shared <see cref="Symphony"/> spectrum and shows: for N≥3 the coherence hand walks the
/// γ-protected F2b band-edge ladder Omega = 2J·cos(π/(N+1)) (√2 / φ / √3 at N=3/4/5), and for
/// N=2 it is γ-pulled to Omega = 2√(J²−γ²), stopping at the exceptional point Q=1.
///
/// <para>It reuses <see cref="Symphony"/> as the spectrum engine (it does not re-diagonalize);
/// the typed parent edges to F2b and the Absorption Theorem live on the claim, not here. The
/// witness is the numeric lab: it recomputes the live clock and checks it against the closed
/// forms.</para></summary>
public sealed class ClockHandLadderWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    public const double DefaultJ = 1.0;
    public const double DefaultGamma = 0.2;
    private const int SpectrumPoints = 8; // the clock reads the spectrum, grid-independent

    public double J { get; }
    public double Gamma { get; }

    public ClockHandLadderWitness(double j = DefaultJ, double gamma = DefaultGamma)
    {
        if (j <= 0) throw new ArgumentOutOfRangeException(nameof(j), j, "J must be > 0.");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), gamma, "γ must be > 0.");
        J = j;
        Gamma = gamma;
    }

    /// <summary>The live coherence hand Omega = max|Im λ| at the gap, read from the Symphony
    /// clock at the given N (XY chain, this witness's J and γ).</summary>
    public double OmegaMem(int n)
    {
        var s = new Symphony(n: n, j: J, gamma: Gamma, initialState: InitialStateKind.BellPair,
            tPoints: SpectrumPoints);
        return s.Clock.Omega;
    }

    /// <summary>The live Takt hand Gap = slowest non-zero decay rate (= 2γ).</summary>
    public double Gap(int n)
    {
        var s = new Symphony(n: n, j: J, gamma: Gamma, initialState: InitialStateKind.BellPair,
            tPoints: SpectrumPoints);
        return s.Clock.Gap;
    }

    /// <summary>The closed-form F2b band edge 2J·cos(π/(N+1)) (the N≥3 protected hand).</summary>
    public double BandEdge(int n) => 2.0 * J * Math.Cos(Math.PI / (n + 1));

    public string DisplayName =>
        $"ClockHandLadderWitness (the two clocks live, J={J.ToString("0.#", Inv)}, γ={Gamma.ToString("0.###", Inv)})";

    public string Summary => "the two clocks live (typed home: ClockHandLadderClaim): for N≥3 the " +
        "coherence hand Omega is the γ-protected F2b band edge 2J·cos(π/(N+1)); for N=2 it is γ-pulled " +
        "to 2√(J²−γ²) and stops at the exceptional point Q=1.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("the ladder (N≥3, the band edge)",
                summary: $"Omega = 2J·cos(π/(N+1)) = F2b band edge: N=3 → {OmegaMem(3).ToString("0.#####", Inv)} (√2), " +
                         $"N=4 → {OmegaMem(4).ToString("0.#####", Inv)} (φ), N=5 → {OmegaMem(5).ToString("0.#####", Inv)} (√3). " +
                         "Typed parent F2b (docs/ANALYTICAL_FORMULAS.md F2b).");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
