using System;
using System.Collections.Generic;
using System.Globalization;
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
/// <para>It reuses <see cref="Symphony"/> as the spectrum engine: each reading builds a Symphony
/// and reads its clock (the per-N spectrum is memoized so repeated renders do not re-diagonalize).
/// The typed parent edges to F2b and the Absorption Theorem live on the claim, not here. The
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

    // Per-N clock cache: one Symphony build (one dense eigendecomposition) per N, shared by
    // OmegaMem(n) and Gap(n). Keeps repeated Children/Summary renders from re-diagonalizing.
    private readonly Dictionary<int, (double Gap, double Omega)> _clockCache = new();

    private (double Gap, double Omega) ClockAt(int n)
    {
        if (_clockCache.TryGetValue(n, out var c)) return c;
        var s = new Symphony(n: n, j: J, gamma: Gamma, initialState: InitialStateKind.BellPair,
            tPoints: SpectrumPoints);
        c = s.Clock;
        _clockCache[n] = c;
        return c;
    }

    /// <summary>The live coherence hand Omega = max|Im λ| at the gap, read from the Symphony
    /// clock at the given N (XY chain, this witness's J and γ).</summary>
    public double OmegaMem(int n) => ClockAt(n).Omega;

    /// <summary>The live Takt hand Gap = slowest non-zero decay rate (= 2γ).</summary>
    public double Gap(int n) => ClockAt(n).Gap;

    /// <summary>The closed-form F2b band edge 2J·cos(π/(N+1)) (the N≥3 protected hand).</summary>
    public double BandEdge(int n) => 2.0 * J * Math.Cos(Math.PI / (n + 1));

    /// <summary>The closed-form N=2 pulled hand 2√(J²−γ²) (vanishes at the exceptional point
    /// Q=1, γ=J). NOT the F2b band edge: a different (population/antisymmetric-coherence) block.
    /// NOTE: the live <see cref="OmegaMem"/>(2) tracks this only for Q ≥ 2/√3 (γ ≤ √3/2·J), the
    /// regime where the pulled coherence sits above the equal-rate {Im=±J} modes that share the
    /// 2γ gap; nearer the EP the live clock's "max |Im| at the gap" tie-breaks to those ±J modes,
    /// so this closed form (not the live clock) is the honest standstill witness near Q=1.</summary>
    public double N2PulledOmega() => Gamma < J ? 2.0 * Math.Sqrt(J * J - Gamma * Gamma) : 0.0;

    /// <summary>The dial angle θ = arctan(Omega/Gap) at the given N (degrees). For N≥3,
    /// θ = arctan(Q·cos(π/(N+1))); for N=2, θ = arctan(√(Q²−1)), zero at the EP Q=1.
    /// NOTE: the method computes arctan(<see cref="OmegaMem"/>(n)/<see cref="Gap"/>(n)) from the
    /// live clock. For N=2 OmegaMem(2) equals the pulled form arctan(√(Q²−1)) only above the
    /// crossover Q=2/√3; near the EP OmegaMem(2) is the ±J band mode (see <see cref="N2PulledOmega"/>),
    /// so the method returns arctan(1/(2γ)) there, not arctan(√(Q²−1)).</summary>
    public double AngleDegrees(int n)
    {
        double gap = Gap(n);
        return gap <= 0 ? 0.0 : Math.Atan(OmegaMem(n) / gap) * 180.0 / Math.PI;
    }

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

            yield return new InspectableNode("the γ-protection (N≥3, the hand holds)",
                summary: $"quadrupling γ leaves the coherence hand unmoved: at N=3, Omega = " +
                         $"{OmegaMem(3).ToString("0.#####", Inv)} at γ={Gamma.ToString("0.###", Inv)}; the Takt hand " +
                         $"Gap = 2γ = {Gap(3).ToString("0.#####", Inv)} tracks γ. The |vac⟩⟨ψ_k| modes are simultaneous " +
                         "eigenoperators of L_D (rate −2γ, the Absorption Theorem) and L_H (frequency E_k), so nothing " +
                         "mixes. Typed parent AbsorptionTheoremClaim (docs/proofs/PROOF_ABSORPTION_THEOREM.md).");

            yield return new InspectableNode("the pull and the exceptional point (N=2)",
                summary: $"at N=2 the coherence hand is a different block: Omega = 2√(J²−γ²) = " +
                         $"{N2PulledOmega().ToString("0.#####", Inv)} (NOT the F2b band edge {BandEdge(2).ToString("0.#####", Inv)}). " +
                         "It stops at the exceptional point Q=1 (γ=J): the two clocks merge into one. " +
                         "The dial angle θ = arctan(√(Q²−1)) → 0 there. See docs/ANALYTICAL_FORMULAS.md F2b corollary. " +
                         "(Live note: OmegaMem(2) tracks this pulled hand only for Q > 2/√3; nearer the EP the live " +
                         "clock's gap shares equal-rate {Im=±J} modes and the closed form is the honest standstill witness.)");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
