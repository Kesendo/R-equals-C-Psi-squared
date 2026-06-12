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
/// <para>It reuses <see cref="Symphony"/> as the spectrum engine: each reading builds a Symphony
/// and reads its clock (the per-N spectrum is memoized so repeated renders do not re-diagonalize).
/// The typed parent edges to F2b and the Absorption Theorem live on the claim, not here. The
/// witness is the numeric lab: it recomputes the live clock and checks it against the closed
/// forms.</para></summary>
public sealed class ClockHandLadderWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // Default operating point: the EP-family coupling J_EP=0.075 and the canonical carrier γ₀=0.05,
    // so Q = J/γ = 1.5 (= EP(2) − ½). This is a deliberately marginal regime: only the N=3 rung of the
    // band-edge ladder still holds (the higher rungs are overdamped), so the witness defaults to the
    // physically loaded point near the exceptional point, not the abstract clean ladder. Pass --J 1
    // (or new ClockHandLadderWitness(j: 1.0)) for the H-competitive regime where all rungs are protected.
    public const double DefaultJ = 0.075;
    public const double DefaultGamma = 0.05;
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

    /// <summary>The dial angle θ = arctan(coherence-hand/Gap) at the given N (degrees). For N≥3,
    /// θ = arctan(Q·cos(π/(N+1))); for N=2, θ = arctan(√(Q²−1)), zero at the EP Q=1. The dial reads the
    /// hand that actually stops at the EP, so for N=2 it uses the closed-form pulled hand
    /// <see cref="N2PulledOmega"/> (the genuine coalescing mode), NOT the raw <see cref="OmegaMem"/>(2):
    /// below the crossover Q=2/√3 the raw max|Im| clock returns the ±J band mode and would give
    /// arctan(1/(2γ)) (e.g. 29.05° rather than the honest 25.84° at γ=0.9), which is not the dial we mean.
    /// For N≥3 the band edge IS the raw clock hand, so OmegaMem(n) is used directly.</summary>
    public double AngleDegrees(int n)
    {
        double gap = Gap(n);
        if (gap <= 0) return 0.0;
        double hand = n == 2 ? N2PulledOmega() : OmegaMem(n);
        return Math.Atan(hand / gap) * 180.0 / Math.PI;
    }

    /// <summary>True when the F2b band edge IS the gap-mode at this N and γ (the protected regime): the
    /// live coherence hand equals the closed-form band edge AND the Takt hand equals 2γ. The band-edge
    /// protection holds only while the |vac⟩⟨ψ_k| sector (decay rate 2γ) is the slowest mode, the
    /// H-competitive regime. At strong dephasing a slower real (non-oscillating) mode takes the gap;
    /// OmegaMem(n) then drops to ~0 and the band edge is no longer the coherence hand, so this returns
    /// false. The H-competitive regime (high Q, J≫γ) is in; the default operating point Q=1.5 is marginal
    /// (only the N=3 rung holds). Only meaningful for n≥3 (at n=2 the gap-mode is the pulled block, not the
    /// band edge; see <see cref="N2PulledOmega"/>).</summary>
    public bool BandEdgeIsTheGapMode(int n)
    {
        const double tol = 1e-6;
        return Math.Abs(OmegaMem(n) - BandEdge(n)) < tol && Math.Abs(Gap(n) - 2.0 * Gamma) < tol;
    }

    public string DisplayName =>
        $"ClockHandLadderWitness (the two clocks live, J={J.ToString("0.###", Inv)}, γ={Gamma.ToString("0.###", Inv)}, Q={(J / Gamma).ToString("0.##", Inv)})";

    public string Summary => "the two clocks live (typed home: ClockHandLadderClaim): for N≥3 the " +
        "coherence hand Omega is the F2b band edge 2J·cos(π/(N+1)), γ-protected in the H-competitive regime " +
        "(high Q); for N=2 it is γ-pulled to 2√(J²−γ²) and stops at the exceptional point Q=1. The default " +
        "operating point is Q=1.5 (the EP family), where only the N=3 rung still holds.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            int protectedCount = new[] { 3, 4, 5 }.Count(BandEdgeIsTheGapMode);
            bool protectedRegime = protectedCount == 3;
            string regimeHead = $"γ={Gamma.ToString("0.###", Inv)}, J={J.ToString("0.###", Inv)} (Q={(J / Gamma).ToString("0.##", Inv)})";

            yield return new InspectableNode("the ladder (N≥3, the band edge)",
                summary: protectedRegime
                    ? $"the band-edge ladder ω_mem/J = √2/φ/√3 (J-independent): live ω_mem at {regimeHead}: " +
                      $"N=3 → {OmegaMem(3).ToString("0.#####", Inv)} (J√2), N=4 → {OmegaMem(4).ToString("0.#####", Inv)} (Jφ), " +
                      $"N=5 → {OmegaMem(5).ToString("0.#####", Inv)} (J√3). The |vac⟩⟨ψ_k| sector (rate 2γ) is the slowest mode. " +
                      "Typed parent F2b (docs/ANALYTICAL_FORMULAS.md F2b)."
                    : $"out of the protected regime at {regimeHead}: the band edge ω_mem/J = √2/φ/√3 is no longer the slowest mode for " +
                      $"{3 - protectedCount} of the three rungs; a real overdamped mode (gap {Gap(3).ToString("0.#####", Inv)}) takes over, so " +
                      $"live ω_mem reads N=3 → {OmegaMem(3).ToString("0.#####", Inv)}, N=4 → {OmegaMem(4).ToString("0.#####", Inv)}, " +
                      $"N=5 → {OmegaMem(5).ToString("0.#####", Inv)} (the higher rungs leave first). The band-edge protection holds only in the " +
                      "H-competitive regime (high Q). Typed parent F2b (docs/ANALYTICAL_FORMULAS.md F2b).");

            yield return new InspectableNode("the γ-protection (N≥3, the hand holds)",
                summary: protectedRegime
                    ? $"in the protected regime the coherence hand sits at the band edge independently of γ: at N=3, ω_mem = " +
                      $"{OmegaMem(3).ToString("0.#####", Inv)} = J√2, the Takt hand Gap = 2γ = {Gap(3).ToString("0.#####", Inv)} tracks γ. " +
                      "The |vac⟩⟨ψ_k| modes are simultaneous eigenoperators of L_D (rate −2γ, the Absorption Theorem) and L_H " +
                      "(frequency E_k), so nothing mixes (γ-independent while all of N=3-5 hold). Typed parent " +
                      "AbsorptionTheoremClaim (docs/proofs/PROOF_ABSORPTION_THEOREM.md)."
                    : $"the band-edge protection is {(protectedCount == 0 ? "exited" : "partially exited")} at {regimeHead}: " +
                      $"{protectedCount}/3 of N=3-5 still sit at the band edge (the higher rungs leave first as a slower real mode takes the " +
                      "gap; see the ladder node for the live per-rung values). The simultaneous-eigenoperator protection (L_D rate −2γ, L_H " +
                      "frequency E_k) holds rung-by-rung only while the |vac⟩⟨ψ_k| sector (rate 2γ) is the slowest mode. Typed parent " +
                      "AbsorptionTheoremClaim (docs/proofs/PROOF_ABSORPTION_THEOREM.md).");

            yield return new InspectableNode("the pull and the exceptional point (N=2)",
                summary: $"at N=2 the coherence hand is a different block: Omega = 2√(J²−γ²) = " +
                         $"{N2PulledOmega().ToString("0.#####", Inv)} (NOT the F2b band edge {BandEdge(2).ToString("0.#####", Inv)}). " +
                         "It stops at the exceptional point Q=1 (γ=J): the two clocks merge into one. " +
                         "The dial angle θ = arctan(√(Q²−1)) → 0 there. See docs/ANALYTICAL_FORMULAS.md F2b corollary. " +
                         "(Live note: OmegaMem(2) tracks this pulled hand only for Q > 2/√3; nearer the EP the live " +
                         "clock's gap shares equal-rate {Im=±J} modes and the closed form is the honest standstill witness.)");

            var ns = new double[] { 2, 3, 4, 5 };
            var omegas = ns.Select(n => OmegaMem((int)n)).ToArray();
            yield return new InspectableNode("the dial (the angle is Q)",
                summary: $"θ = arctan(coherence-hand/Gap): N=3 → {AngleDegrees(3).ToString("0.##", Inv)}°, " +
                         $"N=2 → {AngleDegrees(2).ToString("0.##", Inv)}°. For N≥3, θ = arctan(Q·cos(π/(N+1))) (the band " +
                         "edge); for N=2, θ = arctan(√(Q²−1)) from the pulled hand, zero at the EP. The angle between the " +
                         "hands is Q. (on the curve, the N=2 point is the pulled block, not a ladder rung)",
                payload: new InspectablePayload.Curve("Omega vs N", ns, omegas, "N", "Omega (coherence hand)"));
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
