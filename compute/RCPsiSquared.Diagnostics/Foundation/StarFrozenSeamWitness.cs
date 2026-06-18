using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The star's frozen seam computed live (typed home: <c>StarFrozenSeamClaim</c>). On the star the
/// longest-lived coherence NEVER un-freezes: its survivor is the darkest [H,A]=0 commutant (1,1) coherence
/// (it lies in the ad_H kernel, so it commutes with H and cannot oscillate — frozen, |Im|=0, by construction),
/// and it is the survivor exactly when it undercuts the −2γ floor, i.e. when the structural ceiling
/// g2 = 4/(N−1) ≤ 1 (N ≥ 5). At N=4 (g2 = 4/3 > 1) the commutant mode is brighter than the floor, an
/// oscillating band-edge mode is the slowest, and the star un-freezes — the known (2,2)/K₄ outlier.
///
/// <para>So the star's frozen seam IS the structural ceiling (<see cref="StructuralCeilingWitness"/>),
/// read DYNAMICALLY: the high-Q ceiling g2 ≤ 1 and the all-Q frozenness of the survivor are the same fact.
/// It is the third member of the trichotomy chain (un-freezes at the SE-EP coherence horizon Q*(N)) / ring
/// (frozen (2,2) level crossing) / star (frozen (1,1) commutant). This witness recomputes (a) the exact
/// commutant ceiling g2 = 4/(N−1) (cheap, no 4^N), and (b) the |Im| of the global slowest mode of the full
/// Liouvillian at N=4,5 across Q, showing N=5 frozen at every Q and N=4 un-freezing at high Q. The C# twin
/// of simulations/star_frozen_seam.py (which carries the N=6,7,8 dynamic confirmation).</para></summary>
public sealed class StarFrozenSeamWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Gamma = 0.05;            // the fixed dephasing the witness reads Q = J/γ against
    private const double Frozen = 1e-6;           // |Im| below this = frozen (non-oscillating)

    /// <summary>The (1,1) commutant ceiling g2 = 4/(N−1): the darkest H-commutant coherence in the
    /// single-excitation sector of the star, reused from <see cref="StructuralCeilingWitness.CommutantDarkest"/>.
    /// Because the commutant lives in the ad_H kernel ([H,A]=0), this coherence is FROZEN (|Im|=0) by
    /// construction; the only question is whether it is the survivor (g2 ≤ 1).</summary>
    public static double Ceiling(int n) => StructuralCeilingWitness.CommutantDarkest("star", n, 1, 1) ?? 1.0;

    /// <summary>(g2, max|Im|) of the global slowest non-kernel mode of the full 4^N star Liouvillian at
    /// coupling J = Q·γ. The |Im| is the max over the modes within a tolerance of the slowest decay rate
    /// (so an oscillating mode degenerate with a frozen one is detected). N ≤ 5 only (the full 4^N cost).</summary>
    public static (double G2, double ImMax) FullSlowest(int n, double q)
    {
        int d = 1 << n;
        double scale = q * Gamma;

        // star XY hopping H (bonds (0,i)), real symmetric, amplitude scale per hop
        var H = new double[d, d];
        for (int i = 1; i < n; i++)
            for (int s = 0; s < d; s++)
            {
                int b0 = s & 1, bi = (s >> i) & 1;
                if (b0 != bi) { int s2 = s ^ 1 ^ (1 << i); H[s2, s] += scale; }
            }

        int dim = d * d;
        var L = Matrix<Complex>.Build.Dense(dim, dim);
        var negI = new Complex(0, -1);
        for (int a = 0; a < d; a++)
            for (int c = 0; c < d; c++)
            {
                double h = H[a, c];
                if (h != 0.0) for (int b = 0; b < d; b++) L[a * d + b, c * d + b] += negI * h;
            }
        for (int b = 0; b < d; b++)
            for (int e = 0; e < d; e++)
            {
                double h = H[b, e];
                if (h != 0.0) for (int a = 0; a < d; a++) L[a * d + b, a * d + e] += -negI * h;
            }
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
            {
                int hamm = BitOperations.PopCount((uint)(a ^ b));
                L[a * d + b, a * d + b] += -2.0 * Gamma * hamm;
            }

        var evs = L.Evd().EigenValues.Where(e => e.Real < -1e-9).ToArray();
        if (evs.Length == 0) return (double.NaN, 0.0);
        double slowestRe = evs.Max(e => e.Real);
        double imMax = evs.Where(e => Math.Abs(e.Real - slowestRe) < 1e-6).Max(e => Math.Abs(e.Imaginary));
        return (-slowestRe / (2.0 * Gamma), imMax);
    }

    /// <summary>Does the star's global slowest mode un-freeze anywhere on a Q-sweep? (false for N ≥ 5,
    /// true at N=4.) The max |Im| over the sweep; "frozen at all Q" iff it stays below <c>Frozen</c>.</summary>
    public static double MaxImOverSweep(int n, params double[] qs) =>
        qs.Select(q => FullSlowest(n, q).ImMax).DefaultIfEmpty(0.0).Max();

    private static readonly double[] Sweep = { 1.0, 2.0, 4.0, 8.0, 16.0, 32.0 };

    // ---- IInspectable ----

    public string DisplayName => "StarFrozenSeamWitness (the star survivor never un-freezes; the ceiling read dynamically)";

    public string Summary
    {
        get
        {
            string g4 = Ceiling(4).ToString("0.####", Inv);
            string g6 = Ceiling(6).ToString("0.####", Inv);
            return
                "the star's frozen seam computed live (typed home: StarFrozenSeamClaim): the longest-lived " +
                "coherence is the [H,A]=0 commutant (1,1) coherence (in the ad_H kernel ⟹ frozen, |Im|=0, by " +
                $"construction), the survivor iff it undercuts the −2γ floor, i.e. g2=4/(N−1)≤1 (N≥5). N=4 " +
                $"(g2={g4}>1) un-freezes (the band edge wins, the (2,2)/K₄ outlier); N=6 (g2={g6}<1) frozen. " +
                "The third member of chain(SE-EP un-freeze) / ring(frozen level crossing) / star(frozen commutant): " +
                "the structural ceiling read dynamically. Sector overview: inspect --root ceiling / --root survivor.";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheCeilingLadder();
            yield return TheImSweep();
            yield return TheTrichotomy();
        }
    }

    /// <summary>The exact commutant ceiling g2 = 4/(N−1) across N (cheap, no 4^N), with the frozen/un-freeze
    /// verdict by the g2 ≤ 1 threshold. The (1,1) commutant coherence is frozen for every N (it commutes with
    /// H); g2 decides whether it is the survivor.</summary>
    private InspectableNode TheCeilingLadder()
    {
        var rows = new List<IInspectable>();
        for (int n = 4; n <= 8; n++)
        {
            double g = Ceiling(n);
            double cf = 4.0 / (n - 1);
            bool frozenSeam = g <= 1.0 + 1e-9;
            rows.Add(new InspectableNode($"star N={n}",
                summary: $"g2 = {g.ToString("0.#####", Inv)} (4/(N−1) = {cf.ToString("0.#####", Inv)}, " +
                         $"{(Math.Abs(g - cf) < 1e-6 ? "match" : "MISMATCH")}); the (1,1) commutant " +
                         $"{(frozenSeam ? "undercuts the floor ⟹ it is the FROZEN survivor" : "is brighter than the floor ⟹ the survivor is an oscillating mode (UN-FREEZES)")}"));
        }
        return new InspectableNode("the commutant ceiling g2 = 4/(N−1) (the threshold AND the high-Q darkness)",
            summary: "the darkest [H,A]=0 commutant (1,1) coherence on the star leaf manifold, recomputed live " +
                     "(StructuralCeilingWitness.CommutantDarkest). It is frozen at every N (ad_H kernel); the " +
                     "frozen seam holds iff it is the survivor, i.e. g2 = 4/(N−1) ≤ 1 (N ≥ 5). N=4 (4/3 > 1) is " +
                     "the lone outlier. g2 is also the value the survivor darkness saturates at in the high-Q " +
                     "limit: ⟨n_XY⟩(Q→∞) = min(g2,1) — on g2 below the −2γ floor (N≥6), or pinned to the floor 1 " +
                     "when g2 ≥ 1 (N≤5, where the star un-freezes).",
            children: rows);
    }

    /// <summary>The dynamic |Im|(Q) check on the full Liouvillian: N=5 frozen at every Q, N=4 un-freezing at
    /// high Q (N ≤ 5 only; N=6,7,8 confirmed frozen in simulations/star_frozen_seam.py).</summary>
    private InspectableNode TheImSweep()
    {
        var rows = new List<IInspectable>();
        foreach (int n in new[] { 4, 5 })
        {
            var cells = Sweep.Select(q => $"{FullSlowest(n, q).ImMax:0.00e+0}");
            double mx = MaxImOverSweep(n, Sweep);
            string verdict = mx < Frozen ? "frozen at every Q" : "UN-FREEZES";
            rows.Add(new InspectableNode($"star N={n} |Im| over Q∈{{1..32}}",
                summary: $"[{string.Join(", ", cells)}]; max |Im| = {mx:0.00e+0} → {verdict}"));
        }
        return new InspectableNode("the |Im|(Q) sweep (the dynamic check, full 4^N, N≤5)",
            summary: "the global slowest mode's oscillation frequency across Q. The star N=5 survivor stays " +
                     "frozen (|Im|≈0) at every Q (the marginal g2=1 commutant); the N=4 survivor un-freezes at " +
                     "high Q (the band edge overtakes, the outlier). A single Q does not separate the topologies " +
                     "(below its horizon every survivor is frozen); the |Im|(Q) curve is the signature.",
            children: rows);
    }

    private static InspectableNode TheTrichotomy() =>
        new InspectableNode("the trichotomy (chain / ring / star) and the mechanism",
            summary: "as Q grows (the dephasing weakens), does the survivor acquire a frequency? CHAIN: yes, at " +
                     "the SE-EP coherence horizon Q*(N) (CoherenceHorizonClaim). RING: yes, at the handover where " +
                     "the oscillating band edge overtakes the frozen (2,2) level-crossing seam. STAR: no (N≥5) — " +
                     "its survivor is the (1,1) commutant coherence, which commutes with H and cannot oscillate. " +
                     "The star's freeze is a DIFFERENT route to |Im|=0 than the ring's: a commutant (the mode " +
                     "commutes with H) vs a level crossing (two real eigenvalues coincide). This is the survivor-" +
                     "level reading of PROOF_STRUCTURAL_CEILING §7's 'the star has no coherence horizon'.");

    public InspectablePayload Payload => InspectablePayload.Empty;
}
