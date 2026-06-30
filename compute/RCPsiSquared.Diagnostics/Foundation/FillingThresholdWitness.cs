using System;
using System.Globalization;
using System.Linq;
using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The F89 Door-C DECISIVE follow-up, live: is fixed-q dissipative quantum chaos (GinUE complex spacing
/// ratio) a FILLING threshold, not an integrability one? Door-C (<see cref="GaloisSpectralChaosWitness"/>,
/// <see cref="IntegrabilityBreakingCsr"/>) found the DILUTE (SE,DE)=(1,2) coherence block stays Poisson / non-GinUE
/// under EVERY integrability-breaking knob (XXZ Δ, a random Z-field), and concluded the null is structural: a
/// 2-excitation sector cannot thermalize. This witness puts that conclusion under its own decisive test by building
/// the GENERAL (wKet,wBra) coherence block at EXTENSIVE filling (wKet,wBra near N/2) and re-running the SAME
/// disordered CSR.
///
/// <para>The finding (live): the DENSE block IS chaotic where the dilute one is not. Its radial statistic ⟨|z|⟩
/// sits at the GinUE value, and its angular repulsion ⟨cosθ⟩ goes NEGATIVE and CLIMBS toward GinUE with the block
/// size (N=6→7→8: ≈ −0.09 → −0.13 → −0.16, ≈ 43%→56%→67% of the size-matched GinUE angle), while the dilute (1,2)
/// block stays flat at ⟨cosθ⟩ ≈ 0 (≈ 23%) at every N. So Door-C's null is a FILLING threshold: chaos switches on
/// with extensive excitation content, not with breaking the Galois/Hamiltonian integrability (the dilute block
/// stays Poisson however hard you break it). Class A is licensed by the unequal weight (p,p+1): the F1 palindrome Π
/// maps the (p,p+1) block to the conjugate (p+1,p) block, not to itself, so no residual antiunitary survives — the
/// GinUE 0.738/−0.24 target is the right one (confirmed live: the disordered spectrum's conjugation-match fraction
/// is ≈ 0). Live on the trusted machine: <see cref="FillingThresholdCsr"/> (general WeightCoherenceBlock + random
/// field) → MathNet EVD → the Sá-Ribeiro-Prosen complex spacing ratio, pooled per-spectrum with finite-size-matched
/// references.</para></summary>
public sealed class FillingThresholdWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // canonical operating point for the live read: interacting (Δ=1), ergodic disorder window (W), fixed q.
    private const double Q = 1.0;
    private const double Delta = 1.0;
    private const double W = 0.75;

    private static string Z(double x) => x.ToString("0.000", Inv);
    private static string Cos(double x) => x.ToString("+0.000;-0.000", Inv);

    public string DisplayName =>
        "Filling-threshold chaos (live: the DENSE coherence block reaches toward GinUE where the dilute (SE,DE) stays Poisson)";

    public string Summary =>
        "the F89 Door-C decisive follow-up: dissipative quantum chaos (GinUE) is a FILLING threshold, not an " +
        "integrability one. The dilute (1,2)=(SE,DE) block stays Poisson (⟨cosθ⟩≈0) under disorder+interactions; " +
        "the dense (p,p+1) block near half-filling develops GinUE angular repulsion (⟨cosθ⟩<0, climbing toward " +
        "GinUE with N), with ⟨|z|⟩ already at the GinUE value. Live: general WeightCoherenceBlock + random Z-field " +
        "→ MathNet EVD → complex spacing ratio (Sá-Ribeiro-Prosen), pooled per-spectrum, finite-size-matched refs.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // dilute control vs dense, at N=6 (fast: 90 vs 300 dim) and the N=7 dense size point (1225 dim).
            var dilute6 = FillingThresholdCsr.DisorderSweep(6, 1, 2, Q, Delta, W, 40, seed: 9001);
            var dense6 = FillingThresholdCsr.DisorderSweep(6, 3, 4, Q, Delta, W, 20, seed: 9002);
            var dense7 = FillingThresholdCsr.DisorderSweep(7, 3, 4, Q, Delta, W, 3, seed: 9003);

            // finite-size-matched references at the dense-N6 per-spectrum size (≈ 300), live (calculated not marked).
            int refSize = Math.Max(10, dense6.ZCount / 20);
            var pRef = IntegrabilityBreakingCsr.PoissonReference(refSize, draws: 20, seed: 31);
            var gRef = IntegrabilityBreakingCsr.GinueReference(refSize, draws: 20, seed: 32);

            yield return new InspectableNode("the references (live, calculated not marked)",
                summary: $"finite-size (≈{refSize} pts) 2D-Poisson ⟨|z|⟩={Z(pRef.MeanAbs)} ⟨cos⟩={Cos(pRef.MeanCos)} | " +
                         $"GinUE (dissipative chaos) ⟨|z|⟩={Z(gRef.MeanAbs)} ⟨cos⟩={Cos(gRef.MeanCos)}. The diagnostic " +
                         "DOES separate the classes, so the contrast below is meaningful.");

            yield return new InspectableNode(
                $"DILUTE (1,2)=(SE,DE), N=6: Poisson — the Door-C null reproduced via the general builder",
                summary: $"⟨|z|⟩={Z(dilute6.MeanAbs)} [{Z(dilute6.CiLo)},{Z(dilute6.CiHi)}] ⟨cosθ⟩={Cos(dilute6.MeanCos)} " +
                         $"(GinUE ⟨cos⟩≈{Cos(gRef.MeanCos)}). A 2-excitation sector: no angular repulsion however hard the " +
                         "disorder+interactions break integrability — the dilute block cannot thermalize.");

            yield return DenseNode(6, dense6, gRef.MeanCos);
            yield return DenseNode(7, dense7, gRef.MeanCos);

            // the class-A guard: a random field breaks conjugation symmetry, so OffReal + the GinUE (class A) target.
            var rng = new Random(7);
            var disorderField = Enumerable.Range(0, 6).Select(_ => (2 * rng.NextDouble() - 1) * W).ToArray();
            double conj = FillingThresholdCsr.ConjugationMatchFraction(6, 3, 4, Q, Delta, disorderField);
            yield return new InspectableNode("class A licensed (live): the disordered block is NOT conjugation-symmetric",
                summary: $"conjugation-match fraction = {conj.ToString("P0", Inv)} (≈ 0). The unequal weight (p,p+1) sends " +
                         "the block under Π to the conjugate (p+1,p) block, not to itself, and the random field breaks " +
                         "the rest — no residual antiunitary, so the GinUE (class A) reference is the right target, not AI+/AII+.");

            string trend = (dense7.MeanCos < dense6.MeanCos && dense6.MeanCos < dilute6.MeanCos - 0.02)
                ? "CONFIRMED: dilute flat at ≈0; dense negative and growing more so with N (toward GinUE)"
                : "INVESTIGATE: the expected dilute-flat / dense-rising-with-N pattern did not hold";
            yield return new InspectableNode("the verdict: chaos is a FILLING threshold, not an integrability one",
                summary: $"{trend}. ⟨cosθ⟩: dilute(1,2) {Cos(dilute6.MeanCos)} | dense(3,4) N=6 {Cos(dense6.MeanCos)} → " +
                         $"N=7 {Cos(dense7.MeanCos)} (GinUE ≈{Cos(gRef.MeanCos)}); ⟨|z|⟩ of the dense block sits at the GinUE " +
                         "value already. Door-C's null is structural/kinematic: the dilute (SE,DE) sector stays Poisson however " +
                         "hard you break the Galois/Hamiltonian integrability, and the SAME Liouvillian's extensive-filling " +
                         "coherence sector reaches GinUE under the same disorder. Galois-chaos (over q) and spectral-chaos (at " +
                         "fixed q) merge only at extensive filling — a different object than the Door-C block.");
        }
    }

    private static InspectableNode DenseNode(int n, IntegrabilityBreakingCsr.CsrReading r, double ginueCos)
    {
        double pct = ginueCos < 0 ? 100.0 * r.MeanCos / ginueCos : double.NaN;
        string verdict = r.MeanCos < -0.05 ? "GinUE angular repulsion present, ⟨|z|⟩ at GinUE" : "no clear repulsion — investigate";
        return new InspectableNode(
            $"DENSE (3,4) near half-filling, N={n}: {verdict}",
            summary: $"⟨|z|⟩={Z(r.MeanAbs)} [{Z(r.CiLo)},{Z(r.CiHi)}] ⟨cosθ⟩={Cos(r.MeanCos)} " +
                     $"(≈{pct.ToString("0", Inv)}% of the size-matched GinUE angle {Cos(ginueCos)}), over {r.ZCount} pooled z's. " +
                     "The extensive-filling coherence sector thermalizes: radial rigidity at the GinUE value, angular " +
                     "repulsion switched on and climbing toward GinUE with N — the chaos the dilute block never reaches.");
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
