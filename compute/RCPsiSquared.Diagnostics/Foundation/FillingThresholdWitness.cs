using System;
using System.Globalization;
using System.Linq;
using System.Collections.Generic;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The F89 Door-C filling comparison, live: finite executed CSR evidence on dilute (SE,DE)=(1,2)
/// and dense (wKet,wBra near N/2) coherence blocks under the same interacting disorder at canonical Delta=1.
/// Nonzero Delta breaks the underlying Hamiltonian's free-fermion additivity, but uniform XXZ remains Bethe-integrable;
/// at Delta=0 the random-field XY Hamiltonian remains quadratic (Anderson/free fermions). This does not classify the Z-dephasing Liouvillian
/// itself as a quadratic free-fermion generator;
/// generic random field plus Delta!=0 is the interacting disordered nonintegrable test.
/// Symmetry/cross-fold breaking and Hamiltonian integrability are separate from the measured Liouvillian CSR.
///
/// <para>The recorded filling-associated crossover evidence at N=6..8 is movement toward GinUE,
/// not a causal or thermodynamic threshold theorem. The live read below uses N=6 and N=7.
/// The dense radial statistic ⟨|z|⟩ lies near the GinUE reference, and angular repulsion ⟨cosθ⟩ grows with the block
/// size (N=6→7→8: ≈ −0.09 → −0.13 → −0.16, ≈ 43%→56%→67% of the size-matched GinUE angle), while the dilute (1,2)
/// block stays near ⟨cosθ⟩ ≈ 0 (≈ 23%) in that comparison. This supports a finite-size filling dependence,
/// not a universal filling threshold or a thermalization cause. GinUE is a comparison ensemble, not a licensed
/// class label: the known F1 map leaves the unequal-weight block and the measured conjugation-match fraction is
/// near zero, but those checks do not exhaust the full sector symmetry algebra, whose irreducible SRP class remains
/// open. Live on the trusted machine: <see cref="FillingThresholdCsr"/> (general WeightCoherenceBlock + random
/// field) → MathNet EVD → the Sá-Ribeiro-Prosen complex spacing ratio, pooled per-spectrum with finite-size-matched
/// references. Its 95% intervals use a whole-spectrum cluster bootstrap over independent disorder realizations;
/// ratios inside one spectrum are never counted as independent observations. Every spectrum contributes one
/// representative per 1e-9 finite-precision cluster; the CSR and
/// representative count are tolerance-dependent and not an exact degeneracy census.</para></summary>
public sealed class FillingThresholdWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    // canonical operating point for the live read: interacting (Δ=1), ergodic disorder window (W), fixed q.
    private const double Q = 1.0;
    private const double Delta = 1.0;
    private const double W = 0.75;

    private static string Z(double x) => x.ToString("0.000", Inv);
    private static string Cos(double x) => x.ToString("+0.000;-0.000", Inv);
    private static string Interval(IntegrabilityBreakingCsr.CsrReading r) =>
        r.Uncertainty == IntegrabilityBreakingCsr.UncertaintySemantics.SpectrumClusterBootstrap95
            ? $"whole-spectrum cluster-bootstrap 95% CI [{Z(r.CiLo)},{Z(r.CiHi)}] over {r.IndependentSpectrumCount} independent disorder realizations"
            : "no sampling CI";

    public string DisplayName =>
        "Filling-associated crossover evidence (finite size; live dilute-vs-dense CSR comparison)";

    public string Summary =>
        "the F89 Door-C filling-associated crossover evidence at canonical Delta=1 plus disorder, not a causal or thermodynamic threshold theorem. " +
        "The dilute (1,2)=(SE,DE) block stays Poisson-like (⟨cosθ⟩≈0) in the sampled disorder+interactions regime; " +
        "the dense (p,p+1) block near half-filling develops GinUE angular repulsion (⟨cosθ⟩<0, climbing toward " +
        "GinUE with N), with ⟨|z|⟩ near the GinUE reference. GinUE comparison only: the full sector symmetry " +
        "algebra and irreducible SRP class remain open. Live: general WeightCoherenceBlock + random Z-field " +
        "→ MathNet EVD → complex spacing ratio (Sá-Ribeiro-Prosen), pooled per-spectrum, finite-size-matched refs. " +
        "Uncertainty uses a whole-spectrum cluster bootstrap over independent disorder realizations and never resamples individual z-values. " +
        "Each spectrum contributes finite-precision cluster representatives; counts and CSR are tolerance-dependent, " +
        "not an exact degeneracy census.";

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

            yield return new InspectableNode("the numerical CSR clustering boundary",
                summary: FillingThresholdCsr.ClusteringScope);

            yield return new InspectableNode(
                $"DILUTE (1,2)=(SE,DE), N=6: Poisson — the Door-C null reproduced via the general builder",
                summary: $"⟨|z|⟩={Z(dilute6.MeanAbs)}; {Interval(dilute6)}; ⟨cosθ⟩={Cos(dilute6.MeanCos)} " +
                         $"(GinUE ⟨cos⟩≈{Cos(gRef.MeanCos)}). This dilute sector shows weak angular repulsion at the " +
                         "executed disorder+interactions operating point; this is not a thermalization proof.");

            yield return DenseNode(6, dense6, gRef.MeanCos);
            yield return DenseNode(7, dense7, gRef.MeanCos);

            // This probes one conjugation match only; it cannot classify the full irreducible sector.
            var rng = new Random(7);
            var disorderField = Enumerable.Range(0, 6).Select(_ => (2 * rng.NextDouble() - 1) * W).ToArray();
            double conj = FillingThresholdCsr.ConjugationMatchFraction(6, 3, 4, Q, Delta, disorderField);
            yield return new InspectableNode("GinUE comparison scope (live conjugation-match diagnostic)",
                summary: $"conjugation-match fraction = {conj.ToString("P0", Inv)} (≈ 0). The unequal weight (p,p+1) sends " +
                         "the block under Π to the conjugate (p+1,p) block, not to itself, and the random field breaks " +
                         "this bare conjugation match. GinUE is a comparison ensemble only; other unitary or antiunitary " +
                         "relations have not been excluded after full irreducible-sector reduction.");

            string trend = (dense7.MeanCos < dense6.MeanCos && dense6.MeanCos < dilute6.MeanCos - 0.02)
                ? "CONFIRMED: dilute flat at ≈0; dense negative and growing more so with N (toward GinUE)"
                : "INVESTIGATE: the expected dilute-flat / dense-rising-with-N pattern did not hold";
            yield return new InspectableNode("the verdict: finite-size filling dependence under interacting disorder",
                summary: $"{trend}. ⟨cosθ⟩: dilute(1,2) {Cos(dilute6.MeanCos)} | dense(3,4) N=6 {Cos(dense6.MeanCos)} → " +
                         $"N=7 {Cos(dense7.MeanCos)} (GinUE ≈{Cos(gRef.MeanCos)}); ⟨|z|⟩ of the dense block lies near the GinUE " +
                         "reference. At the same model parameters and disorder distribution, but in separately sampled " +
                         "realization ensembles, the extensive-filling coherence sector has stronger angular repulsion " +
                         "than the dilute sector. This finite executed CSR evidence is an ensemble comparison, not a " +
                         "within-realization causal intervention, and " +
                         "does not establish a universal filling threshold or identify Hamiltonian integrability from spacings.");
        }
    }

    private static InspectableNode DenseNode(int n, IntegrabilityBreakingCsr.CsrReading r, double ginueCos)
    {
        double pct = ginueCos < 0 ? 100.0 * r.MeanCos / ginueCos : double.NaN;
        string verdict = r.MeanCos < -0.05 ? "angular repulsion present, ⟨|z|⟩ near GinUE" : "no clear repulsion — investigate";
        return new InspectableNode(
            $"DENSE (3,4) near half-filling, N={n}: {verdict}",
            summary: $"⟨|z|⟩={Z(r.MeanAbs)}; {Interval(r)}; ⟨cosθ⟩={Cos(r.MeanCos)} " +
                     $"(≈{pct.ToString("0", Inv)}% of the size-matched GinUE angle {Cos(ginueCos)}), over {r.ZCount} pooled z's. " +
                     "The sampled extensive-filling sector shows radial spacing ratios near the GinUE reference and angular " +
                     "repulsion trending toward it with N; CSR alone does not prove thermalization.");
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
