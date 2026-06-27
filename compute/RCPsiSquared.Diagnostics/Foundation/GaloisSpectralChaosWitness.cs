using System;
using System.Collections.Generic;
using System.Globalization;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Galois vs spectral chaos, the sector-resolved RMT test (the door rmt_topology_csr.py left
/// open: "the clean RMT test would be SECTOR-resolved"). The (SE,DE) relaxation rates split into two
/// halves: the AT-locked half (Re/γ ∈ {−2,−6}, free-fermion, radically writable) and the H_B-mixed
/// half (the residue whose chain Galois group is S_n — S_8/18/32/53 — with NO radical closure). The
/// conjecture under test: does Galois-S_n (algebraic chaos over the q=J/γ field) show up as GinUE
/// spectral statistics (geometric chaos at fixed q)?
///
/// <para>The finding is a clean NULL. At accessible N the H_B-mixed half reads Poisson-like /
/// sub-Poisson (⟨|z|⟩ &lt; 0.66, ⟨cos θ⟩ ≈ 0), NOT GinUE (⟨|z|⟩ ≈ 0.738, ⟨cos θ⟩ ≈ −0.24): the
/// "scrambled" half still sits on the integrable frequency lattice. Algebraic chaos (the Galois group,
/// a monodromy statement over q) and spectral chaos (RMT at fixed q) are distinct here; the former does
/// not imply the latter. Live on the trusted machine: the shared SeDeBlockBuilder, MathNet EVD, and the
/// ComplexSpacingRatio diagnostic whose own GinUE reference confirms it CAN see chaos. The where-it-does-
/// live sequel is the q-parametric monodromy (the discriminant/EP loci), not the fixed-q geometry.</para>
///
/// <para>Both sides of the door's comparison are rendered: the H_B-mixed half is the Poisson-like / sub-
/// Poisson residue, and the AT-locked half (rates −2γ/−6γ, free-fermion Bloch frequencies) is the sparse
/// picket-fence — a structured set with low ⟨|z|⟩ (clustering, below GinUE's 0.74), not a chaos cloud.
/// Neither half is GinUE.</para></summary>
public sealed class GaloisSpectralChaosWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double ImTol = 1e-6;
    private const double RateTol = 1e-6;

    // upper-half-plane eigenvalues of L(q) = re + i·q·im, split into the AT-locked half (Re/γ pinned
    // to −2 or −6) and the H_B-mixed half (the spread residue). Upper half mirrors the validated scout.
    private static (List<Complex> at, List<Complex> hb) SectorEigenvalues(int n, string topo, double q)
    {
        var (re, im) = SeDeBlockBuilder.Build(n, topo);
        int d = re.GetLength(0);
        var l = Matrix<Complex>.Build.Dense(d, d, (a, b) => new Complex(re[a, b], q * im[a, b]));
        var vals = l.Evd().EigenValues;
        var at = new List<Complex>();
        var hb = new List<Complex>();
        for (int t = 0; t < d; t++)
        {
            var lam = vals[t];
            if (lam.Imaginary <= ImTol) continue;
            bool locked = Math.Abs(lam.Real + 2) < RateTol || Math.Abs(lam.Real + 6) < RateTol;
            (locked ? at : hb).Add(lam);
        }
        return (at, hb);
    }

    private static (double meanAbs, double meanCos, double avgCount) SectorCsr(
        int n, string topo, double[] qs, bool hbMixed)
    {
        double sa = 0, sc = 0, scount = 0;
        int used = 0;
        foreach (var q in qs)
        {
            var (at, hb) = SectorEigenvalues(n, topo, q);
            var (a, c, cnt) = ComplexSpacingRatio.Of(hbMixed ? hb : at);
            scount += cnt;
            if (!double.IsNaN(a)) { sa += a; sc += c; used++; }
        }
        return used == 0
            ? (double.NaN, double.NaN, scount / qs.Length)
            : (sa / used, sc / used, scount / qs.Length);
    }

    /// <summary>⟨|z|⟩, ⟨cos θ⟩, and avg distinct-point count of the H_B-mixed half over a q-sweep.</summary>
    public static (double meanAbs, double meanCos, double avgCount) HbMixedCsr(int n, string topo, double[] qs)
        => SectorCsr(n, topo, qs, hbMixed: true);

    /// <summary>⟨|z|⟩, ⟨cos θ⟩, and avg distinct-point count of the AT-locked half over a q-sweep.</summary>
    public static (double meanAbs, double meanCos, double avgCount) AtLockedCsr(int n, string topo, double[] qs)
        => SectorCsr(n, topo, qs, hbMixed: false);

    private static double[] QSweep(int count)
    {
        var qs = new double[count];
        for (int i = 0; i < count; i++) qs[i] = 0.3 + i * (3.7 / (count - 1));
        return qs;
    }

    // ---- IInspectable ----

    public string DisplayName =>
        "Galois vs spectral chaos (live: the S_n half is NOT GinUE at fixed q — a clean null)";

    public string Summary =>
        "the sector-resolved RMT test: does the H_B-mixed half (chain Galois S_8/18/32/53, no radical " +
        "closure) read as dissipative quantum chaos (GinUE) at fixed q? It does not — it reads Poisson-like/" +
        "sub-Poisson, still on the integrable frequency lattice. Algebraic chaos (Galois over q) ≠ spectral " +
        "chaos (RMT at q). Live: shared (SE,DE) block → MathNet EVD → complex spacing ratio (Sá-Ribeiro-Prosen).";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var qs = QSweep(8);
            var (pAbs, pCos) = ComplexSpacingRatio.PoissonDiskReference(3000, seed: 11);
            var (gAbs, gCos) = ComplexSpacingRatio.GinueReference(250, seed: 12);

            yield return new InspectableNode("the references (live, calculated not marked)",
                summary: $"2D-Poisson (integrable/fragmented) ⟨|z|⟩={pAbs.ToString("0.000", Inv)} " +
                         $"⟨cos⟩={pCos.ToString("+0.000;-0.000", Inv)} | GinUE (dissipative chaos) " +
                         $"⟨|z|⟩={gAbs.ToString("0.000", Inv)} ⟨cos⟩={gCos.ToString("+0.000;-0.000", Inv)}. " +
                         "The diagnostic DOES separate the two classes — so the null below is meaningful.");

            yield return ChainNode(7, qs, gAbs);
            yield return ChainNode(6, qs, gAbs);

            // the other side of the door's comparison: the AT-locked half should be picket-fence
            // (radically writable, rates pinned to the AT rungs), NOT a 2D chaos cloud.
            yield return AtLockedNode(7, qs);
            yield return AtLockedNode(6, qs);

            var (cAbs, _, cCnt) = HbMixedCsr(7, "complete", qs);
            yield return new InspectableNode("complete K_7 H_B-mixed: solvable ⟹ collapses",
                summary: $"the complete graph's H_B factors are all ≤ quartic (radically writable) and " +
                         $"S_N-degenerate: only {cCnt.ToString("0.#", Inv)} distinct H_B points/q remain — too few " +
                         $"for a CSR (avg ⟨|z|⟩={(double.IsNaN(cAbs) ? "n/a" : cAbs.ToString("0.000", Inv))}). " +
                         "Solvable does not even populate a 2D cloud; it fragments to a handful of levels.");

            yield return new InspectableNode("the verdict: Galois-S_n is NOT spectral GinUE at fixed q",
                summary: "the Abel-Ruffini non-closure of the rates (a statement about the Galois group over the " +
                         "q=J/γ field) leaves NO fingerprint in the fixed-q geometry: the H_B half is Poisson-like, " +
                         "not Ginibre. The two notions of dissipative chaos diverge. Where the Galois structure DOES " +
                         "live spectrally is the q-PARAMETRIC monodromy: the d roots braid as a single S_n orbit, " +
                         "avoided-crossing at the discriminant (the F89 diabolic/EP loci, --root f89octic). The " +
                         "fixed-q point cloud is the wrong object; the q-braid is the right one.");
        }
    }

    private static InspectableNode ChainNode(int n, double[] qs, double ginueAbs)
    {
        var (abs, cos, cnt) = HbMixedCsr(n, "chain", qs);
        string verdict = abs < 0.70 ? "Poisson-like / sub-Poisson, NOT GinUE" : "unexpectedly chaotic — investigate";
        return new InspectableNode($"chain N={n} H_B-mixed (Galois S_{(n == 7 ? "53" : "32")}): {verdict}",
            summary: $"⟨|z|⟩={abs.ToString("0.000", Inv)} ⟨cos θ⟩={cos.ToString("+0.000;-0.000", Inv)} " +
                     $"over {cnt.ToString("0.#", Inv)} distinct points/q (GinUE would be ⟨|z|⟩≈{ginueAbs.ToString("0.000", Inv)}, " +
                     "⟨cos⟩≈−0.24). The half with no radical closure still reads integrable-lattice, not chaos.");
    }

    // The AT-locked half: the two rate-rungs (−2γ, −6γ) carrying free-fermion Bloch frequencies, a
    // sparse STRUCTURED set (only ~10-16 distinct points/q), not a 2D cloud. The discriminator from
    // GinUE is ⟨|z|⟩, not the angle: GinUE chaos needs HIGH ⟨|z|⟩≈0.74 (repulsion spreads NN and NNN
    // apart); this half reads LOW ⟨|z|⟩ (below even 2D-Poisson's 0.66 — clustering). Its ⟨cos θ⟩ can run
    // negative (lattice angular order, plus few-point noise), but that alone is not chaos: genuine GinUE
    // would also demand the high ⟨|z|⟩ this half lacks. NaN-safe (a sparser sector can fall below the
    // CSR's 10-point floor, itself the picket-fence verdict).
    private static InspectableNode AtLockedNode(int n, double[] qs)
    {
        var (abs, cos, cnt) = AtLockedCsr(n, "chain", qs);
        string verdict = double.IsNaN(abs)
            ? "collapses — too few distinct points for a CSR (picket-fence, as expected)"
            : abs < 0.70 ? "picket-fence / sparse structured set, NOT a GinUE chaos cloud" : "unexpectedly high ⟨|z|⟩ — investigate";
        return new InspectableNode($"chain N={n} AT-locked (rates −2γ/−6γ, free-fermion Bloch): {verdict}",
            summary: $"⟨|z|⟩={(double.IsNaN(abs) ? "n/a" : abs.ToString("0.000", Inv))} " +
                     $"⟨cos θ⟩={(double.IsNaN(cos) ? "n/a" : cos.ToString("+0.000;-0.000", Inv))} " +
                     $"over only {cnt.ToString("0.#", Inv)} distinct points/q (vs ~50 for the H_B half). The radically-" +
                     "writable half is two AT rate-rungs carrying free-fermion Bloch frequencies — a structured set, not a " +
                     "cloud: its low ⟨|z|⟩ is clustering (far below GinUE's 0.74), and the negative ⟨cos θ⟩ is lattice " +
                     "angular order plus few-point noise, not the high-⟨|z|⟩ repulsion that genuine GinUE chaos requires.");
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
