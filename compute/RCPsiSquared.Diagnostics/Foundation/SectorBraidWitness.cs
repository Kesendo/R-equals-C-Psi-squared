using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F89PathK;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Multi-Sector Monodromy verdict, live: is the S₈ Galois braid the F89 (1,2) octic carries CONFINED
/// to the (1,2) coherence orbit, or SHARED across the joint-popcount sectors of the Liouvillian L(q)? The answer is
/// <b>N-DEPENDENT</b>, and that N-dependence is the finding.
///
/// <list type="bullet">
///   <item><b>N=4 (degenerate self-fold): CONFINED.</b> The braid-carrying defective √-EP is shared across exactly
///   the D₄ orbit of (1,2) under transpose (p↔q̃) + the F89d cross-fold: {(1,2),(2,1),(2,3),(3,2)}. It does NOT
///   reach the dense half-filled core (2,2). Echoes Door C: the braid lives on the dilute orbit, not the dense
///   core.</item>
///   <item><b>N=5 (first generic): SPREADS.</b> The braid reaches a symmetric 12-sector diamond
///   {(1,2),(1,3),(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(3,3),(3,4),(4,2),(4,3)} — INCLUDING the dense core (2,2). It
///   splits into two cross-fold-conjugate families of 6, each carrying a BYTE-IDENTICAL shared eigenvalue λ (same
///   branch point, not merely the same gap).</item>
/// </list>
///
/// <para><b>The mechanism: free-fermion / AT additivity.</b> The membership is the clean rule
/// <c>diamond = {|bra−ket| = 1, both popcounts ∈ [1, N−1]} ∪ its q̃ ↦ N−q̃ cross-fold</c>. The elementary EP lives
/// in the |bra−ket| = 1 SE-DE coherence "rung"; adding a DIAGONAL mode-spectator (one more excitation on BOTH bra
/// and ket, the same single-particle mode) leaves E_bra − E_ket, and the whole EP (λ, character, gap), invariant —
/// so (1,2) ≡ (2,3) ≡ (3,4) byte-identically (gate-verified over all reference loci). Family B is the F89d
/// cross-fold image: λ_(p, N−q̃)(q̄) = −conj(λ_(p,q̃)(q)) − 2N exactly (<see cref="CrossFoldSimilarityWitness"/>), no
/// separate phenomenon. The λ VALUE is grounded in the Absorption Theorem: L(q) = A + qC with A = −2·diag(n_diff)
/// real and C the ANTI-HERMITIAN free-fermion hopping, so at the real reference loci Re λ = −2·⟨n_diff⟩_v EXACTLY
/// (machine zero), the AT Theorem-2 rate of the coalescing eigenvector. That eigenvector is a q-tuned MIXTURE of
/// the overlap (rate 2) and no-overlap (rate 6) eigenmodes, so ⟨n_diff⟩ is a weighted average between 1 and 3 and
/// the rate lies between the pure −2 and −6 lines (N=5: −4.62; N=4: the −4 midpoint); the defectiveness is the
/// eigenvector coalescence (a Jordan block), SEPARATE from the rate (at N=4 the defective EP sits at ⟨n_diff⟩=2,
/// the midpoint, yet is Jordan). What stays open is a closed form for the mixture ⟨n_diff⟩(q) and general N (the
/// octic is S₈, so no global radical form is expected).
///
/// <para><b>What is live here.</b> The census node runs <see cref="MultiSectorMonodromyCensus.Run"/> at inspect
/// time (via the AT-aware <see cref="SectorEpProbe.ProbeDefectiveAnywhere"/>), so the braid/node maps are recomputed
/// from the raw joint-popcount blocks; the additivity-embedding node re-probes the representative sectors at a real
/// reference locus and shows the byte-identity + cross-fold live. Cost: <c>Run(N)</c> dominates (~15 s at N=4, ~2
/// min at N=5). Default N=4; <c>inspect --root sectorbraid --N 5</c> runs the (slow) generic verdict. The exact
/// 12-set is pinned as the regression test <c>Census_N5_BraidSpreadsBeyondOrbit_ReachesDenseCore</c>. Reading:
/// <c>experiments/F89_MULTI_SECTOR_MONODROMY.md</c>; typed home: the multi-sector monodromy verdict claim. Siblings:
/// <see cref="CrossFoldSimilarityWitness"/> (F89d, <c>--root crossfold</c>) and <c>--root galoismonodromy</c>.</para>
/// </summary>
public sealed class SectorBraidWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private readonly int _n;

    public SectorBraidWitness(int n = 4)
    {
        if (n < 4) throw new ArgumentOutOfRangeException(nameof(n), n,
            "N must be >= 4 (the anchor; the (1,2) octic is the path-(N−1) residual and the D₄ orbit needs N >= 4).");
        _n = n;
    }

    public string DisplayName =>
        $"SectorBraidWitness (is the F89 (1,2) octic's S₈ braid confined to the (1,2) orbit or shared across popcount sectors? live at N={_n})";

    public string Summary =>
        "the S₈ braid is CONFINED to the D₄ orbit {(1,2),(2,1),(2,3),(3,2)} at N=4 (the dense core (2,2) is " +
        "braid-free) but SPREADS at N=5 to a symmetric 12-sector diamond INCLUDING the dense core (2,2), split into " +
        "two cross-fold-conjugate families of 6 sharing a BYTE-IDENTICAL eigenvalue λ (same branch point). The " +
        "N-dependence is the finding. Mechanism: free-fermion / AT additivity — the |bra−ket|=1 SE-DE rung's EP is " +
        "invariant under a diagonal mode-spectator (so (1,2)≡(2,3)≡(3,4) byte-identically), Family B is the F89d " +
        "cross-fold image λ↦−λ̄−2N; membership = {|bra−ket|=1, popcounts∈[1,N−1]} ∪ cross-fold. The λ VALUE is the " +
        "AT Theorem-2 rate of the mixed defective eigenvector: Re λ = −2·⟨n_diff⟩_v exactly at real loci (off the " +
        "INTEGER-quantized AT lines only because ⟨n_diff⟩ is a q-tuned mixture); a closed form for the mixture is open.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // Node 1 — the census verdict, recomputed live from the raw joint-popcount blocks. The clean,
            // group-label-free confined/spread signal: does the braid reach a DIAGONAL core (p = q̃, a dense
            // populational-equal bra/ket coherence)? At N=4 it does not (the core (2,2) is braid-free = CONFINED);
            // at N≥5 it does (the core sectors join = SPREAD). This sidesteps naming the acting group.
            var census = MultiSectorMonodromyCensus.Run(_n);
            var braid = census.DefectiveCoincidenceSectors;
            bool reachesDenseCore = braid.Any(s => s.PCol == s.PRow);
            var coreSectors = braid.Where(s => s.PCol == s.PRow).Select(s => $"({s.PCol},{s.PRow})").ToList();
            string braidList = string.Join(" ", braid.Select(s => $"({s.PCol},{s.PRow})"));
            yield return new InspectableNode(
                displayName: $"N={_n} census verdict: braid in {braid.Count} sector(s) — " +
                             (reachesDenseCore ? "SPREAD (reaches a dense diagonal core)" : "CONFINED (dense core braid-free)"),
                summary: $"defective (braid) sectors: {braidList}. " +
                         (reachesDenseCore
                             ? $"The braid REACHES the dense diagonal core(s) {string.Join(" ", coreSectors)} (p=q̃, populational-equal " +
                               "bra/ket): the generic SPREAD verdict — not confined to the dilute (1,2) orbit."
                             : "No diagonal (p=q̃) core sector carries the braid: the CONFINED verdict — the braid stays on the dilute " +
                               "(1,2) orbit {(1,2),(2,1),(2,3),(3,2)}, the dense half-filled core is braid-free (Door-C echo).") +
                         $" [{census.SharedNodeSectors.Count} node-only (semisimple/AT), {census.ProbedSectors.Count} probed]",
                provenance: NodeProvenance.Live);

            // Node 2 — the additivity + cross-fold, live at a real reference locus: the F89d cross-fold pairing
            // (robust at every N), the dilute node control, and (N≥5 only) the diagonal-spectator byte-identity.
            var loci = ReferenceDefectiveLoci.For(_n);                          // memoized (Run already populated the cache)
            var qStar = loci.OrderBy(q => Math.Abs(q.Imaginary)).ThenBy(q => q.Real).First();   // the "most real" locus
            var r12 = SectorEpProbe.ProbeDefectiveAnywhere(_n, 1, 2, qStar);
            var rFold = SectorEpProbe.ProbeDefectiveAnywhere(_n, 1, _n - 2, Complex.Conjugate(qStar));  // F89d partner at q̄
            var rNode = SectorEpProbe.ProbeDefectiveAnywhere(_n, 1, 1, qStar);  // |Δ|=0 dilute control (node-only)
            Complex foldPredicted = -Complex.Conjugate(r12.DefectiveCenter) - new Complex(2.0 * _n, 0);
            double foldResidual = (rFold.DefectiveCenter - foldPredicted).Magnitude;

            string diagPart;
            if (_n >= 5)
            {
                // Byte-identity across a diagonal-spectator ladder ((1,2)→(2,3)→…) is an N≥5 phenomenon: at N=4 the
                // diamond IS the degenerate self-fold (N−2 = 2), so its 4 members are cross-fold/conjugate-related,
                // not byte-identical. Here (2,3) = (1,2) + one excitation on BOTH bra & ket (the same mode).
                var rDiag = SectorEpProbe.ProbeDefectiveAnywhere(_n, 2, 3, qStar);
                double diagResidual = (r12.DefectiveCenter - rDiag.DefectiveCenter).Magnitude;
                diagPart = $"diagonal spectator (2,3) λ={FmtC(rDiag.DefectiveCenter)} (byte-identity residual " +
                           $"{diagResidual.ToString("E2", Inv)} ⟹ adding one excitation to BOTH bra & ket, same mode, leaves the EP); ";
            }
            else
            {
                diagPart = "(the diagonal-spectator byte-identity is an N≥5 phenomenon; at N=4 the 4-orbit is the " +
                           "degenerate self-fold — same coalescence gap but cross-fold/conjugate-related λ, not byte-identical); ";
            }

            yield return new InspectableNode(
                displayName: $"additivity + cross-fold gate at the real locus q={FmtC(qStar)} (N={_n})",
                summary: $"(1,2) defective λ={FmtC(r12.DefectiveCenter)}; {diagPart}" +
                         $"F89d cross-fold (1,{_n - 2}){(_n == 4 ? "=SELF" : "")} at q̄ λ={FmtC(rFold.DefectiveCenter)} vs predicted " +
                         $"−conj(λ)−2N={FmtC(foldPredicted)} (residual {foldResidual.ToString("E2", Inv)}); dilute control (1,1) " +
                         $"HasDefective={rNode.HasDefective} (a node, not a braid). Free-fermion additivity in action.",
                provenance: NodeProvenance.Live);

            // Node 3 — the membership rule + the mechanism verdict (stored: the structural characterization).
            yield return new InspectableNode(
                displayName: "the membership rule + mechanism (free-fermion / AT additivity)",
                summary: "the braid set = {|bra−ket| = 1, both popcounts ∈ [1, N−1]} ∪ its q̃↦N−q̃ cross-fold image. " +
                         "The |bra−ket|=1 edge sectors that touch popcount 0 or N (empty/full a side) carry no " +
                         "coalescence (the 'neither' set). The elementary EP is a property of the SE-DE coherence rung, " +
                         "made popcount-translation-invariant by diagonal mode-spectators (the same single-particle " +
                         "mode added to both sides, E_bra−E_ket unchanged) and cross-fold-invariant by F89d. So the " +
                         "SAME defective eigenvalue lives in every diamond sector: shared spectral content, a symmetry " +
                         "broader than the naive Klein-four.",
                provenance: NodeProvenance.Stored);

            // Node 4 — the λ VALUE from below (live at a real locus): it IS the Absorption-Theorem rate of the
            // (mixed) coalescing eigenvector. L(q) = A + qC with A = −2·diag(n_diff) real and C the ANTI-HERMITIAN
            // free-fermion hopping (the −i[H,·] part), so v^H C v is purely imaginary; at real q the coherent term
            // contributes only to Im λ, and Re λ = v^H A v = −2·⟨n_diff⟩_v EXACTLY.
            if (Math.Abs(qStar.Imaginary) < 1e-9)
            {
                var (nd, atRes) = AtRateOfDefectiveEigenvector(_n, qStar, r12.DefectiveCenter);
                yield return new InspectableNode(
                    displayName: "the λ VALUE from below: Re λ = −2·⟨n_diff⟩, the AT rate of the mixed defective eigenvector",
                    summary: $"at the real locus q={FmtC(qStar)} the (1,2) defective λ={FmtC(r12.DefectiveCenter)} has " +
                             $"Re λ = −2·⟨n_diff⟩_v with ⟨n_diff⟩_v={nd.ToString("F4", Inv)} (residual {atRes.ToString("E2", Inv)}, " +
                             "machine zero). So the λ VALUE is NOT off the AT structure: it is the AT Theorem-2 rate of a " +
                             "q-tuned MIXTURE of the overlap (rate 2) and no-overlap (rate 6) eigenmodes, so ⟨n_diff⟩ is a " +
                             "weighted average between 1 and 3 and the rate lies between the pure −2 and −6 lines (N=5: " +
                             "−4.62; N=4: the −4 midpoint), exactly on the AT theorem for the mixed eigenvector. The " +
                             "defectiveness is the eigenvector coalescence (a genuine Jordan block), SEPARATE from the rate: " +
                             "at N=4 the defective EP even sits at ⟨n_diff⟩=2 (the midpoint, like the diabolic), yet is Jordan.",
                    provenance: NodeProvenance.Live);
            }

            // Node 4b — the ⟨n_diff⟩(q*) mixture table across ALL real defective loci (Head-1 data, mined
            // from the trusted tools only: the memoized ReferenceDefectiveLoci, ProbeDefectiveAnywhere for the
            // defective cluster λ, and the AT-rate helper — no fresh scan, no reimplementation). This is the
            // per-locus evidence for the closed-form question: at N=4 the self-fold pins EVERY value to N/2 (Re λ
            // = −N is the antiunitary fixed axis of the F89d map λ↦−λ̄−2N); at N≥5 the values are the non-radical
            // P₁₀ mixtures (the octic is S₈, so no per-locus radical form — the numbers do not simplify).
            var realLoci = loci.Where(q => Math.Abs(q.Imaginary) < 1e-9).OrderBy(q => q.Real).ToList();
            var mixtureRows = new List<string>();
            foreach (var q in realLoci)
            {
                var probe = SectorEpProbe.ProbeDefectiveAnywhere(_n, 1, 2, q);
                if (!probe.HasDefective) continue;
                var (nd, res) = AtRateOfDefectiveEigenvector(_n, q, probe.DefectiveCenter);
                mixtureRows.Add($"q*={q.Real.ToString("0.######", Inv)}: λ={FmtC(probe.DefectiveCenter)}, " +
                                $"⟨n_diff⟩={nd.ToString("F6", Inv)} (AT res {res.ToString("E1", Inv)})");
            }
            yield return new InspectableNode(
                displayName: $"the ⟨n_diff⟩ mixture across all {mixtureRows.Count} real defective loci (Head-1 data, N={_n})",
                summary: (mixtureRows.Count == 0 ? "no real defective loci at this N. " : string.Join("  |  ", mixtureRows) + ". ") +
                         (_n == 4
                             ? "At N=4 the self-fold pins every value to ⟨n_diff⟩ = N/2 = 2: the F89d map λ↦−λ̄−2N is an "
                               + "antiunitary symmetry of the block itself, and Re λ = −N is its real-part fixed axis."
                             : "At N≥5 there is no self-fold, so the values are free P₁₀ mixtures; they are non-radical "
                               + "(the octic is S₈), so the OPEN target is the FUNCTION ⟨n_diff⟩(q) from the free-fermion "
                               + "mode geometry, not a per-locus closed form."),
                provenance: NodeProvenance.Live);

            // Node 5 — what remains open (stored: the honest edge, now narrowed).
            yield return new InspectableNode(
                displayName: "what is open: the closed-form mixture (the S₈ wall) and general N",
                summary: "the λ VALUE is grounded (the AT rate of the mixed eigenvector, above) and the additivity " +
                         "embedding is confirmed. What remains: a CLOSED FORM for the overlap mixture ⟨n_diff⟩(q) — " +
                         "predicting it from the free-fermion mode geometry WITHOUT diagonalizing the block — and general " +
                         "N. The defective loci are the P₁₀ roots and the octic is S₈ (non-solvable in radicals), so no " +
                         "global radical form is expected; the per-EP mixture stays numerically determined.",
                provenance: NodeProvenance.Stored);
        }
    }

    /// <summary>The Absorption-Theorem rate of the (1,2) defective coalescing eigenvector at (real) coupling
    /// <paramref name="q"/>: builds the clean L(q) = A + qC (<see cref="WeightCoherenceBlock"/>, A = −2·diag(n_diff)),
    /// takes the right-eigenvector nearest <paramref name="lambdaNear"/>, and returns its ⟨n_diff⟩_v (the overlap
    /// mixture) and the residual |−2·⟨n_diff⟩_v − Re λ| (machine zero at real q, since C is anti-Hermitian).</summary>
    private static (double NDiff, double AtResidual) AtRateOfDefectiveEigenvector(int n, Complex q, Complex lambdaNear)
    {
        var raw = WeightCoherenceBlock.Build(n, 1, 2, q);        // L = A + qC, A = −2·diag(n_diff)
        int d = raw.GetLength(0);
        var m = Matrix<Complex>.Build.DenseOfArray(raw);
        var evd = m.Evd();
        var eigs = evd.EigenValues;
        var vecs = evd.EigenVectors;

        var kets = WeightCoherenceBlock.Configs(n, 1);
        var bras = WeightCoherenceBlock.Configs(n, 2);
        var ndiff = new double[d];
        int idx = 0;
        foreach (var k in kets) foreach (var b in bras) ndiff[idx++] = System.Numerics.BitOperations.PopCount((uint)(k ^ b));

        int col = Enumerable.Range(0, d).OrderBy(i => (eigs[i] - lambdaNear).Magnitude).First();
        double num = 0, den = 0;
        for (int r = 0; r < d; r++)
        {
            var c = vecs[r, col];
            double w = c.Real * c.Real + c.Imaginary * c.Imaginary;
            num += w * ndiff[r]; den += w;
        }
        double nd = num / den;
        return (nd, Math.Abs(-2.0 * nd - eigs[col].Real));
    }

    private static string FmtC(Complex z) =>
        $"({z.Real.ToString("0.####", Inv)},{z.Imaginary.ToString("0.####", Inv)}i)";

    public InspectablePayload Payload => InspectablePayload.Empty;
}
