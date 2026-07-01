using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The Multi-Sector Monodromy Census (Layer 1): the first scientific verdict on whether the S₈ braid the
/// (1,2) octic carries is CONFINED to the dilute (single-/double-excitation) sectors or SHARED across the whole
/// joint-popcount lattice. For a chain length N it iterates EVERY joint-popcount sector (p, q̃) and probes each at
/// the (1,2) octic's reference DEFECTIVE loci (<see cref="ReferenceDefectiveLoci.For"/>) via the AT-aware
/// <see cref="SectorEpProbe.ProbeDefectiveAnywhere"/>, recording where a sector shares a coalescence with the
/// reference braid.
///
/// <para><b>The load-bearing distinction (the design-review's central correction).</b> A shared coalescence is
/// keyed by its CHARACTER, not merely by "two eigenvalues met":
/// <list type="bullet">
///   <item><see cref="EpCharacter.EpKind.Defective"/> ⟹ a Jordan block, eigenvectors coalesced, a braid-carrying
///         √-branch point. A sector reading Defective at a reference locus SHARES THE BRAID (a shared branch point
///         of the S₈ monodromy). These sectors are the <see cref="CensusResult.DefectiveCoincidenceSectors"/> map:
///         the monodromy verdict.</item>
///   <item><see cref="EpCharacter.EpKind.Diabolic"/> ⟹ a semisimple repeated eigenvalue, eigenvectors independent,
///         a SILENT node (the loop about it is the identity). A sector whose near-coalescences are ALL semisimple
///         (no defective cluster) at a reference locus shares a NODE, NOT a braid. These are reported SEPARATELY in
///         <see cref="CensusResult.SharedNodeSectors"/>, never folded into the braid map.</item>
///   <item><see cref="EpCharacter.EpKind.Normal"/> / a size-1 corner ⟹ no near-coalescing (size ≥ 2) cluster at
///         that locus at all; not recorded as an incidence.</item>
/// </list></para>
///
/// <para><b>Two footings, one q-set.</b> <see cref="ReferenceDefectiveLoci"/> extracts its q-values from the (1,2)
/// symmetric-residual octic (the S₈-carrying object), while the census probes every sector's RAW joint-popcount
/// block (<see cref="SectorBlock"/>) at those SAME q-values. The reference set is conjugation-closed, which is what
/// lets the branch-locus-palindrome partner appear: the cross-fold (SE,DE) = (1,2) ↔ (SE, w_{N−2}) = (1, N−2) is an
/// EXACT antiunitary similarity (F89d, <see cref="CrossFoldSimilarityWitness"/>) that preserves the whole Jordan
/// structure, so a (1,2) DEFECTIVE EP at q maps to a (1, N−2) DEFECTIVE EP at q̄; because q̄ is in the reference
/// set, sector (1, N−2) MUST read Defective. That cross-fold coincidence is the pipeline's own correctness gate
/// (<c>MultiSectorMonodromyCensusTests</c>). At N=4 the partner w_{N−2} = w2 = DE, so the partner IS the reference
/// (1,2) block: the census self-coincidence.</para>
///
/// <para><b>Cost.</b> One <see cref="ReferenceDefectiveLoci.For"/> scan (the dominant cost; ~13 s at N=4, ~58 s at
/// N=5 with cell=0.001, 24 loci) plus (#sectors with dim &gt; 1) × (#loci) raw-block eigendecomps.</para>
///
/// <para><b>Validity — AT-aware: the residual √-EP is reached at N=4 AND N ≥ 5.</b> The verdict this instrument
/// delivers is the DEFECTIVE (braid) map. The obstacle it must see through is a PERMANENT exact degeneracy of the
/// RAW blocks: many joint-popcount blocks are permanently degenerate (a gap-0 eigenvalue pair at EVERY q, sitting
/// on the AT dissipative rate lines, an unlifted spatial/AT symmetry of the RAW footing). A global-closest-pair
/// probe locks onto that semisimple pair and reads Diabolic at any q, MASKING the residual braid-carrying √-EP that
/// lives at a DIFFERENT λ. <see cref="Run"/> therefore probes with the AT-aware
/// <see cref="SectorEpProbe.ProbeDefectiveAnywhere"/>, which single-linkage-clusters every near-coalescence and
/// characterizes each on its own contour: the permanent semisimple cluster is skipped, the residual √-EP is
/// reached, and a sector's <c>HasDefective</c> at a reference locus is the braid verdict.</para>
///
/// <para>Why the mask matters by N. At N=4 the (1,2)-orbit blocks (dim 24: (1,2),(2,1),(2,3),(3,2)) are NOT
/// permanently degenerate (generic-q gap ~0.6), so the residual defective EP (gap ~4e-5) IS already the
/// global-closest pair — the AT-aware probe agrees with the old global-closest probe there (a strict improvement,
/// guarded in <c>SectorEpProbeDefectiveAnywhereTests</c>). At N ≥ 5 those same blocks ACQUIRE a permanent
/// degeneracy (dim 50 at N=5, generic-q gap ~5e-15) that MASKS the residual defective EP (gap ~5e-4) from the old
/// probe; <see cref="SectorEpProbe.ProbeDefectiveAnywhere"/> still reaches it, so the N=5 cross-fold coincidence
/// (the F89d fold-image λ = −λ̄ − 2N carries a Defective EP with a gap identical to the (1,2) residual) is now SEEN
/// by the census, not only recoverable by hand. The DIABOLIC (node) map remains a SECONDARY channel: a member
/// flags "carries a (often permanent) semisimple near-coalescence", NOT "a genuine diabolic coincidence with the
/// reference braid".</para></summary>
public static class MultiSectorMonodromyCensus
{
    /// <summary>One census incidence at reference locus <paramref name="ReferenceQ"/> for sector
    /// (<paramref name="PCol"/>, <paramref name="PRow"/>), read by the AT-aware defective-anywhere probe.
    /// <paramref name="HasDefective"/> is true iff at least one near-coalescing cluster of the RAW block carries a
    /// genuine √-EP (a Jordan block) — a SHARED BRAID; false records a near-coalescence that is only semisimple (an
    /// AT node) — a SHARED NODE. <paramref name="MinDefectiveGap"/> is the smallest defective-cluster intra gap
    /// |λ_i − λ_j| (+∞ when <paramref name="HasDefective"/> is false); <paramref name="DefectiveCenter"/> is that
    /// cluster's midpoint λ — the shared defective eigenvalue (default when <paramref name="HasDefective"/> is
    /// false), the datum that tests whether sectors carry the SAME branch point, not merely one at the same gap;
    /// <paramref name="MultiMemberClusterCount"/> is the number of near-coalescing (size ≥ 2) clusters the probe
    /// examined at this locus. A sector with no near-coalescence at a locus is not recorded as an incidence there.</summary>
    public sealed record SectorIncidence(int PCol, int PRow, Complex ReferenceQ, bool HasDefective,
        double MinDefectiveGap, Complex DefectiveCenter, int MultiMemberClusterCount);

    /// <summary>The census for one N: every genuine-coalescence incidence, plus the two disjoint sector maps
    /// (braid-sharing vs node-sharing) and the full list of probed sectors (dim &gt; 1) so a caller can render the
    /// complete localized-vs-shared picture (a sector in neither map read NO coalescence at any reference locus).</summary>
    public sealed class CensusResult
    {
        /// <summary>Chain length.</summary>
        public int N { get; }

        /// <summary>Every (sector × reference-q) probe that read a genuine near-coalescence, with its
        /// <see cref="SectorIncidence.HasDefective"/> braid flag and defective gap. A sector may appear multiple
        /// times (one per locus).</summary>
        public IReadOnlyList<SectorIncidence> Incidences { get; }

        /// <summary>Sectors with ≥ 1 DEFECTIVE incidence at a reference locus: the SHARED-BRAID-LOCUS sectors (the
        /// monodromy verdict). Distinct, ascending by (PCol, PRow).</summary>
        public IReadOnlyList<(int PCol, int PRow)> DefectiveCoincidenceSectors { get; }

        /// <summary>Sectors with ≥ 1 node-only incidence at a reference locus (a near-coalescence whose clusters are
        /// all semisimple, HasDefective = false): a SHARED NODE, NOT a braid (kept separate on purpose). Distinct,
        /// ascending by (PCol, PRow). A sector may appear in BOTH maps if it reads defective at one locus and
        /// node-only at another. CAVEAT (see the class Validity note): this channel is confounded by PERMANENT exact
        /// degeneracies of the RAW blocks (a permanently-degenerate block carries a semisimple near-coalescence at
        /// every q), so a member here is not necessarily a genuine diabolic coincidence with the reference; at N=4
        /// every member is such a permanently-degenerate block.</summary>
        public IReadOnlyList<(int PCol, int PRow)> SharedNodeSectors { get; }

        /// <summary>Every sector actually probed (p + q̃ &gt; 0 and dim &gt; 1; the size-1 corners are skipped).
        /// Distinct, ascending. A probed sector in neither coincidence map read NO coalescence at any locus.</summary>
        public IReadOnlyList<(int PCol, int PRow)> ProbedSectors { get; }

        public CensusResult(int n, IReadOnlyList<SectorIncidence> incidences,
            IReadOnlyList<(int, int)> defectiveCoincidenceSectors, IReadOnlyList<(int, int)> sharedNodeSectors,
            IReadOnlyList<(int, int)> probedSectors)
        {
            N = n;
            Incidences = incidences;
            DefectiveCoincidenceSectors = defectiveCoincidenceSectors;
            SharedNodeSectors = sharedNodeSectors;
            ProbedSectors = probedSectors;
        }
    }

    /// <summary>Run the census at chain length <paramref name="N"/>: pull the (1,2) octic's reference defective loci
    /// (conjugation-closed, diabolic-excluded), then probe every joint-popcount sector (dim &gt; 1, p + q̃ &gt; 0) at
    /// each locus with the AT-aware <see cref="SectorEpProbe.ProbeDefectiveAnywhere"/> and record the defective
    /// (braid) / node-only coincidences. The (1,2) sector self-coincides at its own loci (expected); at N ≥ 5 the
    /// cross-fold partner (1, N−2) coincides Defective at the conjugate locus.</summary>
    /// <param name="N">Chain length (N ≥ 2: the (1,2) octic is the path-(N−1) symmetric residual).</param>
    public static CensusResult Run(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be >= 2 (the (1,2) octic is path-(N-1) >= 1).");

        var loci = ReferenceDefectiveLoci.For(N);
        var decomposition = JointPopcountSectorBuilder.Build(N);

        var incidences = new List<SectorIncidence>();
        var probed = new List<(int PCol, int PRow)>();
        var defectiveSectors = new SortedSet<(int, int)>();
        var nodeSectors = new SortedSet<(int, int)>();

        foreach (var range in decomposition.SectorRanges.OrderBy(r => r.PCol).ThenBy(r => r.PRow))
        {
            int p = range.PCol, qTilde = range.PRow;

            // Skip the vacuum (p + q̃ = 0) and the size-1 corners (C(N,0)² = 1): they carry no eigenvalue pair, so
            // SectorEpProbe guards them (MinGap = +∞, Normal). Enumerate only sectors that can host a coalescence.
            if (p + qTilde == 0) continue;
            if (range.Size <= 1) continue;

            probed.Add((p, qTilde));

            foreach (var q0 in loci)
            {
                var reading = SectorEpProbe.ProbeDefectiveAnywhere(N, p, qTilde, q0);

                // Key the verdict on HasDefective (does ANY near-coalescing cluster carry a genuine √-EP / Jordan
                // block), NOT on the global-closest pair: at N >= 5 a PERMANENT semisimple (AT) degeneracy is the
                // global-closest pair and MASKS the residual braid-carrying EP from the old Probe. A defective
                // cluster ⟹ a SHARED BRAID; a near-coalescence with NO defective cluster ⟹ a SHARED (often
                // permanent/AT) NODE, kept in the separate node map. No near-coalescence here ⟹ not an incidence.
                if (reading.HasDefective)
                {
                    incidences.Add(new SectorIncidence(p, qTilde, q0, true, reading.MinDefectiveGap, reading.DefectiveCenter, reading.MultiMemberClusterCount));
                    defectiveSectors.Add((p, qTilde));
                }
                else if (reading.MultiMemberClusterCount > 0)
                {
                    incidences.Add(new SectorIncidence(p, qTilde, q0, false, reading.MinDefectiveGap, reading.DefectiveCenter, reading.MultiMemberClusterCount));
                    nodeSectors.Add((p, qTilde));
                }
            }
        }

        return new CensusResult(N, incidences, defectiveSectors.ToList(), nodeSectors.ToList(), probed);
    }
}
