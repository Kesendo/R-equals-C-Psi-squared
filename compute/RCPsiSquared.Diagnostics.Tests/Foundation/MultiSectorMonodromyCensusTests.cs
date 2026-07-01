using System.Collections.Generic;
using System.Linq;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Guards <see cref="MultiSectorMonodromyCensus"/>, the Layer-1 census loop: for a chain length N it
/// probes EVERY joint-popcount sector (p, q̃) at the (1,2) octic's reference DEFECTIVE loci
/// (<see cref="ReferenceDefectiveLoci"/>) via the AT-aware <see cref="SectorEpProbe.ProbeDefectiveAnywhere"/> and
/// records where a sector shares a DEFECTIVE (braid-carrying) coalescence, keeping the node-only (silent)
/// coincidences in a separate map.
/// This is the monodromy-verdict instrument: keying on <see cref="EpCharacter.EpKind.Defective"/> (a Jordan
/// block, a shared √-branch) rather than any coalescence is load-bearing, since a shared diabolic is a shared
/// NODE, never a shared braid.
///
/// <para><b>N=4 verdict (validated here).</b> The braid the (1,2) octic carries is SHARED across exactly the
/// Klein-four orbit of the reference under the F89d bit-flip similarities + transpose:
/// {(1,2), (2,1), (2,3), (3,2)} (the SE↔DE coherence blocks and their bra/ket mirrors) all read Defective at all
/// 10 reference loci with an identical coalescence gap. It is CONFINED there: it does NOT reach the half-filled
/// core (2,2) (which reads only a semisimple degeneracy). So the first verdict is "shared but localized to the
/// (1,2) Klein orbit", not "spread to the bulk".</para>
///
/// <para><b>N=5 cross-fold gate — LIVE (AT-aware probe, SLOW_MSM).</b> The gate is real: (1,3) = (SE, w_{N−2}) is
/// the branch-locus-palindrome partner of (1,2), and F89d makes it carry a Defective EP at the conjugate locus.
/// That coincidence IS real (confirmed from below in the Task-5 diagnosis). At N ≥ 5 the RAW (1,2)/(1,3) blocks
/// acquire a PERMANENT exact degeneracy (gap ~5e-15 at every q) that a global-closest-pair probe locks onto,
/// MASKING the residual defective EP (gap ~5e-4); the census now probes with
/// <see cref="SectorEpProbe.ProbeDefectiveAnywhere"/>, which clusters past the permanent pair and reaches the
/// residual √-EP, so <c>Run(5)</c> reports (1,3) Defective and passes the gate.</para></summary>
public class MultiSectorMonodromyCensusTests
{
    [Fact]
    public void Census_N4_ReferenceSectorSelfCoincides_Defective()
    {
        var c = MultiSectorMonodromyCensus.Run(4);

        // The task's core self-coincidence: (1,2) is probed at its own loci, so it trivially shares a braid.
        Assert.Contains(c.Incidences, x => x.PCol == 1 && x.PRow == 2 && x.HasDefective);
        Assert.Contains(c.DefectiveCoincidenceSectors, s => s.PCol == 1 && s.PRow == 2);
    }

    // The scientific verdict as a regression guard: the braid is SHARED across exactly the F89d Klein-four orbit
    // of the reference (SE↔DE and its bra/ket mirrors), and CONFINED there — it does not reach the half-filled
    // core (2,2). This is the first "localized vs shared" output of the census.
    [Fact]
    public void Census_N4_BraidIsSharedAcrossKleinOrbit_NotTheHalfFilledCore()
    {
        var c = MultiSectorMonodromyCensus.Run(4);
        var braid = c.DefectiveCoincidenceSectors.ToHashSet();

        foreach (var s in new[] { (1, 2), (2, 1), (2, 3), (3, 2) })
            Assert.Contains(s, braid);                          // the whole Klein-four orbit shares the braid

        Assert.DoesNotContain((2, 2), braid);                   // the braid does NOT reach the half-filled core
        Assert.Equal(4, braid.Count);                           // and nothing else: exactly the orbit
    }

    // The N=5 SCIENTIFIC VERDICT (and the cross-fold hard gate), re-enabled by the AT-aware probe (commit 50867a2).
    // Two things at once, from one Run(5):
    //   (1) The cross-fold hard gate: (1,3) = (SE, w_{N−2}) is the branch-locus-palindrome partner of (1,2), and
    //       F89d makes it carry a Defective EP at the conjugate locus. The old global-closest probe masked it behind
    //       the permanent AT degeneracy (gap ~5e-15); ProbeDefectiveAnywhere clusters past that and reaches the
    //       residual √-EP (gap ~5e-4). If this gate fails, the census pipeline is wrong.
    //   (2) The verdict N=4 → N=5 CHANGES: unlike the N=4 D₄ orbit (4 sectors, dense core braid-free), at N=5 the
    //       braid SPREADS to a symmetric 12-sector diamond — including the dense half-filled core (2,2) — that
    //       splits into two cross-fold-conjugate families of 6 carrying a byte-identical shared λ-set (same branch
    //       point, not merely the same gap). The braid is NOT confined to the (1,2) orbit at the first generic N.
    [Fact(DisplayName = "N=5 cross-fold gate + braid-spread verdict")]
    [Trait("Category", "SLOW_MSM")]
    public void Census_N5_BraidSpreadsBeyondOrbit_ReachesDenseCore()
    {
        var c = MultiSectorMonodromyCensus.Run(5);

        // (1) Cross-fold hard gate: (1,3), the conjugate-locus partner of (1,2), reads Defective.
        Assert.Contains(c.Incidences, x => x.PCol == 1 && x.PRow == 3 && x.HasDefective);

        // (2) The spread: the dense half-filled core (2,2) — braid-FREE at N=4 — now carries the braid.
        var braid = c.DefectiveCoincidenceSectors.ToHashSet();
        Assert.Contains((2, 2), braid);

        // The full N=5 braid map: the symmetric 12-sector diamond (p,q̃ ∈ 1..4, q̃ in a p-dependent band), strictly
        // larger than the N=4 D₄ orbit {(1,2),(2,1),(2,3),(3,2)}. Pinned as the verdict regression guard.
        var expected = new HashSet<(int, int)>
        {
            (1, 2), (1, 3),
            (2, 1), (2, 2), (2, 3), (2, 4),
            (3, 1), (3, 2), (3, 3), (3, 4),
            (4, 2), (4, 3),
        };
        Assert.Equal(expected, braid);
    }
}
