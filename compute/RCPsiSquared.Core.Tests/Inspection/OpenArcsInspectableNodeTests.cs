using System.Linq;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.OpenArcs;

namespace RCPsiSquared.Core.Tests.Inspection;

public class OpenArcsInspectableNodeTests
{
    [Fact]
    public void Registry_HasFiftySixEntries()
    {
        // The literal was STALE BY TWO before 2026-08-23: 29d09d2 took the registry to 51 and
        // 0144649 to 52, neither updating this line, so this test was RED across two commits.
        // It only runs under a filter here, which is how it stayed unseen. If you add an arc,
        // this number moves in the SAME commit.
        Assert.Equal(56, OpenArcsRegistry.All.Count);   // +mirrorworld_what_is_missing (2026-09-05,
                                                        //   MirrorWorld's extension backlog: thirteen
                                                        //   registry entries that fit the genre and are
                                                        //   not home, each with its typed-carrier status;
                                                        //   the 2026-09-03 F152/F153 lesson with its
                                                        //   N = 35 overflow verified; 19 Formulas.cs
                                                        //   members with no test reference; four private
                                                        //   copies of the GF(p) arithmetic; the engine
                                                        //   survey; inherits item (6) of
                                                        //   f_registry_meets_the_typed_layer).
                                                        // +the_forced_and_the_met (2026-08-30,
                                                        //   F157's gcd seats survive only on a
                                                        //   codimension-1 Delta locus while the
                                                        //   mirror-forced centre survives every
                                                        //   Delta by a proof; which of our
                                                        //   arithmetic laws are forced and which
                                                        //   are met is now a sweep with a queue).
                                                        // +noise_origin_after_candidate_two
                                                        //   (2026-08-29, INCOMPLETENESS_PROOF's
                                                        //   five-candidate elimination lost the two
                                                        //   candidates that carried it; what is left
                                                        //   is the trace identity, which says OPEN
                                                        //   and not EXTERNAL).
                                                        // +one_word_two_objects (2026-08-23, one
                                                        //   word for two objects: the V-shape
                                                        //   orientation inverted between GLOSSARY
                                                        //   and GAMMA_CONTROL, the noise/light
                                                        //   vocabulary split, the 21.5x
                                                        //   attribution; no numeric gate can see
                                                        //   any of them),
                                                        // +polarity_break_dephase_siblings
                                                        //   (2026-08-19, the X- and Y-dephase
                                                        //   siblings of F155; RETIRED the same day,
                                                        //   resolved by PROOF_F155 section (g)),
                                                        // +corner_beat (2026-08-19, PARKED by Tom
                                                        //   before submission; the arc records the
                                                        //   wake-up condition, not a closure),
                                                        // +gamma_book_enforced_nowhere (2026-08-18),
                                                        // +compressed_density_laws (2026-08-15),
                                                        // +f_registry_meets_the_typed_layer (2026-08-13,
                                                        //   a0feaba opened the arc and did not move this
                                                        //   count, so the gate was RED from that commit
                                                        //   until 2026-08-14),
                                                        // +gamma_is_the_sender_not_the_watching (2026-08-08),
                                                        // +sideways_spin_ladder (2026-08-07),
                                                        // +bit_exact_vocabulary (2026-08-06),
                                                        // +f138_converse_failures (2026-08-03),
                                                        // +f138_clause_two_sweep (2026-08-03),
                                                        // +benzene_center_tier_upgrade (2026-08-03),
                                                        // +two_coast_classifier_repair,
                                                        // +ninety_degree_family_bookkeeping,
                                                        // +f91_scope_fences,
                                                        // +f115_valuation_gate_width,
                                                        // +concentrator_amplitude_signs,
                                                        // +f55_cavity_tier_chain (2026-08-01),
                                                        // +site_resolved_vacuum_block (2026-08-02),
                                                        // +unfalsifiable_verification_gates (2026-08-02)
    }

    [Fact]
    public void Registry_OpenAndRetiredCountsPartitionTheRegistry()
    {
        int retired = OpenArcsRegistry.All.Count(a => a.Status == OpenArcStatus.Retired);
        Assert.Equal(OpenArcsRegistry.All.Count, OpenArcsRegistry.OpenCount + retired);
    }

    [Fact]
    public void Registry_ArcsParkedInProseCarryTheTypedParkedReason()
    {
        // The defect this guards: "parked" used to live ONLY as prose inside arcs whose status is
        // Open, so OpenCount reported them as work there is to do. It fails if a future session
        // parks an arc in its prose and leaves the typed field null, which is exactly what had
        // happened to both of the arcs it currently finds.
        var parkedInProse = OpenArcsRegistry.All
            .Where(a => (a.ParkedAt + a.NextStep + a.Origin).Contains("PARKED BY TOM")
                     || (a.ParkedAt + a.NextStep + a.Origin).Contains("ARC PARKED"));

        Assert.NotEmpty(parkedInProse);
        Assert.All(parkedInProse, a =>
        {
            Assert.Equal(OpenArcStatus.Open, a.Status);
            Assert.False(string.IsNullOrWhiteSpace(a.ParkedReason));
        });
    }

    [Fact]
    public void Registry_LiveCountIsOpenMinusParked_AndIsStrictlySmaller()
    {
        Assert.Equal(OpenArcsRegistry.OpenCount - OpenArcsRegistry.ParkedCount,
                     OpenArcsRegistry.LiveCount);
        // Anti-vacuity: if nothing were parked the identity above would hold trivially and this
        // whole distinction would be dead weight.
        Assert.True(OpenArcsRegistry.ParkedCount > 0);
        Assert.True(OpenArcsRegistry.LiveCount < OpenArcsRegistry.OpenCount);
    }

    [Fact]
    public void Registry_OnlyOpenArcsCarryAParkedReason()
    {
        Assert.All(
            OpenArcsRegistry.All.Where(a => a.ParkedReason is not null),
            a => Assert.Equal(OpenArcStatus.Open, a.Status));
    }

    [Fact]
    public void ArcsNode_RendersTheParkedReasonForAParkedArc()
    {
        var parked = OpenArcsRegistry.All.First(a => a.ParkedReason is not null);
        var child = OpenArcsInspectableNode.Build().Children
            .Single(c => c.DisplayName == parked.Name);
        Assert.Contains("PARKED:", child.Summary);
        // The reason must actually reach the reader, not just exist on the record.
        Assert.Contains(parked.ParkedReason!.Substring(0, 40), child.Summary);
    }

    [Fact]
    public void Registry_RetiredArcsCarryAReason()
    {
        Assert.All(
            OpenArcsRegistry.All.Where(a => a.Status == OpenArcStatus.Retired),
            a => Assert.False(string.IsNullOrWhiteSpace(a.RetiredReason)));
    }

    [Fact]
    public void OneDiagonalMirrorGroup_IsRetired_SpunOutToLinearS3()
    {
        var entry = OpenArcsRegistry.Lookup("one_diagonal_mirror_group");
        Assert.NotNull(entry);
        Assert.Equal(OpenArcStatus.Retired, entry!.Status);
        Assert.Contains("linear_s3_mirror_completion", entry.RetiredReason!);
    }

    [Fact]
    public void StrangerDoor_IsRetired_FifthDoorHung()
    {
        var entry = OpenArcsRegistry.Lookup("stranger_door");
        Assert.NotNull(entry);
        Assert.Equal(OpenArcStatus.Retired, entry!.Status);
        Assert.Contains("provenance badge", entry.RetiredReason!);
    }

    [Fact]
    public void Lookup_RoundtripsByName()
    {
        var entry = OpenArcsRegistry.Lookup("witness_coverage");
        Assert.NotNull(entry);
        Assert.Equal("witness_coverage", entry!.Name);
        Assert.Equal("2026-06-11", entry.Opened);
    }

    [Fact]
    public void Lookup_UnknownReturnsNull()
    {
        Assert.Null(OpenArcsRegistry.Lookup("no_such_arc"));
    }

    [Fact]
    public void Build_HasExactlyOneChildPerRegistryEntry()
    {
        var node = OpenArcsInspectableNode.Build();
        Assert.Equal(OpenArcsRegistry.All.Count, node.Children.Count());
    }

    [Fact]
    public void Build_ChildSummaryCarriesParkedAndNextStep()
    {
        var node = OpenArcsInspectableNode.Build();
        // Decoupled from any specific arc name: the exemplar is whichever arc is still Open,
        // so retiring an arc (e.g. birth_canal_surface, 2026-06-29) cannot re-break this test.
        var entry = OpenArcsRegistry.All.First(a => a.Status == OpenArcStatus.Open);
        var child = node.Children.Single(c => c.DisplayName == entry.Name);
        Assert.Equal($"{entry.Opened} parked: {entry.ParkedAt} -> next: {entry.NextStep}", child.Summary);
    }

    [Fact]
    public void Build_RootSummaryReportsOpenAndRetiredCounts()
    {
        var node = OpenArcsInspectableNode.Build();
        int retired = OpenArcsRegistry.All.Count(a => a.Status == OpenArcStatus.Retired);
        Assert.Contains($"{OpenArcsRegistry.OpenCount} open arc(s)", node.Summary);
        Assert.Contains($"{retired} retired", node.Summary);
    }
}
