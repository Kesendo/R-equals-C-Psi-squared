using RCPsiSquared.Core.Confirmations;

namespace RCPsiSquared.Core.Tests.Confirmations;

public class ConfirmationsRegistryTests
{
    [Fact]
    public void Lookup_PalindromeTrichotomy_HasMarrakeshDate()
    {
        var entry = ConfirmationsRegistry.Lookup("palindrome_trichotomy");
        Assert.NotNull(entry);
        Assert.Equal("ibm_marrakesh", entry!.Machine);
        Assert.Equal("2026-04-26", entry.Date);
        Assert.Contains("trotter_n3", entry.PredictedValue);
        Assert.Contains("Δ_soft-truly=-0.722", entry.MeasuredValue);
    }

    [Fact]
    public void Lookup_UnknownName_ReturnsNull()
    {
        Assert.Null(ConfirmationsRegistry.Lookup("does_not_exist"));
    }

    [Fact]
    public void All_HasTwelveEntries()
    {
        Assert.Equal(12, ConfirmationsRegistry.All.Count);
    }

    [Fact]
    public void ListNames_ContainsKeyAnchors()
    {
        var names = ConfirmationsRegistry.ListNames().ToHashSet();
        Assert.Contains("palindrome_trichotomy", names);
        Assert.Contains("f25_cusp_trajectory", names);
        Assert.Contains("chiral_mirror_law", names);
        Assert.Contains("marrakesh_transverse_y_field_detection", names);
        Assert.Contains("lebensader_skeleton_trace_decoupling", names);
        Assert.Contains("f83_pi2_class_signature_marrakesh", names);
        Assert.Contains("f95_angle_steering_kingston_may2026", names);
    }

    [Fact]
    public void Lookup_F95AngleSteering_HasKingstonDateAndBothPairs()
    {
        var entry = ConfirmationsRegistry.Lookup("f95_angle_steering_kingston_may2026");
        Assert.NotNull(entry);
        Assert.Equal("ibm_kingston", entry!.Machine);
        Assert.Equal("2026-05-16", entry.Date);
        Assert.Contains("bxyj5yd4j", entry.JobId);
        Assert.Contains("bzklqwt7f", entry.JobId);
        Assert.Contains("Pair A_mid", entry.MeasuredValue);
        Assert.Contains("Pair B_high", entry.MeasuredValue);
        Assert.Contains("F95AngleAtQuadraticZeroPi2Inheritance", entry.FrameworkPrimitive);
    }

    [Fact]
    public void ByMachine_Marrakesh_FiltersCorrectly()
    {
        var marrakesh = ConfirmationsRegistry.ByMachine("ibm_marrakesh").ToList();
        Assert.True(marrakesh.Count >= 5);
        Assert.All(marrakesh, c => Assert.Equal("ibm_marrakesh", c.Machine));
    }

    [Fact]
    public void ByMachine_Kingston_HasCuspEntries()
    {
        var kingston = ConfirmationsRegistry.ByMachine("ibm_kingston").ToList();
        Assert.True(kingston.Count >= 2);
        Assert.Contains(kingston, c => c.Name == "f25_cusp_trajectory");
        Assert.Contains(kingston, c => c.Name == "f57_kdwell_gamma_invariance");
    }

    [Fact]
    public void ByPath_FrameworkSnapshotsPath_ReturnsBothMarrakeshN3Entries()
    {
        var hits = ConfirmationsRegistry.ByPath(new[] { 0, 1, 2 }).ToList();
        Assert.Contains(hits, c => c.Name == "pi_protected_xiz_yzzy");
        Assert.Contains(hits, c => c.Name == "lebensader_skeleton_trace_decoupling");
    }

    [Fact]
    public void ByPath_SoftBreakPath_ReturnsBothMarrakeshN3Entries()
    {
        var hits = ConfirmationsRegistry.ByPath(new[] { 48, 49, 50 }).ToList();
        Assert.Contains(hits, c => c.Name == "palindrome_trichotomy");
        Assert.Contains(hits, c => c.Name == "marrakesh_transverse_y_field_detection");
    }

    [Fact]
    public void ByPath_F83Path_ReturnsF83Entry()
    {
        var hits = ConfirmationsRegistry.ByPath(new[] { 4, 5, 6 }).ToList();
        Assert.Contains(hits, c => c.Name == "f83_pi2_class_signature_marrakesh");
    }

    [Fact]
    public void ByPath_OrderSensitive_ReversedDoesNotMatch()
    {
        // [50, 49, 48] is the same chain physically but the registered direction
        // is [48, 49, 50]; ByPath uses sequence equality, so the reverse must miss.
        var fwd = ConfirmationsRegistry.ByPath(new[] { 48, 49, 50 }).ToList();
        var rev = ConfirmationsRegistry.ByPath(new[] { 50, 49, 48 }).ToList();
        Assert.NotEmpty(fwd);
        Assert.Empty(rev);
    }

    [Fact]
    public void ByMachineAndPath_RestrictsToBackend()
    {
        var hits = ConfirmationsRegistry
            .ByMachineAndPath("ibm_marrakesh", new[] { 0, 1, 2 }).ToList();
        Assert.NotEmpty(hits);
        Assert.All(hits, c => Assert.Equal("ibm_marrakesh", c.Machine));
        Assert.Empty(ConfirmationsRegistry
            .ByMachineAndPath("ibm_torino", new[] { 0, 1, 2 }));
    }

    [Fact]
    public void ByPathOverlap_AnyQubitMatches()
    {
        // Q49 alone overlaps with both [48,49,50] paths.
        var hits = ConfirmationsRegistry.ByPathOverlap(new[] { 49 }).ToList();
        Assert.Contains(hits, c => c.Name == "palindrome_trichotomy");
        Assert.Contains(hits, c => c.Name == "marrakesh_transverse_y_field_detection");
        // Q1 overlaps with [0,1,2] paths only.
        var hits2 = ConfirmationsRegistry.ByPathOverlap(new[] { 1 }).ToList();
        Assert.Contains(hits2, c => c.Name == "pi_protected_xiz_yzzy");
        Assert.DoesNotContain(hits2, c => c.Name == "palindrome_trichotomy");
    }

    [Fact]
    public void EntriesWithoutDocumentedPath_StayNull()
    {
        // Eight of twelve have paths (five backfilled + the May-5 Kingston entry +
        // the May-8 Kingston Block-CΨ saturation entry + the May-16 F95 angle-steering
        // entry on pair A_mid [82,83]); four remain null (chiral_mirror_law,
        // f57_kdwell_gamma_invariance, bonding_mode_receiver, f25_cusp_trajectory)
        // since their paths are not unambiguously documented.
        int withPath = ConfirmationsRegistry.All.Count(c => c.QubitPath != null);
        int withoutPath = ConfirmationsRegistry.All.Count(c => c.QubitPath == null);
        Assert.Equal(8, withPath);
        Assert.Equal(4, withoutPath);
    }
}
