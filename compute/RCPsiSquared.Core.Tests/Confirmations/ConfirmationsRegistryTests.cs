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
    public void All_HasTwentyEntries()
    {
        // Union discipline with simulations/framework/confirmations.py: both registries
        // hold the same set. Reconciled to 15 on 2026-06-08; ibm_ep_onset_may2026
        // (Kingston EP onset, 2026-05-31) added to both on 2026-06-10 makes 16;
        // f120_moment_tower_kingston_june2026 added to both on 2026-06-11 makes 17;
        // the three Torino calibration-era runs (Feb-Mar 2026) registered 2026-06-18
        // (front_matter_truth arc) make 20.
        Assert.Equal(20, ConfirmationsRegistry.All.Count);
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
        Assert.Contains("gamma_0_marrakesh_calibration", names);
        Assert.Contains("d_zero_sector_trichotomy_marrakesh", names);
        Assert.Contains("ibm_ep_onset_may2026", names);
        Assert.Contains("f120_moment_tower_kingston_june2026", names);
    }

    [Fact]
    public void Lookup_F120MomentTower_HasTheModelTestAndTheNull()
    {
        var entry = ConfirmationsRegistry.Lookup("f120_moment_tower_kingston_june2026");
        Assert.NotNull(entry);
        Assert.Equal("2026-06-11", entry!.Date);
        Assert.Equal("ibm_kingston", entry.Machine);
        Assert.Contains("d8l6c7rqv2lc73863acg", entry.JobId);
        Assert.Contains("d8l6c832d42s73cb16a0", entry.JobId);
        Assert.Contains("d8l6h03nn5bs738rmrug", entry.JobId);
        Assert.Contains("Double null HELD", entry.MeasuredValue);
        Assert.Contains("HOLDS everywhere in-situ", entry.MeasuredValue);
        Assert.Contains("EPOCH ARTIFACT", entry.MeasuredValue);
        Assert.Contains("f120_prep_split_reanalysis", entry.MeasuredValue);
        Assert.Contains("MomentTowerPumpChannelClaim", entry.FrameworkPrimitive);
        Assert.Equal("experiments/F120_MOMENT_TOWER_KINGSTON.md", entry.ExperimentDoc);
        Assert.Equal(new[] { 149, 13, 9 }, entry.QubitPath);
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
    public void Lookup_Gamma0OffTheLever_AnchorsTypedCarrierValue()
    {
        // The first direct hardware read-off of γ₀ from its only lever J, anchoring the
        // typed UniversalCarrierClaim.DefaultGammaZero = 0.05 on Kingston q13-q14.
        var entry = ConfirmationsRegistry.Lookup("gamma0_off_the_lever_kingston_may2026");
        Assert.NotNull(entry);
        Assert.Equal("ibm_kingston", entry!.Machine);
        Assert.Equal("2026-05-29", entry.Date);
        Assert.Contains("DefaultGammaZero = 0.05", entry.PredictedValue);
        Assert.Contains("frequency tracks J", entry.MeasuredValue);
        Assert.Equal(new[] { 13, 14 }, entry.QubitPath);
        // Two confirmations now sit on Kingston q13-q14 (Block-CΨ saturation + this one).
        var onPath = ConfirmationsRegistry.ByPath(new[] { 13, 14 }).Select(c => c.Name).ToList();
        Assert.Contains("gamma0_off_the_lever_kingston_may2026", onPath);
        Assert.Contains("block_cpsi_saturation_kingston_may2026", onPath);
    }

    [Fact]
    public void Lookup_IbmEpOnset_AnchorsEpFieldHardwareTable()
    {
        // The Kingston EP-onset run: the revival lifts off the 1/N floor as Q crosses
        // Q_EP. This entry is the registry anchor for the hard-coded hardware table in
        // Diagnostics/Foundation/EpField.cs (node 5, "the hardware").
        var entry = ConfirmationsRegistry.Lookup("ibm_ep_onset_may2026");
        Assert.NotNull(entry);
        Assert.Equal("ibm_kingston", entry!.Machine);
        Assert.Equal("2026-05-31", entry.Date);
        Assert.Contains("d8dr7dfd0j8c73f4man0", entry.JobId);
        Assert.Contains("d8drjbfd0j8c73f4mobg", entry.JobId);
        Assert.Contains("1/N = 1/3 equipartition floor", entry.PredictedValue);
        Assert.Contains("{0.30, 0.36, 0.34, 0.49, 0.56, 0.70}", entry.MeasuredValue);
        Assert.Contains("ExceptionalPointClock", entry.FrameworkPrimitive);
        Assert.Equal(new[] { 13, 14, 15 }, entry.QubitPath);
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
        // Sixteen of twenty have paths; four remain null (chiral_mirror_law,
        // f57_kdwell_gamma_invariance, bonding_mode_receiver, f25_cusp_trajectory)
        // since their paths are not unambiguously documented. The 2026-06-08
        // reconciliation with the Python registry added gamma_0_marrakesh_calibration
        // and d_zero_sector_trichotomy_marrakesh, both on the April-26 [48,49,50] run.
        // ibm_ep_onset_may2026 (2026-06-10) is documented on Kingston [13,14,15];
        // f120_moment_tower_kingston_june2026 (2026-06-11) on Kingston [149,13,9].
        // The three Torino runs (2026-06-18) carry single-qubit paths: q52 (×2) and q80.
        int withPath = ConfirmationsRegistry.All.Count(c => c.QubitPath != null);
        int withoutPath = ConfirmationsRegistry.All.Count(c => c.QubitPath == null);
        Assert.Equal(16, withPath);
        Assert.Equal(4, withoutPath);
    }
}
