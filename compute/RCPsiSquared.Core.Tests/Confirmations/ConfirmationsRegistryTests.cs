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
    public void PricePairFrameworkPrimitive_NamesItsRateBook()
    {
        // The measured TOTAL coherence rate includes local T1, whereas F1's exact
        // spectral price uses only pure-Z Lindblad coefficients. Mirrored in Python by
        // test_price_pair_framework_primitive_names_its_rate_book
        // (simulations/framework/tests/registry/test_confirmations.py).
        var entry = ConfirmationsRegistry.Lookup("price_pair_locality_marrakesh_july2026");
        Assert.NotNull(entry);
        Assert.Contains("pure-Z structural context", entry!.FrameworkPrimitive);
        Assert.Contains("2·Σγ_Z,j", entry.FrameworkPrimitive);
        Assert.Contains("total fitted transverse rates", entry.FrameworkPrimitive);
        Assert.Contains("including local T1", entry.FrameworkPrimitive);
        Assert.Contains("not an exact F1 center measurement", entry.FrameworkPrimitive);
        Assert.DoesNotContain("γ_F1,j = 1/(2·T2*_j)", entry.FrameworkPrimitive);
        Assert.DoesNotContain("P = 2·Σγ_F1", entry.FrameworkPrimitive);
    }

    [Fact]
    public void All_HasTwentyFourEntries()
    {
        // Union discipline with simulations/framework/confirmations.py: both registries
        // hold the same set. Reconciled to 15 on 2026-06-08; ibm_ep_onset_may2026
        // (Kingston EP onset, 2026-05-31) added to both on 2026-06-10 makes 16;
        // f120_moment_tower_kingston_june2026 added to both on 2026-06-11 makes 17;
        // the three Torino calibration-era runs (Feb-Mar 2026) registered 2026-06-18
        // (front_matter_truth arc) make 20; price_pair_locality_marrakesh_july2026
        // (the four-run F89d-price campaign, 2026-07-04) makes 21;
        // f84_heating_leg_attribution_kingston_july2026 (the two-leg cold-bath
        // attribution, 2026-07-05) makes 22;
        // concentrator_site_contrast_kingston_july2026 (the reload flight's A-sign,
        // registered 2026-07-12 per the pre-registered verdict split) makes 23.
        Assert.Equal(24, ConfirmationsRegistry.All.Count);
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
    public void Lookup_Gamma0OffTheLever_IsFiniteGridResponseBracket_NotCarrierCalibration()
    {
        var entry = ConfirmationsRegistry.Lookup("gamma0_off_the_lever_kingston_may2026");
        Assert.NotNull(entry);
        Assert.Equal("ibm_kingston", entry!.Machine);
        Assert.Equal("2026-05-29", entry.Date);
        Assert.Contains("finite-time transfer maximum", entry.PredictedValue);
        Assert.Contains("brackets the observed overshoot change", entry.MeasuredValue);
        Assert.Contains("does not isolate", entry.Description);
        Assert.DoesNotContain("Confirms the typed", entry.Description);
        Assert.DoesNotContain("hardware-anchored", entry.Description);
        Assert.Equal(new[] { 13, 14 }, entry.QubitPath);
        // Two confirmations now sit on Kingston q13-q14 (Block-CΨ saturation + this one).
        var onPath = ConfirmationsRegistry.ByPath(new[] { 13, 14 }).Select(c => c.Name).ToList();
        Assert.Contains("gamma0_off_the_lever_kingston_may2026", onPath);
        Assert.Contains("block_cpsi_saturation_kingston_may2026", onPath);
    }

    [Fact]
    public void Lookup_IbmEpOnset_RecordsRateBookAwarePopulationHandover()
    {
        // The historical slug is retained, but the result is a population handover in
        // the runner's coherence-rate book, not a critical-damping or EP certificate.
        var entry = ConfirmationsRegistry.Lookup("ibm_ep_onset_may2026");
        Assert.NotNull(entry);
        Assert.Equal("ibm_kingston", entry!.Machine);
        Assert.Equal("2026-05-31", entry.Date);
        Assert.Contains("d8dr7dfd0j8c73f4man0", entry.JobId);
        Assert.Contains("d8drjbfd0j8c73f4mobg", entry.JobId);
        Assert.Contains("1/N = 1/3 reference level", entry.PredictedValue);
        Assert.Contains("0.28 → 0.84", entry.PredictedValue);
        Assert.Contains("0.8417853730254796", entry.MeasuredValue);
        Assert.Contains("{0.2978515625, 0.3623809814453125, 0.34356689453125, 0.4898834228515625, 0.560699462890625, 0.70330810546875}", entry.MeasuredValue);
        Assert.Contains("Q_label", entry.MeasuredValue);
        Assert.Contains("Q_Lindblad = 2 Q_label", entry.MeasuredValue);
        // The walk's spectral transition is placed, not left open: the (1,1)-block EP Q*(3) = √2
        // (Q_label ≈ 0.71) lies below the handover, and the handover moves with the probe time.
        Assert.Contains("Q*(3) = √2", entry.Description);
        Assert.Contains("Q_label ≈ 0.71", entry.Description);
        Assert.Contains("probe-time crossover", entry.Description);
        Assert.Contains("CoherenceHorizonClaim", entry.FrameworkPrimitive);
        Assert.Contains("EpCharacterWitness", entry.FrameworkPrimitive);
        Assert.Contains("K=16 twirl-simulation value", entry.MeasuredValue);
        Assert.DoesNotContain("spectral character remains open", entry.Description);
        Assert.DoesNotContain("exact-twirl", entry.MeasuredValue);
        // The placement, computed on the walk's own (1,1) block (γ = 1, hop Q, coherence decay 4γ):
        // at the two lowest flown points Q_Lindblad = 1 and 2 the slowest non-kernel mode is real and a
        // rotating pair respectively, so the EP lies between them. Both points are far from Q* = √2, where
        // the eigenvalues are well conditioned: rounding leaves |Im| near 1e-15 on the real mode, and the
        // pair at Q = 2 has |Im| = 2 exactly (roots −2 ± 2i of λ²+4λ+8), so the 1e-7 cut separates them by
        // many decades.
        Assert.True(SlowestSeBlockMode(1.0).Imaginary < 1e-7);
        Assert.True(SlowestSeBlockMode(2.0).Imaginary > 1.0);
        Assert.Contains("do not certify convergence", entry.Description);
        Assert.Contains("without a gate/readout/leakage error model", entry.Description);
        Assert.DoesNotContain("floor and the onset are clean", entry.Description, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("critical-damping", entry.PredictedValue, StringComparison.OrdinalIgnoreCase);
        Assert.DoesNotContain("ExceptionalPointClock", entry.FrameworkPrimitive);
        Assert.Equal(new[] { 13, 14, 15 }, entry.QubitPath);
    }

    /// <summary>The N=3 single-excitation (1,1) block, row-major ρ_ij: −iQ[A, ρ] − 4(ρ − diag ρ), and its
    /// slowest non-kernel mode (largest Re among |λ| > 1e-7, the largest |Im| among ties).</summary>
    private static System.Numerics.Complex SlowestSeBlockMode(double q)
    {
        const int n = 3;
        var l = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>.Build.Dense(n * n, n * n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                int r = i * n + j;
                foreach (int k in new[] { i - 1, i + 1 })
                    if (k >= 0 && k < n) l[r, k * n + j] += new System.Numerics.Complex(0, -q);
                foreach (int k in new[] { j - 1, j + 1 })
                    if (k >= 0 && k < n) l[r, i * n + k] += new System.Numerics.Complex(0, q);
                if (i != j) l[r, r] += -4.0;
            }
        var ev = l.Evd().EigenValues.Where(e => e.Magnitude > 1e-7).ToArray();
        double top = ev.Max(e => e.Real);
        var pick = ev.Where(e => Math.Abs(e.Real - top) < 1e-7).OrderByDescending(e => Math.Abs(e.Imaginary)).First();
        return new System.Numerics.Complex(pick.Real, Math.Abs(pick.Imaginary));
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
        Assert.True(kingston.Count >= 3);
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
        // price_pair_locality_marrakesh_july2026 (2026-07-04) is on Marrakesh [93,94,95]
        // (runs 2-4; run 1's [2,3,4] shares the campaign entry), making seventeen with paths.
        // f84_heating_leg_attribution_kingston_july2026 (2026-07-05) is on Kingston
        // [82,83,13], making eighteen with paths.
        // concentrator_site_contrast_kingston_july2026 (flown 2026-07-11, registered
        // 2026-07-12) is on Kingston [109,108,107,106,105], making nineteen with paths.
        // f129_standing_fringe_kingston_july2026 (2026-07-15) is on Kingston
        // [11,12,13,14,15,19,35,34], making twenty with paths.
        int withPath = ConfirmationsRegistry.All.Count(c => c.QubitPath != null);
        int withoutPath = ConfirmationsRegistry.All.Count(c => c.QubitPath == null);
        Assert.Equal(20, withPath);
        Assert.Equal(4, withoutPath);
    }
}
