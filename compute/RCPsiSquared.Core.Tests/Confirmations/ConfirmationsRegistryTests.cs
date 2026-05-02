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
    public void All_HasNineEntries()
    {
        Assert.Equal(9, ConfirmationsRegistry.All.Count);
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
}
