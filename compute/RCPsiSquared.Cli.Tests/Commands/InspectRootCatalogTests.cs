using RCPsiSquared.Cli;
using RCPsiSquared.Cli.Commands;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Cli.Tests.Commands;

public class InspectRootCatalogTests
{
    // The roots that were hardcoded in the InspectCommand switch before the catalog refactor.
    private static readonly string[] PreviouslyHardcodedRoots =
    {
        "f71", "f1", "f87", "pi2", "mirror", "flow", "between",
        "fourmode", "f86", "c2hwhm", "c2cpsi", "c2cpsi-scan",
    };

    [Fact]
    public void Catalog_ContainsEveryPreviouslyHardcodedRoot()
    {
        var names = InspectCommand.Catalog.Select(e => e.Name).ToHashSet();
        foreach (var root in PreviouslyHardcodedRoots)
            Assert.Contains(root, names);
    }

    [Fact]
    public void Catalog_AddsTheWorldRoot()
    {
        var names = InspectCommand.Catalog.Select(e => e.Name).ToHashSet();
        Assert.Contains("world", names);
    }

    [Fact]
    public void Catalog_DefaultRootIsFourmode_WithDepthFour()
    {
        var fourmode = InspectCommand.Catalog.Single(e => e.Name == "fourmode");
        Assert.Equal(4, fourmode.DefaultDepth);
    }

    [Fact]
    public void Catalog_F86KeepsDefaultDepthOne()
    {
        var f86 = InspectCommand.Catalog.Single(e => e.Name == "f86");
        Assert.Equal(1, f86.DefaultDepth);
    }

    [Fact]
    public void Catalog_WorldDefaultDepthIsTwo()
    {
        var world = InspectCommand.Catalog.Single(e => e.Name == "world");
        Assert.Equal(2, world.DefaultDepth);
    }

    [Fact]
    public void Catalog_NamesAreUnique()
    {
        var names = InspectCommand.Catalog.Select(e => e.Name).ToList();
        Assert.Equal(names.Count, names.Distinct().Count());
    }

    [Fact]
    public void Catalog_EveryEntryHasADescription()
    {
        Assert.All(InspectCommand.Catalog, e => Assert.False(string.IsNullOrWhiteSpace(e.Description)));
    }

    [Fact]
    public void Catalog_AddsTheQuditRoot_WhichDoesNotRequireN()
    {
        var qudit = InspectCommand.Catalog.Single(e => e.Name == "qudit");
        Assert.Equal("F121 qudit partial palindrome, recomputed live", qudit.Description);
        Assert.False(qudit.RequiresN);
    }

    [Fact]
    public void Catalog_QuditFactory_BuildsTheLiveWitness_DefaultD3N2()
    {
        var qudit = InspectCommand.Catalog.Single(e => e.Name == "qudit");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = qudit.Factory(ctx);
        Assert.Contains("d=3", root.Summary);
        Assert.Contains("ceiling met", root.Summary);
    }

    [Fact]
    public void Catalog_WorldAndQudit_DoNotRequireN()
    {
        Assert.False(InspectCommand.Catalog.Single(e => e.Name == "world").RequiresN);
        Assert.False(InspectCommand.Catalog.Single(e => e.Name == "qudit").RequiresN);
        Assert.True(InspectCommand.Catalog.Single(e => e.Name == "fourmode").RequiresN);
    }

    [Fact]
    public void Catalog_AddsTheGlossaryRoot_WhichDoesNotRequireN()
    {
        var glossary = InspectCommand.Catalog.Single(e => e.Name == "glossary");
        Assert.False(glossary.RequiresN);
        Assert.False(string.IsNullOrWhiteSpace(glossary.Description));
    }

    // The load-bearing house terms the glossary root must define for a stranger.
    private static readonly string[] ExpectedGlossaryTerms =
    {
        "claim", "tier", "confirmation", "witness", "arc", "F-number",
        "palindrome / Π", "light / lens", "truly / soft / hard", "Q and γ", "--N and --n",
    };

    [Fact]
    public void Glossary_DefinesEveryLoadBearingTerm_WithNonEmptyDefinitions()
    {
        var glossary = InspectCommand.Catalog.Single(e => e.Name == "glossary");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = glossary.Factory(ctx);

        var children = root.Children.ToList();
        var names = children.Select(c => c.DisplayName).ToList();
        foreach (var term in ExpectedGlossaryTerms)
            Assert.Contains(term, names);

        // Every term carries a real, non-trivial definition (a stranger must be able to read it).
        Assert.All(children, c => Assert.False(string.IsNullOrWhiteSpace(c.Summary)));
        Assert.All(children, c => Assert.True(c.Summary.Length >= 40));
    }

    [Fact]
    public void Glossary_RootSummary_PointsAStrangerOnward()
    {
        var glossary = InspectCommand.Catalog.Single(e => e.Name == "glossary");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = glossary.Factory(ctx);
        Assert.Contains("--root world", root.Summary);
    }

    [Fact]
    public void Catalog_HasEnvelopeRoot_WithDescription()
    {
        var envelope = InspectCommand.Catalog.Single(e => e.Name == "envelope");
        Assert.False(string.IsNullOrWhiteSpace(envelope.Description));
    }

    [Fact]
    public void Catalog_EnvelopeFactory_BuildsTheLiveWitness()
    {
        var envelope = InspectCommand.Catalog.Single(e => e.Name == "envelope");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 3,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = envelope.Factory(ctx);
        Assert.IsType<EnvelopeTheoremWitness>(root);
        Assert.Contains("Envelope Theorem", root.Summary);
    }

    [Fact]
    public void Catalog_HasEpCharacterRoot_NFree()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "epcharacter");
        Assert.False(entry.RequiresN);
        Assert.Contains("artifact-free", entry.Description);
        Assert.Contains("DEFECTIVE", entry.Description);
    }

    [Fact]
    public void Catalog_EpCharacterFactory_BuildsTheLiveWitness_GatePassesDefectiveHorizon()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "epcharacter");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = entry.Factory(ctx);
        Assert.IsType<EpCharacterWitness>(root);
        Assert.Contains("GATE PASSED", root.Summary);
        Assert.Contains("Defective", root.Summary);
    }

    [Fact]
    public void Catalog_HasF89OcticRoot_NFree()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "f89octic");
        Assert.False(entry.RequiresN);
        Assert.Contains("DIABOLIC", entry.Description);
        Assert.Contains("epcharacter", entry.Description);
    }

    [Fact]
    public void Catalog_F89OcticFactory_BuildsTheLiveWitness_ReadsDiabolic()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "f89octic");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = entry.Factory(ctx);
        Assert.IsType<F89OcticCharacterWitness>(root);
        Assert.Contains("Diabolic", root.Summary);
    }

    [Fact]
    public void Catalog_HasF89GaloisRoot_NFree()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "f89galois");
        Assert.False(entry.RequiresN);
        Assert.Contains("S_d", entry.Description);
    }

    [Fact]
    public void Catalog_F89GaloisFactory_BuildsThePathKWitness_ReadsSd()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "f89galois");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = entry.Factory(ctx);
        // f89galois is now the path-3..6 composite F89PathKGaloisWitness (F89OcticGaloisWitness is its path-3 child).
        Assert.IsType<F89PathKGaloisWitness>(root);
        Assert.Contains("S_d", root.Summary);
    }

    [Fact]
    public void Catalog_HasMonodromyMirrorRoot_NFree()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "monodromymirror");
        Assert.False(entry.RequiresN);
        Assert.Contains("Z(S_8)=1", entry.Description);
        Assert.Contains("galoismonodromy", entry.Description);
    }

    [Fact]
    public void Catalog_MonodromyMirrorFactory_BuildsTheLiveWitness()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "monodromymirror");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = entry.Factory(ctx);
        Assert.IsType<RCPsiSquared.Diagnostics.Foundation.MonodromyMirrorWitness>(root);
    }

    [Fact]
    public void Catalog_HasDiabolicParityRoot_NFree()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "diabolicparity");
        Assert.False(entry.RequiresN);
        Assert.Contains("dimension-mismatch", entry.Description);
        Assert.Contains("sector-swap", entry.Description);
    }

    [Fact]
    public void Catalog_DiabolicParityFactory_BuildsTheLiveWitness_StatesParityVerdict()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "diabolicparity");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = entry.Factory(ctx);
        Assert.IsType<DiabolicReflectionParityWitness>(root);
        Assert.Contains("odd", root.Summary, StringComparison.OrdinalIgnoreCase);
    }

    [Fact]
    public void Catalog_HasDecoderRoot_NFree()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "decoder");
        Assert.False(entry.RequiresN);
        Assert.Contains("Q-factor", entry.Description);
    }

    [Fact]
    public void DecoderRoot_HonorsOptionalN()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "decoder");
        Assert.False(entry.RequiresN);
        Assert.True(entry.HonorsOptionalN);
    }

    [Fact]
    public void QuditRoot_IgnoresOptionalN_WarningStaysCorrect()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "qudit");
        Assert.False(entry.RequiresN);
        Assert.False(entry.HonorsOptionalN);
    }

    [Fact]
    public void Catalog_HasBlockSpectrumRoot_NFree_HonorsOptionalN()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "blockspectrum");
        Assert.False(entry.RequiresN);
        Assert.True(entry.HonorsOptionalN);
        Assert.Contains("block_spectrum_n9", entry.Description);
    }

    [Fact]
    public void Catalog_BlockSpectrumFactory_BuildsTheLiveWitness_DefaultN6()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "blockspectrum");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = entry.Factory(ctx);
        Assert.IsType<BlockSpectrumWitness>(root);
        Assert.Contains("N=6", root.DisplayName);
    }

    [Fact]
    public void Catalog_HasSectorBraidRoot_NFree_HonorsOptionalN()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "sectorbraid");
        Assert.False(entry.RequiresN);
        Assert.True(entry.HonorsOptionalN);
        Assert.Contains("Multi-Sector Monodromy", entry.Description);
        Assert.Contains("N-DEPENDENT", entry.Description);
    }

    [Fact]
    public void Catalog_SectorBraidFactory_BuildsTheLiveWitness_DefaultN4()
    {
        var entry = InspectCommand.Catalog.Single(e => e.Name == "sectorbraid");
        var ctx = new InspectRootContext(new ArgParser(Array.Empty<string>()), N: 1,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = entry.Factory(ctx);   // no --N flag ⟹ the witness defaults to N=4 (the fast anchor)
        Assert.IsType<SectorBraidWitness>(root);
        Assert.Contains("N=4", root.DisplayName);
        Assert.Contains("CONFINED", root.Summary);   // the Summary states the full N-dependent verdict (static, no census run)
        Assert.Contains("SPREAD", root.Summary);
    }

    [Fact]
    public void SymphonyFactory_TempoRatio_GrowsTheClockMovement()
    {
        var symphony = InspectCommand.Catalog.Single(e => e.Name == "symphony");
        var ctx = new InspectRootContext(new ArgParser(new[] { "--N", "3", "--tempo-ratio", "20" }), N: 3,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var root = symphony.Factory(ctx);
        var labels = root.Children.Select(c => c.DisplayName).ToList();
        Assert.Contains("movement: the clock", labels);
    }

    [Fact]
    public void SymphonyFactory_Calibrate_GrowsTheSeamMovement()
    {
        var symphony = InspectCommand.Catalog.Single(e => e.Name == "symphony");

        // With --calibrate, the seam movement appears (N=3 defaults J=1, γ=0.1 ⟹ Q=10, in the protected regime).
        var ctxWith = new InspectRootContext(new ArgParser(new[] { "--N", "3", "--calibrate" }), N: 3,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var rootWith = symphony.Factory(ctxWith);
        var labelsWith = rootWith.Children.Select(c => c.DisplayName).ToList();
        Assert.Contains("movement: the seam", labelsWith);

        // Without --calibrate, no seam movement is grown.
        var ctxWithout = new InspectRootContext(new ArgParser(new[] { "--N", "3" }), N: 3,
            WithQSweep: false, WithMeasured: false, QGridPoints: null);
        var rootWithout = symphony.Factory(ctxWithout);
        var labelsWithout = rootWithout.Children.Select(c => c.DisplayName).ToList();
        Assert.DoesNotContain("movement: the seam", labelsWithout);
    }
}
