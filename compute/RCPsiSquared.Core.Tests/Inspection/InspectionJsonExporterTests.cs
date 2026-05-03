using System.Text.Json;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.Decomposition.Views;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Resonance;

namespace RCPsiSquared.Core.Tests.Inspection;

public class InspectionJsonExporterTests
{
    [Fact]
    public void ToJson_FourModeEffective_HasExpectedRootShape()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);
        var json = InspectionJsonExporter.ToJson(eff);
        using var doc = JsonDocument.Parse(json);
        var root = doc.RootElement;

        Assert.True(root.TryGetProperty("displayName", out var name));
        Assert.Contains("FourModeEffective", name.GetString());
        Assert.True(root.TryGetProperty("summary", out _));
        Assert.True(root.TryGetProperty("payload", out var payload));
        Assert.Equal("none", payload.GetProperty("kind").GetString());
        Assert.True(root.TryGetProperty("children", out var children));
        Assert.True(children.GetArrayLength() >= 5,
            $"FourModeEffective should have ≥5 structural children; got {children.GetArrayLength()}");
    }

    [Fact]
    public void ToJson_DiagonalRates_HasMatrixPayloadWithFourModeLabels()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);
        var json = InspectionJsonExporter.ToJson(eff.DEffView);
        using var doc = JsonDocument.Parse(json);
        var root = doc.RootElement;
        var payload = root.GetProperty("payload");

        Assert.Equal("matrix", payload.GetProperty("kind").GetString());
        var rowLabels = payload.GetProperty("rowLabels");
        Assert.Equal(4, rowLabels.GetArrayLength());
        Assert.Equal("|c_1⟩", rowLabels[0].GetString());
        Assert.Equal("|v_0⟩", rowLabels[3].GetString());

        var matrix = payload.GetProperty("matrixValues");
        Assert.Equal(4, matrix.GetArrayLength());
        Assert.Equal(-0.1, matrix[0][0].GetProperty("re").GetDouble(), 10);
        Assert.Equal(-0.3, matrix[1][1].GetProperty("re").GetDouble(), 10);
    }

    [Fact]
    public void ToJson_PerBondCoupling_PreservesBondClassMetadata()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);
        var json = InspectionJsonExporter.ToJson(eff.MhPerBondViews[0]);
        using var doc = JsonDocument.Parse(json);
        var root = doc.RootElement;
        var children = root.GetProperty("children");

        // First two children are bond index (real 0) and bond class (Endpoint summary).
        Assert.Equal("bond index", children[0].GetProperty("displayName").GetString());
        Assert.Equal(0, children[0].GetProperty("payload").GetProperty("realValue").GetDouble());
        Assert.Equal("bond class", children[1].GetProperty("displayName").GetString());
        Assert.Equal("Endpoint", children[1].GetProperty("summary").GetString());
    }

    [Fact]
    public void ToJson_LEffSweep_Exports3DTensor_QTimes4x4()
    {
        // Q × 4 × 4 = 3D structure. Add EVD per Q (eigenvalues 4-vec + eigenvectors 4×4) → 4D.
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);
        double[] qGrid = { 0.5, 1.0, 1.5, 2.0 };
        var sweep = new LEffSweepView(eff, qGrid);
        var json = InspectionJsonExporter.ToJson(sweep);
        using var doc = JsonDocument.Parse(json);
        var root = doc.RootElement;

        // sweep children: Q-at-peak-Im (real), max-Im-curve (curve), snapshots (group of N)
        var children = root.GetProperty("children");
        Assert.Equal(3, children.GetArrayLength());

        // max-Im-curve is a curve payload with X.length == qGrid.length and Y matching
        var curveNode = children[1];
        var curvePayload = curveNode.GetProperty("payload");
        Assert.Equal("curve", curvePayload.GetProperty("kind").GetString());
        Assert.Equal(qGrid.Length, curvePayload.GetProperty("x").GetArrayLength());
        Assert.Equal(qGrid.Length, curvePayload.GetProperty("y").GetArrayLength());

        // snapshots group has Q snapshots, each with L_eff (4×4) + eigenvalues (4-vec) + eigenvectors (4×4)
        var snapshots = children[2].GetProperty("children");
        Assert.Equal(qGrid.Length, snapshots.GetArrayLength());
        var firstSnap = snapshots[0];
        var snapChildren = firstSnap.GetProperty("children");
        // first child: L_eff (BlockMatrixIn4Mode) — payload kind matrix 4×4
        var lEffPayload = snapChildren[0].GetProperty("payload");
        Assert.Equal("matrix", lEffPayload.GetProperty("kind").GetString());
        Assert.Equal(4, lEffPayload.GetProperty("matrixValues").GetArrayLength());
        // second child: eigenvalues — payload kind vector 4-component
        var evalsPayload = snapChildren[1].GetProperty("payload");
        Assert.Equal("vector", evalsPayload.GetProperty("kind").GetString());
        Assert.Equal(4, evalsPayload.GetProperty("vectorValues").GetArrayLength());
        // third child: eigenvectors — payload kind matrix 4×4
        var evecsPayload = snapChildren[2].GetProperty("payload");
        Assert.Equal("matrix", evecsPayload.GetProperty("kind").GetString());
        Assert.Equal(4, evecsPayload.GetProperty("matrixValues").GetArrayLength());
    }

    [Fact]
    public void WriteToFile_RoundtripsThroughDisk()
    {
        var block = new CoherenceBlock(N: 5, n: 1, gammaZero: 0.05);
        var eff = FourModeEffective.Build(block);
        string tempPath = Path.GetTempFileName();
        try
        {
            InspectionJsonExporter.WriteToFile(eff, tempPath);
            string content = File.ReadAllText(tempPath);
            Assert.True(content.Length > 100, "exported file should have content");
            using var doc = JsonDocument.Parse(content);
            Assert.Contains("FourModeEffective",
                doc.RootElement.GetProperty("displayName").GetString());
        }
        finally
        {
            if (File.Exists(tempPath)) File.Delete(tempPath);
        }
    }
}
