using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public sealed class RouteBA2N6AtlasManifestTests
{
    private const string SourceSha256 =
        "c4b19015a750e56a55b0c8fdba3d05534c0c133bc5355a74eb89b77020532d47";

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_ATLAS")]
    public void BuildAtlasManifest_BindsTheExecutedInventoryCharacterAndPartnerOrbits()
    {
        var inventory = RouteBA2N6Inventory.LoadDefault();
        var manifest = new RouteBA2N6CharacterClassifier(inventory).BuildAtlasManifest();

        Assert.Equal(1, manifest.SchemaVersion);
        Assert.Equal(6, manifest.N);
        Assert.Equal("simulations/results/route_b_a2_n6.json", manifest.SourceArtifact);
        Assert.Equal(SourceSha256, manifest.SourceSha256);
        Assert.Equal("exact t-box midpoint converted to display double; not an algebraic root",
            manifest.CoordinateMeaning);
        Assert.Equal(266, manifest.Loci.Count);
        Assert.Equal(96, manifest.Orbits.Count);
        Assert.Equal(59, manifest.Orbits.Count(orbit => orbit.MemberIds.Count == 2));
        Assert.Equal(37, manifest.Orbits.Count(orbit => orbit.MemberIds.Count == 4));
        Assert.All(manifest.Loci, locus => Assert.Equal(3, locus.Contours.Count));
        Assert.Equal(133, manifest.Loci.Count(locus => locus.Parity == "E"));
        Assert.Equal(133, manifest.Loci.Count(locus => locus.Parity == "O"));
        Assert.Equal(118, manifest.Loci.Count(locus => locus.CharacterSource == "HermitianAxis"));
        Assert.Equal(148, manifest.Loci.Count(locus => locus.CharacterSource == "EpCharacterStable"));
        Assert.All(manifest.Loci, locus =>
        {
            Assert.Equal("Diabolic", locus.Verdict);
            Assert.Equal(2, locus.AlgebraicMultiplicity);
            Assert.Equal(2, locus.GeometricMultiplicity);
        });
        Assert.Equal(266, manifest.Reconciliation.ConsumedLoci);
        Assert.Equal(0, manifest.Reconciliation.UnresolvedLoci);

        var byId = manifest.Loci.ToDictionary(locus => locus.Id, StringComparer.Ordinal);
        Assert.All(manifest.Orbits, orbit =>
        {
            var members = orbit.MemberIds.ToHashSet(StringComparer.Ordinal);
            Assert.All(members, id =>
            {
                Assert.True(byId.ContainsKey(id));
                Assert.Contains(byId[id].ConjugationPartnerId, members);
                Assert.Contains(byId[id].ParityPartnerId, members);
            });
        });
        Assert.Equal(manifest.Loci.Select(locus => locus.Id).Order(StringComparer.Ordinal),
            manifest.Orbits.SelectMany(orbit => orbit.MemberIds).Order(StringComparer.Ordinal));
    }
}
