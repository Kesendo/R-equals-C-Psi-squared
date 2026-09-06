using System.Diagnostics;
using System.Globalization;
using System.Numerics;
using System.Text.Json.Nodes;
using System.Text.Json;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public sealed class RouteBA2InventoryTests(ITestOutputHelper output)
{
    [Fact]
    [Trait("Category", "ROUTE_B_A2_RECONCILE")]
    public void ClassifyAll_ReconcilesEveryExactN5LocusWithLocalCharacter()
    {
        var watch = Stopwatch.StartNew();
        var inventory = RouteBA2Inventory.LoadDefault();
        var readings = new RouteBA2CharacterClassifier(inventory, new ArtifactExactRankFallback(inventory)).ClassifyAll();
        Assert.Equal(29, inventory.Roots.Count);
        Assert.Empty(inventory.ExactRankCertificates);
        Assert.Equal(58, readings.Count);
        Assert.Equal(58, readings.Select(r => r.LocusId).Distinct().Count());
        Assert.Equal(24, readings.Count(r => r.Source == A2CharacterSource.HermitianAxis));
        Assert.Equal(34, readings.Count(r => r.Source is A2CharacterSource.EpCharacter or A2CharacterSource.ExactRank));
        foreach (var root in inventory.Roots)
        foreach (var locus in root.QLoci)
        {
            var reading = Assert.Single(readings, r => r.LocusId == locus.Id);
            Assert.Equal(root.RootKind == A2RootKind.NegativeReal
                ? A2CharacterSource.HermitianAxis : A2CharacterSource.EpCharacter, reading.Source);
            Assert.Equal(2, reading.Algebraic);
            Assert.True(reading.Kind is EpCharacter.EpKind.Diabolic or EpCharacter.EpKind.Defective);
            Assert.Equal(EpCharacter.EpKind.Diabolic, reading.Kind);
            Assert.Equal(2, reading.Geometric);
            Assert.True(reading.IsolationMargin > 0);
            if (root.RootKind == A2RootKind.NegativeReal)
            {
                Assert.True(reading.FullBlockHermiticityResidual.HasValue);
                Assert.True(reading.FullBlockHermiticityResidual < 1e-12, reading.ToString());
            }
            else Assert.Null(reading.FullBlockHermiticityResidual);
        }
        double maximum = readings.Where(r => r.Source == A2CharacterSource.HermitianAxis)
            .Max(r => r.FullBlockHermiticityResidual!.Value);
        output.WriteLine($"Roots=29; unique loci=58; HermitianAxis=24; stable EpCharacter=34; ExactRank=0; " +
            $"Diabolic=58; alg=geo=2; maximum full-sector Hermiticity residual={maximum:G17}; seconds={watch.Elapsed.TotalSeconds:F3}");
    }

    [Theory]
    [InlineData(A2Parity.Even, 26)]
    [InlineData(A2Parity.Odd, 24)]
    [Trait("Category", "ROUTE_B_A2_RECONCILE")]
    public void FullSectorHermiticity_RejectsRealQOnTheSameMatrixPath(A2Parity parity, int dimension)
    {
        var inventory = RouteBA2Inventory.LoadDefault();
        var root = inventory.Roots.First(r => r.Parity == parity && r.RootKind == A2RootKind.NegativeReal);
        var q = root.QLoci[0].QPhysicalCSharpSeed;
        var block = RouteBA2CharacterClassifier.FullParityBlock(5, parity, q);
        Assert.Equal(dimension, block.RowCount);
        Assert.True(RouteBA2CharacterClassifier.HermiticityResidual(block) < 1e-12);
        var realQBlock = RouteBA2CharacterClassifier.FullParityBlock(5, parity, new Complex(q.Imaginary, 0));
        double control = RouteBA2CharacterClassifier.HermiticityResidual(realQBlock);
        Assert.True(control > 1e-2, $"Real-q control failed to break Hermiticity: {control:G17}");
        output.WriteLine($"{parity}: full dimension={dimension}; real-q Hermiticity residual={control:G17}");
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_FALLBACK")]
    public void ForcedAmbiguity_WithoutCertificate_DoesNotLoosenTolerances()
    {
        var inventory = RouteBA2Inventory.LoadDefault();
        IA2ExactRankFallback fallback = new NoExactRankFallback();
        Assert.Empty(inventory.ExactRankCertificates);
        Assert.False(new ArtifactExactRankFallback(inventory).TryGet("missing", out _));
        var classifier = RouteBA2CharacterClassifier.ForcedAmbiguousForTest(inventory, fallback);
        var root = inventory.Roots.First(r => r.RootKind == A2RootKind.Nonreal);
        var error = Assert.Throws<A2CharacterUncertifiedException>(() => classifier.Classify(root, root.QLoci[0]));
        Assert.Contains(root.QLoci[0].Id, error.Message);
        Assert.Contains("exact rank certificate absent", error.Message);
        Assert.Contains("tolerances not loosened", error.Message);
        Assert.Equal(3, error.Readings.Count);
        Assert.All(error.Readings, r => Assert.Equal(EpCharacter.EpKind.NearEp, r.Kind));
    }

    [Theory]
    [InlineData(1, EpCharacter.EpKind.Defective)]
    [InlineData(2, EpCharacter.EpKind.Diabolic)]
    [Trait("Category", "ROUTE_B_A2_FALLBACK")]
    public void ForcedAmbiguity_ConsumesSyntheticCertificate(int nullity, EpCharacter.EpKind kind)
    {
        var inventory = LoadWithCertificates(SyntheticCertificate(nullity));
        var fallback = new ArtifactExactRankFallback(inventory);
        var certificate = Assert.Single(inventory.ExactRankCertificates);
        Assert.True(fallback.TryGet(certificate.LocusId, out var found));
        Assert.Equal(certificate, found);
        Assert.False(fallback.TryGet(certificate.LocusId.ToLowerInvariant(), out _));
        var root = inventory.Roots.Single(r => r.QLoci.Any(q => q.Id == certificate.LocusId));
        var locus = root.QLoci.Single(q => q.Id == certificate.LocusId);
        var classifier = RouteBA2CharacterClassifier.ForcedAmbiguousForTest(inventory, fallback);
        var result = classifier.Classify(root, locus);
        Assert.Equal(kind, result.Kind);
        Assert.Equal(2, result.Algebraic);
        Assert.Equal(nullity, result.Geometric);
        Assert.Equal(A2CharacterSource.ExactRank, result.Source);
        Assert.Null(result.RelativeDeparture);
        Assert.Null(result.Radius);
        Assert.Null(result.IsolationMargin);
        Assert.Contains(locus.Id, classifier.FindAmbiguousNonHermitianLocusIds());
    }

    [Theory]
    [InlineData("duplicate")]
    [InlineData("unknown-locus")]
    [InlineData("hermitian-locus")]
    [InlineData("dimension")]
    [InlineData("multiplicity")]
    [InlineData("rank")]
    [InlineData("criterion")]
    [InlineData("deleted-row")]
    [InlineData("deleted-column")]
    [InlineData("missing-remainder")]
    [InlineData("zero-remainder")]
    [InlineData("missing-digest")]
    [InlineData("bad-digest")]
    [InlineData("unexpected-minor")]
    [InlineData("transport")]
    [InlineData("missing-field")]
    [InlineData("unknown-field")]
    [Trait("Category", "ROUTE_B_A2_FALLBACK")]
    public void Load_RejectsMalformedExactRankCertificate(string mutation)
    {
        var certificate = SyntheticCertificate(mutation is "missing-digest" or "bad-digest" or "unexpected-minor" ? 2 : 1);
        var node = JsonSerializer.SerializeToNode(certificate, new JsonSerializerOptions { PropertyNamingPolicy = JsonNamingPolicy.CamelCase })!;
        switch (mutation)
        {
            case "unknown-locus": node["locusId"] = "unknown"; break;
            case "hermitian-locus": node["locusId"] = RouteBA2Inventory.LoadDefault().Roots.First(r => r.RootKind == A2RootKind.NegativeReal).QLoci[0].Id; break;
            case "dimension": node["matrixDimension"] = 2; node["exactRank"] = 1; break;
            case "multiplicity": node["algebraicMultiplicity"] = 3; break;
            case "rank": node["exactRank"] = 0; break;
            case "criterion": node["criterion"] = "floating-rank"; break;
            case "deleted-row": node["deletedRow"] = -1; break;
            case "deleted-column": node["deletedColumn"] = certificate.MatrixDimension; break;
            case "missing-remainder": node["remainderPolynomial"] = null; break;
            case "zero-remainder": node["remainderPolynomial"] = "0"; break;
            case "missing-digest": node["zeroRemaindersSha256"] = null; break;
            case "bad-digest": node["zeroRemaindersSha256"] = new string('z', 64); break;
            case "unexpected-minor": node["deletedRow"] = 0; break;
            case "transport": node["transport"] = "unknown automorphism"; break;
            case "missing-field": node.AsObject().Remove("transport"); break;
            case "unknown-field": node["unverified"] = true; break;
        }
        var document = JsonNode.Parse(File.ReadAllText(ArtifactPath()))!;
        document["exactRankCertificates"] = mutation == "duplicate"
            ? new JsonArray(node, node.DeepClone()) : new JsonArray(node);
        Assert.Throws<InvalidDataException>(() => LoadDocument(document));
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_FALLBACK")]
    public void InsufficientIsolation_ConsultsFallbackAndRejectsWrongCertificateId()
    {
        const string rootId = "N5-E-A2-W-012";
        var originalRoot = RouteBA2Inventory.LoadDefault().Roots.Single(r => r.Id == rootId);
        var certificate = SyntheticCertificate(1) with { LocusId = originalRoot.QLoci[0].Id };
        var document = JsonNode.Parse(File.ReadAllText(ArtifactPath()))!;
        var rootDto = document["sectors"]![0]!["a2Roots"]!.AsArray()
            .Single(r => r!["id"]!.GetValue<string>() == rootId)!;
        rootDto["lambdaSeed"]!["real"] = originalRoot.LambdaBox.Midpoint.Real.ToString("R", CultureInfo.InvariantCulture);
        rootDto["lambdaSeed"]!["imag"] = originalRoot.LambdaBox.Midpoint.Imaginary.ToString("R", CultureInfo.InvariantCulture);
        var inventory = LoadDocument(document);
        var root = inventory.Roots.Single(r => r.Id == rootId);
        var error = Assert.Throws<A2CharacterUncertifiedException>(() =>
            new RouteBA2CharacterClassifier(inventory).Classify(root, root.QLoci[0]));
        Assert.Contains("Pair is not isolated", error.Message);
        var result = new RouteBA2CharacterClassifier(inventory, new SuppliedCertificate(certificate))
            .Classify(root, root.QLoci[0]);
        Assert.Equal(A2CharacterSource.ExactRank, result.Source);
        Assert.Equal(EpCharacter.EpKind.Defective, result.Kind);
        Assert.Null(result.Radius);
        var wrongId = certificate with { LocusId = originalRoot.QLoci[1].Id };
        Assert.Throws<InvalidDataException>(() =>
            new RouteBA2CharacterClassifier(inventory, new SuppliedCertificate(wrongId)).Classify(root, root.QLoci[0]));
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_FALLBACK")]
    public void InjectedMalformedCertificate_FailsClosed()
    {
        var inventory = RouteBA2Inventory.LoadDefault();
        var certificate = SyntheticCertificate(1) with { AlgebraicMultiplicity = 3 };
        var root = inventory.Roots.Single(r => r.QLoci.Any(q => q.Id == certificate.LocusId));
        var classifier = RouteBA2CharacterClassifier.ForcedAmbiguousForTest(inventory, new SuppliedCertificate(certificate));
        Assert.Throws<InvalidDataException>(() => classifier.Classify(root, root.QLoci[0]));
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_AMBIGUITY_MANIFEST")]
    public void AmbiguityManifest_RealReadingsAreStableAtAll34NonHermitianLoci()
    {
        var inventory = RouteBA2Inventory.LoadDefault();
        var classifier = new RouteBA2CharacterClassifier(inventory, new RejectFallbackUse());
        var ids = classifier.FindAmbiguousNonHermitianLocusIds();
        string? path = Environment.GetEnvironmentVariable("ROUTE_B_A2_AMBIGUOUS_OUT");
        if (!string.IsNullOrWhiteSpace(path)) File.WriteAllLines(path, ids);
        output.WriteLine($"Non-Hermitian loci=34; ambiguous={ids.Count}; stable={34 - ids.Count}");
        Assert.Empty(ids);
        var readings = inventory.Roots.Where(r => r.RootKind != A2RootKind.NegativeReal)
            .SelectMany(r => r.QLoci.Select(q => classifier.Classify(r, q))).ToArray();
        Assert.Equal(34, readings.Length);
        Assert.All(readings, r => Assert.Equal(A2CharacterSource.EpCharacter, r.Source));
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_AMBIGUITY_MANIFEST")]
    public void AmbiguityManifest_ForcedNearEpReportsAll34SortedIds()
    {
        var inventory = RouteBA2Inventory.LoadDefault();
        var expected = inventory.Roots.Where(r => r.RootKind != A2RootKind.NegativeReal)
            .SelectMany(r => r.QLoci).Select(q => q.Id).Order(StringComparer.Ordinal).ToArray();
        Assert.Equal(34, expected.Length);
        Assert.Equal(2, inventory.Roots.Where(r => r.RootKind == A2RootKind.PositiveReal).Sum(r => r.QLoci.Count));
        Assert.Equal(expected, RouteBA2CharacterClassifier.ForcedAmbiguousForTest(inventory, new NoExactRankFallback())
            .FindAmbiguousNonHermitianLocusIds());
        output.WriteLine("Forced ambiguity control: 34 sorted IDs, including both positive anchors and 32 nonreal loci.");
    }

    // Synthetic records test consumption; they are not algebraic evidence for the fixture locus.
    private static A2ExactRankCertificate SyntheticCertificate(int nullity)
    {
        var root = RouteBA2Inventory.LoadDefault().Roots.First(r => r.RootKind == A2RootKind.Nonreal);
        const int dimension = 26;
        return new(root.QLoci[0].Id, dimension, dimension - nullity, 2,
            nullity == 1 ? "nonzero-codimension-one-minor" : "all-codimension-one-minors-zero",
            nullity == 1 ? 0 : null, nullity == 1 ? 0 : null, nullity == 1 ? "1" : null,
            nullity == 2 ? new string('a', 64) : null, null);
    }

    private sealed class SuppliedCertificate(A2ExactRankCertificate certificate) : IA2ExactRankFallback
    {
        public bool TryGet(string locusId, out A2ExactRankCertificate result) { result = certificate; return true; }
    }

    private sealed class RejectFallbackUse : IA2ExactRankFallback
    {
        public bool TryGet(string locusId, out A2ExactRankCertificate result) =>
            throw new InvalidDataException($"A stable reading must not consult a fallback: {locusId}.");
    }

    private static RouteBA2Inventory LoadWithCertificates(params A2ExactRankCertificate[] certificates)
    {
        var document = JsonNode.Parse(File.ReadAllText(ArtifactPath()))!;
        document["exactRankCertificates"] = JsonSerializer.SerializeToNode(certificates,
            new JsonSerializerOptions { PropertyNamingPolicy = JsonNamingPolicy.CamelCase });
        return LoadDocument(document);
    }

    private static RouteBA2Inventory LoadDocument(JsonNode document)
    {
        string path = Path.Combine(Path.GetTempPath(), $"route-b-a2-{Guid.NewGuid():N}.json");
        try { File.WriteAllText(path, document.ToJsonString()); return RouteBA2Inventory.Load(path); }
        finally { File.Delete(path); }
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_SCHEMA")]
    public void LoadDefault_PreservesExactInventoryAndLocusCounts()
    {
        var inventory = RouteBA2Inventory.LoadDefault();
        Assert.Equal(1, inventory.SchemaVersion);
        Assert.Equal(5, inventory.N);
        Assert.Equal(29, inventory.Roots.Count);
        Assert.Equal(16, inventory.Roots.Count(r => r.Parity == A2Parity.Even));
        Assert.Equal(13, inventory.Roots.Count(r => r.Parity == A2Parity.Odd));
        Assert.Equal(12, inventory.Roots.Count(r => r.RootKind == A2RootKind.NegativeReal));
        Assert.Equal(1, inventory.Roots.Count(r => r.RootKind == A2RootKind.PositiveReal));
        Assert.Equal(16, inventory.Roots.Count(r => r.RootKind == A2RootKind.Nonreal));
        Assert.Equal(29, inventory.Roots.Select(r => r.Id).Distinct().Count());
        var loci = inventory.Roots.SelectMany(r => r.QLoci).ToArray();
        Assert.Equal(58, loci.Length);
        Assert.Equal(58, loci.Select(q => q.Id).Distinct().Count());
        Assert.All(inventory.Roots, r => Assert.Equal(2, r.QLoci.Count));
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_SCHEMA")]
    public void Load_UsesExportedSeedsRatherThanMidpointsOfPropagatedLambdaBoxes()
    {
        var inventory = RouteBA2Inventory.LoadDefault();
        var document = JsonNode.Parse(File.ReadAllText(ArtifactPath()))!;
        foreach (var sector in document["sectors"]!.AsArray())
        foreach (var dto in sector!["a2Roots"]!.AsArray())
        {
            var root = inventory.Roots.Single(r => r.Id == dto!["id"]!.GetValue<string>());
            Assert.Equal(ParseSeed(dto!["lambdaSeed"]!), root.LambdaSeed);
            foreach (var locusDto in dto["qLoci"]!.AsArray())
            {
                var locus = root.QLoci.Single(q => q.Id == locusDto!["id"]!.GetValue<string>());
                Assert.Equal(ParseSeed(locusDto!["qPhysicalCSharpSeed"]!), locus.QPhysicalCSharpSeed);
            }
        }
        static Complex ParseSeed(JsonNode seed) => new(
            double.Parse(seed["real"]!.GetValue<string>(), CultureInfo.InvariantCulture),
            double.Parse(seed["imag"]!.GetValue<string>(), CultureInfo.InvariantCulture));
    }

    [Theory]
    [InlineData("duplicate-root")]
    [InlineData("duplicate-locus")]
    [InlineData("denominator")]
    [InlineData("schema")]
    [InlineData("n")]
    [InlineData("count")]
    [InlineData("q-half")]
    [InlineData("q-seed")]
    [InlineData("lambda-seed")]
    [InlineData("parity")]
    [InlineData("sector-parity")]
    [InlineData("inverted-box")]
    [InlineData("convention")]
    [Trait("Category", "ROUTE_B_A2_SCHEMA")]
    public void Load_RejectsCorruptedArtifact(string mutation)
    {
        var document = JsonNode.Parse(File.ReadAllText(ArtifactPath()))!;
        var sector = document["sectors"]![0]!;
        var roots = sector["a2Roots"]!.AsArray();
        var root = roots[0]!;
        var locus = root["qLoci"]![0]!;
        switch (mutation)
        {
            case "duplicate-root": roots[1]!["id"] = root["id"]!.GetValue<string>(); break;
            case "duplicate-locus": root["qLoci"]![1]!["id"] = locus["id"]!.GetValue<string>(); break;
            case "denominator": root["wBox"]!["reLo"]!["denominator"] = "0"; break;
            case "schema": document["schemaVersion"] = 2; break;
            case "n": document["n"] = 6; break;
            case "count": roots.RemoveAt(0); break;
            case "q-half":
                var endpoint = locus["qPhysicalCSharpBox"]!["imHi"]!;
                endpoint["numerator"] = (BigInteger.Parse(endpoint["numerator"]!.GetValue<string>(), CultureInfo.InvariantCulture) + 1)
                    .ToString(CultureInfo.InvariantCulture);
                break;
            case "q-seed": locus["qPhysicalCSharpSeed"]!["real"] = "100"; break;
            case "lambda-seed": root["lambdaSeed"]!["real"] = "100"; break;
            case "parity": root["parity"] = "O"; break;
            case "sector-parity": sector["parity"] = "O"; break;
            case "inverted-box": root["wBox"]!["reLo"]!["numerator"] = "100"; break;
            case "convention": document["conventions"]!["qPhysicalCSharp"] = "qUnitHop"; break;
        }
        string path = Path.Combine(Path.GetTempPath(), $"route-b-a2-{Guid.NewGuid():N}.json");
        try
        {
            File.WriteAllText(path, document.ToJsonString());
            Assert.Throws<InvalidDataException>(() => RouteBA2Inventory.Load(path));
        }
        finally { File.Delete(path); }
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_ALL")]
    public void Classify_RejectsAnInBoxLambdaSeedThatDoesNotIsolateThePair()
    {
        const string rootId = "N5-E-A2-W-012";
        var inventory = RouteBA2Inventory.LoadDefault();
        var originalRoot = inventory.Roots.Single(r => r.Id == rootId);
        var document = JsonNode.Parse(File.ReadAllText(ArtifactPath()))!;
        var rootDto = document["sectors"]![0]!["a2Roots"]!.AsArray()
            .Single(r => r!["id"]!.GetValue<string>() == rootId)!;
        rootDto["lambdaSeed"]!["real"] = originalRoot.LambdaBox.Midpoint.Real.ToString("R", CultureInfo.InvariantCulture);
        rootDto["lambdaSeed"]!["imag"] = originalRoot.LambdaBox.Midpoint.Imaginary.ToString("R", CultureInfo.InvariantCulture);
        string path = Path.Combine(Path.GetTempPath(), $"route-b-a2-{Guid.NewGuid():N}.json");
        try
        {
            File.WriteAllText(path, document.ToJsonString());
            var mutated = RouteBA2Inventory.Load(path);
            var root = mutated.Roots.Single(r => r.Id == rootId);
            var classifier = new RouteBA2CharacterClassifier(mutated);
            var error = Assert.Throws<A2CharacterUncertifiedException>(() => classifier.Classify(root, root.QLoci[0]));
            Assert.Equal(rootId + "-Q-plus", error.LocusId);
            Assert.Contains("Pair is not isolated", error.Message);
            Assert.Empty(error.Readings);
        }
        finally { File.Delete(path); }
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_ALL")]
    public void ClassifyNonreal_All32LociHaveStableTwoDimensionalCharacter()
    {
        var watch = Stopwatch.StartNew();
        var inventory = RouteBA2Inventory.LoadDefault();
        var classifier = new RouteBA2CharacterClassifier(inventory);
        var readings = new List<A2CharacterReading>();
        var failures = new List<string>();
        foreach (var root in inventory.Roots.Where(r => r.RootKind == A2RootKind.Nonreal))
        foreach (var locus in root.QLoci)
        {
            try
            {
                var reading = classifier.Classify(root, locus);
                readings.Add(reading);
                output.WriteLine(reading.ToString());
                Assert.Equal(A2CharacterSource.EpCharacter, reading.Source);
                Assert.Equal(2, reading.Algebraic);
                Assert.True(reading.IsolationMargin > 0);
            }
            catch (A2CharacterUncertifiedException error)
            {
                failures.Add(error.Message);
                output.WriteLine(error.Message);
            }
        }
        output.WriteLine($"Stable={readings.Count}; uncertified={failures.Count}; seconds={watch.Elapsed.TotalSeconds:F3}; " +
            $"minimum margin={(readings.Count > 0 ? readings.Min(r => r.IsolationMargin) : null)}");
        Assert.True(failures.Count == 0, string.Join(Environment.NewLine, failures));
        Assert.Equal(32, readings.Count);
        Assert.Equal(32, readings.Select(r => r.LocusId).Distinct().Count());
        Assert.Equal(readings, classifier.ClassifyNonreal());
    }

    private static string ArtifactPath()
    {
        foreach (string start in new[] { AppContext.BaseDirectory, Directory.GetCurrentDirectory() })
        for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
        {
            string path = Path.Combine(directory.FullName, "simulations", "results", "route_b_a2_n5.json");
            if (File.Exists(path)) return path;
        }
        throw new FileNotFoundException("The Route B A2 inventory artifact was not found.");
    }
}
