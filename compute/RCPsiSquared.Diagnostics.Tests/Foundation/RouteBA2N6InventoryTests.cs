using System.Numerics;
using System.Text.Json.Nodes;
using System.Globalization;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit.Abstractions;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

public sealed class RouteBA2N6InventoryTests(ITestOutputHelper output)
{
    private const double N6MinimumIsolationMargin = 0.06589058580248003;

    private static void AssertFixedN6Anchors(RouteBA2N6ReconciliationReport report, bool perParity = false)
    {
        // Fixed anchors from the executed schema-3 N6 reconciliation, not a
        // second classifier run. Self-consistent source/verdict drift must fail.
        int total = perParity ? 133 : 266;
        Assert.Equal(total, report.TotalLoci);
        Assert.Equal(total, report.ConsumedLoci);
        Assert.Equal(0, report.UnresolvedLoci);
        Assert.Equal(total, report.ByAlgebraicSource["ExactAlgebraic"]);
        Assert.Equal(perParity ? 59 : 118, report.ByCharacterSource["HermitianAxis"]);
        Assert.Equal(perParity ? 74 : 148, report.ByCharacterSource["EpCharacterStable"]);
        Assert.Equal(0, report.ByCharacterSource["ExactRankExecuted"]);
        Assert.Equal(total, report.ByVerdict[EpCharacter.EpKind.Diabolic]);
        Assert.Equal(0, report.ByVerdict.TryGetValue(EpCharacter.EpKind.Defective, out int defective) ? defective : 0);
        if (!perParity)
        {
            // Absolute 1e-12 allows only eigensolver rounding noise in this
            // deterministic middle-contour margin; the +1e-8 control fails.
            Assert.InRange(Math.Abs(report.MinimumIsolationMargin - N6MinimumIsolationMargin), 0.0, 1e-12);
        }
        Assert.Equal(report.TotalLoci, report.ConsumedLoci);
        Assert.Equal(report.TotalLoci, report.ByAlgebraicSource.Values.Sum());
        Assert.Equal(report.TotalLoci, report.ByCharacterSource.Values.Sum());
        Assert.Equal(report.TotalLoci, report.ByVerdict.Values.Sum());
    }

    [Theory]
    [InlineData("total", false)]
    [InlineData("total", true)]
    [InlineData("unresolved", false)]
    [InlineData("algebraic", false)]
    [InlineData("algebraic", true)]
    [InlineData("hermitian-numerical", false)]
    [InlineData("hermitian-numerical", true)]
    [InlineData("exact-rank", false)]
    [InlineData("exact-rank", true)]
    [InlineData("verdict", false)]
    [InlineData("verdict", true)]
    [InlineData("margin", false)]
    [Trait("Category", "ROUTE_B_A2_N6_RECONCILE")]
    public void FixedN6Anchors_RejectSelfConsistentDrift(string mutation, bool perParity)
    {
        // Synthetic report fixtures test the fixed regression oracle, not physics.
        // Counts are shifted together so aggregate agreement cannot expose the drift.
        int total = perParity ? 133 : 266;
        var algebraic = new Dictionary<string, int> { ["ExactAlgebraic"] = total };
        var sources = new Dictionary<string, int>
        {
            ["HermitianAxis"] = perParity ? 59 : 118,
            ["EpCharacterStable"] = perParity ? 74 : 148, ["ExactRankExecuted"] = 0
        };
        var verdicts = new Dictionary<EpCharacter.EpKind, int> { [EpCharacter.EpKind.Diabolic] = total };
        var report = new RouteBA2N6ReconciliationReport(total, total, 0,
            algebraic, sources, verdicts, N6MinimumIsolationMargin);
        AssertFixedN6Anchors(report, perParity);
        switch (mutation)
        {
            case "total":
                report = report with { TotalLoci = total + 1, ConsumedLoci = total + 1 };
                algebraic["ExactAlgebraic"]++; sources["EpCharacterStable"]++; verdicts[EpCharacter.EpKind.Diabolic]++;
                break;
            case "unresolved": report = report with { UnresolvedLoci = 1 }; break;
            case "algebraic": algebraic["ExactAlgebraic"]--; algebraic["WrongAlgebraic"] = 1; break;
            case "hermitian-numerical": sources["HermitianAxis"]--; sources["EpCharacterStable"]++; break;
            case "exact-rank": sources["EpCharacterStable"]--; sources["ExactRankExecuted"]++; break;
            case "verdict": verdicts[EpCharacter.EpKind.Diabolic]--; verdicts[EpCharacter.EpKind.Defective] = 1; break;
            case "margin": report = report with { MinimumIsolationMargin = N6MinimumIsolationMargin + 1e-8 }; break;
            default: throw new ArgumentOutOfRangeException(nameof(mutation));
        }
        Assert.NotNull(Record.Exception(() => AssertFixedN6Anchors(report, perParity)));
    }

    [Theory]
    [InlineData(0.25, "0.25")]
    [InlineData(double.PositiveInfinity, "\"not-applicable\"")]
    [Trait("Category", "ROUTE_B_A2_N6_RECONCILE")]
    public void ReconciliationReport_FormatsAbsentNumericalMarginWithoutJsonNonfinite(double margin, string expected)
    {
        // DTO output control, not an inventory: finite data must remain numeric,
        // while the all-structural-character sentinel must serialize as not-applicable.
        var report = new RouteBA2N6ReconciliationReport(0, 0, 0,
            new Dictionary<string, int>(), new Dictionary<string, int>(),
            new Dictionary<EpCharacter.EpKind, int>(), margin);
        string? serialized = null;
        var error = Record.Exception(() => serialized = System.Text.Json.JsonSerializer.Serialize(report));
        Assert.True(error == null, $"N6 reconciliation output must encode absent numerical margin as not-applicable: {error?.Message}");
        Assert.Equal(expected, JsonNode.Parse(serialized!)!["MinimumIsolationMargin"]!.ToJsonString());
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_RECONCILE")]
    public void ReconcileAll_ConsumesCertifiedDegreesAndReportsExecutedProvenance()
    {
        var inventory = RouteBA2N6Inventory.LoadDefault();
        var classifier = new RouteBA2N6CharacterClassifier(inventory);
        RouteBA2N6ReconciliationReport report = classifier.ReconcileAll();
        AssertFixedN6Anchors(report);
        var json = System.Text.Json.JsonSerializer.SerializeToNode(report)!.AsObject();
        int total = inventory.A2Degrees.E + inventory.A2Degrees.O;
        Assert.Equal(total, json["TotalLoci"]!.GetValue<int>());
        Assert.Equal(total, json["ConsumedLoci"]!.GetValue<int>());
        Assert.Equal(0, json["UnresolvedLoci"]!.GetValue<int>());
        Assert.Equal(total, json["ByAlgebraicSource"]!["ExactAlgebraic"]!.GetValue<int>());
        // Repeated local execution checks aggregation, not independent numerical truth. The controls below
        // remove/duplicate/mislabel these returned readings through the same reconciliation door.
        var readings = inventory.Loci.OrderBy(l => l.Id, StringComparer.Ordinal).Select(classifier.Classify).ToArray();
        var sources = json["ByCharacterSource"]!.AsObject();
        Assert.Equal(readings.Count(r => r.Source == A2CharacterSource.HermitianAxis), sources["HermitianAxis"]!.GetValue<int>());
        Assert.Equal(readings.Count(r => r.Source == A2CharacterSource.EpCharacter), sources["EpCharacterStable"]!.GetValue<int>());
        Assert.Equal(readings.Count(r => r.Source == A2CharacterSource.ExactRank), sources["ExactRankExecuted"]!.GetValue<int>());
        Assert.Equal(total, sources.Sum(kv => kv.Value!.GetValue<int>()));
        foreach (var group in readings.GroupBy(r => r.Kind))
            Assert.Equal(group.Count(), json["ByVerdict"]![group.Key.ToString()]!.GetValue<int>());
        Assert.Equal(total, json["ByVerdict"]!.AsObject().Sum(kv => kv.Value!.GetValue<int>()));
        Assert.Equal(readings.Where(r => r.Source == A2CharacterSource.EpCharacter)
            .Min(r => r.IsolationMargin!.Value), json["MinimumIsolationMargin"]!.GetValue<double>());
        Assert.DoesNotContain(sources.Where(kv => kv.Value!.GetValue<int>() > 0),
            kv => kv.Key.Contains("exact", StringComparison.OrdinalIgnoreCase));
        output.WriteLine(System.Text.Json.JsonSerializer.Serialize(report));
        foreach (var parity in new[] { A2Parity.Even, A2Parity.Odd })
        {
            var ids = inventory.Loci.Where(l => l.Parity == parity).Select(l => l.Id).ToHashSet();
            var subset = readings.Where(r => ids.Contains(r.LocusId)).ToArray();
            AssertFixedN6Anchors(new RouteBA2N6ReconciliationReport(subset.Length, subset.Length, 0,
                new Dictionary<string, int> { ["ExactAlgebraic"] = parity == A2Parity.Even ? inventory.A2Degrees.E : inventory.A2Degrees.O },
                new Dictionary<string, int>
                {
                    ["HermitianAxis"] = subset.Count(r => r.Source == A2CharacterSource.HermitianAxis),
                    ["EpCharacterStable"] = subset.Count(r => r.Source == A2CharacterSource.EpCharacter),
                    ["ExactRankExecuted"] = subset.Count(r => r.Source == A2CharacterSource.ExactRank)
                }, subset.GroupBy(r => r.Kind).ToDictionary(g => g.Key, g => g.Count()), double.PositiveInfinity), perParity: true);
            Assert.Equal(parity == A2Parity.Even ? inventory.A2Degrees.E : inventory.A2Degrees.O, subset.Length);
            output.WriteLine($"{parity}: loci={subset.Length}; HermitianAxis={subset.Count(r => r.Source == A2CharacterSource.HermitianAxis)}; " +
                $"EpCharacterStable={subset.Count(r => r.Source == A2CharacterSource.EpCharacter)}; " +
                $"Diabolic={subset.Count(r => r.Kind == EpCharacter.EpKind.Diabolic)}; Defective={subset.Count(r => r.Kind == EpCharacter.EpKind.Defective)}");
        }
    }

    [Theory]
    [InlineData("missing")]
    [InlineData("duplicate")]
    [InlineData("foreign-id")]
    [InlineData("unresolved")]
    [InlineData("wrong-algebraic")]
    [InlineData("wrong-geometric")]
    [InlineData("exact-label")]
    [InlineData("hermitian-label")]
    [InlineData("margin")]
    [InlineData("conjugation-partner")]
    [InlineData("parity-partner")]
    [Trait("Category", "ROUTE_B_A2_N6_RECONCILE")]
    public void Reconciliation_RejectsWrongConsumptionAndProvenance(string mutation)
    {
        var inventory = RouteBA2N6Inventory.LoadDefault();
        var classifier = new RouteBA2N6CharacterClassifier(inventory);
        // Synthetic readings exercise bookkeeping only; no claim of numerical character.
        var readings = inventory.Loci.OrderBy(l => l.Id, StringComparer.Ordinal).Select(l =>
            new A2CharacterReading(l.Id, EpCharacter.EpKind.Diabolic, 2, 2,
                RouteBA2N6CharacterClassifier.IsHermitianEligible(l) ? null : 0,
                0.5, 0.5, RouteBA2N6CharacterClassifier.IsHermitianEligible(l)
                    ? A2CharacterSource.HermitianAxis : A2CharacterSource.EpCharacter)
            { FullBlockHermiticityResidual = RouteBA2N6CharacterClassifier.IsHermitianEligible(l) ? 0 : null }).ToList();
        int index = readings.FindIndex(r => r.Source == A2CharacterSource.EpCharacter);
        switch (mutation)
        {
            case "missing": readings.RemoveAt(index); break;
            case "duplicate": readings.Add(readings[index]); break;
            case "foreign-id": readings[index] = readings[index] with { LocusId = "N6-UNKNOWN" }; break;
            case "unresolved": readings[index] = readings[index] with { Kind = EpCharacter.EpKind.NearEp }; break;
            case "wrong-algebraic": readings[index] = readings[index] with { Algebraic = 3 }; break;
            case "wrong-geometric": readings[index] = readings[index] with { Geometric = 1 }; break;
            case "exact-label": readings[index] = readings[index] with { Source = A2CharacterSource.ExactRank }; break;
            case "hermitian-label": readings[index] = readings[index] with { Source = A2CharacterSource.HermitianAxis }; break;
            case "margin": readings[index] = readings[index] with { IsolationMargin = double.NaN }; break;
            case "conjugation-partner":
            case "parity-partner":
                var loci = inventory.Loci.ToArray();
                loci[0] = mutation == "conjugation-partner"
                    ? loci[0] with { ConjugationPartnerId = loci[1].Id }
                    : loci[0] with { ParityPartnerId = loci[0].Id };
                // Bypass the loader only to prove ReconcileAll rechecks the loaded maps.
                var ctor = typeof(RouteBA2N6Inventory).GetConstructors(System.Reflection.BindingFlags.Instance
                    | System.Reflection.BindingFlags.NonPublic).Single();
                inventory = (RouteBA2N6Inventory)ctor.Invoke([inventory.SchemaVersion, inventory.N,
                    inventory.A2Degrees, inventory.LayerIdentity, loci, inventory.ExactRankCertificates]);
                classifier = new RouteBA2N6CharacterClassifier(inventory);
                break;
        }
        var error = Assert.Throws<InvalidOperationException>(() => classifier.ReconcileReadings(readings));
        Assert.Contains("N6 reconciliation", error.Message);
    }

    [Theory]
    [InlineData(EpCharacter.EpKind.Diabolic, 2, 0.0)]
    [InlineData(EpCharacter.EpKind.Defective, 1, 0.2)]
    [Trait("Category", "ROUTE_B_A2_N6_CHARACTER")]
    public void Kernel_RequiresTheFixedThreeRadii(EpCharacter.EpKind kind, int geo, double departure)
    {
        var radii = new List<double>();
        var readings = RouteBA2CharacterKernel.Read("control", [0, 0.001, 1, 3], 0, radius =>
        {
            radii.Add(radius);
            return Character(kind, geo, departure);
        });
        Assert.Equal(new[] { 0.25075, 0.5005, 0.75025 }, radii);
        Assert.Equal(3, readings.Count);
        Assert.All(readings, r =>
        {
            Assert.Equal(kind, r.Kind);
            Assert.Equal(2, r.Algebraic);
            Assert.Equal(geo, r.Geometric);
            Assert.True(r.IsolationMargin > 0);
        });
    }

    [Theory]
    [InlineData("Normal")]
    [InlineData("NearEp")]
    [InlineData("algebraic")]
    [InlineData("geometric")]
    [InlineData("departure-nan")]
    [InlineData("norm-infinity")]
    [InlineData("insufficient-margin")]
    [InlineData("AT-overlap")]
    [InlineData("spectrum-nan")]
    [InlineData("radius-disagreement")]
    [InlineData("diabolic-boundary")]
    [InlineData("defective-boundary")]
    [Trait("Category", "ROUTE_B_A2_N6_CHARACTER")]
    public void Kernel_RejectsWrongInputsThroughTheProductionDoor(string mutation)
    {
        Complex[] spectrum = [0, 0.001, 1, 3];
        var character = Character(EpCharacter.EpKind.Diabolic, 2, 0);
        switch (mutation)
        {
            case "Normal": character = character with { Kind = EpCharacter.EpKind.Normal }; break;
            case "NearEp": character = character with { Kind = EpCharacter.EpKind.NearEp }; break;
            case "algebraic": character = character with { Algebraic = 3 }; break;
            case "geometric": character = character with { Geometric = 1 }; break;
            case "departure-nan": character = character with { Departure = double.NaN }; break;
            case "norm-infinity": character = character with { CompressionNorm = double.PositiveInfinity }; break;
            case "insufficient-margin": spectrum = [0, 0.01, 1, 3]; break;
            case "AT-overlap": spectrum = [0, 0, 0, 3]; break;
            case "spectrum-nan": spectrum = [0, 0, double.NaN, 3]; break;
            case "diabolic-boundary": character = character with { Departure = 1e-6 }; break;
            case "defective-boundary": character = Character(EpCharacter.EpKind.Defective, 1, 5e-2); break;
        }
        int calls = 0;
        var error = Assert.Throws<A2CharacterUncertifiedException>(() =>
            RouteBA2CharacterKernel.Read("control", spectrum, 0, _ =>
                mutation == "radius-disagreement" && ++calls == 2
                    ? Character(EpCharacter.EpKind.Defective, 1, 0.2) : character));
        Assert.Equal("control", error.LocusId);
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_CHARACTER")]
    public void Kernel_UsesTheSuppliedLambdaRatherThanTheClosestPairAnywhere()
    {
        Assert.Throws<A2CharacterUncertifiedException>(() => RouteBA2CharacterKernel.Read(
            "wrong-seed", [0, 0, 5, 9], 5, _ => Character(EpCharacter.EpKind.Diabolic, 2, 0)));
    }

    [Theory]
    [InlineData(A2Parity.Even)]
    [InlineData(A2Parity.Odd)]
    [Trait("Category", "ROUTE_B_A2_N6_CHARACTER")]
    public void ExactRealT_ExecutesFull45DimensionalHermiticity(A2Parity parity)
    {
        var inventory = RouteBA2N6Inventory.LoadDefault();
        var locus = inventory.Loci.First(l => l.Parity == parity
            && l.TBox.ImLo.Numerator.IsZero && l.TBox.ImHi.Numerator.IsZero);
        var classifier = new RouteBA2N6CharacterClassifier(inventory);
        var readings = classifier.ClassifyRadii(locus);
        Assert.Equal(3, readings.Count);
        Assert.All(readings, r =>
        {
            Assert.Equal(A2CharacterSource.HermitianAxis, r.Source);
            Assert.Equal(EpCharacter.EpKind.Diabolic, r.Kind);
            Assert.Equal(2, r.Algebraic);
            Assert.Equal(2, r.Geometric);
            Assert.Equal(0.0, r.FullBlockHermiticityResidual);
            Assert.True(r.IsolationMargin > 0);
        });
        double t = double.Parse(locus.TSeed.Real, CultureInfo.InvariantCulture);
        var realQ = RouteBA2CharacterClassifier.FullParityBlock(6, parity, new Complex(t, 0));
        Assert.Equal(45, realQ.RowCount);
        Assert.True(RouteBA2CharacterClassifier.HermiticityResidual(realQ) > 1e-2);
        Assert.False(RouteBA2N6CharacterClassifier.IsHermitianEligible(locus with
        { TBox = new ExactComplexBox(locus.TBox.ReLo, locus.TBox.ReHi,
            locus.TBox.ImLo, new ExactRational(1, BigInteger.Pow(10, 400))) }));
    }

    [Theory]
    [InlineData(EpCharacter.EpKind.Normal)]
    [InlineData(EpCharacter.EpKind.NearEp)]
    [Trait("Category", "ROUTE_B_A2_N6_CHARACTER")]
    public void N6Dispatcher_RejectsForcedAmbiguity(EpCharacter.EpKind kind)
    {
        var inventory = RouteBA2N6Inventory.LoadDefault();
        var locus = inventory.Loci.First(l => !RouteBA2N6CharacterClassifier.IsHermitianEligible(l));
        var classifier = new RouteBA2N6CharacterClassifier(inventory, (_, _, _) => Character(kind, 1, 0));
        var error = Assert.Throws<A2CharacterUncertifiedException>(() => { classifier.Classify(locus); });
        Assert.Equal(locus.Id, error.LocusId);
        Assert.Equal(3, error.Readings.Count);
        Assert.All(error.Readings, r => Assert.Equal(kind, r.Kind));
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "ROUTE_B_A2_N6_CHARACTER")]
    public void ManifestOutput_ResolvesFromRepositoryAndCreatesMissingDirectory(bool absolute)
    {
        // Git is an independent root anchor for the expected destination; the writer uses
        // the source-project layout. The test executes under VSTest's output-directory cwd.
        var root = new DirectoryInfo(AppContext.BaseDirectory);
        while (root != null && !Directory.Exists(Path.Combine(root.FullName, ".git"))
            && !File.Exists(Path.Combine(root.FullName, ".git"))) root = root.Parent;
        Assert.NotNull(root);
        string relative = Path.Combine("simulations", "results", $"manifest-path-{Guid.NewGuid():N}", "ids.txt");
        string expected = absolute
            ? Path.Combine(Path.GetTempPath(), $"manifest-path-{Guid.NewGuid():N}", "ids.txt")
            : Path.Combine(root.FullName, relative);
        string directory = Path.GetDirectoryName(expected)!;
        Assert.False(Directory.Exists(directory));
        try
        {
            string written = WriteAmbiguityManifest(absolute ? expected : relative, ["N6-E-CONTROL"]);
            Assert.Equal(expected, written);
            Assert.Equal(new[] { "N6-E-CONTROL" }, File.ReadAllLines(expected));
            Assert.True(Directory.Exists(directory));
        }
        finally
        {
            if (File.Exists(expected)) File.Delete(expected);
            if (Directory.Exists(directory)) Directory.Delete(directory);
        }
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_AMBIGUITY")]
    public void RealManifest_AccountsForEveryLocusBeforeAnyFallback()
    {
        var inventory = RouteBA2N6Inventory.LoadDefault();
        var classifier = new RouteBA2N6CharacterClassifier(inventory);
        var accepted = new List<A2CharacterReading>();
        var ambiguous = new List<A2CharacterUncertifiedException>();
        foreach (var locus in inventory.Loci.OrderBy(l => l.Id, StringComparer.Ordinal))
        {
            try { accepted.Add(classifier.Classify(locus)); }
            catch (A2CharacterUncertifiedException error) { ambiguous.Add(error); }
        }
        string[] ids = ambiguous.Select(e => e.LocusId).ToArray();
        string? path = Environment.GetEnvironmentVariable("ROUTE_B_A2_N6_AMBIGUOUS_OUT");
        if (!string.IsNullOrWhiteSpace(path)) WriteAmbiguityManifest(path, ids);
        Assert.Equal(inventory.Loci.Select(l => l.Id).Order(StringComparer.Ordinal),
            accepted.Select(r => r.LocusId).Concat(ids).Order(StringComparer.Ordinal));
        Assert.Equal(ids.Distinct().Count(), ids.Length);
        foreach (var error in ambiguous) output.WriteLine(error.Message);
        output.WriteLine($"Loci={inventory.Loci.Count}; accepted={accepted.Count}; ambiguous={ids.Length}; " +
            $"HermitianAxis={accepted.Count(r => r.Source == A2CharacterSource.HermitianAxis)}; " +
            $"EpCharacter={accepted.Count(r => r.Source == A2CharacterSource.EpCharacter)}; " +
            $"Diabolic={accepted.Count(r => r.Kind == EpCharacter.EpKind.Diabolic)}; " +
            $"Defective={accepted.Count(r => r.Kind == EpCharacter.EpKind.Defective)}");
        Assert.Empty(ambiguous);
        Assert.Empty(inventory.ExactRankCertificates);
    }

    private static string WriteAmbiguityManifest(string path, IEnumerable<string> ids)
    {
        if (!Path.IsPathFullyQualified(path))
        {
            string? repository = null;
            foreach (string start in new[] { AppContext.BaseDirectory, Directory.GetCurrentDirectory() })
            {
                for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
                {
                    // The artifact is copied into bin; the source project identifies the repo root.
                    string project = Path.Combine(directory.FullName, "compute", "RCPsiSquared.Diagnostics.Tests",
                        "RCPsiSquared.Diagnostics.Tests.csproj");
                    if (!File.Exists(project)) continue;
                    repository = directory.FullName;
                    break;
                }
                if (repository != null) break;
            }
            path = Path.GetFullPath(path, repository
                ?? throw new DirectoryNotFoundException("Cannot locate the repository for the N=6 ambiguity manifest."));
        }
        Directory.CreateDirectory(Path.GetDirectoryName(path)!);
        File.WriteAllLines(path, ids);
        return path;
    }

    private static EpCharacter.Reading Character(EpCharacter.EpKind kind, int geometric, double departure) =>
        new(kind, 1, 2, 2, geometric, departure, 1, [], 0);

    private const string ExpectedSourcePencilDigest =
        "cf1549c54a116132373e481d0ce7a7ea03409c07f75e3f6ab5b9fd0f738dc916";

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void LoadDefault_ExposesTheExactN6DirectTInventory()
    {
        RouteBA2N6Inventory inventory = RouteBA2N6Inventory.LoadDefault();

        Assert.Equal(3, inventory.SchemaVersion);
        Assert.Equal(6, inventory.N);
        Assert.Equal(new N6A2Degrees(133, 133), inventory.A2Degrees);
        Assert.Equal(266, inventory.Loci.Count);
        Assert.Equal(inventory.A2Degrees.E + inventory.A2Degrees.O, inventory.Loci.Count);
        Assert.Equal(133, inventory.Loci.Count(x => x.Parity == A2Parity.Even));
        Assert.Equal(133, inventory.Loci.Count(x => x.Parity == A2Parity.Odd));
        Assert.All(inventory.Loci, x => Assert.Equal(2, x.AlgebraicMultiplicity));
        Assert.All(inventory.Loci,
            locus => Assert.Equal(locus.LambdaClearedBox.Half(), locus.LambdaPhysicalBox));
        Assert.All(inventory.Loci,
            locus => Assert.Equal(locus.TBox.RotateTToPhysicalQ(), locus.QPhysicalCSharpBox));
        Assert.Empty(inventory.ExactRankCertificates);
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void LoadDefault_ExposesStrictIdentityAndTotalPartnerInvolutions()
    {
        var inventory = RouteBA2N6Inventory.LoadDefault();
        Assert.Equal(926, inventory.LayerIdentity.DiscriminantDegree);
        Assert.Equal(536, inventory.LayerIdentity.Valuation);
        Assert.Equal(124, inventory.LayerIdentity.A1Degree);
        Assert.Equal(133, inventory.LayerIdentity.A2Degree);
        Assert.True(BigInteger.Parse(inventory.LayerIdentity.Constant) > 0);
        Assert.True(BigInteger.Parse(inventory.LayerIdentity.ProofModulus)
            > 2 * BigInteger.Parse(inventory.LayerIdentity.ProofBound));

        var byId = inventory.Loci.ToDictionary(x => x.Id, StringComparer.Ordinal);
        foreach (var locus in inventory.Loci)
        {
            var conjugate = byId[locus.ConjugationPartnerId];
            var parity = byId[locus.ParityPartnerId];
            Assert.Equal(locus.Id, conjugate.ConjugationPartnerId);
            Assert.Equal(locus.Id, parity.ParityPartnerId);
            Assert.Equal(locus.Parity, conjugate.Parity);
            Assert.NotEqual(locus.Parity, parity.Parity);
        }
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_SCHEMA")]
    public void ExistingN5Loader_RemainsSchemaOneAndUnchanged()
    {
        var inventory = RouteBA2Inventory.LoadDefault();
        Assert.Equal(1, inventory.SchemaVersion);
        Assert.Equal(5, inventory.N);
        Assert.Equal(29, inventory.Roots.Count);
        Assert.Equal(58, inventory.Roots.Sum(root => root.QLoci.Count));
    }

    [Theory]
    [InlineData("schema-2")]
    [InlineData("schema-4")]
    [InlineData("missing-source-digest")]
    [InlineData("wrong-source-digest")]
    [InlineData("uppercase-source-digest")]
    [InlineData("n")]
    [InlineData("model")]
    [InlineData("parameter-convention")]
    [InlineData("eigenvalue-convention")]
    [InlineData("coefficient-convention")]
    [InlineData("degree-e")]
    [InlineData("degree-o")]
    [InlineData("layer-degree")]
    [InlineData("layer-valuation")]
    [InlineData("layer-constant")]
    [InlineData("layer-proof-bound")]
    [InlineData("missing-locus")]
    [InlineData("duplicate-id")]
    [InlineData("unstable-id")]
    [InlineData("parity")]
    [InlineData("multiplicity")]
    [InlineData("zero-denominator")]
    [InlineData("negative-denominator")]
    [InlineData("unreduced-rational")]
    [InlineData("inverted-box")]
    [InlineData("missing-conjugation-partner")]
    [InlineData("conjugation-not-involution")]
    [InlineData("missing-parity-partner")]
    [InlineData("parity-not-involution")]
    [InlineData("wrong-q-rotation")]
    [InlineData("wrong-lambda-half")]
    [InlineData("nonempty-rank-certificates")]
    [InlineData("n5-sectors-lift")]
    [InlineData("unknown-top-field")]
    [InlineData("unknown-convention-field")]
    [InlineData("unknown-locus-field")]
    [InlineData("unknown-box-field")]
    [InlineData("unknown-rational-field")]
    [InlineData("unknown-seed-field")]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsStructuralCorruption(string mutation)
    {
        JsonObject document = LoadDocumentNode();
        JsonArray loci = document["loci"]!.AsArray();
        JsonObject first = loci[0]!.AsObject();
        JsonObject second = loci[1]!.AsObject();
        JsonObject nonreal = loci.Select(x => x!.AsObject()).First(x => !IsReal(x));

        switch (mutation)
        {
            case "schema-2": document["schemaVersion"] = 2; break;
            case "schema-4": document["schemaVersion"] = 4; break;
            case "missing-source-digest": document.Remove("sourcePencilDigest"); break;
            case "wrong-source-digest": document["sourcePencilDigest"] = new string('0', 64); break;
            case "uppercase-source-digest": document["sourcePencilDigest"] = ExpectedSourcePencilDigest.ToUpperInvariant(); break;
            case "n": document["n"] = 5; break;
            case "model": document["model"] = "closed-chain"; break;
            case "parameter-convention": document["conventions"]!["parameter"] = "q=t"; break;
            case "eigenvalue-convention": document["conventions"]!["eigenvalue"] = "Lambda=lambda"; break;
            case "coefficient-convention": document["conventions"]!["coefficientOrder"] = "highest-first"; break;
            case "degree-e": document["a2Degrees"]!["E"] = 132; break;
            case "degree-o": document["a2Degrees"]!["O"] = 134; break;
            case "layer-degree": document["layerIdentity"]!["a2Degree"] = 132; break;
            case "layer-valuation": document["layerIdentity"]!["valuation"] = 535; break;
            case "layer-constant": document["layerIdentity"]!["constant"] = "0"; break;
            case "layer-proof-bound":
                document["layerIdentity"]!["proofBound"] = document["layerIdentity"]!["proofModulus"]!.DeepClone();
                break;
            case "missing-locus": loci.RemoveAt(0); break;
            case "duplicate-id": second["id"] = first["id"]!.GetValue<string>(); break;
            case "unstable-id": first["id"] = "N6-E-A2-T-999"; break;
            case "parity": first["parity"] = "Even"; break;
            case "multiplicity": first["algebraicMultiplicity"] = 1; break;
            case "zero-denominator": Rational(first, "tBox", "real", "lower")["denominator"] = "0"; break;
            case "negative-denominator": Rational(first, "tBox", "real", "lower")["denominator"] = "-1"; break;
            case "unreduced-rational":
                Rational(first, "tBox", "real", "lower")["numerator"] = "2";
                Rational(first, "tBox", "real", "lower")["denominator"] = "4";
                break;
            case "inverted-box":
                first["tBox"]!["real"]!["lower"] = first["tBox"]!["real"]!["upper"]!.DeepClone();
                Rational(first, "tBox", "real", "lower")["numerator"] = "100";
                Rational(first, "tBox", "real", "lower")["denominator"] = "1";
                break;
            case "missing-conjugation-partner": nonreal["conjugationPartnerId"] = "missing"; break;
            case "conjugation-not-involution": nonreal["conjugationPartnerId"] = first["id"]!.GetValue<string>(); break;
            case "missing-parity-partner": first["parityPartnerId"] = "missing"; break;
            case "parity-not-involution": first["parityPartnerId"] = "N6-O-A2-T-001"; break;
            case "wrong-q-rotation": first["qPhysicalCSharpBox"] = first["tBox"]!.DeepClone(); break;
            case "wrong-lambda-half": first["lambdaPhysicalBox"] = first["lambdaClearedBox"]!.DeepClone(); break;
            case "nonempty-rank-certificates": document["exactRankCertificates"] = new JsonArray(new JsonObject()); break;
            case "n5-sectors-lift": document["sectors"] = new JsonArray(); break;
            case "unknown-top-field": document["unchecked"] = true; break;
            case "unknown-convention-field": document["conventions"]!["unchecked"] = true; break;
            case "unknown-locus-field": first["unchecked"] = true; break;
            case "unknown-box-field": first["tBox"]!["unchecked"] = true; break;
            case "unknown-rational-field": Rational(first, "tBox", "real", "lower")["unchecked"] = true; break;
            case "unknown-seed-field": first["tSeed"]!["unchecked"] = true; break;
            default: throw new InvalidOperationException(mutation);
        }

        Assert.Throws<InvalidDataException>(() => LoadMutated(document));
    }

    [Theory]
    [InlineData("")]
    [InlineData("NaN")]
    [InlineData("Infinity")]
    [InlineData("-Infinity")]
    [InlineData("1,25")]
    [InlineData("1.2.5")]
    [InlineData(" 1")]
    [InlineData("1 ")]
    [InlineData("0x10")]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsNonDecimalSeedText(string text)
    {
        JsonObject document = LoadDocumentNode();
        document["loci"]![0]!["tSeed"]!["real"] = text;
        Assert.Throws<InvalidDataException>(() => LoadMutated(document));
    }

    [Theory]
    [InlineData("", "schemaVersion")]
    [InlineData("", "sourcePencilDigest")]
    [InlineData("\"loci\": [", "id")]
    [InlineData("\"tBox\": {", "real")]
    [InlineData("\"tBox\": {", "numerator")]
    [InlineData("\"tSeed\": {", "real")]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsDuplicateRawJsonPropertiesAtEveryDtoDepth(string anchor, string property)
    {
        string duplicate = DuplicateFirstProperty(File.ReadAllText(ArtifactPath()), anchor, property);
        var error = Assert.Throws<InvalidDataException>(() => LoadRaw(duplicate));
        Assert.Contains("Duplicate JSON property", error.Message);
    }

    [Theory]
    [InlineData("positive-wrong-constant")]
    [InlineData("different-valid-modulus-bound")]
    [InlineData("single-digit-corruption")]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_BindsTheCompleteCommittedLayerIdentity(string mutation)
    {
        JsonObject document = LoadDocumentNode();
        JsonObject identity = document["layerIdentity"]!.AsObject();
        switch (mutation)
        {
            case "positive-wrong-constant":
                identity["constant"] = "1";
                break;
            case "different-valid-modulus-bound":
                identity["proofModulus"] = "101";
                identity["proofBound"] = "50";
                break;
            case "single-digit-corruption":
                string value = identity["constant"]!.GetValue<string>();
                identity["constant"] = value[..^1] + (value[^1] == '9' ? '8' : (char)(value[^1] + 1));
                break;
            default:
                throw new InvalidOperationException(mutation);
        }

        var error = Assert.Throws<InvalidDataException>(() => LoadMutated(document));
        Assert.Contains("canonical N=6 LayerIdentity", error.Message);
    }

    [Theory]
    [InlineData("t-outside")]
    [InlineData("q-outside")]
    [InlineData("lambda-cleared-outside")]
    [InlineData("lambda-physical-outside")]
    [InlineData("q-seed-transform")]
    [InlineData("lambda-seed-transform")]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsSeedContainmentAndTransformCorruption(string mutation)
    {
        JsonObject document = LoadDocumentNode();
        JsonObject locus = document["loci"]!.AsArray().Select(x => x!.AsObject()).First(x => !IsReal(x));
        switch (mutation)
        {
            case "t-outside": locus["tSeed"]!["real"] = "1000000"; break;
            case "q-outside": locus["qPhysicalCSharpSeed"]!["real"] = "1000000"; break;
            case "lambda-cleared-outside": locus["lambdaClearedSeed"]!["real"] = "1000000"; break;
            case "lambda-physical-outside": locus["lambdaPhysicalSeed"]!["real"] = "1000000"; break;
            case "q-seed-transform": PerturbLastDigit(locus["qPhysicalCSharpSeed"]!, "real"); break;
            case "lambda-seed-transform": PerturbLastDigit(locus["lambdaPhysicalSeed"]!, "real"); break;
            default: throw new InvalidOperationException(mutation);
        }
        Assert.Throws<InvalidDataException>(() => LoadMutated(document));
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsACoherentParityPairWhoseTBoxContainsZero()
    {
        JsonObject document = LoadDocumentNode();
        JsonArray loci = document["loci"]!.AsArray();
        JsonObject even = loci.Select(x => x!.AsObject()).First(IsReal);
        JsonObject odd = Partner(loci, even, "parityPartnerId");
        SetTPoint(even, 0, 0);
        SetTPoint(odd, 0, 0);

        var error = Assert.Throws<InvalidDataException>(() => LoadMutated(document));
        Assert.Contains("must exclude t=0", error.Message);
    }

    [Theory]
    [InlineData(false)]
    [InlineData(true)]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsACoherentOrbitThatTouchesOrStraddlesTheRealAxis(bool straddles)
    {
        JsonObject document = LoadDocumentNode();
        JsonArray loci = document["loci"]!.AsArray();
        JsonObject upper = loci.Select(x => x!.AsObject()).First(IsStrictUpper);
        JsonObject lower = Partner(loci, upper, "conjugationPartnerId");
        JsonObject oddOfUpper = Partner(loci, upper, "parityPartnerId");
        JsonObject oddOfLower = Partner(loci, lower, "parityPartnerId");
        JsonObject widened = upper["tBox"]!.DeepClone().AsObject();
        widened["imag"]!["lower"] = straddles
            ? NegateRational(widened["imag"]!["upper"]!)
            : RationalNode(0);
        SetTBox(upper, widened);
        SetTBox(lower, ConjugateBox(widened));
        SetTBox(oddOfUpper, NegateBox(widened));
        SetTBox(oddOfLower, NegateBox(ConjugateBox(widened)));

        var error = Assert.Throws<InvalidDataException>(() => LoadMutated(document));
        Assert.Contains("root-box form", error.Message);
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsACoherentInventoryWithTheWrongRealUpperLowerDistribution()
    {
        JsonObject document = LoadDocumentNode();
        JsonArray loci = document["loci"]!.AsArray();
        JsonObject upper = loci.Select(x => x!.AsObject()).First(IsStrictUpper);
        JsonObject lower = Partner(loci, upper, "conjugationPartnerId");
        JsonObject oddOfUpper = Partner(loci, upper, "parityPartnerId");
        JsonObject oddOfLower = Partner(loci, lower, "parityPartnerId");

        SetCompletePoint(upper, 1000, 0, 0, 0);
        SetCompletePoint(lower, 1001, 0, 0, 0);
        SetCompletePoint(oddOfUpper, -1000, 0, 0, 0);
        SetCompletePoint(oddOfLower, -1001, 0, 0, 0);
        upper["conjugationPartnerId"] = upper["id"]!.GetValue<string>();
        lower["conjugationPartnerId"] = lower["id"]!.GetValue<string>();
        oddOfUpper["conjugationPartnerId"] = oddOfUpper["id"]!.GetValue<string>();
        oddOfLower["conjugationPartnerId"] = oddOfLower["id"]!.GetValue<string>();

        var error = Assert.Throws<InvalidDataException>(() => LoadMutated(document));
        Assert.Contains("59 real, 37 upper and 37 lower", error.Message);
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsCompleteDisjointOrbitSwapThatOnlyBreaksExactTBoxOrdering()
    {
        JsonObject document = LoadDocumentNode();
        JsonArray loci = document["loci"]!.AsArray();
        JsonObject[] evenReal = loci.Select(x => x!.AsObject())
            .Where(x => x["parity"]!.GetValue<string>() == "E" && IsReal(x)).Take(2).ToArray();
        JsonObject firstOdd = Partner(loci, evenReal[0], "parityPartnerId");
        JsonObject secondOdd = Partner(loci, evenReal[1], "parityPartnerId");
        SwapLocalGeometry(evenReal[0], evenReal[1]);
        SwapLocalGeometry(firstOdd, secondOdd);

        var error = Assert.Throws<InvalidDataException>(() => LoadMutated(document));
        Assert.Contains("lexicographic exact t-box order", error.Message);
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsSameParityContactWhileExactOrderingAndPartnerMapsRemainValid()
    {
        JsonObject document = LoadDocumentNode();
        JsonArray loci = document["loci"]!.AsArray();
        JsonObject[] evenReal = loci.Select(x => x!.AsObject())
            .Where(x => x["parity"]!.GetValue<string>() == "E" && IsReal(x)).Take(2).ToArray();
        JsonObject left = evenReal[0], right = evenReal[1];
        JsonObject leftOdd = Partner(loci, left, "parityPartnerId");
        JsonObject widened = left["tBox"]!.DeepClone().AsObject();
        widened["real"]!["upper"] = right["tBox"]!["real"]!["lower"]!.DeepClone();
        SetTBox(left, widened);
        SetTBox(leftOdd, NegateBox(widened));

        Assert.True(CompareJsonBoxes(left["tBox"]!.AsObject(), right["tBox"]!.AsObject()) < 0);
        Assert.True(JsonNode.DeepEquals(RotateBox(left["tBox"]!.AsObject()),
            left["qPhysicalCSharpBox"]));
        Assert.True(JsonNode.DeepEquals(NegateBox(left["tBox"]!.AsObject()), leftOdd["tBox"]));
        Assert.Null(Record.Exception(() => AssertPartnerMapsWithoutCrossBoxChecks(document)));
        var error = Assert.Throws<InvalidDataException>(() => LoadMutated(document));
        Assert.Contains("Same-parity t boxes overlap", error.Message);
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsLocallyCoherentConjugationBoxMismatchAtTheCrossPartnerGate()
    {
        JsonObject document = LoadDocumentNode();
        JsonArray loci = document["loci"]!.AsArray();
        JsonObject upper = loci.Select(x => x!.AsObject()).First(IsStrictUpper);
        JsonObject oddPartner = Partner(loci, upper, "parityPartnerId");
        JsonObject expanded = ExpandImagUpperEndpoint(upper["tBox"]!.AsObject());
        SetTBox(upper, expanded);
        SetTBox(oddPartner, NegateBox(expanded));

        Assert.Null(Record.Exception(() => AssertPartnerMapsWithoutCrossBoxChecks(document)));
        var error = Assert.Throws<InvalidDataException>(() => LoadMutated(document));
        Assert.Contains("Conjugation partner boxes disagree", error.Message);
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsLocallyCoherentParityBoxMismatchAtTheCrossPartnerGate()
    {
        JsonObject document = LoadDocumentNode();
        JsonArray loci = document["loci"]!.AsArray();
        JsonObject evenReal = loci.Select(x => x!.AsObject())
            .First(x => x["parity"]!.GetValue<string>() == "E" && IsReal(x));
        SetTBox(evenReal, ExpandRealLowerEndpoint(evenReal["tBox"]!.AsObject()));

        Assert.Null(Record.Exception(() => AssertPartnerMapsWithoutCrossBoxChecks(document)));
        var error = Assert.Throws<InvalidDataException>(() => LoadMutated(document));
        Assert.Contains("Parity transport boxes disagree", error.Message);
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsLocallyCoherentConjugationLambdaBoxMismatchAtTheCrossPartnerGate()
    {
        JsonObject document = LoadDocumentNode();
        JsonArray loci = document["loci"]!.AsArray();
        JsonObject upper = loci.Select(x => x!.AsObject()).First(IsStrictUpper);
        JsonObject oddPartner = Partner(loci, upper, "parityPartnerId");
        JsonObject originalT = upper["tBox"]!.DeepClone().AsObject();
        JsonObject expandedLambda = ExpandRealLowerEndpoint(upper["lambdaClearedBox"]!.AsObject());
        SetLambdaBox(upper, expandedLambda);
        SetLambdaBox(oddPartner, expandedLambda);

        Assert.True(JsonNode.DeepEquals(originalT, upper["tBox"]));
        Assert.True(JsonNode.DeepEquals(HalfBox(expandedLambda), upper["lambdaPhysicalBox"]));
        Assert.True(JsonNode.DeepEquals(expandedLambda, oddPartner["lambdaClearedBox"]));
        Assert.Null(Record.Exception(() => AssertPartnerMapsWithoutCrossBoxChecks(document)));
        var error = Assert.Throws<InvalidDataException>(() => LoadMutated(document));
        Assert.Contains("Conjugation partner boxes disagree", error.Message);
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_RejectsLocallyCoherentParityLambdaBoxMismatchAtTheCrossPartnerGate()
    {
        JsonObject document = LoadDocumentNode();
        JsonArray loci = document["loci"]!.AsArray();
        JsonObject evenReal = loci.Select(x => x!.AsObject())
            .First(x => x["parity"]!.GetValue<string>() == "E" && IsReal(x));
        JsonObject originalT = evenReal["tBox"]!.DeepClone().AsObject();
        JsonObject expandedLambda = ExpandRealLowerEndpoint(evenReal["lambdaClearedBox"]!.AsObject());
        SetLambdaBox(evenReal, expandedLambda);

        Assert.True(JsonNode.DeepEquals(originalT, evenReal["tBox"]));
        Assert.True(JsonNode.DeepEquals(HalfBox(expandedLambda), evenReal["lambdaPhysicalBox"]));
        Assert.Null(Record.Exception(() => AssertPartnerMapsWithoutCrossBoxChecks(document)));
        var error = Assert.Throws<InvalidDataException>(() => LoadMutated(document));
        Assert.Contains("Parity transport boxes disagree", error.Message);
    }

    [Fact]
    [Trait("Category", "ROUTE_B_A2_N6_SCHEMA")]
    public void Load_AllowsIndependentContainedPartnerSeedApproximants()
    {
        var inventory = RouteBA2N6Inventory.LoadDefault();
        var byId = inventory.Loci.ToDictionary(x => x.Id, StringComparer.Ordinal);
        var locus = inventory.Loci.First(x => x.TSeed.Real !=
            NegateDecimalText(byId[x.ParityPartnerId].TSeed.Real));
        var partner = byId[locus.ParityPartnerId];

        Assert.NotEqual(locus.TSeed.Real, NegateDecimalText(partner.TSeed.Real));
        Assert.Equal(locus.TBox.Negate(), partner.TBox);
    }

    private static void CopyLocalGeometry(JsonObject source, JsonObject target)
    {
        foreach (string name in new[] { "tBox", "qPhysicalCSharpBox", "lambdaClearedBox", "lambdaPhysicalBox",
                     "tSeed", "qPhysicalCSharpSeed", "lambdaClearedSeed", "lambdaPhysicalSeed" })
            target[name] = source[name]!.DeepClone();
    }

    private static void SwapLocalGeometry(JsonObject left, JsonObject right)
    {
        foreach (string name in new[] { "tBox", "qPhysicalCSharpBox", "lambdaClearedBox", "lambdaPhysicalBox",
                     "tSeed", "qPhysicalCSharpSeed", "lambdaClearedSeed", "lambdaPhysicalSeed" })
        {
            JsonNode saved = left[name]!.DeepClone();
            left[name] = right[name]!.DeepClone();
            right[name] = saved;
        }
    }

    private static void SetCompletePoint(JsonObject locus, int tReal, int tImag, int lambdaReal, int lambdaImag)
    {
        SetTPoint(locus, tReal, tImag);
        locus["lambdaClearedBox"] = PointBox(lambdaReal, lambdaImag);
        locus["lambdaPhysicalBox"] = PointBox(lambdaReal / 2, lambdaImag / 2);
        locus["lambdaClearedSeed"] = SeedNode(lambdaReal, lambdaImag);
        locus["lambdaPhysicalSeed"] = SeedNode(lambdaReal / 2, lambdaImag / 2);
    }

    private static void SetTPoint(JsonObject locus, int real, int imag)
    {
        SetTBox(locus, PointBox(real, imag));
        locus["tSeed"] = SeedNode(real, imag);
        locus["qPhysicalCSharpSeed"] = SeedNode(imag, -real);
    }

    private static void SetTBox(JsonObject locus, JsonObject box)
    {
        locus["tBox"] = box.DeepClone();
        locus["qPhysicalCSharpBox"] = RotateBox(box);
    }

    private static void SetLambdaBox(JsonObject locus, JsonObject cleared)
    {
        locus["lambdaClearedBox"] = cleared.DeepClone();
        locus["lambdaPhysicalBox"] = HalfBox(cleared);
    }

    private static JsonObject ExpandRealLowerEndpoint(JsonObject box)
    {
        JsonObject expanded = box.DeepClone().AsObject();
        JsonObject lower = expanded["real"]!["lower"]!.AsObject();
        BigInteger numerator = BigInteger.Parse(lower["numerator"]!.GetValue<string>());
        lower["numerator"] = (numerator - 2).ToString();
        return expanded;
    }

    private static JsonObject ExpandImagUpperEndpoint(JsonObject box)
    {
        JsonObject expanded = box.DeepClone().AsObject();
        JsonObject upper = expanded["imag"]!["upper"]!.AsObject();
        BigInteger numerator = BigInteger.Parse(upper["numerator"]!.GetValue<string>());
        upper["numerator"] = (numerator + 2).ToString();
        return expanded;
    }

    private static JsonObject PointBox(int real, int imag) => new()
    {
        ["real"] = new JsonObject { ["lower"] = RationalNode(real), ["upper"] = RationalNode(real) },
        ["imag"] = new JsonObject { ["lower"] = RationalNode(imag), ["upper"] = RationalNode(imag) }
    };

    private static JsonObject SeedNode(int real, int imag) => new()
    {
        ["real"] = real.ToString(),
        ["imag"] = imag.ToString()
    };

    private static JsonObject RationalNode(int value) => new()
    {
        ["numerator"] = value.ToString(),
        ["denominator"] = "1"
    };

    private static JsonObject RotateBox(JsonObject box) => new()
    {
        ["real"] = box["imag"]!.DeepClone(),
        ["imag"] = new JsonObject
        {
            ["lower"] = NegateRational(box["real"]!["upper"]!),
            ["upper"] = NegateRational(box["real"]!["lower"]!)
        }
    };

    private static JsonObject HalfBox(JsonObject box) => new()
    {
        ["real"] = new JsonObject
        {
            ["lower"] = HalfRational(box["real"]!["lower"]!),
            ["upper"] = HalfRational(box["real"]!["upper"]!)
        },
        ["imag"] = new JsonObject
        {
            ["lower"] = HalfRational(box["imag"]!["lower"]!),
            ["upper"] = HalfRational(box["imag"]!["upper"]!)
        }
    };

    private static JsonObject HalfRational(JsonNode rational)
    {
        BigInteger numerator = BigInteger.Parse(rational["numerator"]!.GetValue<string>());
        BigInteger denominator = 2 * BigInteger.Parse(rational["denominator"]!.GetValue<string>());
        BigInteger gcd = BigInteger.GreatestCommonDivisor(BigInteger.Abs(numerator), denominator);
        return new JsonObject
        {
            ["numerator"] = (numerator / gcd).ToString(),
            ["denominator"] = (denominator / gcd).ToString()
        };
    }

    private static JsonObject ConjugateBox(JsonObject box) => new()
    {
        ["real"] = box["real"]!.DeepClone(),
        ["imag"] = new JsonObject
        {
            ["lower"] = NegateRational(box["imag"]!["upper"]!),
            ["upper"] = NegateRational(box["imag"]!["lower"]!)
        }
    };

    private static JsonObject NegateBox(JsonObject box) => new()
    {
        ["real"] = new JsonObject
        {
            ["lower"] = NegateRational(box["real"]!["upper"]!),
            ["upper"] = NegateRational(box["real"]!["lower"]!)
        },
        ["imag"] = new JsonObject
        {
            ["lower"] = NegateRational(box["imag"]!["upper"]!),
            ["upper"] = NegateRational(box["imag"]!["lower"]!)
        }
    };

    private static JsonObject NegateRational(JsonNode rational)
    {
        var result = rational.DeepClone().AsObject();
        BigInteger numerator = BigInteger.Parse(result["numerator"]!.GetValue<string>());
        result["numerator"] = (-numerator).ToString();
        return result;
    }

    private static bool IsStrictUpper(JsonObject locus) =>
        BigInteger.Parse(Rational(locus, "tBox", "imag", "lower")["numerator"]!.GetValue<string>()) > 0;

    private static void PerturbLastDigit(JsonNode seed, string component)
    {
        string value = seed[component]!.GetValue<string>();
        char last = value[^1];
        seed[component] = value[..^1] + (last == '9' ? '8' : (char)(last + 1));
    }

    private static string NegateDecimalText(string value) =>
        value.StartsWith('-') ? value[1..] : "-" + value;

    private static string DuplicateFirstProperty(string json, string anchor, string property)
    {
        int searchStart = string.IsNullOrEmpty(anchor)
            ? 0 : json.IndexOf(anchor, StringComparison.Ordinal);
        Assert.True(searchStart >= 0, $"Anchor not found: {anchor}");
        int propertyStart = json.IndexOf($"\"{property}\"", searchStart, StringComparison.Ordinal);
        Assert.True(propertyStart >= 0, $"Property not found after anchor: {property}");
        int colon = json.IndexOf(':', propertyStart);
        int valueStart = colon + 1;
        while (char.IsWhiteSpace(json[valueStart])) valueStart++;
        int valueEnd = JsonValueEnd(json, valueStart);
        string propertyText = json[propertyStart..valueEnd];
        return json.Insert(propertyStart, propertyText + ", ");
    }

    private static int JsonValueEnd(string json, int start)
    {
        if (json[start] == '"')
        {
            bool escaped = false;
            for (int index = start + 1; index < json.Length; index++)
            {
                if (!escaped && json[index] == '"') return index + 1;
                escaped = !escaped && json[index] == '\\';
                if (json[index] != '\\') escaped = false;
            }
        }
        if (json[start] is '{' or '[')
        {
            char open = json[start], close = open == '{' ? '}' : ']';
            int depth = 0;
            bool inString = false, escaped = false;
            for (int index = start; index < json.Length; index++)
            {
                char current = json[index];
                if (inString)
                {
                    if (!escaped && current == '"') inString = false;
                    escaped = !escaped && current == '\\';
                    if (current != '\\') escaped = false;
                    continue;
                }
                if (current == '"') { inString = true; continue; }
                if (current == open) depth++;
                if (current == close && --depth == 0) return index + 1;
            }
        }
        int delimiter = json.IndexOfAny([',', '}', ']'], start);
        return delimiter < 0 ? json.Length : delimiter;
    }

    private static JsonObject Partner(JsonArray loci, JsonObject locus, string property)
    {
        string id = locus[property]!.GetValue<string>();
        return loci.Select(x => x!.AsObject()).Single(x => x["id"]!.GetValue<string>() == id);
    }

    private static void AssertPartnerMapsWithoutCrossBoxChecks(JsonObject document)
    {
        JsonObject[] loci = document["loci"]!.AsArray().Select(x => x!.AsObject()).ToArray();
        var byId = loci.ToDictionary(x => x["id"]!.GetValue<string>(), StringComparer.Ordinal);
        foreach (JsonObject locus in loci)
        {
            string id = locus["id"]!.GetValue<string>();
            string parity = locus["parity"]!.GetValue<string>();
            JsonObject conjugate = byId[locus["conjugationPartnerId"]!.GetValue<string>()];
            JsonObject parityPartner = byId[locus["parityPartnerId"]!.GetValue<string>()];
            Assert.Equal(id, conjugate["conjugationPartnerId"]!.GetValue<string>());
            Assert.Equal(parity, conjugate["parity"]!.GetValue<string>());
            Assert.Equal(id, parityPartner["parityPartnerId"]!.GetValue<string>());
            Assert.NotEqual(parity, parityPartner["parity"]!.GetValue<string>());
        }
    }

    private static int CompareJsonBoxes(JsonObject left, JsonObject right)
    {
        foreach ((string component, string endpoint) in new[]
                 { ("real", "lower"), ("real", "upper"), ("imag", "lower"), ("imag", "upper") })
        {
            int comparison = CompareRationals(left[component]![endpoint]!, right[component]![endpoint]!);
            if (comparison != 0) return comparison;
        }
        return 0;
    }

    private static int CompareRationals(JsonNode left, JsonNode right)
    {
        BigInteger leftNumerator = BigInteger.Parse(left["numerator"]!.GetValue<string>());
        BigInteger leftDenominator = BigInteger.Parse(left["denominator"]!.GetValue<string>());
        BigInteger rightNumerator = BigInteger.Parse(right["numerator"]!.GetValue<string>());
        BigInteger rightDenominator = BigInteger.Parse(right["denominator"]!.GetValue<string>());
        return (leftNumerator * rightDenominator).CompareTo(rightNumerator * leftDenominator);
    }

    private static bool IsReal(JsonObject locus) =>
        Rational(locus, "tBox", "imag", "lower")["numerator"]!.GetValue<string>() == "0"
        && Rational(locus, "tBox", "imag", "upper")["numerator"]!.GetValue<string>() == "0";

    private static JsonObject Rational(JsonObject locus, string box, string component, string endpoint) =>
        locus[box]![component]![endpoint]!.AsObject();

    private static JsonObject LoadDocumentNode() => JsonNode.Parse(File.ReadAllText(ArtifactPath()))!.AsObject();

    private static RouteBA2N6Inventory LoadMutated(JsonNode document)
    {
        string path = Path.Combine(Path.GetTempPath(), $"route-b-a2-n6-{Guid.NewGuid():N}.json");
        try
        {
            File.WriteAllText(path, document.ToJsonString());
            return RouteBA2N6Inventory.Load(path);
        }
        finally { File.Delete(path); }
    }

    private static RouteBA2N6Inventory LoadRaw(string json)
    {
        string path = Path.Combine(Path.GetTempPath(), $"route-b-a2-n6-raw-{Guid.NewGuid():N}.json");
        try
        {
            File.WriteAllText(path, json);
            return RouteBA2N6Inventory.Load(path);
        }
        finally { File.Delete(path); }
    }

    private static string ArtifactPath()
    {
        foreach (string start in new[] { AppContext.BaseDirectory, Directory.GetCurrentDirectory() })
        for (var directory = new DirectoryInfo(start); directory != null; directory = directory.Parent)
        {
            string path = Path.Combine(directory.FullName, "simulations", "results", "route_b_a2_n6.json");
            if (File.Exists(path)) return path;
        }
        throw new FileNotFoundException("The N=6 Route B A2 inventory artifact was not found.");
    }
}
