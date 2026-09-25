using System.Text.RegularExpressions;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Knowledge;
using RCPsiSquared.Runtime.ObjectManager;
using RCPsiSquared.Runtime.PolarityArchitecture;

namespace RCPsiSquared.Runtime.Tests.PolarityArchitecture;

public class AbsorptionTheoremClaimRegistrationTests
{
    private static ClaimRegistryBuilder BuildBaseRegistry() =>
        new ClaimRegistryBuilder()
            .RegisterPi2Family()
            .RegisterPi2DyadicLadder();

    [Fact]
    public void RegisterAbsorptionTheoremClaim_AddsClaim()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .Build();

        Assert.True(registry.Contains<AbsorptionTheoremClaim>());
    }

    [Fact]
    public void RegisterAbsorptionTheoremClaim_TierIsTier1Derived()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .Build();

        Assert.Equal(Tier.Tier1Derived,
            registry.Get<AbsorptionTheoremClaim>().Tier);
    }

    [Fact]
    public void RegisterAbsorptionTheoremClaim_AncestorsContainPi2DyadicLadder()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .Build();

        var ancestors = registry.AncestorsOf<AbsorptionTheoremClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(Pi2DyadicLadderClaim), ancestors);
    }

    [Fact]
    public void RegisterAbsorptionTheoremClaim_DissipatorCoefficientIsTwo()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .Build();

        // Exact route: Term(0) = 2^(1−0) is exactly 2.0.
        Assert.True(registry.Get<AbsorptionTheoremClaim>().DissipatorCoefficient == 2.0);
    }

    [Theory]
    [InlineData(0.0625, 0.125)]   // dyadic inputs: 2γ₀ is exact, so the comparison is exact
    [InlineData(1.0, 2.0)]
    public void RegisterAbsorptionTheoremClaim_SingleDisagreementCellCostAcrossRegistry(double gammaZero, double expected)
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .Build();

        Assert.Equal(expected, registry.Get<AbsorptionTheoremClaim>().SingleDisagreementCellCost(gammaZero));
    }

    [Fact]
    public void RegisterAbsorptionTheoremClaim_F66BecomesDescendantWhenWired()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        var descendants = registry.DescendantsOf<AbsorptionTheoremClaim>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(F66PoleModesPi2Inheritance), descendants);
    }

    /// <summary>The self-naming rule, read off the child alone: the F-number of
    /// <c>docs/ANALYTICAL_FORMULAS.md</c> a direct child names itself by, or null. A child counts when its
    /// class name opens with the number (F25CPsiBellPlusPi2Inheritance), when its Name opens with it
    /// ("F153 pinning criterion: ..."), or when the title of its Name, the text before the first colon,
    /// closes with it in parentheses ("The record parity law (F135): ..."). A number cited anywhere else
    /// does not count. Where more than one form applies they must agree.</summary>
    private static string? SelfNamedFNumber(Type type, string name)
    {
        var found = new List<string>();
        foreach (var match in new[] { ClassNamePrefix.Match(type.Name), NamePrefix.Match(name), TitleSuffix.Match(name) })
            if (match.Success) found.Add(match.Groups["f"].Value);
        var distinct = found.Distinct(StringComparer.Ordinal).ToList();
        Assert.True(distinct.Count <= 1, $"{type.Name} names itself by more than one number: {string.Join(", ", distinct)}");
        return distinct.Count == 1 ? distinct[0] : null;
    }

    private static readonly Regex ClassNamePrefix = new(@"^(?<f>F\d+[a-z]?)(?=[A-Z])", RegexOptions.CultureInvariant);
    private static readonly Regex NamePrefix = new(@"^(?<f>F\d+[a-z]?)(?![0-9A-Za-z])", RegexOptions.CultureInvariant);
    private static readonly Regex TitleSuffix = new(@"^[^:]*\((?<f>F\d+[a-z]?)\)\s*:", RegexOptions.CultureInvariant);

    /// <summary>The claim's F-numbered list, pinned by equality against the default registry: the
    /// self-naming rule runs on every direct child the registry returns, and the numbers it finds must be
    /// exactly <see cref="AbsorptionTheoremClaim.FNumberedDirectChildren"/>, each with its entry heading
    /// in <c>docs/ANALYTICAL_FORMULAS.md</c>. The two counts the claim's doc comment states, thirty-four
    /// direct children and fourteen that name themselves by a number, are read off the same
    /// registry.</summary>
    [Fact]
    public void DefaultRegistry_SelfNamedDirectChildren_EqualTheClaimsList()
    {
        var registry = KnowledgeRegistryFactory.BuildDefault();
        var direct = registry.EdgesFrom<AbsorptionTheoremClaim>().Select(e => e.Child).Distinct().ToList();

        var labelled = new List<(string Child, string Number)>();
        foreach (var type in direct)
        {
            var number = SelfNamedFNumber(type, registry.Get(type).Name);
            if (number is not null) labelled.Add((type.Name, number));
        }

        var found = labelled.Select(x => x.Number).Distinct().OrderBy(s => s, StringComparer.Ordinal).ToList();
        var claimed = AbsorptionTheoremClaim.FNumberedDirectChildren.OrderBy(s => s, StringComparer.Ordinal).ToList();
        Assert.True(found.SequenceEqual(claimed),
            "the self-named direct children changed: the registry gives " +
            $"[{string.Join(", ", labelled.Select(x => x.Child + " " + x.Number))}], the claim lists " +
            $"[{string.Join(", ", claimed)}]. Update AbsorptionTheoremClaim.FNumberedDirectChildren, its doc " +
            "comment and the Preface of docs/proofs/PROOF_ABSORPTION_THEOREM.md.");
        Assert.Equal(34, direct.Count);
        Assert.Equal(14, labelled.Count);

        // Each number is a registry entry, not only a label: its heading "### Fnn." exists.
        var root = RepoRootLocator.Find();
        Assert.NotNull(root);
        var formulaRegistry = File.ReadAllText(Path.Combine(root!, "docs", "ANALYTICAL_FORMULAS.md"));
        foreach (var number in found)
            Assert.True(Regex.IsMatch(formulaRegistry, $@"^### {Regex.Escape(number)}\.", RegexOptions.Multiline),
                $"{number} has no entry heading '### {number}.' in docs/ANALYTICAL_FORMULAS.md");

        // The control the direct-edge reading must fail on: F3 reaches the theorem only through F50,
        // so it is a descendant and not a direct child.
        var descendants = registry.DescendantsOf<AbsorptionTheoremClaim>()
            .Select(c => c.GetType()).ToHashSet();
        Assert.Contains(typeof(F3DecayRateBoundsPi2Inheritance), descendants);
        Assert.DoesNotContain(typeof(F3DecayRateBoundsPi2Inheritance), direct);
    }

    [Fact]
    public void RegisterAbsorptionTheoremClaim_F66AncestorsContainAbsorptionTheorem()
    {
        var registry = BuildBaseRegistry()
            .RegisterAbsorptionTheoremClaim()
            .RegisterF66PoleModesPi2Inheritance()
            .Build();

        var ancestors = registry.AncestorsOf<F66PoleModesPi2Inheritance>()
            .Select(c => c.GetType()).ToHashSet();

        Assert.Contains(typeof(AbsorptionTheoremClaim), ancestors);
    }
}
