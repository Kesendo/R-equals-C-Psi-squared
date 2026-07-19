using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below gates for <see cref="RecordLawWitness"/> (F135 + F136, the record laws
/// live). The battery itself is two-sided (each case compares the closed-form classifier
/// against an independent full-state path, and the catalogue mixes luminous and dark
/// configurations), so "all pass" is not a tautology; these tests additionally pin the
/// canonical priced numbers and a handful of classifier verdicts directly.</summary>
public class RecordLawWitnessTests
{
    private const double G = 0.05;                        // the canonical dephasing rate
    private const double PointerPriced = 0.7680396679759238;
    private const double BellPriced = 0.6241464301142645;

    [Fact]
    public void Battery_AllCasesPass()
    {
        var witness = new RecordLawWitness();
        Assert.True(witness.Cases.Count >= 18);
        Assert.Equal(witness.Cases.Count, witness.PassCount);
    }

    [Fact]
    public void Battery_MixesLuminousAndDarkCases()
    {
        // the two-sidedness guard: the battery must contain exact-zero verdicts (dark) AND
        // full-bit verdicts, so a broken engine cannot pass by uniformity.
        var witness = new RecordLawWitness();
        Assert.Contains(witness.Cases, c => c.Expected.StartsWith("Dark"));
        Assert.Contains(witness.Cases, c => c.Expected.Contains("I=1.000000"));
    }

    [Fact]
    public void CanonicalPrices_MatchTheClosedForms()
    {
        var witness = new RecordLawWitness(G);
        Assert.Equal(PointerPriced, witness.PointerPricedBits, 12);
        Assert.Equal(BellPriced, witness.BellPricedBits, 12);
    }

    [Fact]
    public void Classifier_UniformCatalogue()
    {
        var chain3 = new (int, int, double)[] { (0, 1, 1.0), (1, 2, 1.0) };
        var triangle = new (int, int, double)[] { (0, 1, 1.0), (0, 2, 1.0), (1, 2, 1.0) };
        var square = new (int, int, double)[] { (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0) };
        var pentagon = new (int, int, double)[] { (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 4, 1.0), (4, 0, 1.0) };

        var leaf = RecordLawWitness.Classify(chain3, 1, 0);
        Assert.Equal(RecordLawWitness.RecordFamily.Pointer, leaf.Family);
        Assert.Equal("ZY", leaf.Channel);
        Assert.Equal(+1, leaf.Sign);

        var swap = RecordLawWitness.Classify(chain3, 0, 1);
        Assert.Equal(RecordLawWitness.RecordFamily.RoleSwap, swap.Family);
        Assert.Equal("YZ", swap.Channel);

        var tri = RecordLawWitness.Classify(triangle, 0, 1);
        Assert.Equal(RecordLawWitness.RecordFamily.Bell, tri.Family);
        Assert.Equal("YY", tri.Channel);
        Assert.Equal(+1, tri.Sign);

        var plaq = RecordLawWitness.Classify(square, 0, 2);
        Assert.Equal(RecordLawWitness.RecordFamily.Bell, plaq.Family);
        Assert.Equal("XX", plaq.Channel);

        Assert.Equal(RecordLawWitness.RecordFamily.Dark, RecordLawWitness.Classify(pentagon, 0, 1).Family);
        Assert.Equal(RecordLawWitness.RecordFamily.Dark, RecordLawWitness.Classify(pentagon, 0, 2).Family);
    }

    [Fact]
    public void Classifier_SignWalk_K21()
    {
        // K₂,₁ at r = 1, 3, 5: YY = +1, −1, +1 (the signed-coherence corollary made observable).
        foreach ((double r, int sign) in new[] { (1.0, +1), (3.0, -1), (5.0, +1) })
        {
            var bonds = new (int, int, double)[] { (0, 2, 1.0), (1, 2, r) };
            var reading = RecordLawWitness.Classify(bonds, 0, 1);
            Assert.Equal(RecordLawWitness.RecordFamily.Bell, reading.Family);
            Assert.Equal("YY", reading.Channel);
            Assert.Equal(sign, reading.Sign);
        }
    }

    [Fact]
    public void Classifier_PricesFollowTheTwoSpecies()
    {
        var square = new (int, int, double)[] { (0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (3, 0, 1.0) };
        var chain3 = new (int, int, double)[] { (0, 1, 1.0), (1, 2, 1.0) };

        var bell = RecordLawWitness.Classify(square, 0, 2, gammaS: G, gammaJ: G);
        Assert.Equal(BellPriced, bell.Bits, 12);                    // both sites pay

        var pointer = RecordLawWitness.Classify(chain3, 1, 0, gammaS: 0.3, gammaJ: G);
        Assert.Equal(PointerPriced, pointer.Bits, 12);              // γ_S exactly invisible
    }

    [Fact]
    public void Ctor_RejectsNegativeGamma()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new RecordLawWitness(-0.1));
    }
}
