using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class Pi2Z4KleinDistinctionClaimTests
{
    private readonly ITestOutputHelper _out;

    public Pi2Z4KleinDistinctionClaimTests(ITestOutputHelper output) => _out = output;

    [Fact]
    public void Tier_IsTier1Derived()
    {
        var c = new Pi2Z4KleinDistinctionClaim();
        Assert.Equal(Tier.Tier1Derived, c.Tier);
    }

    [Fact]
    public void OrderOfBothGroups_IsFour()
    {
        // Both Z₄ and Klein Z₂ × Z₂ have exactly 4 elements; this coincidence motivated
        // Tom's "are they the same partition?" question.
        Assert.Equal(4, Pi2Z4KleinDistinctionClaim.OrderOfBothGroups);
    }

    [Fact]
    public void Z4MaxElementOrder_IsFour_KleinMaxIsTwo()
    {
        // The structural distinction: Z₄ has an order-4 element (i), Klein does not.
        Assert.Equal(4, Pi2Z4KleinDistinctionClaim.Z4MaxElementOrder);
        Assert.Equal(2, Pi2Z4KleinDistinctionClaim.KleinMaxNonIdentityElementOrder);
    }

    [Fact]
    public void AreIsomorphic_IsFalse()
    {
        // Tom's hypothesis answered: Z₄ and Klein Z₂ × Z₂ are not isomorphic. The
        // four Z₄ sectors and the four Klein cells are different partitions of the
        // Pauli space, even though both have cardinality 4^(N−1).
        var c = new Pi2Z4KleinDistinctionClaim();
        Assert.False(c.AreIsomorphic);
    }

    [Theory]
    [InlineData(0, 0)]    // 1² = 1
    [InlineData(1, 1)]    // i² = −1
    [InlineData(2, 0)]    // (−1)² = 1
    [InlineData(3, 1)]    // (−i)² = −1
    [InlineData(4, 0)]    // i^4 = 1 → cyclic wrap
    [InlineData(5, 1)]    // i^5 = i → squared = −1
    [InlineData(-1, 1)]   // i^{−1} = −i → squared = −1
    [InlineData(-2, 0)]   // i^{−2} = −1 → squared = 1
    public void SquareInZ4ToZ2_FoldsZ4OntoOneZ2Factor(int z4Power, int expectedZ2)
    {
        // The Z₄ → Z₂ squaring map. {1, −1} (Z₄ subgroup of order 2) maps to 1 in Z₂;
        // {i, −i} (the order-4 elements) map to −1. This is exactly Π² as a homomorphism
        // from the cyclic Π Z₄ down to the involutive Π²_Z = squaring image.
        var c = new Pi2Z4KleinDistinctionClaim();
        Assert.Equal(expectedZ2, c.SquareInZ4ToZ2(z4Power));
    }

    [Theory]
    [InlineData(1, 1)]      // 4^0 = 1 (one sector at N=1, trivial)
    [InlineData(2, 4)]      // 4^1 = 4 (each sector has 4 elements at N=2)
    [InlineData(3, 16)]     // 4^2 = 16 (the ERROR_CORRECTION_PALINDROME N=3 case)
    [InlineData(4, 64)]
    [InlineData(5, 256)]
    public void SectorDimension_EqualsFourToTheNMinusOne(int N, int expected)
    {
        // Both Z₄ and Klein decompose 4^N Pauli space into 4 sectors of equal size.
        // At N=3: 16-dimensional sectors. The cardinality coincidence is a Lagrange
        // theorem consequence of any order-4 group acting regularly, NOT evidence that
        // the two partitions are the same.
        var c = new Pi2Z4KleinDistinctionClaim();
        Assert.Equal(expected, c.SectorDimension(N));
    }

    [Fact]
    public void SectorDimension_NLessThanOne_Throws()
    {
        var c = new Pi2Z4KleinDistinctionClaim();
        Assert.Throws<ArgumentOutOfRangeException>(() => c.SectorDimension(0));
    }

    [Fact]
    public void Anchor_References_Z4_AndKleinFour_AndErrorCorrectionPalindrome()
    {
        var c = new Pi2Z4KleinDistinctionClaim();
        Assert.Contains("Pi2I4MemoryLoopClaim", c.Anchor);
        Assert.Contains("Pi2KnowledgeBaseClaims.cs", c.Anchor);
        Assert.Contains("ERROR_CORRECTION_PALINDROME.md", c.Anchor);
        Assert.Contains("ANALYTICAL_FORMULAS.md", c.Anchor);
    }

    [Fact]
    public void Reconnaissance_EmitsDistinctionTable()
    {
        // Documents the answer to Tom's question as a queryable Schicht-1 table:
        // Z₄ and Klein are both order 4, both partition 4^N into 4-sector structures,
        // but they are not isomorphic. Connection: Z₄ → Z₂ via squaring → one Klein axis.
        var c = new Pi2Z4KleinDistinctionClaim();
        _out.WriteLine("");
        _out.WriteLine("    structure   |  order  |  generator order  |  axes covered");
        _out.WriteLine("    ------------+---------+-------------------+---------------");
        _out.WriteLine($"    Z₄ (Π)     |    {Pi2Z4KleinDistinctionClaim.OrderOfBothGroups}    |        {Pi2Z4KleinDistinctionClaim.Z4MaxElementOrder}          |  Z only (one axis, fine)");
        _out.WriteLine($"    Z₂ × Z₂    |    {Pi2Z4KleinDistinctionClaim.OrderOfBothGroups}    |        {Pi2Z4KleinDistinctionClaim.KleinMaxNonIdentityElementOrder}          |  Z and X (two axes, coarse)");
        _out.WriteLine("");
        _out.WriteLine("    sector dim by N:");
        for (int N = 1; N <= 5; N++)
            _out.WriteLine($"      N={N}:  4^(N−1) = {c.SectorDimension(N)}");
        _out.WriteLine("");
        _out.WriteLine($"    AreIsomorphic = {c.AreIsomorphic}");
        _out.WriteLine("    bridge Z₄ → Z₂: squaring (1↦1, i↦−1, −1↦1, −i↦−1)");
    }
}
