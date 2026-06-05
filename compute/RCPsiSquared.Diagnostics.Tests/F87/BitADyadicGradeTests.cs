using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.F87;

/// <summary>The graded bit_a (<see cref="BitADyadicGrade"/>) is the (1 + x)-adic valuation read as a
/// tower of dyadic moments, with the Klein bit_a as its bottom rung. These tests tie it to the GF(2)[x]
/// valuation (<see cref="WindowedObstructionScan.ValuationAtOnePlusX"/>), to the Klein bit_a, and to the
/// F115 valuation classes c_v = 2^(k-1-v), and check the dyadic-distance reading.</summary>
public class BitADyadicGradeTests
{
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void Grade_EqualsGF2Valuation_AndBitAIsTheBottomRung(int k)
    {
        for (ulong mask = 1; mask < (1UL << k); mask++)
        {
            Assert.Equal(WindowedObstructionScan.ValuationAtOnePlusX(mask), BitADyadicGrade.Grade(mask, k));
            int bitA = BitOperations.PopCount(mask) & 1;            // Klein bit_a = #(X+Y) mod 2
            Assert.Equal(bitA, BitADyadicGrade.BitA(mask, k));
            Assert.Equal(bitA, BitADyadicGrade.Moment(mask, 0, k)); // moment_0 = total parity = bit_a
        }
    }

    [Theory]
    [InlineData(4)]
    [InlineData(5)]
    [InlineData(6)]
    public void GradeClasses_ReproduceTheF115ValuationClasses(int k)
    {
        var classCount = new Dictionary<int, int>();
        for (ulong mask = 1; mask < (1UL << k); mask++)
            if ((BitOperations.PopCount(mask) & 1) == 0)            // the bit_a = 0 cell (even popcount)
            {
                int v = BitADyadicGrade.Grade(mask, k);
                classCount[v] = classCount.GetValueOrDefault(v) + 1;
            }
        for (int v = 1; v <= k - 1; v++)
            Assert.Equal(1 << (k - 1 - v), classCount[v]);          // c_v = 2^(k-1-v)
    }

    [Fact]
    public void DyadicPairs_GradeAsTheirDistance()
    {
        // A pair at distance 2^m is (1 + x)^(2^m), of grade 2^m.
        Assert.Equal(1, BitADyadicGrade.Grade(0b11, 6));     // adjacent pair, distance 1
        Assert.Equal(2, BitADyadicGrade.Grade(0b101, 6));    // distance-2 pair
        Assert.Equal(4, BitADyadicGrade.Grade(0b10001, 6));  // distance-4 pair
    }
}
