using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="OperatorSchmidtRank"/>: operators written by hand on the Pauli strings,
/// their ranks known without the routine.</summary>
public class OperatorSchmidtRankTests
{
    // letters 0 = I, 1 = X, 2 = Y, 3 = Z; site 0 the most significant base-4 digit
    static BigRational[] Op(int n, params (string Letters, long C)[] terms)
    {
        var v = Enumerable.Repeat(BigRational.Zero, (int)Math.Pow(4, n)).ToArray();
        foreach (var (letters, c) in terms)
        {
            int s = 0;
            foreach (var ch in letters) s = s * 4 + "IXYZ".IndexOf(ch);
            v[s] = v[s] + new BigRational(c, 1);
        }
        return v;
    }

    [Fact]
    public void A_Product_Has_Rank_One_Across_Every_Split()
    {
        // (I + X) ⊗ Z ⊗ (Y + 2Z), multiplied out
        var v = Op(3, ("IZY", 1), ("IZZ", 2), ("XZY", 1), ("XZZ", 2));
        foreach (var left in new[] { new[] { 0 }, new[] { 1 }, new[] { 2 }, new[] { 0, 2 } })
            Assert.Equal(1, OperatorSchmidtRank.Across(v, 3, left));
    }

    [Fact]
    public void An_Entangled_Pair_Is_Seen_Only_Across_Its_Own_Cut()
    {
        // (XX + ZZ) on sites 0, 2 and Y on site 1: rank 2 where the cut separates 0 from 2, 1 where it does not
        var v = Op(3, ("XYX", 1), ("ZYZ", 1));
        Assert.Equal(2, OperatorSchmidtRank.Across(v, 3, new[] { 0 }));
        Assert.Equal(1, OperatorSchmidtRank.Across(v, 3, new[] { 1 }));
        Assert.Equal(2, OperatorSchmidtRank.Across(v, 3, new[] { 2 }));
        Assert.Equal(1, OperatorSchmidtRank.Across(v, 3, new[] { 0, 2 }));
    }

    [Fact]
    public void The_Site_Order_Is_Read_As_Written()
    {
        // (XX + ZZ) on sites 0, 1 and Y on site 2: not symmetric under 0 <-> 2, so a reversed digit order turns it red
        var v = Op(3, ("XXY", 1), ("ZZY", 1));
        Assert.Equal(2, OperatorSchmidtRank.Across(v, 3, new[] { 0 }));
        Assert.Equal(2, OperatorSchmidtRank.Across(v, 3, new[] { 1 }));
        Assert.Equal(1, OperatorSchmidtRank.Across(v, 3, new[] { 2 }));
    }

    [Fact]
    public void The_Swap_Has_Rank_Four()
    {
        // 2·SWAP = II + XX + YY + ZZ
        var v = Op(2, ("II", 1), ("XX", 1), ("YY", 1), ("ZZ", 1));
        Assert.Equal(4, OperatorSchmidtRank.AcrossCut(v, 2, 1));
    }
}
