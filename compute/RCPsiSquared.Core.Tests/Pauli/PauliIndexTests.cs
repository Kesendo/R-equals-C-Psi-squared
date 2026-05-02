using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Tests.Pauli;

public class PauliIndexTests
{
    [Theory]
    [InlineData(2)]
    [InlineData(3)]
    [InlineData(4)]
    public void FromFlat_ToFlat_RoundTrip(int N)
    {
        long d2 = 1L << (2 * N);
        for (long k = 0; k < d2; k++)
        {
            var letters = PauliIndex.FromFlat(k, N);
            Assert.Equal(N, letters.Length);
            long back = PauliIndex.ToFlat(letters);
            Assert.Equal(k, back);
        }
    }

    [Fact]
    public void TotalBitA_CountsXAndY()
    {
        // String XYZI: X (bit_a=1), Y (bit_a=1), Z (bit_a=0), I (bit_a=0) → total 2.
        var letters = PauliLabel.Parse("XYZI");
        Assert.Equal(2, PauliIndex.TotalBitA(letters));
    }

    [Fact]
    public void TotalBitBParity_CountsYAndZ_Mod2()
    {
        // String XYZI: X (bit_b=0), Y (bit_b=1), Z (bit_b=1), I (bit_b=0) → 2 mod 2 = 0.
        var letters = PauliLabel.Parse("XYZI");
        Assert.Equal(0, PauliIndex.TotalBitBParity(letters));

        // String XYY: X (bit_b=0), Y (bit_b=1), Y (bit_b=1) → 2 mod 2 = 0.
        Assert.Equal(0, PauliIndex.TotalBitBParity(PauliLabel.Parse("XYY")));

        // String XZ: X (bit_b=0), Z (bit_b=1) → 1 mod 2 = 1.
        Assert.Equal(1, PauliIndex.TotalBitBParity(PauliLabel.Parse("XZ")));
    }

    [Fact]
    public void Label_ParseFormat_RoundTrip()
    {
        var letters = PauliLabel.Parse("IXYZ");
        Assert.Equal("IXYZ", PauliLabel.Format(letters));
    }
}
