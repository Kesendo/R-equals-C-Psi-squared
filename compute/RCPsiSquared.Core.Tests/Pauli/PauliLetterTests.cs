using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Tests.Pauli;

public class PauliLetterTests
{
    [Theory]
    [InlineData(PauliLetter.I, 0, 0)]
    [InlineData(PauliLetter.X, 1, 0)]
    [InlineData(PauliLetter.Z, 0, 1)]
    [InlineData(PauliLetter.Y, 1, 1)]
    public void BitA_BitB_MatchFrameworkConvention(PauliLetter letter, int expectedBitA, int expectedBitB)
    {
        Assert.Equal(expectedBitA, letter.BitA());
        Assert.Equal(expectedBitB, letter.BitB());
    }

    [Theory]
    [InlineData(PauliLetter.I, 'I')]
    [InlineData(PauliLetter.X, 'X')]
    [InlineData(PauliLetter.Y, 'Y')]
    [InlineData(PauliLetter.Z, 'Z')]
    public void Symbol_Roundtrip(PauliLetter letter, char symbol)
    {
        Assert.Equal(symbol, letter.Symbol());
        Assert.Equal(letter, PauliLetterExtensions.FromSymbol(symbol));
    }

    [Theory]
    [InlineData((int)PauliLetter.I, 0)]
    [InlineData((int)PauliLetter.X, 1)]
    [InlineData((int)PauliLetter.Z, 2)]
    [InlineData((int)PauliLetter.Y, 3)]
    public void EnumValues_FollowFlatEncoding(int actual, int expected)
    {
        // Flat encoding: a + 2*b. I=(0,0)=0, X=(1,0)=1, Z=(0,1)=2, Y=(1,1)=3.
        Assert.Equal(expected, actual);
    }
}
