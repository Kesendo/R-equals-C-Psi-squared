using MirrorWorld;

namespace MirrorWorldTests;

public class NeuralPalindromeTests
{
    [Fact]
    public void Constructed_Complex_Pair_Has_Exact_Dyadic_Residual()
    {
        double[,] j = { { -0.5, -0.25 }, { 0.25, -0.25 } };
        Assert.Equal(0.0, NeuralPalindrome.MaxResidual(j, [1, 0], 0.375));

        // Changing one coupling must reach the same residual door and break the identity.
        j[0, 1] = -0.125;
        Assert.Equal(0.125, NeuralPalindrome.MaxResidual(j, [1, 0], 0.375));
    }

    [Fact]
    public void Fixed_Seat_Diagonal_Can_Break_The_Identity()
    {
        double[,] j = { { -0.5, 0, 0 }, { 0, -0.25, 0 }, { 0, 0, -0.5 } };
        Assert.Equal(0.25, NeuralPalindrome.MaxResidual(j, [1, 0, 2], 0.375));
    }

    [Fact]
    public void Involution_Check_Distinguishes_A_Three_Cycle()
    {
        Assert.True(NeuralPalindrome.IsInvolution([1, 0, 3, 2]));
        Assert.True(NeuralPalindrome.IsInvolution([1, 0, 2]));
        Assert.False(NeuralPalindrome.IsInvolution([1, 2, 0]));
    }

    [Fact]
    public void Invalid_Permutations_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => NeuralPalindrome.IsInvolution(null!));
        Assert.Throws<ArgumentException>(() => NeuralPalindrome.IsInvolution([0, 0]));
        Assert.Throws<ArgumentException>(() => NeuralPalindrome.IsInvolution([-1, 0]));
        Assert.Throws<ArgumentException>(() => NeuralPalindrome.IsInvolution([1, 2]));
    }

    [Fact]
    public void Residual_Rejects_Invalid_Shapes_And_NonInvolutions()
    {
        Assert.Throws<ArgumentNullException>(() => NeuralPalindrome.MaxResidual(null!, [], 0));
        Assert.Throws<ArgumentException>(() => NeuralPalindrome.MaxResidual(new double[2, 3], [1, 0], 0));
        Assert.Throws<ArgumentException>(() => NeuralPalindrome.MaxResidual(new double[2, 2], [0], 0));
        Assert.Throws<ArgumentException>(() => NeuralPalindrome.MaxResidual(new double[3, 3], [1, 2, 0], 0));
        Assert.Throws<ArgumentException>(() => NeuralPalindrome.MaxResidual(new double[2, 2], [0, 0], 0));
    }

    [Theory]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    [InlineData(double.NegativeInfinity)]
    public void Residual_Rejects_Nonfinite_Entries_And_Centre(double value)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => NeuralPalindrome.MaxResidual(new double[1, 1], [0], value));
        Assert.Throws<ArgumentException>(() => NeuralPalindrome.MaxResidual(new double[,] { { value } }, [0], 0));
    }

    [Fact]
    public void Centre_And_PairSum_Have_The_Expected_Sign_And_Factor()
    {
        Assert.Equal(0.375, NeuralPalindrome.Centre(2, 4));
        Assert.Equal(-0.75, Formulas.F37_NeuralPairSum(2, 4));
        Assert.Equal(-2 * NeuralPalindrome.Centre(2, 4), Formulas.F37_NeuralPairSum(2, 4));
        // Decimal 0.3 is not binary exact; this pins the biological-parameter reading.
        Assert.Equal(-0.3, Formulas.F37_NeuralPairSum(5, 10), 15);
    }

    [Theory]
    [InlineData(0)]
    [InlineData(-1)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    [InlineData(double.NegativeInfinity)]
    public void Time_Constants_Must_Be_Finite_And_Positive(double tau)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => NeuralPalindrome.Centre(tau, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => NeuralPalindrome.Centre(2, tau));
        Assert.Throws<ArgumentOutOfRangeException>(() => Formulas.F37_NeuralPairSum(tau, 4));
        Assert.Throws<ArgumentOutOfRangeException>(() => Formulas.F37_NeuralPairSum(2, tau));
    }
}
