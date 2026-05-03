using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class Pi2ProjectionTests
{
    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void Projectors_SumToIdentity_AndAreOrthogonal(int N)
    {
        var pPlus = Pi2Projection.SymmetricProjector(N);
        var pMinus = Pi2Projection.AntisymmetricProjector(N);
        long d2 = 1L << (2 * N);
        var I = ComplexMatrix.Build.SparseIdentity((int)d2);

        Assert.True((pPlus + pMinus - I).FrobeniusNorm() < 1e-12,
            $"P_+ + P_- ≠ I at N={N}");
        Assert.True((pPlus * pMinus).FrobeniusNorm() < 1e-12,
            $"P_+ · P_- ≠ 0 at N={N}");
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void Projectors_AreIdempotent(int N)
    {
        var pPlus = Pi2Projection.SymmetricProjector(N);
        var pMinus = Pi2Projection.AntisymmetricProjector(N);
        Assert.True((pPlus * pPlus - pPlus).FrobeniusNorm() < 1e-12);
        Assert.True((pMinus * pMinus - pMinus).FrobeniusNorm() < 1e-12);
    }

    [Theory]
    [InlineData(1, 2L, 2L)]    // I, X are even (bit_b=0); Y, Z are odd (bit_b=1)
    [InlineData(2, 8L, 8L)]    // half/half: 4^2 / 2 = 8 each
    [InlineData(3, 32L, 32L)]  // 4^3 / 2 = 32 each
    [InlineData(4, 128L, 128L)]
    public void Counts_AreBalancedHalfHalf_ForZDephasing(int N, long expectedEven, long expectedOdd)
    {
        var (even, odd) = Pi2Projection.Counts(N, PauliLetter.Z);
        Assert.Equal(expectedEven, even);
        Assert.Equal(expectedOdd, odd);
    }

    [Fact]
    public void Counts_RotateUnderDephaseLetter_ButRemainHalfHalf()
    {
        // Each Klein cell has 2 of 4 letters; the sum-mod-2 over N letters partitions strings
        // into balanced half/half regardless of which Klein index we count.
        foreach (var letter in new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z })
        {
            var (even, odd) = Pi2Projection.Counts(N: 3, dephaseLetter: letter);
            Assert.Equal(32L, even);
            Assert.Equal(32L, odd);
        }
    }

    [Fact]
    public void Split_PauliEvenString_GoesEntirelyToEven()
    {
        // X⊗X is bit_b=0+0=0 → Π²-even under Z-dephasing.
        int N = 2;
        var XX = PauliString.Build(new[] { PauliLetter.X, PauliLetter.X });
        var (even, odd) = Pi2Projection.Split(XX, N);
        Assert.True((even - XX).FrobeniusNorm() < 1e-12);
        Assert.True(odd.FrobeniusNorm() < 1e-12);
    }

    [Fact]
    public void Split_PauliOddString_GoesEntirelyToOdd()
    {
        // X⊗Y is bit_b=0+1=1 → Π²-odd under Z-dephasing.
        int N = 2;
        var XY = PauliString.Build(new[] { PauliLetter.X, PauliLetter.Y });
        var (even, odd) = Pi2Projection.Split(XY, N);
        Assert.True(even.FrobeniusNorm() < 1e-12);
        Assert.True((odd - XY).FrobeniusNorm() < 1e-12);
    }

    [Fact]
    public void Split_RoundTrip_EvenPlusOddEqualsInput()
    {
        // Build a mixed-class operator: H = X⊗X (Π²-even) + X⊗Y (Π²-odd).
        // Decomposition must round-trip exactly.
        int N = 2;
        var XX = PauliString.Build(new[] { PauliLetter.X, PauliLetter.X });
        var XY = PauliString.Build(new[] { PauliLetter.X, PauliLetter.Y });
        var H = XX + 2.0 * XY;

        var (even, odd) = Pi2Projection.Split(H, N);
        Assert.True((even + odd - H).FrobeniusNorm() < 1e-12);
        // Each piece must match its input fragment.
        Assert.True((even - XX).FrobeniusNorm() < 1e-12);
        Assert.True((odd - 2.0 * XY).FrobeniusNorm() < 1e-12);
    }

    [Fact]
    public void Split_FrobeniusOrthogonality_HoldsBetweenEvenAndOdd()
    {
        // The Π²-even and Π²-odd Pauli subspaces are mutually Frobenius-orthogonal because
        // any cross term Tr(σ_α† σ_β) vanishes for distinct α, β (Pauli orthogonality).
        int N = 3;
        var XYZ = PauliString.Build(new[] { PauliLetter.X, PauliLetter.Y, PauliLetter.Z });
        var XXX = PauliString.Build(new[] { PauliLetter.X, PauliLetter.X, PauliLetter.X });
        var H = XXX + 0.7 * XYZ;

        var (even, odd) = Pi2Projection.Split(H, N);
        Complex hsInner = (even.ConjugateTranspose() * odd).Trace();
        Assert.True(hsInner.Magnitude < 1e-12,
            $"Π²-even and Π²-odd should be Frobenius-orthogonal; got tr(even† · odd) = {hsInner}");
    }

    [Theory]
    [InlineData(1)]
    [InlineData(2)]
    [InlineData(3)]
    public void SymmetricProjector_AppliedToVecForm_MatchesPi2EigenSign(int N)
    {
        // Internal consistency: P_+² in Pauli-vector space picks exactly the strings
        // with eigenvalue +1.
        var pPlus = Pi2Projection.SymmetricProjector(N);
        long d2 = 1L << (2 * N);
        for (long k = 0; k < d2; k++)
        {
            int eig = PiOperator.SquaredEigenvalue(PauliIndex.FromFlat(k, N));
            double diag = pPlus[(int)k, (int)k].Real;
            Assert.Equal(eig == +1 ? 1.0 : 0.0, diag, precision: 12);
        }
    }
}
