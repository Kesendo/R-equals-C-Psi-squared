using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using Xunit;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.Lindblad;

/// <summary>The physical (no-jump) generator builder. The end-to-end reading, where the generator
/// meets F155's closed form, is <c>PhysicalGeneratorPolarityBreakWitnessTests</c> in
/// RCPsiSquared.Diagnostics.Tests; what is pinned here is the shape and the guards.</summary>
public class NoJumpGeneratorTests
{
    private static ComplexMatrix Pauli(params PauliLetter[] letters) =>
        PauliString.Build(letters);

    private static ComplexMatrix Zero(int d) => Matrix<Complex>.Build.Dense(d, d);

    [Fact]
    public void WithBZero_ItIsExactlyTheOrdinaryCommutator()
    {
        // The B = 0 face must reproduce the commutator the dephasing dissipator builds, entry for
        // entry: same Kronecker order, same sign. If the two ever drift apart, every polarity
        // number read through this path moves with them.
        var a = 0.7 * Pauli(PauliLetter.Z, PauliLetter.I) + 1.3 * Pauli(PauliLetter.X, PauliLetter.X);
        var mine = NoJumpGenerator.Build(a, Zero(4));
        var theirs = PauliDephasingDissipator.Build(a, new[] { 0.0, 0.0 });

        Assert.Equal(0.0, (mine - theirs).FrobeniusNorm());
    }

    [Fact]
    public void WithAZero_TheAntiHermitianHalfLiftsWithoutASign()
    {
        // B enters as B⊗I + I⊗Bᵀ, the transpose-paired ANTI-commutator lift, and nothing else.
        var b = Pauli(PauliLetter.Z, PauliLetter.I);
        var identity = Matrix<Complex>.Build.DenseIdentity(4);
        var expected = b.KroneckerProduct(identity) + identity.KroneckerProduct(b.Transpose());

        Assert.Equal(0.0, (NoJumpGenerator.Build(Zero(4), b) - expected).FrobeniusNorm());
    }

    [Fact]
    public void TheGeneratorIsLinearInEachHalfSeparately()
    {
        var a = Pauli(PauliLetter.X, PauliLetter.I);
        var b = Pauli(PauliLetter.Z, PauliLetter.I);

        var split = NoJumpGenerator.Build(a, Zero(4)) + NoJumpGenerator.Build(Zero(4), b);
        Assert.Equal(0.0, (NoJumpGenerator.Build(a, b) - split).FrobeniusNorm());
    }

    [Fact]
    public void ItRefusesANonHermitianHalf()
    {
        // The whole point of the A + iB split is that both halves are Hermitian; a caller who has
        // not split H would otherwise get a silently meaningless generator.
        var notHermitian = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
        {
            { Complex.Zero, Complex.One },
            { Complex.Zero, Complex.Zero },
        });

        var ex = Assert.Throws<ArgumentException>(
            () => NoJumpGenerator.Build(notHermitian, Zero(2)));
        Assert.Contains("not exactly Hermitian", ex.Message);
    }

    [Fact]
    public void ItRefusesMismatchedHalves() =>
        Assert.Throws<ArgumentException>(() => NoJumpGenerator.Build(Zero(4), Zero(2)));

    [Fact]
    public void ItRefusesANonSquareMatrix() =>
        Assert.Throws<ArgumentException>(
            () => NoJumpGenerator.Build(Matrix<Complex>.Build.Dense(2, 3), Zero(2)));

    [Fact]
    public void ItRefusesADimensionThatIsNotAPowerOfTwo() =>
        Assert.Throws<ArgumentException>(() => NoJumpGenerator.Build(Zero(3), Zero(3)));

    [Fact]
    public void ItRefusesNull() =>
        Assert.Throws<ArgumentNullException>(() => NoJumpGenerator.Build(null!, Zero(2)));
}
