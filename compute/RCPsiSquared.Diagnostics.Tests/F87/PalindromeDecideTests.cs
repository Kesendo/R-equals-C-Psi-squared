using System.Numerics;
using System.Collections.Generic;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.F87;

namespace RCPsiSquared.Diagnostics.Tests.F87;

public class PalindromeDecideTests
{
    private static PauliTerm T(params PauliLetter[] ls) => new(ls, Complex.One);
    private const PauliLetter I = PauliLetter.I, X = PauliLetter.X, Y = PauliLetter.Y, Z = PauliLetter.Z;

    [Fact]
    public void CertifyHard_True_OnDiagonalCellHardOddCyclePair()
    {
        // XXZ + XZX: Klein (0,1), Mixed, the §7 hard-odd-cycle witness (BipartiteChirality StandardSet).
        Assert.True(PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation(new[] { T(X, X, Z), T(X, Z, X) }));
    }

    [Fact]
    public void CertifyHard_False_OnDiagonalCellSoftPair()
    {
        // XXZ + ZXX: Klein (0,1), Mixed, soft (bipartite); equal valuation => not hard.
        Assert.False(PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation(new[] { T(X, X, Z), T(Z, X, X) }));
    }

    [Fact]
    public void CertifyHard_False_OutOfScope()
    {
        Assert.False(PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation(new[] { T(X, X, Z) }));               // 1 term
        Assert.False(PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation(new[] { T(X, X, X), T(X, X, Y), T(Y, X, X) })); // 3 terms
        Assert.False(PalindromeSoftCertifier.CertifyHardByDiagonalCellValuation(new[] { T(Z), T(Z) }));               // pure-diagonal (na<2), longitudinal field
    }
}
