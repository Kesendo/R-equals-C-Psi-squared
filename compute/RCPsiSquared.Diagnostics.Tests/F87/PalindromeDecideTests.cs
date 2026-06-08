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

    [Fact]
    public void Decide_Hard_OnDiagonalCellHardPair_ExhibitsObstruction()
    {
        var d = PalindromeSoftCertifier.Decide(new[] { T(X, X, Z), T(X, Z, X) }, n: 4);
        Assert.Equal(PalindromeSoftCertifier.Decision.Hard, d.Verdict);
        Assert.Equal(PalindromeSoftCertifier.HardStrategy.DiagonalCellValuation, d.HardStrategy);
        Assert.Contains("valuation", d.Reason);          // the exhibited (1+x)-valuation obstruction
    }

    [Fact]
    public void Decide_Soft_OnXYModel()
    {
        var d = PalindromeSoftCertifier.Decide(new[] { T(X, X), T(Y, Y) }, n: 4);
        Assert.Equal(PalindromeSoftCertifier.Decision.Soft, d.Verdict);
        Assert.NotEqual(PalindromeSoftCertifier.SoftStrategy.None, d.SoftStrategy);
    }

    [Fact]
    public void Decide_Undetermined_OutOfScope()
    {
        // A 3-term frustrated hard set: soft strategies decline AND it is out of the F115 pair scope.
        var d = PalindromeSoftCertifier.Decide(new[] { T(X, X, X), T(X, X, Y), T(Y, X, X) }, n: 4);
        Assert.Equal(PalindromeSoftCertifier.Decision.Undetermined, d.Verdict);
    }
}
