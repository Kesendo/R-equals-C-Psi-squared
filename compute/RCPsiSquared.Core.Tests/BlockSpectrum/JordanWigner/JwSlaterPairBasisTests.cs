using RCPsiSquared.Core.BlockSpectrum.JordanWigner;
using RCPsiSquared.Core.Knowledge;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.BlockSpectrum.JordanWigner;

/// <summary>Witnesses that <see cref="JwSlaterPairBasis"/> builds an orthonormal U mapping
/// the (p_c, p_r) sector computational basis to JW Slater-pair basis, diagonalises the
/// Hamiltonian-only Liouvillian L_H, and reproduces the dispersion eigenvalues
/// −i·(Σε(L) − Σε(K)).</summary>
public class JwSlaterPairBasisTests
{
    private readonly ITestOutputHelper _out;
    public JwSlaterPairBasisTests(ITestOutputHelper output) => _out = output;

    [Theory]
    [InlineData(4, 1, 1)]  // (p_c, p_r) = (1, 1): single-excitation/single-excitation block, dim 16
    [InlineData(4, 1, 2)]  // c=2 F86 coherence block, dim 24
    [InlineData(4, 2, 2)]  // dim 36
    [InlineData(5, 1, 1)]  // dim 25
    [InlineData(5, 1, 2)]  // dim 50 = c=2 N=5
    [InlineData(5, 2, 2)]  // dim 100
    [InlineData(5, 2, 3)]  // dim 100
    [InlineData(6, 3, 3)]  // dim 400
    public void Build_AtSmallN_PassesAllWitnesses(int N, int pCol, int pRow)
    {
        var basis = JwSlaterPairBasis.Build(N, pCol, pRow);
        _out.WriteLine($"N={N}, (p_c,p_r)=({pCol},{pRow}), dim={basis.U.RowCount}: " +
            $"orth={basis.OrthonormalityResidual:G3}, " +
            $"H-diag={basis.MhTotalDiagonalityResidual:G3}, " +
            $"eig-match={basis.MhTotalEigenvalueMatchResidual:G3}");
        Assert.True(basis.OrthonormalityResidual < JwSlaterPairBasis.Tolerance);
        Assert.True(basis.MhTotalDiagonalityResidual < JwSlaterPairBasis.Tolerance);
        Assert.True(basis.MhTotalEigenvalueMatchResidual < JwSlaterPairBasis.Tolerance);
        Assert.Equal(Tier.Tier1Derived, basis.Tier);
    }

    [Fact]
    public void Build_LargeBlock_N7_p3p4_StillPassesWitnesses()
    {
        // N=7 (p_c=3, p_r=4): C(7,3)*C(7,4) = 35*35 = 1225 — the largest block at N=7.
        var basis = JwSlaterPairBasis.Build(N: 7, pCol: 3, pRow: 4);
        _out.WriteLine($"N=7, (p_c,p_r)=(3,4), dim={basis.U.RowCount}: " +
            $"orth={basis.OrthonormalityResidual:G3}, " +
            $"H-diag={basis.MhTotalDiagonalityResidual:G3}, " +
            $"eig-match={basis.MhTotalEigenvalueMatchResidual:G3}");
        Assert.True(basis.OrthonormalityResidual < JwSlaterPairBasis.Tolerance);
        Assert.True(basis.MhTotalDiagonalityResidual < JwSlaterPairBasis.Tolerance);
        Assert.True(basis.MhTotalEigenvalueMatchResidual < JwSlaterPairBasis.Tolerance);
    }
}
