using System.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Core.Tests.Pauli;

public class PauliHamiltonianTests
{
    [Fact]
    public void XYChain_N2_HasExpectedTerms()
    {
        var H = PauliHamiltonian.XYChain(N: 2, J: 1.0);
        // 1 bond × 2 terms (XX, YY) = 2 PauliTerms, each with coefficient 0.5.
        Assert.Equal(2, H.Terms.Count);
        Assert.Equal("XX", H.Terms[0].Label);
        Assert.Equal("YY", H.Terms[1].Label);
        Assert.Equal(0.5, H.Terms[0].Coefficient.Real, 12);
        Assert.Equal(0.5, H.Terms[1].Coefficient.Real, 12);
    }

    [Fact]
    public void XYChain_N2_MatrixIsHermitian()
    {
        var H = PauliHamiltonian.XYChain(N: 2, J: 1.0).ToMatrix();
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
            {
                Assert.Equal(H[i, j].Real, H[j, i].Real, 12);
                Assert.Equal(H[i, j].Imaginary, -H[j, i].Imaginary, 12);
            }
    }

    [Fact]
    public void HeisenbergChain_N2_HasThreeTermsPerBond()
    {
        var H = PauliHamiltonian.HeisenbergChain(N: 2, J: 1.0);
        Assert.Equal(3, H.Terms.Count);
        var labels = H.Terms.Select(t => t.Label).ToList();
        Assert.Contains("XX", labels);
        Assert.Contains("YY", labels);
        Assert.Contains("ZZ", labels);
    }

    [Theory]
    [InlineData("II", 0, 0, 0, 0)]   // k_body=0, no Y
    [InlineData("XY", 2, 1, 1, 1)]   // k_body=2, n_y=1, π2_parity=1, y_parity=1
    [InlineData("YYY", 3, 3, 1, 1)]  // k_body=3, n_y=3, π2_parity=3 mod 2=1, y_parity=3 mod 2=1
    [InlineData("ZZ", 2, 0, 0, 0)]   // k_body=2, n_y=0, π2_parity=2 mod 2=0
    public void PauliTerm_StructuralProperties(string label, int kBody, int ny, int pi2, int yParity)
    {
        var term = new PauliTerm(PauliLabel.Parse(label), Complex.One);
        Assert.Equal(kBody, term.KBody);
        Assert.Equal(ny, term.Ny);
        Assert.Equal(pi2, term.Pi2Parity);
        Assert.Equal(yParity, term.YParity);
    }
}
