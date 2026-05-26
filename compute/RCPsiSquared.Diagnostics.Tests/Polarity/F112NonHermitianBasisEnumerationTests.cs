using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Diagnostics.Polarity;

namespace RCPsiSquared.Diagnostics.Tests.Polarity;

public class F112NonHermitianBasisEnumerationTests
{
    [Fact]
    public void BuildLHInPauliBasis_PauliStringXI_AtN2_HasExpectedShape()
    {
        // XI is the Pauli string "X ⊗ I" at N=2; its L_H = -i[XI, ·] in Pauli basis
        // must be a 16x16 complex matrix (4^N x 4^N at N=2).
        var letters = new[] { PauliLetter.X, PauliLetter.I };
        var H = PauliString.Build(letters);

        var L = F112NonHermitianBasisEnumeration.BuildLHInPauliBasis(H, N: 2);

        Assert.Equal(16, L.RowCount);
        Assert.Equal(16, L.ColumnCount);
    }
}
