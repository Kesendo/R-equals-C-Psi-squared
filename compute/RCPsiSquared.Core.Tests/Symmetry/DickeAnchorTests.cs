using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Core.Tests.Symmetry;

public class DickeAnchorTests
{
    [Theory]
    [InlineData(3, 1, DickeAnchor.Mirror)]          // 2n+1=N, odd
    [InlineData(5, 2, DickeAnchor.Mirror)]
    [InlineData(7, 3, DickeAnchor.Mirror)]
    [InlineData(4, 1, DickeAnchor.KIntermediate)]   // N=4, n=N/2-1
    [InlineData(4, 2, DickeAnchor.KIntermediate)]   // N=4, n=N/2
    [InlineData(6, 2, DickeAnchor.KIntermediate)]
    [InlineData(6, 3, DickeAnchor.KIntermediate)]
    [InlineData(8, 3, DickeAnchor.KIntermediate)]
    [InlineData(8, 4, DickeAnchor.KIntermediate)]
    [InlineData(3, 0, DickeAnchor.Generic)]
    [InlineData(4, 0, DickeAnchor.Generic)]
    [InlineData(5, 0, DickeAnchor.Generic)]
    [InlineData(6, 0, DickeAnchor.Generic)]
    [InlineData(7, 0, DickeAnchor.Generic)]
    public void Classify_ThreeAnchors(int N, int n, DickeAnchor expected)
    {
        Assert.Equal(expected, DickeAnchorExtensions.Classify(N, n));
    }

    [Theory]
    [InlineData(DickeAnchor.Mirror, 0.0)]
    [InlineData(DickeAnchor.KIntermediate, 3.0 / 8.0)]
    [InlineData(DickeAnchor.Generic, 0.5)]
    public void AlphaTotal_MatchesClosedForm(DickeAnchor anchor, double expected)
    {
        Assert.Equal(expected, anchor.AlphaTotal(), 12);
    }

    [Theory]
    [InlineData(DickeAnchor.Mirror, 1.0)]
    [InlineData(DickeAnchor.KIntermediate, 0.5)]
    [InlineData(DickeAnchor.Generic, 0.0)]
    public void Gamma_MatchesXTotalOverlap(DickeAnchor anchor, double expected)
    {
        Assert.Equal(expected, anchor.Gamma(), 12);
    }

    [Theory]
    [InlineData(DickeAnchor.Mirror)]
    [InlineData(DickeAnchor.KIntermediate)]
    [InlineData(DickeAnchor.Generic)]
    public void ClosedForm_AlphaEqualsOneMinusGammaSquaredOverTwo(DickeAnchor anchor)
    {
        double gamma = anchor.Gamma();
        double expectedAlpha = (1 - gamma * gamma) / 2;
        Assert.Equal(expectedAlpha, anchor.AlphaTotal(), 12);
    }

    [Theory]
    [InlineData(3, 1, 0.0)]                   // mirror
    [InlineData(4, 1, 3.0 / 8.0)]             // K-intermediate
    [InlineData(4, 2, 3.0 / 8.0)]
    [InlineData(6, 2, 3.0 / 8.0)]
    [InlineData(6, 3, 3.0 / 8.0)]
    [InlineData(7, 0, 0.5)]                   // generic
    public void Pi2OddTotalDickeSuperposition_DispatchesViaDickeAnchor(int N, int n, double expected)
    {
        Assert.Equal(expected, PopcountCoherencePi2Odd.Pi2OddTotalDickeSuperposition(N, n), 12);
        Assert.Equal(expected, DickeAnchorExtensions.Classify(N, n).AlphaTotal(), 12);
    }
}
