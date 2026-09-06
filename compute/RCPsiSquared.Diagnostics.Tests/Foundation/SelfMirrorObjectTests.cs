using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The centre-line object inherits its x/y/z frame from the system and keeps the two
/// different spectral orbits separate: the linear F1 map fixes only λ = −σ, while F1 followed by
/// conjugation fixes the whole line Re λ = −σ.</summary>
public class SelfMirrorObjectTests
{
    private static SelfMirrorObject Build(
        int n,
        double gamma = 0.1,
        HamiltonianType hamiltonianType = HamiltonianType.XY)
    {
        var h = new ChainSystem(n, 1.0, gamma, hamiltonianType, TopologyKind.Chain).BuildHamiltonian();
        var channels = Enumerable.Range(0, n).Select(l => new ChannelRate($"q{l}", gamma)).ToList();
        return new SelfMirrorObject(new MirrorSystem(n, h, channels));
    }

    [Fact]
    public void Object_InheritsFrame_AndSitsAtMinusSigma()
    {
        var obj = Build(4, gamma: 0.1);
        Assert.Equal(4, obj.N);              // N inherited from the system, not owned
        Assert.Equal(0.4, obj.Sigma, 10);    // σ = Σγ = Nγ inherited
        Assert.Equal(-0.4, obj.Center, 10);  // the object sits at Re λ = −σ
    }

    [Theory]
    [InlineData(2, true)]
    [InlineData(4, true)]
    [InlineData(3, false)]
    [InlineData(5, false)]
    public void CompositeFixedLine_PopulatedIffEvenN(int n, bool populated)
    {
        int count = Build(n).CompositeFixedLineCount;
        if (populated)
            Assert.True(count > 0, $"N={n} (even): the k=N/2 self-mirror sector should be populated, got {count}");
        else
            Assert.Equal(0, count);          // odd N: half-integer w_XY = N/2, the sector is empty
    }

    [Fact]
    public void N2_DistinguishesCompositeFixedLine_FromLinearF1FixedPoint()
    {
        var obj = Build(2, hamiltonianType: HamiltonianType.Heisenberg);

        Assert.Equal(10, obj.CompositeFixedLineCount);
        Assert.Equal(4, obj.LinearF1FixedPointCount);
        Assert.True(obj.CompositeFixedLineCount > obj.LinearF1FixedPointCount);
    }

    [Fact]
    public void LiveStrings_NameCompositeMap_AndFenceItFromLinearF1()
    {
        var obj = Build(2);
        string rendered = string.Join("\n", new[] { obj.DisplayName, obj.Summary }
            .Concat(obj.Children.Select(c => $"{c.DisplayName}\n{c.Summary}")));

        Assert.Contains("λ ↦ −2σ − conj(λ)", rendered);
        Assert.Contains("linear F1", rendered);
        Assert.Contains("λ = −σ", rendered);
        Assert.DoesNotContain("each its own mirror under λ ↦ −2σ − λ", rendered);
    }
}
