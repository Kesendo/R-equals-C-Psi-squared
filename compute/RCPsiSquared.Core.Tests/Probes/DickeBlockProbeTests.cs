using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Probes;

namespace RCPsiSquared.Core.Tests.Probes;

public class DickeBlockProbeTests
{
    [Theory]
    [InlineData(5, 1)]   // c=2
    [InlineData(7, 3)]   // c=4
    [InlineData(8, 4)]   // c=4
    public void AllEntriesEqual_OneOverTwoSqrtMpMq(int N, int n)
    {
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var probe = DickeBlockProbe.Build(block);

        double expected = 1.0 / (2.0 * Math.Sqrt((double)block.Basis.Mp * block.Basis.Mq));
        Assert.Equal(block.Basis.MTotal, probe.Count);
        for (int i = 0; i < probe.Count; i++)
        {
            Assert.Equal(expected, probe[i].Real, 12);
            Assert.Equal(0.0, probe[i].Imaginary, 12);
        }
    }

    [Theory]
    [InlineData(5, 1, 0.25)]   // 1/(2·√(5·10)) · √(50) = 1/(2·√(50))·√(50) = 1/2; |v|² = sum of (1/2√50)² · 50 = 50/(4·50) = 0.25
    [InlineData(7, 3, 0.25)]
    public void NormSquared_IsOneQuarter(int N, int n, double expected)
    {
        // ||probe||² = MTotal · (1/(2·√(Mp·Mq)))² = (Mp·Mq) / (4·Mp·Mq) = 1/4
        var block = new CoherenceBlock(N, n, gammaZero: 0.05);
        var probe = DickeBlockProbe.Build(block);
        double normSq = 0;
        for (int i = 0; i < probe.Count; i++)
            normSq += probe[i].Real * probe[i].Real + probe[i].Imaginary * probe[i].Imaginary;
        Assert.Equal(expected, normSq, 12);
    }
}
