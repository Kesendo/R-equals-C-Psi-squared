using RCPsiSquared.Core.Lindblad;

namespace RCPsiSquared.Core.Tests.Lindblad;

public class HardwareDissipatorsTests
{
    [Fact]
    public void All_HasFiveCanonicalChannels()
    {
        var names = HardwareDissipators.All.Select(d => d.Name).ToHashSet();
        Assert.Contains("T1", names);
        Assert.Contains("T1pump", names);
        Assert.Contains("Tphi", names);
        Assert.Contains("Xnoise", names);
        Assert.Contains("Ynoise", names);
    }

    [Theory]
    [MemberData(nameof(AllDissipators))]
    public void Tabulated_C1C2_MatchClosedFormFromPauliDecomposition(HardwareDissipator d)
    {
        var (c1, c2) = DissipatorClosedForms.C1C2FromPauli(d.Alpha, d.Beta, d.Delta);
        Assert.Equal(d.C1, c1, 10);
        Assert.Equal(d.C2, c2, 10);
    }

    public static IEnumerable<object[]> AllDissipators() =>
        HardwareDissipators.All.Select(d => new object[] { d });
}
