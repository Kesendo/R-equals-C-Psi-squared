using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="TrianglePalindromeLocusWitness"/> (section "The triangle's palindromic
/// locus at N = 3" of <c>docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md</c>): the rows the proof and
/// <c>EndCountTests</c> pin, read on the witness's own exact route, and the grid with a planted mutation.</summary>
public class TrianglePalindromeLocusWitnessTests
{
    static BigRational[] T(string s) => s.Split(',').Select(p =>
    {
        var q = p.Split('/');
        return new BigRational(long.Parse(q[0]), q.Length > 1 ? long.Parse(q[1]) : 1);
    }).ToArray();

    [Theory]
    [InlineData("ZXY", "3,6,-2", "0,1,1", 1, 1, true, 'C')]     // the reciprocal family's pinned row
    [InlineData("ZXY", "3,6,-2", "0,0,0", 10, 10, true, 'A')]   // at h = 0 it is coloured
    [InlineData("ZXY", "3,6,-3", "0,1,1", 1, 0, false, null)]   // off the locus it breaks
    [InlineData("ZXY", "2,2,7", "0,1,1", 1, 1, true, 'B')]      // the swap branch
    [InlineData("XXX", "1,2,3", "1,2,3", 1, 1, true, 'A')]      // a colouring
    public void Pinned_Rows(string word, string j, string h, int near, int far, bool pal, char? cls)
    {
        var r = TrianglePalindromeLocusWitness.ReadRow(word, T(j), T(h));
        Assert.Equal((near, far, pal, cls), (r.Near, r.Far, r.Palindrome, r.InClass));
        Assert.True(r.Agrees);
    }

    [Fact]
    public void The_Grid_Agrees_And_Every_Class_Occurs()
    {
        var g = TrianglePalindromeLocusWitness.ReadGrid();
        Assert.Equal(27 * 8 * 6, g.Rows);
        Assert.Equal(0, g.Disagreements);
        Assert.True(g.A > 0 && g.B > 0 && g.C > 0, $"A {g.A}, B {g.B}, C {g.C}");
    }

    [Theory]
    [InlineData('B')]
    [InlineData('C')]
    public void Leaving_A_Class_Out_Breaks_The_Grid(char without)
    {
        Assert.True(TrianglePalindromeLocusWitness.ReadGrid(without).Disagreements > 0);
    }

    [Fact]
    public void A_Disconnected_Row_Is_Refused()
    {
        Assert.Throws<ArgumentException>(() => new TrianglePalindromeLocusWitness("ZXY", "0,0,3", "0,1,1"));
    }
}
