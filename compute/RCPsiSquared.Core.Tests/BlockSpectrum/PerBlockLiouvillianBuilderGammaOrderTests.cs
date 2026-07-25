using System.Numerics;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using Xunit;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>Which index does <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/>'s
/// <c>gammaPerSite</c> run over: the SITE (site 0 = leftmost Kronecker factor, the convention of
/// <see cref="PauliString.SiteOp"/> and therefore of
/// <see cref="PauliDephasingDissipator.BuildZ"/>) or the BIT of the basis index (the block
/// builder reads <c>(diff >> l) &amp; 1</c>, and bit l is site N−1−l)?
///
/// <para>Under a uniform profile the two coincide, which is why every existing caller and the
/// native-memory parity test are blind to the difference; a non-uniform profile separates them.
/// The block builder's own summary says it "matches" <c>PauliDephasingDissipator.BuildZ</c>, so
/// the two must agree entry for entry on a common set of Liouville-space flat indices, with the
/// SAME array handed to both.</para></summary>
public class PerBlockLiouvillianBuilderGammaOrderTests
{
    /// <summary>Every flat index r·d + c of the joint-popcount sector (p_row, p_col).</summary>
    private static int[] SectorFlatIndices(int n, int pRow, int pCol)
    {
        int d = 1 << n;
        var flat = new List<int>();
        for (int r = 0; r < d; r++)
        {
            if (System.Numerics.BitOperations.PopCount((uint)r) != pRow) continue;
            for (int c = 0; c < d; c++)
                if (System.Numerics.BitOperations.PopCount((uint)c) == pCol)
                    flat.Add(r * d + c);
        }
        return flat.ToArray();
    }

    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void BuildBlockZ_WithNonUniformGamma_MatchesFullDissipatorEntryForEntry(int n)
    {
        // A deliberately non-palindromic profile: reversing it changes the array, so a
        // site/bit mix-up cannot hide.
        var gamma = new double[n];
        for (int l = 0; l < n; l++) gamma[l] = 0.05 + 0.11 * l;

        var h = PauliHamiltonian.XYChain(n, 0.75).ToMatrix();
        var full = PauliDephasingDissipator.BuildZ(h, gamma);
        var flat = SectorFlatIndices(n, 1, 1);
        var block = PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, flat);

        double worst = 0.0;
        for (int i = 0; i < flat.Length; i++)
            for (int j = 0; j < flat.Length; j++)
                worst = System.Math.Max(worst, (block[i, j] - full[flat[i], flat[j]]).Magnitude);

        Assert.True(worst < 1e-12,
            $"N={n}: block builder and full dissipator differ by {worst:0.0e0} under a non-uniform " +
            "gamma profile; the two disagree on what gammaPerSite[l] indexes");
    }

    /// <summary>The control that makes the test above two-sided: the same comparison against the
    /// REVERSED profile must FAIL, otherwise the profile was effectively palindromic and the test
    /// could not have detected a mix-up in the first place.</summary>
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    public void BuildBlockZ_AgainstReversedGamma_Differs(int n)
    {
        var gamma = new double[n];
        for (int l = 0; l < n; l++) gamma[l] = 0.05 + 0.11 * l;
        var reversed = (double[])gamma.Clone();
        System.Array.Reverse(reversed);

        var h = PauliHamiltonian.XYChain(n, 0.75).ToMatrix();
        var full = PauliDephasingDissipator.BuildZ(h, reversed);
        var flat = SectorFlatIndices(n, 1, 1);
        var block = PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, flat);

        double worst = 0.0;
        for (int i = 0; i < flat.Length; i++)
            for (int j = 0; j < flat.Length; j++)
                worst = System.Math.Max(worst, (block[i, j] - full[flat[i], flat[j]]).Magnitude);

        Assert.True(worst > 1e-6,
            $"N={n}: the reversed profile is indistinguishable (max diff {worst:0.0e0}); the probe " +
            "profile is too symmetric to separate the site and bit conventions");
    }
}
