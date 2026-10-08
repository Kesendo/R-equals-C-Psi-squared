using RCPsiSquared.Diagnostics.Foundation;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>From-below pins for <see cref="DepolarizingAttainmentWitness"/>, the F5 live lab
/// (registry F5, gate <c>simulations/f5_depolarizing_attainment.py</c>).
///
/// <para>The witness rebuilds the Hamiltonians as Gaussian-integer matrices, finds the commuting
/// global strings exactly, applies the dissipator to them exactly, and reads the spectrum with its
/// own eigensolver. The expected verdicts are the gate's CLASS verdicts and the registry's (the
/// generic and ladder Hamiltonians here are the witness's own seeded instances of those classes,
/// which the gate never ran row by row); the float thresholds are separations read against stated
/// models, with the measured values in the comments.</para>
///
/// <para><b>Two-sided throughout.</b> Every attained row sits beside a row where the bound is
/// strict, because a witness that had lost the ability to report "not attained" would pass a
/// one-sided suite unchanged. The sharpest pair is the Heisenberg chain under one smallest letter
/// everywhere against the SAME chain with the smallest letter differing between sites: the fast
/// span is one string in both, and only the uniform one commutes.</para></summary>
public class DepolarizingAttainmentWitnessTests
{
    private const string Rates = "3,7,11,5";

    // the registry's witness per chain: Z^N where the excitation number or its parity is
    // conserved, X^N for Ising in an X field; other uniform strings may commute as well
    // (Heisenberg with all three), so the test asks for membership, and for isingx for the
    // list to be exactly X^N
    [Theory]
    [InlineData("heisenberg", "ZZZZ", false)]
    [InlineData("xy", "ZZZZ", false)]
    [InlineData("xx", "ZZZZ", false)]
    [InlineData("ising", "ZZZZ", false)]
    [InlineData("xxz2", "ZZZZ", false)]
    [InlineData("dm", "ZZZZ", false)]
    [InlineData("heisdm", "ZZZZ", false)]
    [InlineData("isingx", "XXXX", true)]
    public void MeasuredChains_AttainTheBound_WithAnExactCertificate(string model, string certificate4, bool only)
    {
        foreach (int n in new[] { 2, 3, 4 })
        {
            var r = new DepolarizingAttainmentWitness(n, model, Rates).Read();
            Assert.Contains(certificate4[..n], r.Certificates);
            if (only) Assert.Single(r.Certificates);
            Assert.True(r.CertificateEigenvalueExact);
            Assert.True(r.CriterionDimension >= 1);
            Assert.True(r.BoundThirtieths > 0);
            // the float side: the shortfall equals the bound to the eigensolver's rounding, read as
            // the ratio of |shortfall - bound| to eps * ||L||_F. The model is a law, not a number: on
            // an attained row the edge eigenvalue lies on the boundary of L's numerical range, where
            // L acts normally, so it is semisimple with condition number 1 and eps ||L|| is its scale.
            // Measured 0.0 to 1.0 under MathNet on these rows at N = 2..4 (a LAPACK replication read
            // up to 3); the strict rows read 1e12 and more, the gate below
            Assert.True(r.ShortfallDeviationOverModel < 10.0,
                $"{model} N={n}: deviation/model = {r.ShortfallDeviationOverModel}");
            // the criterion's count is a separation: nearest uncounted |1 - s| against the worst
            // counted (floored at eps), at least six decades as the gate asks (measured 14 to 15),
            // and the worst counted value within 1e3 eps of 1 (measured 0 to 20 eps), so that a
            // count at 1e-7 could not pass as a count; the grouping of ad_H's eigenvalue differences
            // must leave decades between its gap and its inner spread (measured 13 to 16)
            Assert.True(r.CriterionSeparationDecades >= 6.0,
                $"{model} N={n}: criterion separation {r.CriterionSeparationDecades} decades");
            Assert.True(r.CriterionWorstCountedOverEps < 1e3,
                $"{model} N={n}: worst counted {r.CriterionWorstCountedOverEps} eps from 1");
            Assert.True(r.GroupingDecades >= 6.0, $"{model} N={n}: grouping {r.GroupingDecades} decades");
        }
    }

    [Fact]
    public void GenericChain_DoesNotAttain_NoStringNoEigenvectorAndAStrictShortfall()
    {
        foreach (int n in new[] { 2, 3, 4 })
        {
            var r = new DepolarizingAttainmentWitness(n, "generic", Rates).Read();
            Assert.True(r.Certificate is null, $"generic N={n}: certificate {r.Certificate}");
            Assert.True(r.CriterionDimension == 0,
                $"generic N={n}: criterion dim {r.CriterionDimension}, max s {r.CriterionMaxSingular}, ratio {r.ShortfallOverBound}");
            // the strict side read against the SAME rounding model as the attained rows: the
            // excess over the bound is the physics, decades above eps ||L||_F (measured 1e12 and
            // more; the gate reads 1.04 to 1.24 times the bound on these instances)
            Assert.True(r.ShortfallDeviationOverModel > 1e6,
                $"generic N={n}: deviation/model = {r.ShortfallDeviationOverModel}, ratio {r.ShortfallOverBound}");
            // and the largest singular value stays decades below the count (1 - s measured 0.019,
            // 0.11 and 0.23 against the 1e-6 cut)
            Assert.True(1 - r.CriterionMaxSingular > 1e-3, $"generic N={n}: max s = {r.CriterionMaxSingular}");
        }
    }

    [Fact]
    public void Ladder_AttainsWithoutAnyCommutingString()
    {
        // |0..0><1..1| is the witness: both are eigenvectors of H at distinct energies, and no
        // global string without an identity letter commutes (at N = 2 the ladder commutes with ZZ,
        // so the gate and this test start at N = 3)
        foreach (int n in new[] { 3, 4 })
        {
            var r = new DepolarizingAttainmentWitness(n, "ladder", Rates).Read();
            Assert.Null(r.Certificate);
            Assert.True(r.CriterionDimension >= 1);
            Assert.True(r.ShortfallDeviationOverModel < 10.0,
                $"ladder N={n}: deviation/model = {r.ShortfallDeviationOverModel}");
            Assert.True(r.CriterionSeparationDecades >= 6.0, $"ladder N={n}: separation {r.CriterionSeparationDecades}");
        }
    }

    [Fact]
    public void PauliChannel_OneSmallestLetterEverywhere_AttainsThroughThatString()
    {
        // X the strictly smallest rate on every site: the fast span is X^N, which commutes with
        // Heisenberg; the bound is 2 * sum of the smallest rates
        var r = new DepolarizingAttainmentWitness(3, "heisenberg", pauliRates: "1:5:7,2:9:4,1:3:8").Read();
        Assert.Equal("XXX", r.Certificate);
        Assert.True(r.CertificateEigenvalueExact);
        Assert.Equal(2 * (1 + 2 + 1) * 3, r.BoundThirtieths);        // tenths -> thirtieths
        Assert.True(r.CriterionDimension >= 1);
        Assert.True(r.ShortfallDeviationOverModel < 10.0, $"deviation/model = {r.ShortfallDeviationOverModel}");
    }

    [Fact]
    public void PauliChannel_SmallestLettersDiffer_DoesNotAttain()
    {
        // X smallest on site 0, Z on site 1, X on site 2: the fast span is the one mixed string XZX,
        // and Heisenberg commutes with a string without identity letters only when its letters agree
        var r = new DepolarizingAttainmentWitness(3, "heisenberg", pauliRates: "1:5:7,9:4:2,1:3:8").Read();
        Assert.Null(r.Certificate);
        Assert.True(r.CriterionDimension == 0,
            $"criterion dim {r.CriterionDimension}, separation {r.CriterionSeparationDecades}, ratio {r.ShortfallOverBound}");
        Assert.True(r.ShortfallDeviationOverModel > 1e6, $"deviation/model = {r.ShortfallDeviationOverModel}");
        Assert.True(1 - r.CriterionMaxSingular > 1e-3, $"max s = {r.CriterionMaxSingular}");
    }

    [Fact]
    public void PauliChannel_ADarkSiteAndATie_AttainThroughTheUniformString()
    {
        // site 0 dark (every rate zero): its fast letters are all four, the identity among them, and
        // X^N lies in the span and commutes; the bound counts the dark site at zero
        var dark = new DepolarizingAttainmentWitness(3, "heisenberg", pauliRates: "0:0:0,2:9:4,1:3:8").Read();
        Assert.Contains("XXX", dark.Certificates);
        Assert.True(dark.CertificateEigenvalueExact);
        Assert.Equal(2 * (0 + 2 + 1) * 3, dark.BoundThirtieths);
        Assert.True(dark.CriterionDimension >= 1);
        // X and Y tied smallest at site 0, X smallest elsewhere: the span {X, Y} x X x X holds X^N
        var tie = new DepolarizingAttainmentWitness(3, "heisenberg", pauliRates: "1:1:7,1:5:7,1:3:8").Read();
        Assert.Contains("XXX", tie.Certificates);
        Assert.True(tie.CertificateEigenvalueExact);
        // X and Z dephasing with Y dark on every site: the fastest letter is Y everywhere, the
        // bound is 0, and Y^N commutes with Heisenberg
        var ydark = new DepolarizingAttainmentWitness(3, "heisenberg", pauliRates: "3:0:5,2:0:4,1:0:8").Read();
        Assert.Equal(new[] { "YYY" }, ydark.Certificates);
        Assert.Equal(0, ydark.BoundThirtieths);
        Assert.True(ydark.CertificateEigenvalueExact);
        // and XXZ at Delta = 2 with one smallest letter everywhere attains as well (gate Stage D)
        var xxz = new DepolarizingAttainmentWitness(3, "xxz2", pauliRates: "7:5:1,9:4:2,3:8:1").Read();
        Assert.Contains("ZZZ", xxz.Certificates);
        Assert.True(xxz.CertificateEigenvalueExact);
    }

    [Fact]
    public void PauliChannel_ZDephasingOnly_BoundIsZeroAndAttained()
    {
        // every site's smallest rate is zero (X and Y dark), so the bound is 0 and the fast span is
        // {X, Y} per site; X^N commutes with Heisenberg, so the shortfall is exactly 0
        var r = new DepolarizingAttainmentWitness(3, "heisenberg", pauliRates: "0:0:5,0:0:3,0:0:7").Read();
        Assert.Equal(0, r.BoundThirtieths);
        Assert.NotNull(r.Certificate);
        Assert.True(r.CertificateEigenvalueExact);
        Assert.True(r.CriterionDimension >= 1);
    }

    [Fact]
    public void TheFastestRateIsReadOffTheDissipatorExactly_OnEveryRow()
    {
        // 2 sigma - (fastest decay rate among all 4^N strings) == the bound, in thirtieths, on an
        // isotropic row and on an anisotropic one
        foreach (var w in new[]
                 {
                     new DepolarizingAttainmentWitness(3, "heisenberg", Rates),
                     new DepolarizingAttainmentWitness(3, "heisenberg", pauliRates: "1:5:7,9:4:2,1:3:8"),
                     new DepolarizingAttainmentWitness(3, "generic", pauliRates: "0:0:5,0:0:3,0:0:7"),
                 })
        {
            var r = w.Read();
            Assert.Equal(r.BoundThirtieths, 2 * r.SigmaThirtieths - r.FastestRateThirtieths);
        }
    }

    [Fact]
    public void Guards()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new DepolarizingAttainmentWitness(5, "heisenberg"));
        Assert.Throws<ArgumentException>(() => new DepolarizingAttainmentWitness(3, "nosuchmodel"));
        Assert.Throws<ArgumentException>(() => new DepolarizingAttainmentWitness(3, "heisenberg", "3,7"));
        Assert.Throws<ArgumentException>(() => new DepolarizingAttainmentWitness(3, "heisenberg", "0,7,11"));
        Assert.Throws<ArgumentException>(() => new DepolarizingAttainmentWitness(3, "heisenberg", "a,7,11"));
        Assert.Throws<ArgumentException>(() => new DepolarizingAttainmentWitness(3, "heisenberg", pauliRates: "0:0:0,0:0:0,0:0:0"));
        Assert.Throws<ArgumentException>(() => new DepolarizingAttainmentWitness(3, "heisenberg", pauliRates: "1:2,1:2:3,1:2:3"));
    }
}
