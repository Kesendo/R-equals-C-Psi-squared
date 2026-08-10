using System;
using RCPsiSquared.Core.F89PathK;
using Xunit;
using Xunit.Abstractions;

namespace RCPsiSquared.Core.Tests.F89PathK;

/// <summary>Scout for the sideways_spin_ladder arc's Sturm closing route: the even-N disc layers of
/// the (1,2) residual through the D-only certificate path (<see cref="FoldResultantCertificate"/>,
/// even-N guard lifted 2026-08-10 — the D-only path never touches the corner block).
///
/// <para>N=4 splits into a CONTROL and a derived twin: the R-even layers are sourced from the
/// committed closed form, disc(F₈) = const·q²⁴·(3q⁴+q²−1)²·P₂₀(q) with P₂₀(q) = P₁₀(q²)
/// (<c>F89Path3OcticGaloisClaim</c>, which scopes itself to the R-even sector: deg 52, v_q 24,
/// factorisation (24, 8, 20) ⇒ layers [20, 4]); the R-odd row is NOT separately committed — at even
/// N the R-odd residual is the exact complex CONJUGATE of the R-even one (measured; the F89 doc's
/// sector-swap conjugacy), and conjugation fixes every layer degree, so the identical reading is
/// asserted as the derived twin, not quoted as a second source.</para>
///
/// <para>N=6 is the first data for the certificate design (residual degree 32 = the F_32 oracle's
/// symmetric-sector degree; R-odd again the conjugate). Gate: the layer aggregation identity
/// deg D − v_q = Σ (m+1)·deg(layer_m), the same consistency the N=7 disc-multiplicity test pins,
/// plus layers.Length ≤ 2 (no multiplicity-3 layer, the β-exotic point read at even N).</para></summary>
public class EvenNDiscLayerScoutTests
{
    private readonly ITestOutputHelper _out;
    public EvenNDiscLayerScoutTests(ITestOutputHelper o) => _out = o;

    [Fact(DisplayName = "N=4: disc layers reproduce the octic factorization [20, 4] (R-even sourced, R-odd the conjugate twin)")]
    [Trait("Category", "EVEN_DISC_SCOUT")]
    public void N4_Control_KnownOcticFactorization()
    {
        foreach (bool rOdd in new[] { false, true })
        {
            var (degD, vD, layers) = FoldResultantCertificate.DiscLayersAtNthPrime(4, rOdd, 0);
            _out.WriteLine($"N=4 {(rOdd ? "R-odd (conjugate twin)" : "R-even (control)")}: deg_q D = {degD}, v_q = {vD}, layers [{string.Join(", ", layers)}]");
            Assert.Equal(52, degD);
            Assert.Equal(24, vD);
            Assert.Equal(new[] { 20, 4 }, layers);
        }
    }

    [Fact(DisplayName = "N=4 Re/Im control: the self-fold disc is REAL, Im D ≡ 0 mod p")]
    [Trait("Category", "EVEN_DISC_SCOUT")]
    public void N4_ReImControl_DiscIsReal()
    {
        var (p, _, degRe, vRe, degIm, _, imZero, _) =
            FoldResultantCertificate.DiscReImGcdAtNthPrime(4, rOdd: false, nth: 0);
        _out.WriteLine($"N=4 R-even: p={p}, deg Re = {degRe} (v_q {vRe}), deg Im = {degIm}, ImIsZero = {imZero}");
        Assert.True(imZero, "the N=4 disc must be real (self-fold antiunitary, spec B3)");
        Assert.Equal(52, degRe);
        Assert.Equal(24, vRe);
    }

    [Fact(DisplayName = "N=6 scout: gcd(Re D, Im D) mod p = 1 — the closing certificate's live question, pinned")]
    [Trait("Category", "SLOW_EVEN_DISC")]
    public void N6_Scout_ReImGcd()
    {
        var rows = new System.Collections.Generic.List<(int DegRe, int VRe, int DegIm, int VIm, int GcdDeg)>();
        foreach (bool rOdd in new[] { false, true })
        {
            var (p, _, degRe, vRe, degIm, vIm, imZero, gcdDeg) =
                FoldResultantCertificate.DiscReImGcdAtNthPrime(6, rOdd, 0);
            _out.WriteLine($"N=6 {(rOdd ? "R-odd " : "R-even")}: p={p}, deg Re = {degRe} (v_q {vRe}), "
                           + $"deg Im = {degIm} (v_q {vIm}), gcd deg = {gcdDeg}");
            Assert.False(imZero, "the N=6 disc must be genuinely complex (disc-reality re-gate, 13-order split)");
            // PINNED, observed 2026-08-10: gcd degree 0 at the first split prime. Mod-p only — the
            // window-free certificate additionally needs the degree- and valuation-preservation
            // halves (the method's scope note); a real-q coalescence appearing at N=6 would flip
            // this to a nonzero gcd degree and fail here first.
            Assert.Equal(0, gcdDeg);
            rows.Add((degRe, vRe, degIm, vIm, gcdDeg));
        }
        // the even-N conjugation makes R-odd the exact conjugate of R-even (Re shared, Im negated):
        // one measurement, one forced twin — gate the agreement rather than presenting two sources.
        Assert.Equal(rows[0], rows[1]);
    }

    [Fact(DisplayName = "N=6 scout: disc layers of the F_32 residual, first prime, both parities (layer identity gated)")]
    [Trait("Category", "SLOW_EVEN_DISC")]
    public void N6_Scout_DiscLayers()
    {
        foreach (bool rOdd in new[] { false, true })
        {
            var (degD, vD, layers) = FoldResultantCertificate.DiscLayersAtNthPrime(6, rOdd, 0);
            _out.WriteLine($"N=6 {(rOdd ? "R-odd " : "R-even")}: deg_q D = {degD}, v_q = {vD}, layers [{string.Join(", ", layers)}]");
            int aggregate = 0;
            for (int m = 0; m < layers.Length; m++) aggregate += (m + 1) * layers[m];
            Assert.True(degD - vD == aggregate,
                $"layer aggregation identity broken: {degD} − {vD} != Σ (m+1)·deg = {aggregate}");
            Assert.True(layers.Length <= 2,
                $"a multiplicity-≥3 disc layer at N=6 ({layers.Length} layers) — the β-exotic reading fails at even N");
        }
    }
}
