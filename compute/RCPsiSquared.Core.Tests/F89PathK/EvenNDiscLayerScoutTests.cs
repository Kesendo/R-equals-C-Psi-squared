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
