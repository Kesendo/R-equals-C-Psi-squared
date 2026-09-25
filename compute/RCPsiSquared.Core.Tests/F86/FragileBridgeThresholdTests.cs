using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>The FRAGILE_BRIDGE threshold (F40's value at J_b = 1.0), recomputed from the C# builders: two 2-qubit Heisenberg
/// chains (J = 1) joined by a Heisenberg bridge J_b, chain A dephasing at +γ and chain B at −γ
/// (Σγ = 0). That is a 4-site chain with bonds (1, J_b, 1) and the γ profile (γ, γ, −γ, −γ).
/// <see cref="PauliHamiltonian.HeisenbergChain(int, IReadOnlyList{double})"/> writes J/4 per bond,
/// so the Pauli-convention couplings of <c>simulations/fragile_bridge_ep_signature.py</c> enter as
/// 4·(1, J_b, 1); <see cref="PauliDephasingDissipator.BuildZ"/> uses the producer's row-major vec
/// convention and takes the negative rates of the gain chain as they are.
///
/// <para>At the two couplings tested, J_b = 1.0 and 1.9, the first popcount block of L to go
/// unstable is a half of the (2,2) block, which the spin flip and S = G∘R (G the chain reflection
/// on both sides of ρ, R the one-sided flip ρ ↦ ρ·X^⊗N of F118) map to themselves, so the full L holds a single Jordan pair at the collision. Where the first block
/// has partners under the spin flip and S, the full L holds one Jordan pair per block of that
/// orbit: two at J_b = 0.1, 0.5, 1.2, 1.4, 5 and 10 among the producer's samples, four in the
/// (0,1) window near 1.46, with as many singular values at the rounding level, which the gates
/// below would reject by design. At five exact couplings (3/4, √5/2, 4/3 and two cubic roots)
/// the threshold is zero with no EP at all; neither tested coupling is one of them.</para>
///
/// <para>There is no exact route here (an eigensolver on a non-normal matrix), so every gate is a
/// stated law, not a bare number:</para>
/// <list type="bullet">
/// <item><b>Location.</b> Newton on the real γ axis for Re f, f = (λ_a − λ_b)² of the pair nearest the
/// running midpoint. f is analytic in γ with a simple zero at the collision, and along the real
/// axis it is real (negative below, positive above). By the √-law the pair gap is ∝ √|γ − γ*|, so a
/// gap a hundred times below the gap at γ*(1 + 10⁻⁶) puts the located point within 10⁻¹⁰ of the
/// collision in relative γ. The located γ* is also compared with the producer's, located on the
/// same sign change inside the pair's own block, a second route with its own builder.</item>
/// <item><b>Jordan rank.</b> The smallest singular value of L(γ*) − λ*·I sits at the rounding level
/// ε·‖L‖₂ of a backward-stable SVD. The second one exceeds the pair gap by three orders of
/// magnitude: a semisimple double eigenvalue with two independent eigenvectors would put it at the
/// pair's own scale. So the double eigenvalue carries one eigenvector: nullity 1 at algebraic
/// multiplicity 2, isolated from the third eigenvalue by three orders of magnitude. (At J_b = 1.0
/// the full-L second and third singular values, both 2.9·10⁻³, come from the (1,2) block and its
/// flip copy (3,2), whose nearest eigenvalues sit 3.07·10⁻³ away; inside the pair's own half of
/// the (2,2) block the second is 0.176, so the gate holds with room.)</item>
/// <item><b>√-law.</b> Half the pair's split over √|δ|, δ = γ/γ* − 1, is the same at δ = 10⁻⁴ and
/// 10⁻⁶, below γ* (split along Im, both eigenvalues on the imaginary axis) and above (split along
/// Re, a mirror pair at one Im). The next Puiseux term changes it by O(δ): 3.3·10⁻⁴ relative at
/// δ = 10⁻³ for J_b = 1.9 in the producer's output, so about 3·10⁻⁵ at the δ = 10⁻⁴ read here; a
/// transversal crossing or a diabolic touch (split ∝ δ) would change the ratio by a factor of 10
/// between the two decades.</item>
/// </list></summary>
public class FragileBridgeThresholdTests
{
    private const double Eps = 2.220446049250313e-16;

    /// <summary>L(γ) = L0 + γ·D. The Hamiltonian part has a purely imaginary diagonal and the
    /// dissipator a real one, so the difference of the two builds is the gain-loss dissipator D
    /// exactly.</summary>
    private static (ComplexMatrix L0, ComplexMatrix D) Parts(double jBridge)
    {
        var h = PauliHamiltonian.HeisenbergChain(4, new[] { 4.0, 4.0 * jBridge, 4.0 }).ToMatrix();
        var l0 = PauliDephasingDissipator.BuildZ(h, new[] { 0.0, 0.0, 0.0, 0.0 });
        var l1 = PauliDephasingDissipator.BuildZ(h, new[] { 1.0, 1.0, -1.0, -1.0 });
        return (l0, l1 - l0);
    }

    private static ComplexMatrix At(ComplexMatrix l0, ComplexMatrix d, double gamma) =>
        l0 + new Complex(gamma, 0.0) * d;

    private static (Complex A, Complex B, Complex[] Spectrum) PairNear(ComplexMatrix l, Complex anchor)
    {
        var spectrum = l.Evd().EigenValues.ToArray();
        var near = spectrum.OrderBy(z => (z - anchor).Magnitude).Take(2).ToArray();
        return (near[0], near[1], spectrum);
    }

    private static double SquaredGapReal(Complex a, Complex b) => ((a - b) * (a - b)).Real;

    /// <summary>Newton on the real γ axis for Re (λ_a − λ_b)², central-difference slope.</summary>
    private static (double GammaStar, Complex LambdaStar) Locate(
        ComplexMatrix l0, ComplexMatrix d, double gammaSeed, Complex anchor)
    {
        const double h = 1e-7;
        double g = gammaSeed;
        Complex lam = anchor;
        for (int iteration = 0; iteration < 40; iteration++)
        {
            var (a, b, _) = PairNear(At(l0, d, g), lam);
            lam = 0.5 * (a + b);
            var (ap, bp, _) = PairNear(At(l0, d, g + h), lam);
            var (am, bm, _) = PairNear(At(l0, d, g - h), lam);
            double slope = (SquaredGapReal(ap, bp) - SquaredGapReal(am, bm)) / (2.0 * h);
            double step = SquaredGapReal(a, b) / slope;
            g -= step;
            if (Math.Abs(step) <= 1e-14 * g) break;
        }
        var (a2, b2, _) = PairNear(At(l0, d, g), lam);
        return (g, 0.5 * (a2 + b2));
    }

    [Theory]
    [InlineData(1.0, 0.1873101, 2.6515, 0.187310108345)]
    [InlineData(1.9, 0.4058524, 6.9835, 0.405853184723)]
    public void Threshold_IsADefectiveEp2OnTheRealGammaAxis(
        double jBridge, double gammaSeed, double omegaSeed, double gammaStarProducer)
    {
        var (l0, d) = Parts(jBridge);
        var (gStar, lamStar) = Locate(l0, d, gammaSeed, new Complex(0.0, -omegaSeed));
        var lStar = At(l0, d, gStar);
        var (a, b, spectrum) = PairNear(lStar, lamStar);
        double gap = (a - b).Magnitude;
        double third = spectrum.Select(z => (z - lamStar).Magnitude).OrderBy(x => x).Skip(2).First();
        var sv = (lStar - lamStar * Matrix<Complex>.Build.DenseIdentity(lStar.RowCount))
            .Svd(computeVectors: false).S.Select(s => s.Real).OrderBy(s => s).ToArray();
        double norm = sv[^1];
        var (pa, pb, _) = PairNear(At(l0, d, gStar * (1.0 + 1e-6)), lamStar);
        double gapAtMicro = (pa - pb).Magnitude;

        // Location: on the collision, and in agreement with the producer's sign-change locator.
        // Both locators find the zero of the same analytic Re f, whose rounding noise moves that zero
        // by its noise over |f'| (8e-16 at J_b = 1.0, 2.6e-15 at 1.9, rms), budgeted at 1e-13; the
        // producer's γ* is quoted to 12 decimals, which adds half a unit of the last one.
        Assert.True(gap <= 1e-2 * gapAtMicro,
            $"J_b={jBridge}: pair gap {gap:E2} at γ*={gStar:R} vs {gapAtMicro:E2} at γ*(1+1e-6)");
        Assert.True(Math.Abs(gStar - gammaStarProducer) <= 0.5e-12 + 1e-13,
            $"J_b={jBridge}: γ*={gStar:R} vs the producer's {gammaStarProducer:R}");
        // The collision sits on the imaginary axis: Σγ = 0 centres the palindrome at 0.
        Assert.True(Math.Abs(lamStar.Real) <= 1e3 * Eps * norm, $"Re λ* = {lamStar.Real:E2}");

        // Jordan rank: one null direction at the rounding level, none at the pair's scale, and the
        // pair isolated from the rest of the spectrum.
        Assert.True(sv[0] <= 10.0 * Eps * norm, $"J_b={jBridge}: σ_min = {sv[0]:E2}, ε‖L‖ = {Eps * norm:E2}");
        Assert.True(sv[1] >= 1e3 * gap, $"J_b={jBridge}: σ_2 = {sv[1]:E2} vs pair gap {gap:E2}");
        Assert.True(third >= 1e3 * gap, $"J_b={jBridge}: third eigenvalue at {third:E2} vs pair gap {gap:E2}");

        // The √-law, one Puiseux branch through the collision.
        var halfSplit = new Dictionary<(int Side, double Delta), double>();
        foreach (double delta in new[] { 1e-4, 1e-6 })
        {
            foreach (int side in new[] { -1, 1 })
            {
                var (u, v, _) = PairNear(At(l0, d, gStar * (1.0 + side * delta)), lamStar);
                double split = (u - v).Magnitude;
                if (side < 0)
                {
                    // Below: two distinct eigenvalues on the imaginary axis.
                    Assert.True(Math.Max(Math.Abs(u.Real), Math.Abs(v.Real)) <= 1e-6 * split,
                        $"J_b={jBridge}, δ=−{delta}: pair ({u}, {v}) not on the axis");
                }
                else
                {
                    // Above: a mirror pair ±Re at one Im.
                    Assert.True(Math.Abs(u.Imaginary - v.Imaginary) <= 1e-6 * split,
                        $"J_b={jBridge}, δ=+{delta}: pair ({u}, {v}) not a mirror pair");
                }
                halfSplit[(side, delta)] = split / 2.0 / Math.Sqrt(delta);
            }
        }
        foreach (int side in new[] { -1, 1 })
        {
            double decadeRatio = halfSplit[(side, 1e-4)] / halfSplit[(side, 1e-6)];
            Assert.InRange(decadeRatio, 0.99, 1.01);
        }
        Assert.InRange(halfSplit[(-1, 1e-6)] / halfSplit[(1, 1e-6)], 0.99, 1.01);
    }

    [Fact]
    public void LocalGlobalEpLink_CarriesTheRecomputedThreshold()
    {
        var (l0, d) = Parts(1.0);
        var (gStar, lamStar) = Locate(l0, d, 0.1873101, new Complex(0.0, -2.6515));
        // The constants carry 12 and 9 decimals: half a unit of the last one each, plus the locator's
        // rounding noise (γ*: 8e-16 rms; the pair's midpoint: 6e-15 rms), budgeted at 1e-13.
        Assert.True(Math.Abs(gStar - LocalGlobalEpLink.FragileBridgeThresholdGamma) <= 0.5e-12 + 1e-13,
            $"γ* = {gStar:R} vs the claim's {LocalGlobalEpLink.FragileBridgeThresholdGamma:R}");
        Assert.True(Math.Abs(lamStar.Imaginary + LocalGlobalEpLink.FragileBridgeThresholdOmega) <= 0.5e-9 + 1e-13,
            $"Im λ* = {lamStar.Imaginary:R} vs the claim's −{LocalGlobalEpLink.FragileBridgeThresholdOmega:R}");
    }
}
