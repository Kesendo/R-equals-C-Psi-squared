using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Resonance;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.F86;

/// <summary>F86 Petermann-factor analysis at deep post-EP. The two-regime K-curve structure
/// at c=2 N≥7 (orbit-1 saturation plateau at Q ≈ 7-8) suggests a near-EP region where
/// eigenvectors collapse onto each other. Petermann factor K_n = 1/|⟨l_n|r_n⟩|² (unit-norm
/// left/right eigvecs) diverges at exact EP and is large near EPs. If the saturation
/// plateau is a near-EP region, max K_n should spike around Q ≈ 7-8.
///
/// <para>Connection to <see href="hypotheses/FRAGILE_BRIDGE.md">FRAGILE_BRIDGE</see>: there
/// the Hopf bifurcation in the complex γ plane signals an EP via Petermann K = 403 above
/// γ_crit. Same diagnostic, different parameter axis (Q vs γ).</para>
/// </summary>
public class F86PetermannProbe(ITestOutputHelper output)
{
    [Fact]
    public void Probe_PetermannFactor_C2N7_VsQ()
    {
        var block = new CoherenceBlock(N: 7, n: 1, gammaZero: 0.05);
        var decomp = block.Decomposition;

        // Q points: pre-EP, around canonical Q_EP, plateau region, deep post-EP.
        var qValues = new[]
        {
            0.3, 0.5, 0.71, 1.0, 1.5, 2.0, 3.0, 5.0,
            6.0, 6.5, 7.0, 7.2, 7.5, 8.0, 8.5, 9.0,
            10.0, 12.0, 15.0, 20.0,
        };

        output.WriteLine($"c=2 N=7, block-L dim {decomp.D.RowCount}");
        output.WriteLine("Q       | max K_n | mode index | top-3 K_n values");
        output.WriteLine("--------+---------+------------+--------------------");

        foreach (double q in qValues)
        {
            double j = q * block.GammaZero;
            var L = decomp.D + (Complex)j * decomp.MhTotal;

            var evd = L.Evd();
            var R = evd.EigenVectors;
            var Rinv = R.Inverse();

            int dim = R.RowCount;
            var k_n = new double[dim];
            for (int n = 0; n < dim; n++)
            {
                double rNorm2 = 0;
                double lNorm2 = 0;
                for (int i = 0; i < dim; i++)
                {
                    double rMag = R[i, n].Magnitude;
                    rNorm2 += rMag * rMag;
                    double lMag = Rinv[n, i].Magnitude;
                    lNorm2 += lMag * lMag;
                }
                k_n[n] = rNorm2 * lNorm2;  // ⟨l_n|r_n⟩ = 1 in biorthogonal, so K_n = ||r_n||² · ||l_n||²
            }

            int iMax = 0;
            for (int n = 1; n < dim; n++) if (k_n[n] > k_n[iMax]) iMax = n;

            var top3 = k_n
                .Select((k, idx) => (K: k, Index: idx))
                .OrderByDescending(x => x.K)
                .Take(3)
                .ToArray();
            string top3Str = string.Join(", ", top3.Select(x => $"K[{x.Index}]={x.K:F1}"));

            output.WriteLine($"{q,6:F2}  | {k_n[iMax],7:F1} | {iMax,10} | {top3Str}");
        }
    }
}
