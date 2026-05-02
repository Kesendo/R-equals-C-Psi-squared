using System.Numerics;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Polarity;

/// <summary>Polarity-axis diagnostic: read the +0 / 0 / −0 structure of a state.
///
/// <para>The polarity layer (<c>docs/THE_POLARITY_LAYER.md</c>) has three locations along
/// the d=0 axis: +0 (X = +1), 0 (X = 0, boundary), −0 (X = −1). Reflection at the ends;
/// 0 is the boundary that gets crossed. The qubit is the natural window: ⟨X_i⟩ per site
/// is the coordinate along the −0 → 0 → +0 axis.</para>
///
/// <para>No specific F-formula entry; the d=0 / polarity reading is its own structural
/// layer alongside F4 (sector projectors). Cross-references:
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F4 (kernel modes), <c>docs/THE_POLARITY_LAYER.md</c>
/// (the +0/0/−0 layered reading), <c>docs/PRIMORDIAL_QUBIT.md</c> §9 (Inside-Observability).</para>
///
/// <para>Verkettung von <see cref="PauliString.SitePaulis"/> + Bloch-component readout per
/// site — no new primitive, composition only.</para>
/// </summary>
public sealed record PolarityResult(
    double[] PolarityAxis,
    double[] OffAxis,
    double[] DistanceToPlusZero,
    double[] DistanceToMinusZero,
    bool[] OnBoundary,
    double AggregatePolarity,
    double OnAxisFraction,
    double[][] SiteBlochs);

public static class PolarityDiagnostic
{
    public static PolarityResult FromDensityMatrix(ComplexMatrix rho)
    {
        int d = rho.RowCount;
        int N = (int)Math.Round(Math.Log2(d));
        if ((1 << N) != d) throw new ArgumentException($"rho dim {d} not 2^N");
        var sitePaulis = PauliString.SitePaulis(N);

        var siteBlochs = new double[N][];
        var polarityAxis = new double[N];
        var offAxis = new double[N];
        var distPlus = new double[N];
        var distMinus = new double[N];
        var onBoundary = new bool[N];

        for (int i = 0; i < N; i++)
        {
            double x = (sitePaulis[i].X * rho).Trace().Real;
            double y = (sitePaulis[i].Y * rho).Trace().Real;
            double z = (sitePaulis[i].Z * rho).Trace().Real;
            siteBlochs[i] = new[] { x, y, z };
            polarityAxis[i] = x;
            offAxis[i] = Math.Sqrt(y * y + z * z);
            distPlus[i] = 1.0 - x;
            distMinus[i] = 1.0 + x;
            onBoundary[i] = Math.Abs(x) < 1e-10;
        }

        double aggPolarity = polarityAxis.Average(Math.Abs);
        double onAxisFraction = 1.0 - offAxis.Average(o => o * o);

        return new PolarityResult(
            PolarityAxis: polarityAxis,
            OffAxis: offAxis,
            DistanceToPlusZero: distPlus,
            DistanceToMinusZero: distMinus,
            OnBoundary: onBoundary,
            AggregatePolarity: aggPolarity,
            OnAxisFraction: onAxisFraction,
            SiteBlochs: siteBlochs);
    }
}
