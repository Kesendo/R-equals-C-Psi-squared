using System;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Numerics;
using M = RCPsiSquared.Core.Numerics.GaussianRationalMatrix;

namespace RCPsiSquared.Core.Quadrupolar;

/// <summary>The rank-2 tensors T₂₀ = (3I_z² − I(I + 1))/√6, T₂,±₁ = ∓(I_zI_± + I_±I_z)/2,
/// T₂,±₂ = I_±²/2. Exact: √6·T₂₀ (diagonal, the same in both bases) and T₂q, q = ±1, ±2, in the
/// gauge of <see cref="SpinSpace"/>, where the image of T₂q† is (−1)^q·T₂,−q. Float: T₂q in the
/// physical basis, q = 0, ±1, ±2.</summary>
public static class QuadrupolarTensors
{
    public static GaussianRational[,] Root6T20(SpinSpace s)
    {
        var iz = s.GaugedIz();
        return M.Subtract(M.Scale(M.Multiply(iz, iz), 3), M.Scale(M.Identity(s.Levels), s.SpinSquared));
    }

    public static GaussianRational[,] Gauged(SpinSpace s, int q)
    {
        var iz = s.GaugedIz();
        var half = new BigRational(1, 2);
        return q switch
        {
            1 => M.Scale(M.Add(M.Multiply(iz, s.GaugedIplus()), M.Multiply(s.GaugedIplus(), iz)), -half),
            -1 => M.Scale(M.Add(M.Multiply(iz, s.GaugedIminus()), M.Multiply(s.GaugedIminus(), iz)), half),
            2 => M.Scale(M.Multiply(s.GaugedIplus(), s.GaugedIplus()), half),
            -2 => M.Scale(M.Multiply(s.GaugedIminus(), s.GaugedIminus()), half),
            _ => throw new ArgumentOutOfRangeException(nameof(q), q, "q = ±1, ±2 (T₂₀ is Root6T20/√6)."),
        };
    }

    public static Matrix<Complex> Physical(SpinSpace s, int q)
    {
        var iz = s.Iz();
        var ip = s.Iplus();
        var im = s.Iminus();
        return q switch
        {
            0 => (3 * iz * iz - SpinSpace.ToDouble(s.SpinSquared) * Matrix<Complex>.Build.DenseIdentity(s.Levels)) / Math.Sqrt(6),
            1 => -(iz * ip + ip * iz) / 2,
            -1 => (iz * im + im * iz) / 2,
            2 => ip * ip / 2,
            -2 => im * im / 2,
            _ => throw new ArgumentOutOfRangeException(nameof(q)),
        };
    }

    /// <summary>(−1)^q.</summary>
    public static int Parity(int q) => (q % 2 == 0) ? 1 : -1;
}
