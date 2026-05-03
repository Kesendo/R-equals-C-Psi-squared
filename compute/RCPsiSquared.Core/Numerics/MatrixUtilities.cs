using System.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Element-wise scans on small complex matrices used by the F86 inspection wrappers.
/// MathNet's <c>InfinityNorm</c> measures max-row-sum, not max-element-magnitude — these
/// helpers fill the gap.
/// </summary>
public static class MatrixUtilities
{
    /// <summary>Maximum |M[i, j]| across all entries.</summary>
    public static double MaxElementMagnitude(ComplexMatrix m)
    {
        double max = 0;
        for (int i = 0; i < m.RowCount; i++)
            for (int j = 0; j < m.ColumnCount; j++)
            {
                double mag = m[i, j].Magnitude;
                if (mag > max) max = mag;
            }
        return max;
    }

    /// <summary>Maximum |M[i, j]| over off-diagonal entries (i ≠ j) — used to verify that
    /// a matrix is diagonal in a chosen basis (residual ~ 0 ⇔ diagonal).</summary>
    public static double MaxOffDiagonalMagnitude(ComplexMatrix m)
    {
        double max = 0;
        for (int i = 0; i < m.RowCount; i++)
            for (int j = 0; j < m.ColumnCount; j++)
            {
                if (i == j) continue;
                double mag = m[i, j].Magnitude;
                if (mag > max) max = mag;
            }
        return max;
    }
}
