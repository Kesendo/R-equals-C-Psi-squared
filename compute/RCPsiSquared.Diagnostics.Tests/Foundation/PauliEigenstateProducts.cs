using System.Numerics;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>Tensor-product state-vector helpers for X/Y/Z-basis eigenstates, used by
/// the Foundation/ exploration and unit tests. Big-endian site convention matches
/// <c>RCPsiSquared.Core.States.PolarityState</c>.
///
/// <para>For 'X' axis: sign +1 = |+⟩, −1 = |−⟩. For 'Y': +1 = |+i⟩, −1 = |−i⟩.
/// For 'Z': +1 = |0⟩, −1 = |1⟩ (computational basis).</para>
///
/// <para>Test-internal because the test-data ergonomics use <c>char</c> axes
/// (compatible with xUnit <c>InlineData</c> arrays); promotion to a Core/States
/// primitive would prefer <see cref="RCPsiSquared.Core.Pauli.PauliLetter"/> typing.</para>
/// </summary>
internal static class PauliEigenstateProducts
{
    /// <summary>Y-basis tensor product: each site in <c>|+i⟩</c> (sign +1) or <c>|−i⟩</c> (sign −1).</summary>
    public static ComplexVector YBasis(int N, IReadOnlyList<int> signs)
    {
        int d = 1 << N;
        var vec = ComplexVector.Build.Dense(d);
        double norm = 1.0 / Math.Sqrt(d);
        for (int idx = 0; idx < d; idx++)
        {
            Complex amp = Complex.One;
            for (int k = 0; k < N; k++)
            {
                int bit = (idx >> (N - 1 - k)) & 1;
                if (bit == 1)
                    amp *= signs[k] == +1 ? Complex.ImaginaryOne : -Complex.ImaginaryOne;
            }
            vec[idx] = amp * norm;
        }
        return vec;
    }

    /// <summary>General X/Y/Z-basis tensor product. <paramref name="axes"/>[k] ∈ {'X','Y','Z'};
    /// <paramref name="signs"/>[k] ∈ {+1, −1}. Throws on unknown axis.</summary>
    public static ComplexVector General(int N, IReadOnlyList<char> axes, IReadOnlyList<int> signs)
    {
        int d = 1 << N;
        var vec = ComplexVector.Build.Dense(d);
        double sqrt2 = Math.Sqrt(2.0);
        for (int idx = 0; idx < d; idx++)
        {
            Complex amp = Complex.One;
            bool zero = false;
            for (int k = 0; k < N; k++)
            {
                int bit = (idx >> (N - 1 - k)) & 1;
                switch (axes[k])
                {
                    case 'X':
                        amp *= (bit == 0 ? 1.0 : (double)signs[k]) / sqrt2;
                        break;
                    case 'Y':
                        if (bit == 0) amp *= 1.0 / sqrt2;
                        else amp *= (signs[k] == +1 ? Complex.ImaginaryOne : -Complex.ImaginaryOne) / sqrt2;
                        break;
                    case 'Z':
                        if ((signs[k] == +1 && bit == 1) || (signs[k] == -1 && bit == 0))
                            zero = true;
                        break;
                    default:
                        throw new ArgumentException($"unknown axis '{axes[k]}'; expected X/Y/Z");
                }
                if (zero) break;
            }
            vec[idx] = zero ? Complex.Zero : amp;
        }
        return vec;
    }

    /// <summary>X/Y-only mix (subset of <see cref="General"/>): each site in σ_x or σ_y eigenstate.</summary>
    public static ComplexVector MixedXY(int N, IReadOnlyList<char> axes, IReadOnlyList<int> signs) =>
        General(N, axes, signs);
}
