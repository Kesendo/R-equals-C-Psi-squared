using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Lindblad;

/// <summary>Z-dephasing + T1 amplitude damping.
///
/// L(ρ) = −i[H, ρ] + Σ_l γ_l (Z_l ρ Z_l − ρ)
///                 + Σ_l γ^{T1}_l (σ⁻_l ρ σ⁺_l − ½ {σ⁺_l σ⁻_l, ρ})
///
/// σ⁻ = (X − iY)/2 = [[0, 1], [0, 0]] is the lowering operator (|1⟩→|0⟩). With γ^{T1}_l = 0
/// reduces to <see cref="PauliDephasingDissipator.BuildZ"/>. T1 introduces palindrome-breaking;
/// used to measure how Π-protected count shifts under amplitude damping.
/// </summary>
public static class T1Dissipator
{
    private static readonly ComplexMatrix _sigmaMinusSingleSite = Matrix<Complex>.Build.DenseOfArray(new Complex[,]
    {
        { Complex.Zero, Complex.One },
        { Complex.Zero, Complex.Zero },
    });

    public static ComplexMatrix Build(ComplexMatrix H, IReadOnlyList<double> gammaZPerSite, IReadOnlyList<double> gammaT1PerSite)
    {
        int d = H.RowCount;
        int N = (int)Math.Round(Math.Log2(d));
        if (gammaZPerSite.Count != N)
            throw new ArgumentException($"gammaZ list has length {gammaZPerSite.Count}, expected N={N}");
        if (gammaT1PerSite.Count != N)
            throw new ArgumentException($"gammaT1 list has length {gammaT1PerSite.Count}, expected N={N}");

        var cOps = new List<ComplexMatrix>();
        for (int l = 0; l < N; l++)
        {
            if (gammaZPerSite[l] != 0)
                cOps.Add((Complex)Math.Sqrt(gammaZPerSite[l]) * PauliString.SiteOp(N, l, PauliLetter.Z));
            if (gammaT1PerSite[l] != 0)
                cOps.Add((Complex)Math.Sqrt(gammaT1PerSite[l]) * SigmaMinusOnSite(N, l));
        }
        return LindbladianBuilder.Build(H, cOps);
    }

    private static ComplexMatrix SigmaMinusOnSite(int N, int site)
    {
        var I2 = Matrix<Complex>.Build.DenseIdentity(2);
        ComplexMatrix result = (site == 0) ? _sigmaMinusSingleSite : I2;
        for (int i = 1; i < N; i++)
            result = result.KroneckerProduct(i == site ? _sigmaMinusSingleSite : I2);
        return result;
    }
}
