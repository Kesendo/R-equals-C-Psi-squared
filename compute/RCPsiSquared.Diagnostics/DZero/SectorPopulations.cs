using System.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.DZero;

/// <summary>Sector populations p_n = Σ_{|k|=n} ρ_kk for n ∈ {0, 1, …, N}.
///
/// <para>The natural observables on the d=0 axis. For uniform XY/Heisenberg + Z-dephasing,
/// p_n = Tr(P_n · ρ) where P_n is the n-excitation sector projector — conserved by H and
/// by Z-dephasing per <b>F4</b> in <c>docs/ANALYTICAL_FORMULAS.md</c>. Drift in p_n(t)
/// is a fingerprint of non-{Z, XY/Heisenberg} noise: T1 σ⁻ shifts mass to lower n,
/// transverse fields h_x·X, h_y·Y mix sectors with Δn = ±1.</para>
///
/// <para>Hardware-trivial measurement: Z-basis tomography, group shots by Hamming weight.
/// Used in the Marrakesh d_zero retrospective
/// (<see cref="Core.Confirmations.ConfirmationsRegistry"/> entry
/// "d_zero_sector_trichotomy_marrakesh") to discriminate F87 truly/soft/hard from ⟨n⟩ alone
/// at 1/9 the Pauli-measurement cost of full ⟨X₀Z₂⟩ tomography.</para>
/// </summary>
public sealed record SectorPopulationResult(
    double[] P,
    double MeanN,
    double VarN,
    double Entropy,
    int N);

public static class SectorPopulations
{
    public static SectorPopulationResult FromDensityMatrix(ComplexMatrix rho)
    {
        int d = rho.RowCount;
        int N = (int)Math.Round(Math.Log2(d));
        if ((1 << N) != d) throw new ArgumentException($"rho dim {d} not 2^N");
        var p = new double[N + 1];
        for (int k = 0; k < d; k++)
        {
            int n_k = System.Numerics.BitOperations.PopCount((uint)k);
            p[n_k] += rho[k, k].Real;
        }
        return Aggregate(p, N);
    }

    public static SectorPopulationResult FromStateVector(ComplexVector psi)
    {
        int d = psi.Count;
        int N = (int)Math.Round(Math.Log2(d));
        if ((1 << N) != d) throw new ArgumentException($"psi length {d} not 2^N");
        var p = new double[N + 1];
        for (int k = 0; k < d; k++)
        {
            int n_k = System.Numerics.BitOperations.PopCount((uint)k);
            p[n_k] += psi[k].Magnitude * psi[k].Magnitude;
        }
        return Aggregate(p, N);
    }

    private static SectorPopulationResult Aggregate(double[] p, int N)
    {
        double mean = 0;
        for (int n = 0; n <= N; n++) mean += n * p[n];
        double var = 0;
        for (int n = 0; n <= N; n++) var += (n - mean) * (n - mean) * p[n];
        double entropy = 0;
        foreach (var pn in p)
            if (pn > 0) entropy -= pn * Math.Log(pn);
        return new SectorPopulationResult(p, mean, var, entropy, N);
    }
}
