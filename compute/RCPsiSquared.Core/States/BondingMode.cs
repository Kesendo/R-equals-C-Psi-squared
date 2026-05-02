using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.States;

/// <summary>F65 single-excitation bonding-mode states |ψ_k⟩.
///
/// ψ_k(j) = √(2/(N+1)) · sin(π·k·(j+1)/(N+1)) is the j-th site amplitude of the
/// k-th single-excitation eigenmode of the open uniform XX chain. The state vector
/// embeds these amplitudes at single-excitation kets |1_j⟩ — site j → bit position
/// N−1−j in the 2^N index (big-endian: site 0 = leftmost Kronecker factor).
///
/// HANDSHAKE_ALGEBRA.md: selecting a value of k IS the per-side handshake. Two callers
/// who independently construct <c>BondingMode.Build(N, k)</c> with the same k have
/// agreed on the receiver — no exchange step needed.
/// </summary>
public static class BondingMode
{
    /// <summary>F65 single-excitation bonding mode |ψ_k⟩ for k ∈ [1, N], length-2^N normalised.</summary>
    public static ComplexVector Build(int N, int k)
    {
        if (k < 1 || k > N) throw new ArgumentOutOfRangeException(nameof(k), $"k={k} outside [1, N={N}]");
        int d = 1 << N;
        var psi = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d);
        double norm = Math.Sqrt(2.0 / (N + 1));
        for (int j = 0; j < N; j++)
            psi[1 << (N - 1 - j)] = new Complex(norm * Math.Sin(Math.PI * k * (j + 1) / (N + 1)), 0);
        return psi;
    }

    /// <summary>Bonding-mode pair state (|vac⟩ + |ψ_k⟩) / √2 — the canonical PTF / handshake
    /// initial state. Single-qubit purity trajectories from this state are the painter's
    /// canvas in PERSPECTIVAL_TIME_FIELD.md.</summary>
    public static ComplexVector PairState(int N, int k)
    {
        if (k < 1 || k > N) throw new ArgumentOutOfRangeException(nameof(k), $"k={k} outside [1, N={N}]");
        int d = 1 << N;
        var psi = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(d);
        psi[0] = Complex.One;
        double norm = Math.Sqrt(2.0 / (N + 1));
        for (int j = 0; j < N; j++)
            psi[1 << (N - 1 - j)] += new Complex(norm * Math.Sin(Math.PI * k * (j + 1) / (N + 1)), 0);
        var n = (Complex)Math.Sqrt((psi.ConjugateDotProduct(psi)).Real);
        return psi / n;
    }
}
