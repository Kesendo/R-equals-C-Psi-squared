using System.Numerics;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>Per-qubit Bloch-vector reading: the direct ±0.5 pair structure of a density
/// matrix. For each qubit k, computes the Bloch vector
/// <c>r_k = (Tr(σ_x^k · ρ), Tr(σ_y^k · ρ), Tr(σ_z^k · ρ))</c>; the reduced single-qubit
/// density matrix ρ_k = Tr_{not k}(ρ) has eigenvalues 1/2 ± |r_k|/2.
///
/// <para>This is the most direct reading of the trio's ±0.5 anchor at the state level:
/// for any product state with each qubit in a Pauli-axis eigenstate, |r_k| = 1, so the
/// eigenvalues land exactly on {1, 0} = {1/2 + 0.5, 1/2 − 0.5}, the structural pair.
/// The dominant axis ('X', 'Y', 'Z') tells us which Bloch axis carries the ±0.5; the
/// dominant sign tells us which direction.</para>
///
/// <para>Companion to <see cref="MemoryAxisRho"/>: that diagnostic gives Frobenius-norm
/// partition (static / memory / Π²-odd-fraction-within-memory); this one gives the
/// signed-pair reading per qubit per axis. Same trio, different lens.</para>
/// </summary>
public sealed record QubitBlochReading(int QubitIndex, double Rx, double Ry, double Rz)
{
    /// <summary>‖r_k‖₂ = sqrt(rx² + ry² + rz²); equals 1 for any pure single-qubit
    /// reduced density matrix, &lt; 1 if the qubit is entangled with the rest.</summary>
    public double RMagnitude => Math.Sqrt(Rx * Rx + Ry * Ry + Rz * Rz);

    /// <summary>Eigenvalue deviation from 1/2: ρ_k has eigenvalues 1/2 ± EigenDeviation.
    /// Equals 0.5 for any pure single-qubit state aligned with a Pauli axis (the
    /// structural memory-axis pair); equals 0 at the maximally mixed state.</summary>
    public double EigenDeviation => RMagnitude / 2.0;

    /// <summary>Dominant Bloch axis: 'X', 'Y', 'Z', or 'I' (maximally mixed). Under the
    /// framework's standard Z-dephasing convention the axis also reports the qubit's
    /// Π²-class: 'X' is Π²-even, 'Y' and 'Z' are Π²-odd. The bit_a/bit_b roles shift
    /// under X-dephasing; see F88 in <c>docs/ANALYTICAL_FORMULAS.md</c>.</summary>
    public char DominantAxis
    {
        get
        {
            const double tol = 1e-9;
            if (RMagnitude < tol) return 'I';
            double absX = Math.Abs(Rx), absY = Math.Abs(Ry), absZ = Math.Abs(Rz);
            if (absX >= absY && absX >= absZ) return 'X';
            if (absY >= absZ) return 'Y';
            return 'Z';
        }
    }

    /// <summary>Sign of the dominant Bloch component (+1 or −1; 0 if maximally mixed).
    /// For pure-state product qubits: +1 / −1 corresponds to the eigenstate orientation
    /// (e.g., 'Y' + +1 = |+i⟩, 'Y' + −1 = |−i⟩, 'Z' + +1 = |0⟩, 'Z' + −1 = |1⟩).</summary>
    public int DominantSign => DominantAxis switch
    {
        'X' => Math.Sign(Rx),
        'Y' => Math.Sign(Ry),
        'Z' => Math.Sign(Rz),
        _ => 0,
    };
}

public sealed record BlochAxisReadingResult(int N, IReadOnlyList<QubitBlochReading> Qubits);

/// <summary>Compute the per-qubit Bloch reading of a density matrix. Each qubit's
/// reading is independent; this is the per-site projection, not a full multi-qubit
/// state characterisation. For product states, the per-qubit readings together
/// determine the state up to a global phase; for entangled states they capture
/// only the marginal information.</summary>
public static class BlochAxisReading
{
    public static BlochAxisReadingResult Compute(ComplexMatrix rho, int N)
    {
        int d = 1 << N;
        if (rho.RowCount != d || rho.ColumnCount != d)
            throw new ArgumentException(
                $"expected {d}×{d} matrix for N={N}; got {rho.RowCount}×{rho.ColumnCount}");

        var sitePaulis = PauliString.SitePaulis(N);
        var qubits = new List<QubitBlochReading>(N);
        for (int k = 0; k < N; k++)
        {
            var (sigmaX, sigmaY, sigmaZ) = sitePaulis[k];
            double rx = (sigmaX * rho).Trace().Real;
            double ry = (sigmaY * rho).Trace().Real;
            double rz = (sigmaZ * rho).Trace().Real;
            qubits.Add(new QubitBlochReading(k, rx, ry, rz));
        }
        return new BlochAxisReadingResult(N, qubits);
    }
}
