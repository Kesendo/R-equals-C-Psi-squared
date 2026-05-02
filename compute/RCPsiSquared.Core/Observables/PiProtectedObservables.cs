using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Observables;

/// <summary>One Pauli-string observable's classification: Π-protected (always zero under L)
/// or active (nonzero, with a dominant eigenvalue).</summary>
public sealed record PiProtectedEntry(
    long PauliIndex,
    string PauliLabel,
    double MaxClusterContribution,
    Complex? DominantEigenvalue);

/// <summary>Aggregate of <see cref="PiProtectedObservables.Compute"/> output: Π-protected vs
/// active Pauli observables, plus the underlying L-eigenvalues and cluster count.</summary>
public sealed record PiProtectedResult(
    IReadOnlyList<PiProtectedEntry> Protected,
    IReadOnlyList<PiProtectedEntry> Active,
    IReadOnlyList<Complex> Eigenvalues,
    int NumClusters);

/// <summary>Identify Pauli-string observables σ_α whose expectation ⟨σ_α(t)⟩ stays exactly zero
/// for all t under L = −i[H, ·] + Σ_l γ_l (Z_l ρ Z_l − ρ).
///
/// Under the spectral decomposition L_pauli = V · diag(λ) · V⁻¹, the time-evolved expectation
/// is a sum of exponentials with rates λ_k (eigenvalues of L_pauli):
///
///   ⟨σ_α(t)⟩ = 2^N · Σ_λ S_λ(α) · exp(λ t)
///
/// where S_λ(α) = Σ_{k: λ_k = λ} V[α, k] · c[k] sums right-eigenvector components within each
/// degenerate cluster. σ_α is Π-protected iff S_λ(α) = 0 for every cluster — strictly weaker
/// than "each V[α,k]·c[k] vanishes" (degenerate-cluster cancellations are real, e.g.
/// ⟨X_0 I Z_2⟩ = 0 for Heisenberg chain on |+−+⟩ via SU(2)-multiplet cancellation).
///
/// Hardware-confirmed at Marrakesh (April 26, 2026) on YZ+ZY soft Hamiltonian, EQ-030.
/// See <see cref="Confirmations.ConfirmationsRegistry"/> entry "pi_protected_xiz_yzzy".
/// </summary>
public static class PiProtectedObservables
{
    public static PiProtectedResult Compute(ComplexMatrix H, IReadOnlyList<double> gammaPerSite,
        ComplexMatrix rho0, int N, double threshold = 1e-9, double clusterTolerance = 1e-8)
    {
        var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);
        var transform = PauliBasis.VecToPauliBasisTransform(N);
        var lPauli = (transform.ConjugateTranspose() * L * transform) / Math.Pow(2, N);

        var evd = lPauli.Evd();
        var evals = evd.EigenValues;
        var V = evd.EigenVectors;
        var Vinv = V.Inverse();

        var rhoPauli = PauliBasis.ToPauliVector(rho0, N);
        var c = Vinv * rhoPauli;

        long d2 = 1L << (2 * N);
        var evalsArr = evals.ToArray();
        var clusters = ClusterEigenvalues(evalsArr, clusterTolerance);
        var protectedList = new List<PiProtectedEntry>();
        var activeList = new List<PiProtectedEntry>();

        for (long alpha = 1; alpha < d2; alpha++)
        {
            double maxS = 0;
            int dominantIdx = 0;
            Complex dominantS = Complex.Zero;
            foreach (var cluster in clusters)
            {
                Complex S = Complex.Zero;
                foreach (int k in cluster) S += V[(int)alpha, k] * c[k];
                double absS = S.Magnitude;
                if (absS > maxS) { maxS = absS; dominantIdx = cluster[0]; dominantS = S; }
            }

            string label = PauliLabel.Format(PauliIndex.FromFlat(alpha, N));
            if (maxS < threshold)
                protectedList.Add(new PiProtectedEntry(alpha, label, maxS, null));
            else
                activeList.Add(new PiProtectedEntry(alpha, label, maxS, evals[dominantIdx]));
        }

        return new PiProtectedResult(protectedList, activeList, evals.ToArray(), clusters.Count);
    }

    private static List<List<int>> ClusterEigenvalues(IReadOnlyList<Complex> evals, double tolerance)
    {
        int n = evals.Count;
        var used = new bool[n];
        var clusters = new List<List<int>>();
        for (int i = 0; i < n; i++)
        {
            if (used[i]) continue;
            var cluster = new List<int> { i };
            used[i] = true;
            for (int j = i + 1; j < n; j++)
            {
                if (!used[j] && (evals[j] - evals[i]).Magnitude < tolerance)
                {
                    cluster.Add(j);
                    used[j] = true;
                }
            }
            clusters.Add(cluster);
        }
        return clusters;
    }

}
