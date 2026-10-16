using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Numerics;

/// <summary>Phase rigidity of each eigenmode of a (generally non-normal) complex matrix L. For right
/// eigenvectors |R_i⟩ (L|R_i⟩ = λ_i|R_i⟩) and left eigenvectors ⟨L_i| (⟨L_i|L = λ_i⟨L_i|),
/// r_i = |⟨L_i|R_i⟩| / (‖L_i‖·‖R_i‖) with the Hermitian inner product. r_i = 1 for a mode of a normal matrix.
/// The rigidity tends to zero at an exceptional point (the coalescing left/right eigenvectors become orthogonal,
/// the Petermann factor 1/r_i² diverges). Computed in Petermann form r_i = 1/(‖R⁻¹_row_i‖·‖R_col_i‖)
/// from a single Evd(L) (R = eigenvector matrix): near an EP the matching R⁻¹ row diverges so r → 0.
/// For a simple isolated eigenvalue, 0 < r_i ≤ 1 is the usual per-mode Petermann rigidity.
/// This construction avoids the separate-eigensolver left/right matching failure that can fake r → 0.
/// At an exact repeated semisimple eigenvalue, however, the individual values depend on the eigenbasis chosen inside its invariant subspace; use <see cref="ComputeClusters"/> (the spectral-projector reading) or a Jordan diagnostic there.</summary>
public static class PhaseRigidity
{
    /// <summary>One eigenmode: its eigenvalue, its phase rigidity, and its right eigenvector.</summary>
    public readonly record struct Mode(Complex Lambda, double Rigidity, ComplexVector Right);

    /// <summary>Per-eigenvalue phase rigidity of L, returned in L's Evd eigenvalue order.</summary>
    public static IReadOnlyList<Mode> Compute(Matrix<Complex> L)
    {
        var evd = L.Evd();
        var lam = evd.EigenValues;
        var r = evd.EigenVectors;     // columns = right eigenvectors
        var rigidities = RigiditiesFromRightEigenvectors(r);
        int n = lam.Count;

        var modes = new List<Mode>(n);
        for (int i = 0; i < n; i++)
        {
            var right = r.Column(i);
            modes.Add(new Mode(lam[i], rigidities[i], right));
        }
        return modes;
    }

    /// <summary>One eigenvalue cluster: its eigenvalue (the first member's), how many eigenvalues lie within
    /// the clustering tolerance of it, and the basis-free rigidity 1/‖P‖₂ of its spectral projector.</summary>
    public readonly record struct Cluster(Complex Lambda, int Size, double Rigidity);

    /// <summary>Phase rigidity per eigenvalue cluster, free of the eigenbasis. Eigenvalues within
    /// <paramref name="tolerance"/> of a cluster's first member join it; its spectral projector is
    /// P = R[:, c]·R⁻¹[c, :] and its rigidity is 1/‖P‖₂. For a simple eigenvalue ‖P‖₂ = ‖R_col‖·‖R⁻¹_row‖, so
    /// the value equals the per-mode rigidity of <see cref="Compute"/>; at an exact semisimple degeneracy,
    /// where the per-mode values depend on the basis the eigensolver returned, P and hence 1/‖P‖₂ do not.
    /// Clusters are returned in L's Evd eigenvalue order of their first members.</summary>
    public static IReadOnlyList<Cluster> ComputeClusters(Matrix<Complex> L, double tolerance)
    {
        var evd = L.Evd();
        var lam = evd.EigenValues;
        var r = evd.EigenVectors;
        var inverse = r.Inverse();
        int n = lam.Count;
        var assigned = new bool[n];
        var clusters = new List<Cluster>();
        for (int i = 0; i < n; i++)
        {
            if (assigned[i]) continue;
            var members = new List<int>();
            for (int j = i; j < n; j++)
            {
                if (!assigned[j] && (lam[j] - lam[i]).Magnitude < tolerance)
                {
                    assigned[j] = true;
                    members.Add(j);
                }
            }
            var columns = Matrix<Complex>.Build.DenseOfColumnVectors(members.Select(r.Column));
            var rows = Matrix<Complex>.Build.DenseOfRowVectors(members.Select(inverse.Row));
            // ‖P‖₂² = λ_max(P†P) = λ_max((C†C)(RR†)) for P = C·R: a k×k problem (k = cluster size)
            // in place of an SVD of the full n×n projector.
            var gram = columns.ConjugateTranspose() * columns * (rows * rows.ConjugateTranspose());
            double normSquared = gram.Evd().EigenValues.Max(z => z.Real);
            // Guard: if the cluster's eigenvectors are numerically dependent (an exceptional pair closer than the
            // tolerance, |Q − Q*| ≲ 1e-16 on a parameter grid), R⁻¹ is not trustworthy and the norm can come
            // out non-finite; report r = 0 then rather than NaN. The value there is rounding, not a reading.
            double rigidity = double.IsFinite(normSquared) && normSquared > 0 ? 1.0 / Math.Sqrt(normSquared) : 0.0;
            clusters.Add(new Cluster(lam[i], members.Count, rigidity));
        }
        return clusters;
    }

    /// <summary>Petermann rigidities for one chosen right-eigenvector matrix. At an exact repeated
    /// eigenvalue these are properties of that chosen eigenbasis, not invariants of each mode.</summary>
    internal static double[] RigiditiesFromRightEigenvectors(Matrix<Complex> rightEigenvectors)
    {
        var inverse = rightEigenvectors.Inverse(); // rows = biorthogonal left duals
        int n = rightEigenvectors.ColumnCount;
        var rigidities = new double[n];
        for (int i = 0; i < n; i++)
        {
            rigidities[i] = 1.0 /
                (inverse.Row(i).L2Norm() * rightEigenvectors.Column(i).L2Norm());
        }
        return rigidities;
    }
}
