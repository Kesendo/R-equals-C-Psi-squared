using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The bridge between Schicht 1+2 (algebraic Klein-cell decomposition of Pauli
/// operator space) and the dynamical spectrum of the Liouvillian. For each right
/// eigenmode v of L with eigenvalue λ, this computes the **normalised** mass that v
/// carries in each of the 4 Klein-cells (Pp, Pm, Mp, Mm) — the four masses sum to 1.
///
/// <para>Both Π²-axes commute with L for Heisenberg + Z-dephasing:
/// <list type="bullet">
///   <item><b>Π²_Z:</b> [L, Π²_Z] = 0 via F63 (w_YZ parity conservation; consequence of
///         F1's palindrome identity Π·L·Π⁻¹ = −L − 2σ·I squared).</item>
///   <item><b>Π²_X:</b> [L, Π²_X] = 0 via F61 (n_XY parity conservation; the orthogonal
///         axis bit_a-parity).</item>
/// </list>
/// Together these make the (Π²_Z, Π²_X) Klein-cell membership a conserved quantity, so
/// the L-eigenspace for any non-degenerate eigenvalue lives entirely in one Klein cell.
/// For degenerate eigenvalues spanning both subspaces of either axis, a generic basis
/// (such as the one MathNet's Evd returns) can mix across them — the *eigenspace*
/// respects (Π²_Z, Π²_X), but individual basis vectors may not. The Klein-cell masses
/// below describe the chosen basis vector's distribution; aggregated over a degenerate
/// subspace they recover the conserved Klein structure.</para>
///
/// <para>Use this to read which Klein-cells the dynamics populates for a given Hamiltonian
/// across the slow / fast spectrum: do slow modes (Re(λ) ≈ 0) populate the same Klein-cells
/// as the Hamiltonian's bilinears, or do they spread differently? The probe shows where
/// Π² breaking surfaces dynamically.</para>
/// </summary>
public static class Pi2KleinSpectralView
{
    /// <summary>Compute the per-eigenmode Klein-cell masses for a given Liouvillian in vec form.
    /// Returns one <see cref="KleinSpectralMode"/> per eigenmode (4^N total), unsorted.</summary>
    public static IReadOnlyList<KleinSpectralMode> ComputeFor(ComplexMatrix LVec, int N)
    {
        long d2 = 1L << (2 * N);
        if (LVec.RowCount != d2 || LVec.ColumnCount != d2)
            throw new ArgumentException(
                $"expected {d2}×{d2} Liouvillian for N={N}; got {LVec.RowCount}×{LVec.ColumnCount}");

        var T = PauliBasis.VecToPauliBasisTransform(N);
        var Tdag = T.ConjugateTranspose();
        double inv = 1.0 / (1 << N);

        var evd = LVec.Evd();
        var modes = new List<KleinSpectralMode>((int)d2);

        for (int i = 0; i < (int)d2; i++)
        {
            var lambda = evd.EigenValues[i];
            var vVec = evd.EigenVectors.Column(i);
            var vPauli = (Tdag * vVec) * inv;

            double pp = 0, pm = 0, mp = 0, mm = 0;
            for (long k = 0; k < d2; k++)
            {
                double mag = vPauli[(int)k].Magnitude;
                double m = mag * mag;
                if (m < 1e-18) continue;
                var letters = PauliIndex.FromFlat(k, N);
                int eigZ = PiOperator.SquaredEigenvalue(letters, PauliLetter.Z);
                int eigX = PiOperator.SquaredEigenvalue(letters, PauliLetter.X);
                if (eigZ == +1 && eigX == +1) pp += m;
                else if (eigZ == +1 && eigX == -1) pm += m;
                else if (eigZ == -1 && eigX == +1) mp += m;
                else mm += m;
            }
            // Normalise so the four cell masses sum to 1 — the eigenvector's relative
            // distribution across cells is what's structurally interesting; absolute norm
            // depends on MathNet's Evd convention for non-Hermitian L.
            double totalRaw = pp + pm + mp + mm;
            if (totalRaw > 1e-15)
            {
                pp /= totalRaw; pm /= totalRaw; mp /= totalRaw; mm /= totalRaw;
            }
            modes.Add(new KleinSpectralMode(lambda, pp, pm, mp, mm));
        }
        return modes;
    }
}

/// <summary>One Liouvillian eigenmode read through the Klein lens. <c>Eigenvalue</c> is the
/// complex spectral value; the four mass fields are |coefficient|² sums over the Pauli
/// strings in each Klein-cell. For a normalized eigenvector the four masses sum to ~1.</summary>
public sealed record KleinSpectralMode(
    Complex Eigenvalue,
    double MassPp, double MassPm, double MassMp, double MassMm)
{
    public double TotalMass => MassPp + MassPm + MassMp + MassMm;

    /// <summary>Π²_Z parity inferred from the mass concentration. +1 when (Pp + Pm) carries
    /// the mass, −1 when (Mp + Mm), 0 if mass is split between Π²_Z eigenspaces (which
    /// indicates a numerical issue or a degenerate eigenvalue).</summary>
    public int Pi2ZParity
    {
        get
        {
            double zPlus = MassPp + MassPm;
            double zMinus = MassMp + MassMm;
            if (zPlus > 0.99 * TotalMass) return +1;
            if (zMinus > 0.99 * TotalMass) return -1;
            return 0;
        }
    }

    /// <summary>The Klein-cell that holds the largest share of mass, as a label.</summary>
    public string DominantCell
    {
        get
        {
            double max = MassPp;
            string label = "Pp";
            if (MassPm > max) { max = MassPm; label = "Pm"; }
            if (MassMp > max) { max = MassMp; label = "Mp"; }
            if (MassMm > max) { max = MassMm; label = "Mm"; }
            return label;
        }
    }
}
