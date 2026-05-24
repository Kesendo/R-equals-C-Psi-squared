using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Stage 3b extension of <see cref="Pi2KleinSpectralView"/> from the
/// 4-cell Klein-Vierergruppe (bit_a, bit_b) to the full 8-cell Z₂³ decomposition
/// (bit_a, bit_b, y_par), as typed by <see cref="KleinEightCellClaim"/>.
///
/// <para>For each right eigenmode v of L with eigenvalue λ, this computes the
/// normalised mass v carries in each of the 8 cells. Unlike the 4-cell view,
/// where (Π²_Z, Π²_X) commute with L exactly (F61 + F63) and the mass
/// concentrates in a single cell per eigenvalue, the y_par axis is NOT a
/// conserved quantity of L (Pauli-string commutators don't preserve #Y mod 2 in
/// general). The 8-cell view is therefore a passive diagnostic: it shows how each
/// eigenmode distributes across the y_par axis, but the masses can span multiple
/// y_par cells within the same (bit_a, bit_b) Klein cell.</para>
///
/// <para>Use to answer: for a Hamiltonian H built from pairs in a specific
/// (Klein, y_par) cell (per F103/F105/F106 enumerations), do the slow modes
/// (Re(λ) ≈ 0) also concentrate in that y_par cell, or do they spread? The
/// diagnostic surfaces y_par-axis dynamical content that the 4-cell view
/// averages over.</para>
/// </summary>
public static class Pi2KleinEightSpectralView
{
    /// <summary>Compute the per-eigenmode 8-cell masses for a given Liouvillian
    /// in vec form. Returns one <see cref="KleinEightSpectralMode"/> per eigenmode
    /// (4^N total), unsorted.</summary>
    public static IReadOnlyList<KleinEightSpectralMode> ComputeFor(ComplexMatrix LVec, int N)
    {
        long d2 = 1L << (2 * N);
        if (LVec.RowCount != d2 || LVec.ColumnCount != d2)
            throw new ArgumentException(
                $"expected {d2}×{d2} Liouvillian for N={N}; got {LVec.RowCount}×{LVec.ColumnCount}");

        var T = PauliBasis.VecToPauliBasisTransform(N);
        var Tdag = T.ConjugateTranspose();
        double inv = 1.0 / (1 << N);

        var evd = LVec.Evd();
        var modes = new List<KleinEightSpectralMode>((int)d2);

        for (int i = 0; i < (int)d2; i++)
        {
            var lambda = evd.EigenValues[i];
            var vVec = evd.EigenVectors.Column(i);
            var vPauli = (Tdag * vVec) * inv;

            // 8 cells indexed by (bit_a, bit_b, y_par) in {0,1}^3.
            // Index encoding: cell = bit_a * 4 + bit_b * 2 + y_par.
            var cellMass = new double[8];
            for (long k = 0; k < d2; k++)
            {
                double mag = vPauli[(int)k].Magnitude;
                double m = mag * mag;
                if (m < 1e-18) continue;
                var letters = PauliIndex.FromFlat(k, N);
                int bitB = PiOperator.SquaredEigenvalue(letters, PauliLetter.Z) == +1 ? 0 : 1;
                int bitA = PiOperator.SquaredEigenvalue(letters, PauliLetter.X) == +1 ? 0 : 1;
                int yPar = CountY(letters) & 1;
                cellMass[bitA * 4 + bitB * 2 + yPar] += m;
            }

            // Normalise so the eight cell masses sum to 1 (eigenmode's relative
            // distribution across cells; absolute norm depends on MathNet's Evd
            // convention for non-Hermitian L).
            double totalRaw = 0;
            for (int c = 0; c < 8; c++) totalRaw += cellMass[c];
            if (totalRaw > 1e-15)
                for (int c = 0; c < 8; c++) cellMass[c] /= totalRaw;

            modes.Add(new KleinEightSpectralMode(
                lambda,
                Mass000: cellMass[0], Mass001: cellMass[1],
                Mass010: cellMass[2], Mass011: cellMass[3],
                Mass100: cellMass[4], Mass101: cellMass[5],
                Mass110: cellMass[6], Mass111: cellMass[7]));
        }
        return modes;
    }

    private static int CountY(IReadOnlyList<PauliLetter> letters)
    {
        int n = 0;
        foreach (var l in letters) if (l == PauliLetter.Y) n++;
        return n;
    }
}

/// <summary>One Liouvillian eigenmode read through the 8-cell Z₂³ lens. Mass
/// fields are |coefficient|² sums over Pauli strings in each cell, indexed by
/// (bit_a, bit_b, y_par) in {0,1}³. For a normalized eigenvector the 8 masses
/// sum to ~1.
///
/// <para>Naming convention: <c>Mass{bit_a}{bit_b}{y_par}</c>. E.g.,
/// <c>Mass001</c> = mass in (bit_a=0, bit_b=0, y_par=1) = Mother / odd-Y cell
/// per <see cref="KleinEightCellClaim"/>'s cell labels.</para></summary>
public sealed record KleinEightSpectralMode(
    Complex Eigenvalue,
    double Mass000, double Mass001,
    double Mass010, double Mass011,
    double Mass100, double Mass101,
    double Mass110, double Mass111)
{
    public double TotalMass =>
        Mass000 + Mass001 + Mass010 + Mass011 +
        Mass100 + Mass101 + Mass110 + Mass111;

    /// <summary>The 4-cell Klein mass at (bit_a, bit_b), summing both y_par
    /// values. Recovers the 4-cell <see cref="Pi2KleinSpectralView"/> mass for
    /// this eigenmode.</summary>
    public double KleinMass(int bitA, int bitB)
    {
        if (bitA != 0 && bitA != 1) throw new ArgumentOutOfRangeException(nameof(bitA));
        if (bitB != 0 && bitB != 1) throw new ArgumentOutOfRangeException(nameof(bitB));
        int idx0 = bitA * 4 + bitB * 2 + 0;
        int idx1 = bitA * 4 + bitB * 2 + 1;
        return CellMass(idx0) + CellMass(idx1);
    }

    /// <summary>Total y_par=0 mass (across all 4 Klein cells).</summary>
    public double YParityZeroMass => Mass000 + Mass010 + Mass100 + Mass110;

    /// <summary>Total y_par=1 mass (across all 4 Klein cells).</summary>
    public double YParityOneMass => Mass001 + Mass011 + Mass101 + Mass111;

    /// <summary>The y_par axis split inferred from the mass concentration. +1 when
    /// y_par=0 dominates (&gt; 99% of mass), −1 when y_par=1 dominates, 0 if mass
    /// is split (which is expected and not pathological since y_par is not
    /// L-conserved).</summary>
    public int YParityDominant
    {
        get
        {
            double y0 = YParityZeroMass;
            double y1 = YParityOneMass;
            if (y0 > 0.99 * TotalMass) return +1;
            if (y1 > 0.99 * TotalMass) return -1;
            return 0;
        }
    }

    /// <summary>The 8-cell with the largest mass share, as a (bit_a, bit_b, y_par) triple.</summary>
    public (int BitA, int BitB, int YPar) DominantCell
    {
        get
        {
            int bestIdx = 0;
            double bestMass = Mass000;
            for (int c = 1; c < 8; c++)
            {
                double m = CellMass(c);
                if (m > bestMass) { bestMass = m; bestIdx = c; }
            }
            return (bestIdx / 4, (bestIdx / 2) & 1, bestIdx & 1);
        }
    }

    /// <summary>Per-cell mass label per <see cref="KleinEightCellClaim"/>'s naming
    /// convention. Useful for human-readable diagnostic output.</summary>
    public string DominantCellLabel
    {
        get
        {
            var (a, b, y) = DominantCell;
            return (a, b, y) switch
            {
                (0, 0, 0) => "Mother/no-Y",
                (0, 0, 1) => "Mother/odd-Y",
                (0, 1, 0) => "Z-Klein/no-Y",
                (0, 1, 1) => "Z-Klein/odd-Y",
                (1, 0, 0) => "X-Klein/no-Y",
                (1, 0, 1) => "X-Klein/odd-Y",
                (1, 1, 0) => "Y-Klein/no-Y (paradox)",
                (1, 1, 1) => "Y-Klein/odd-Y (canonical)",
                _ => throw new InvalidOperationException($"unreachable: ({a},{b},{y})"),
            };
        }
    }

    private double CellMass(int idx) => idx switch
    {
        0 => Mass000, 1 => Mass001, 2 => Mass010, 3 => Mass011,
        4 => Mass100, 5 => Mass101, 6 => Mass110, 7 => Mass111,
        _ => throw new ArgumentOutOfRangeException(nameof(idx)),
    };
}
