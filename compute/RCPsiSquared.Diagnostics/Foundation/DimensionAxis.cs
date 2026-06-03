using System.Numerics;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>A one-parameter dimension: a θ-grid plus the Hamiltonian H(θ) whose Liouvillian
/// the Object Manager sweeps. Each θ is a single position in the dimension's in-between;
/// sweeping θ reads the marks (the Liouvillian eigenvalues, the fixed contract) against the
/// in-between (how the slow eigen-subspace rotates, the content). See
/// <see cref="DimensionSweep"/> for the read.
///
/// <para>The first and verified axis is <see cref="Crossover"/>: the bond-angle crossover
/// H(θ) = Σ_bonds (cos θ · X_iZ_{i+1} + sin θ · Y_iZ_{i+1}). Its Liouvillian is a similarity
/// transform of L(0), L(θ) = Ad_{R_z(θ)}·L(0)·Ad_{R_z(θ)}⁻¹, so the eigenvalue multiset is
/// exactly invariant across θ (the marks do not move) while the slow subspace rotates with
/// R_z(θ) (the in-between). θ = 0 is pure XZ, θ = π/4 is the T-gate symmetric crossover, and
/// θ = π/2 is pure YZ.</para>
/// </summary>
public sealed record DimensionAxis(
    string Name,
    IReadOnlyList<double> Theta,
    int N,
    IReadOnlyList<double> GammaPerSite,
    Func<double, ComplexMatrix> Hamiltonian,
    IReadOnlyList<int> LitSites)
{
    /// <summary>The verified crossover bond-angle axis. θ-grid is linearly spaced on [0, π/2]
    /// with <paramref name="thetaPoints"/> points (so θ_k = (π/2)·k/(thetaPoints−1)); the
    /// dephasing is uniform Z at rate <paramref name="gamma"/> on every one of the
    /// <paramref name="N"/> sites; the Hamiltonian sums the chain bonds
    /// cos θ · X_iZ_{i+1} + sin θ · Y_iZ_{i+1} for i = 0 .. N−2. The lit sites (the X/Y carriers
    /// the turn rotates, the first site of every bond) are {0 .. N−2}; the last site is pure
    /// shadow Z.</summary>
    public static DimensionAxis Crossover(int N, double gamma, int thetaPoints = 25)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"crossover needs at least one bond; got N={N}");
        if (thetaPoints < 2) throw new ArgumentOutOfRangeException(nameof(thetaPoints), $"need at least two θ points; got {thetaPoints}");

        var theta = LinSpace(0.0, Math.PI / 2.0, thetaPoints);
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();

        ComplexMatrix Hamiltonian(double t)
        {
            double c = Math.Cos(t);
            double s = Math.Sin(t);
            ComplexMatrix? H = null;
            for (int i = 0; i < N - 1; i++)
            {
                var xz = PauliString.SiteOp(N, i, PauliLetter.X) * PauliString.SiteOp(N, i + 1, PauliLetter.Z);
                var yz = PauliString.SiteOp(N, i, PauliLetter.Y) * PauliString.SiteOp(N, i + 1, PauliLetter.Z);
                var bond = (Complex)c * xz + (Complex)s * yz;
                H = H is null ? bond : H + bond;
            }
            return H!; // N ≥ 2 guarantees at least one bond
        }

        return new DimensionAxis(
            Name: "crossover",
            Theta: theta,
            N: N,
            GammaPerSite: gammaPerSite,
            Hamiltonian: Hamiltonian,
            LitSites: Enumerable.Range(0, N - 1).ToArray());
    }

    /// <summary>The PTF J-defect axis. The parameter is δJ (stored in <see cref="Theta"/>), swept on
    /// [0, <paramref name="deltaJMax"/>] with <paramref name="points"/> points; the Hamiltonian is the
    /// uniform XY chain H = Σ_b (J_b/2)(X_bX_{b+1} + Y_bY_{b+1}) with one bond detuned, J_{defectBond}
    /// = 1 + δJ. Dephasing is uniform Z at rate <paramref name="gamma"/> on all N sites.
    ///
    /// <para>This is the contrast to <see cref="Crossover"/>. The crossover is an EXACT similarity
    /// (the eigenvalues are frozen, the in-between is a pure rotation). The J-defect is Π-invariant
    /// (Π H Π⁻¹ = −H), so it keeps the spectrum palindromic (the contract holds, palindrome residual
    /// ~10⁻¹⁵) but it is NOT a similarity: the spectrum genuinely moves, the slow-mode eigenvalues are
    /// protected only to first order (shift O(δJ²)), and the slow eigenVECTORS mix. The in-between
    /// here carries more than a rotation, the eigenvector mixing PTF computed. LitSites is empty: the
    /// defect is not a per-site rotation, so the crossover's lit/shadow split does not apply.</para>
    ///
    /// <para>The defect bond is a chain bond index b ∈ [0, N−2] (the bond on sites b, b+1); the
    /// canonical PTF defect is bond 0 (sites 0, 1).</para></summary>
    public static DimensionAxis JDefect(int N, double gamma, int defectBond,
        double deltaJMax = 0.1, int points = 25)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), $"J-defect needs at least one bond; got N={N}");
        if (defectBond < 0 || defectBond >= N - 1)
            throw new ArgumentOutOfRangeException(nameof(defectBond), $"bond index must be in [0, {N - 2}]; got {defectBond}");
        if (points < 2) throw new ArgumentOutOfRangeException(nameof(points), $"need at least two δJ points; got {points}");

        var theta = LinSpace(0.0, deltaJMax, points); // the parameter here is δJ, not an angle
        var gammaPerSite = Enumerable.Repeat(gamma, N).ToArray();

        ComplexMatrix Hamiltonian(double deltaJ)
        {
            var bondJ = new double[N - 1];
            for (int b = 0; b < N - 1; b++) bondJ[b] = 1.0;
            bondJ[defectBond] = 1.0 + deltaJ;
            return PauliHamiltonian.XYChain(N, bondJ).ToMatrix();
        }

        return new DimensionAxis(
            Name: "jdefect",
            Theta: theta,
            N: N,
            GammaPerSite: gammaPerSite,
            Hamiltonian: Hamiltonian,
            LitSites: Array.Empty<int>());
    }

    /// <summary>Linearly spaced grid of <paramref name="count"/> points on [<paramref name="lo"/>,
    /// <paramref name="hi"/>], endpoints included. The k-th point is lo + (hi−lo)·k/(count−1),
    /// so the last point is exactly <paramref name="hi"/> (no float drift at the endpoints).</summary>
    private static double[] LinSpace(double lo, double hi, int count)
    {
        var grid = new double[count];
        for (int k = 0; k < count; k++)
            grid[k] = lo + (hi - lo) * k / (count - 1);
        return grid;
    }
}
