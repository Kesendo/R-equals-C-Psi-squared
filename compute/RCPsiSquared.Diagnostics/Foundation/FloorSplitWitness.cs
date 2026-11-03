using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Probes;
using RCPsiSquared.Core.Resonance;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F86's K split by the free-fermion floor. In F86's c = 2 block, the (1,2) block of the
/// uniform open XY chain under uniform Z-dephasing γ₀, the Jordan-Wigner family u_j = P_1·c_j·P_2
/// (<see cref="FreeFermionFloorWitness"/>) spans a subspace W that every single bond maps onto itself
/// (checked here in integers, bond by bond) and on which the dephasing is the scalar −2γ₀. So W is
/// invariant under L, under L† and under every ∂L/∂J_b, independent of J, and its orthogonal
/// projector P commutes with the whole run: ρ(t) = Pρ(t) + (1−P)ρ(t), and the same for ∂ρ/∂J_b.
/// F86's observable K_b(Q, t) = 2·Re⟨ρ(t)|S_kernel|∂ρ/∂J_b⟩ is quadratic, so it splits into three
/// terms: floor² (Pρ against P∂ρ), rest² ((1−P)ρ against (1−P)∂ρ), and the interference between the
/// two, read at each bond's canonical Q_peak and envelope time t*.
///
/// <para>Routes. (1) The three terms add up to the production scanner's envelope at the same Q
/// (<see cref="ResonanceScan"/>), two eigen routes on one quantity. (2) The floor² term is read twice:
/// projected in the full block, and from the single-particle closed form, where the readout of site l
/// sees only u_l, through the integer s_l = Σ entries of u_l, and the floor evolves as
/// −2γ₀ + iJσ·A on the chain's adjacency A, whose modes are the sines sin(πk(j+1)/(N+1)) at
/// 2cos(πk/(N+1)); no block and no eigensolver enter that route. Its time dependence is
/// t·e^(−4γ₀t) times a function of J·t alone.</para>
///
/// <para>Readings, at the canonical peak of each bond: the shares of the three terms in K (signed, K
/// itself signed there), so the endpoint bond's floor share and the interference can be followed
/// across N; and the bonds into the middle site at odd N, where the floor term is silent at every J and t:
/// s_l = N−1−2l is odd under the chain reflection, the floor dynamics keep it odd so x_m(t) = 0, a bond
/// into m feeds only x_(m∓1) onto the even vector e_m so ∂x stays even, and Σ_l s_l²·x̄_l·∂x_l is a sum
/// of an odd function. The interference term is not silent there. Docs: <c>docs/STERILE_BIRTHCANAL_AND_THE_JUNCTION.md</c>
/// (the floor), <c>docs/proofs/PROOF_CHAIN_GAP_DOMINANCE.md</c> (the family), <c>docs/F86_VALUES_INVENTORY.md</c>.</para></summary>
public sealed class FloorSplitWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    public const int DefaultN = 6;
    public const double DefaultGammaZero = 0.05;
    public const int MinN = 4, MaxN = 8;

    public int N { get; }
    public double GammaZero { get; }
    public IReadOnlyList<FloorSplitReading> Bonds { get; }

    public FloorSplitWitness(int n = DefaultN, double gammaZero = DefaultGammaZero)
    {
        if (n < MinN || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n), n, $"N must be in [{MinN}, {MaxN}].");
        if (!(gammaZero > 0) || double.IsInfinity(gammaZero))
            throw new ArgumentOutOfRangeException(nameof(gammaZero), gammaZero, "γ₀ must be positive and finite.");
        MathNetSetup.EnsureInitialized();
        N = n;
        GammaZero = gammaZero;
        Bonds = Read();
    }

    /// <summary>Every bond closes the family, the three terms meet the scanner, and the two floor
    /// routes meet.</summary>
    public bool AllRoutesMeet
    {
        get
        {
            double scale = Bonds.Max(b => Math.Abs(b.K));
            return Bonds.All(b => b.FamilyClosedUnderThisBond
                && Math.Abs(Math.Abs(b.KFloor + b.KRest + b.KCross) - b.ScanEnvelope) <= 1e-9 * b.ScanEnvelope
                && Math.Abs(b.KFloor - b.KFloorClosedForm) <= 1e-10 * scale);
        }
    }

    private IReadOnlyList<FloorSplitReading> Read()
    {
        var block = new CoherenceBlock(N, 1, GammaZero);
        var basis = block.Basis;
        var decomp = block.Decomposition;
        int m = basis.MTotal, nb = N - 1;
        var family = FreeFermionFloorWitness.JordanWignerFamily(basis);
        int cells = family[0].Count;                       // C(N−1, 1), the same for every u_j
        var s = family.Select(u => (double)u.Values.Sum()).ToArray();

        // per-bond closure in integers: K_b = −i·M_b maps u_j onto Σ_j' T_b[j', j]·u_j'
        var closed = new bool[nb];
        var bondSign = new int[nb];
        var owner = new Dictionary<int, int>();
        for (int j = 0; j < N; j++) foreach (int c in family[j].Keys) owner[c] = j;
        for (int b = 0; b < nb; b++)
            (closed[b], bondSign[b]) = ClosureUnderBond(decomp.MhPerBond[b], family, owner, b);
        if (closed.All(x => x) && bondSign.Distinct().Count() != 1)
            throw new InvalidOperationException("the bonds act on the family with different signs; the closed form assumes one.");

        // orthonormal basis of W (disjoint supports, equal norms)
        var bw = ComplexMatrix.Build.Dense(m, N);
        double inv = 1.0 / Math.Sqrt(cells);
        for (int j = 0; j < N; j++) foreach (var (c, v) in family[j]) bw[c, j] = v * inv;
        var bwH = bw.ConjugateTranspose();

        var probe = DickeBlockProbe.Build(block);
        var sKernel = SpatialSumKernel.Build(block);

        // the scanner: canonical peak per bond on its default grid, then its envelope at exactly that Q
        var scan = new ResonanceScan(block);
        var curve = scan.ComputeKCurve();
        var qPeaks = Enumerable.Range(0, nb).Select(b => curve.CanonicalPeakAtBond(b).QPeak).ToArray();
        var atPeaks = scan.ComputeKCurve(qPeaks);

        var readings = new List<FloorSplitReading>(nb);
        for (int b = 0; b < nb; b++)
        {
            double q = qPeaks[b];
            double tStar = atPeaks.TAtPeakByBond[b, b];
            double scanK = atPeaks.KByBond[b, b];
            var (rho, drho) = Propagate(decomp, probe, q * GammaZero, b, tStar);
            var rhoW = bw * (bwH * rho); var drhoW = bw * (bwH * drho);
            var rhoR = rho - rhoW; var drhoR = drho - drhoW;
            double k = TwoRe(rho, sKernel * drho);
            double kWW = TwoRe(rhoW, sKernel * drhoW);
            double kRR = TwoRe(rhoR, sKernel * drhoR);
            double kX = TwoRe(rhoW, sKernel * drhoR) + TwoRe(rhoR, sKernel * drhoW);
            double kClosed = FloorClosedForm(N, b, q * GammaZero, GammaZero, tStar, s, probe[0].Real, cells, bondSign[0]);
            readings.Add(new FloorSplitReading(b, q, tStar, k, kWW, kRR, kX, kClosed, scanK, closed[b]));
        }
        return readings;
    }

    private static double TwoRe(ComplexVector a, ComplexVector sb) => 2.0 * a.ConjugateDotProduct(sb).Real;

    /// <summary>ρ(t) and ∂ρ/∂J_b at uniform J by the eigen-Duhamel route.</summary>
    private static (ComplexVector Rho, ComplexVector DRho) Propagate(BlockLDecomposition decomp, ComplexVector probe,
        double j, int bond, double t)
    {
        var l = decomp.AssembleUniform(j);
        var evd = l.Evd();
        var r = evd.EigenVectors;
        var rInv = r.Inverse();
        var lam = evd.EigenValues;
        int dim = lam.Count;
        var c0 = rInv * probe;
        var e = ComplexVector.Build.Dense(dim, i => Complex.Exp(lam[i] * t));
        var rho = r * e.PointwiseMultiply(c0);
        var xb = rInv * decomp.MhPerBond[bond] * r;
        var f = ComplexVector.Build.Dense(dim);
        for (int row = 0; row < dim; row++)
        {
            Complex acc = Complex.Zero;
            for (int col = 0; col < dim; col++)
            {
                Complex d = lam[col] - lam[row];
                Complex w = d.Magnitude > 1e-10 ? (e[col] - e[row]) / d : t * e[row];
                acc += xb[row, col] * w * c0[col];
            }
            f[row] = acc;
        }
        return (rho, r * f);
    }

    /// <summary>The floor² term from the single-particle problem alone: x(t) = e^{t(−2γ₀ + iJσA)}·x(0)
    /// with x_j(0) = probe·s_j/C, the readout amplitude of site l is a_l = s_l·x_l, ∂a/∂J_b by Duhamel
    /// in the sine modes of A, and K = 4·Re Σ_l conj(a_l)·∂a_l (S_kernel = 2 Σ_l w_l w_lᵀ).</summary>
    internal static double FloorClosedForm(int n, int bond, double j, double gamma, double t,
        double[] s, double probeValue, int cells, int sigma)
    {
        var phi = new double[n, n];   // phi[site, k-1]
        var lam = new double[n];
        for (int k = 1; k <= n; k++)
        {
            double th = Math.PI * k / (n + 1);
            lam[k - 1] = 2.0 * Math.Cos(th);
            for (int site = 0; site < n; site++)
                phi[site, k - 1] = Math.Sqrt(2.0 / (n + 1)) * Math.Sin(th * (site + 1));
        }
        var shat = new double[n];
        for (int k = 0; k < n; k++) for (int site = 0; site < n; site++) shat[k] += phi[site, k] * s[site];
        // A_b in the mode basis
        var abk = new double[n, n];
        for (int k = 0; k < n; k++)
            for (int kk = 0; kk < n; kk++)
                abk[k, kk] = phi[bond, k] * phi[bond + 1, kk] + phi[bond + 1, k] * phi[bond, kk];
        // generator eigenvalues μ_k = −2γ + iJσλ_k; x̂(t) = e^{μt}·x̂(0), x̂(0) = (probe/C)·ŝ
        var mu = lam.Select(l => new Complex(-2.0 * gamma, j * sigma * l)).ToArray();
        var ex = mu.Select(z => Complex.Exp(z * t)).ToArray();
        double pre = probeValue / cells;
        var xHat = new Complex[n];
        var dxHat = new Complex[n];
        for (int k = 0; k < n; k++)
        {
            xHat[k] = ex[k] * pre * shat[k];
            Complex acc = Complex.Zero;
            for (int kk = 0; kk < n; kk++)
            {
                Complex d = mu[kk] - mu[k];
                Complex w = d.Magnitude > 1e-12 ? (ex[kk] - ex[k]) / d : t * ex[k];
                acc += new Complex(0, sigma) * abk[k, kk] * w * pre * shat[kk];
            }
            dxHat[k] = acc;
        }
        double kSum = 0.0;
        for (int site = 0; site < n; site++)
        {
            Complex a = Complex.Zero, da = Complex.Zero;
            for (int k = 0; k < n; k++) { a += phi[site, k] * xHat[k]; da += phi[site, k] * dxHat[k]; }
            kSum += (Complex.Conjugate(s[site] * a) * (s[site] * da)).Real;
        }
        return 4.0 * kSum;
    }

    /// <summary>Whether the single bond's integer hopping K_b = −i·M_b maps every u_j into the span of
    /// the family, and the sign with which it does (the coefficient on its neighbour across the bond).</summary>
    private static (bool Closed, int Sign) ClosureUnderBond(ComplexMatrix mh, Dictionary<int, int>[] family,
        Dictionary<int, int> owner, int bond)
    {
        int n = family.Length, m = mh.RowCount;
        int sign = 0;
        for (int j = 0; j < n; j++)
        {
            var image = new Dictionary<int, long>();
            foreach (var (cell, val) in family[j])
                for (int row = 0; row < m; row++)
                {
                    Complex z = mh[row, cell];
                    if (z == Complex.Zero) continue;
                    double kr = z.Imaginary;                       // K = −i·M: entry ±i becomes ±1
                    long kv = (long)Math.Round(kr);
                    if (z.Real != 0.0 || kr != kv) return (false, 0);
                    image[row] = image.GetValueOrDefault(row) + kv * val;
                }
            var coef = new Dictionary<int, long>();
            foreach (var (cell, v) in image)
            {
                if (v == 0) continue;
                if (!owner.TryGetValue(cell, out int jj)) return (false, 0);
                long c = v * family[jj][cell];
                if (coef.TryGetValue(jj, out long prev) && prev != c) return (false, 0);
                coef[jj] = c;
            }
            foreach (var (jj, c) in coef)
                foreach (var (cell, v) in family[jj])
                    if (image.GetValueOrDefault(cell) != c * v) return (false, 0);
            foreach (var (jj, c) in coef)
            {
                bool neighbourAcross = (j == bond && jj == bond + 1) || (j == bond + 1 && jj == bond);
                if (!neighbourAcross) return (false, 0);
                if (sign == 0) sign = (int)c; else if (sign != c) return (false, 0);
            }
            if ((j == bond || j == bond + 1) && coef.Count != 1) return (false, 0);
            if (j != bond && j != bond + 1 && coef.Count != 0) return (false, 0);
        }
        return (sign is 1 or -1, sign);
    }

    // ---- Object Manager ----
    public string DisplayName =>
        $"FloorSplitWitness (N={N}, c = 2 block, γ₀={GammaZero.ToString("0.###", Inv)})";

    public string Summary =>
        $"F86's K at each bond's canonical Q_peak, split into floor², rest² and their interference; " +
        $"endpoint floor share {Bonds[0].FloorShare.ToString("0.000", Inv)}, interference {Bonds[0].CrossShare.ToString("0.000", Inv)}";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("the gate (three routes on one K)",
                summary: AllRoutesMeet
                    ? "routes meet: every bond closes the family in integers, floor²+rest²+interference equals ResonanceScan's " +
                      "envelope at the same Q, and the floor term projected in the block equals its single-particle closed form"
                    : "routes do NOT meet: " + string.Join("; ", Bonds.Select(b =>
                        $"bond {b.Bond}: closed {b.FamilyClosedUnderThisBond}, sum {(b.KFloor + b.KRest + b.KCross).ToString("E3", Inv)} " +
                        $"vs scan {b.ScanEnvelope.ToString("E3", Inv)}, floor {b.KFloor.ToString("E3", Inv)} vs closed {b.KFloorClosedForm.ToString("E3", Inv)}")));
            foreach (var b in Bonds)
                yield return new InspectableNode(
                    $"bond {b.Bond} (Q_peak {b.QPeak.ToString("0.0000", Inv)}, t* {b.TStar.ToString("0.000", Inv)})",
                    summary: $"K {b.K.ToString("0.000000", Inv)} = floor² {b.KFloor.ToString("0.000000", Inv)} ({b.FloorShare.ToString("0.000", Inv)})" +
                             $" + rest² {b.KRest.ToString("0.000000", Inv)} ({b.RestShare.ToString("0.000", Inv)})" +
                             $" + interference {b.KCross.ToString("0.000000", Inv)} ({b.CrossShare.ToString("0.000", Inv)})");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}

/// <summary>One bond's split at its canonical Q_peak and envelope time t*: K (signed) and its three
/// terms; <paramref name="KFloorClosedForm"/> the floor term from the single-particle route;
/// <paramref name="ScanEnvelope"/> ResonanceScan's |K| envelope at the same Q.</summary>
public sealed record FloorSplitReading(
    int Bond,
    double QPeak,
    double TStar,
    double K,
    double KFloor,
    double KRest,
    double KCross,
    double KFloorClosedForm,
    double ScanEnvelope,
    bool FamilyClosedUnderThisBond)
{
    public double FloorShare => KFloor / K;
    public double RestShare => KRest / K;
    public double CrossShare => KCross / K;
}
