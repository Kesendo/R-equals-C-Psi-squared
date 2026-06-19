using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>One defect bond's reading of the survivor's diffusion mode: the first-order rate shift
/// |dRe(b)| under a delta-J bond defect, the squared density gradient grad^2 = (n(j)-n(j+1))^2 of the
/// slow mode across that bond, and their ratio (bond-independent iff the diffusion-Rayleigh law holds).</summary>
public readonly record struct GradientBond(int SiteA, int SiteB, double RateShift, double GradSq, double Ratio);

/// <summary>(D) THE CLOSURE FUNCTIONAL (felt_time_dimensions arc, D follow-up), graduated to a live
/// witness (<c>inspect --root gradient</c>): the C# landing of <c>simulations/_felt_time_amplitude_law.py</c>.
///
/// <para>The PTF painter closure Sum f(b) (step B / <see cref="StoneSurvivorClosureWitness"/>) reads the
/// survivor's first-order RATE shift dRe(b). This witness gives that rate shift its EXACT bond functional:
/// the slow survivor is a DENSITY / diffusion mode n(j) (not a current -- its hopping content is zero), so
/// a delta-J defect on bond b=(j,j+1) perturbs the LOCAL diffusion coefficient D_b, and the first-order
/// rate shift is the diffusion Rayleigh-quotient derivative dRe(b) = dlambda/dD_b ~ (n(j)-n(j+1))^2 -- the
/// SQUARED density GRADIENT ("amplitude^2"). It vanishes at the no-flux (reflecting) chain ends (the
/// gradient -> 0 there) and peaks in the interior, mirror-symmetric, Q-invariant (the lowest diffusion
/// harmonic k_min is Q-fixed).</para>
///
/// <para>dRe(b) is the biorthogonal first-order shift of the slowest (p,p) block mode: with R the block
/// eigenvector matrix, dRe(b) = Re[(R^-1 row)_s . V_b . (R col)_s], where V_b = BondPerturbation (the
/// commutator superoperator of the unit-coupling XY bond, restricted to the sector). The density n(j) is
/// the occupation profile of the embedded slow mode. The earlier single-particle phi*phi candidate used
/// the WRONG standing wave (single-particle, not the multi-magnon density mode): right power, wrong wave.
/// Gate-first verified N=4..7 (slope dRe vs |grad| = 2.0..2.2, dRe/grad^2 bond-independent).</para></summary>
public sealed class SurvivorDiffusionGradientWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Qh = 0.5;          // H = Qh*(XX+YY) (carbon J=1; ChainSystem J = 2*Qh)
    private const double KernelTol = 1e-7;

    public int N { get; }
    public double Q { get; }

    public SurvivorDiffusionGradientWitness(int n = 4, double q = 1.5)
    {
        if (n < 4 || n > 5)
            throw new ArgumentOutOfRangeException(nameof(n), n, "N in 4..5 (interior survivor needs N>=4; the per-bond full-V build is cheap only through N=5)");
        if (q <= 0) throw new ArgumentOutOfRangeException(nameof(q), q, "Q must be > 0");
        N = n; Q = q;
    }

    private IReadOnlyList<GradientBond>? _bonds;
    private int _survivorP;
    private double[]? _density;
    private double _offDiag;

    /// <summary>Per-bond readout (all chain bonds, in order): |dRe(b)|, grad^2(b), and the ratio.</summary>
    public IReadOnlyList<GradientBond> Bonds { get { _bonds ??= Compute(); return _bonds; } }

    /// <summary>The survivor's half-filling sector p (the (p,p) block carrying the global-slowest mode).</summary>
    public int SurvivorP { get { _ = Bonds; return _survivorP; } }

    /// <summary>The slow diffusion mode's per-site density profile n(j).</summary>
    public IReadOnlyList<double> Density { get { _ = Bonds; return _density!; } }

    /// <summary>Density-mode character: ‖M − diag(M)‖ / ‖M‖ of the embedded slow mode (0 = a pure density
    /// mode, the diffusion-Rayleigh premise). ~0 in the strong-dephasing limit (Q -> 0); grows with Q as
    /// the survivor picks up off-diagonal coherence dressing, and the law (slope, CV) degrades with it.</summary>
    public double OffDiagonalWeight { get { _ = Bonds; return _offDiag; } }

    /// <summary>Coefficient of variation of |dRe(b)|/grad^2(b) over the bonds: ~0 means the rate shift is
    /// LINEAR in grad^2 with a bond-independent constant (the diffusion-Rayleigh law holds).</summary>
    public double RatioCv
    {
        get
        {
            var r = Bonds.Where(b => !double.IsNaN(b.Ratio)).Select(b => b.Ratio).ToArray();
            if (r.Length < 2) return double.NaN;
            double mean = r.Average();
            double var = r.Select(x => (x - mean) * (x - mean)).Sum() / r.Length;
            return Math.Sqrt(var) / Math.Max(Math.Abs(mean), 1e-15);
        }
    }

    /// <summary>Log-log slope of |dRe(b)| vs |grad(b)|: ~2 confirms dRe ~ (gradient)^2.</summary>
    public double PowerSlope
    {
        get
        {
            var pts = Bonds.Where(b => b.GradSq > 1e-12)
                .Select(b => (x: 0.5 * Math.Log(b.GradSq), y: Math.Log(Math.Max(b.RateShift, 1e-30)))).ToArray();
            if (pts.Length < 2) return double.NaN;
            double mx = pts.Average(p => p.x), my = pts.Average(p => p.y);
            double cov = pts.Sum(p => (p.x - mx) * (p.y - my));
            double vx = pts.Sum(p => (p.x - mx) * (p.x - mx));
            return vx > 1e-15 ? cov / vx : double.NaN;
        }
    }

    /// <summary>The law holds: rate shift is amplitude^2 in the density gradient (slope ~2, ratio constant).</summary>
    public bool LawHolds => RatioCv < 0.15 && Math.Abs(PowerSlope - 2.0) < 0.4;

    private IReadOnlyList<GradientBond> Compute()
    {
        int d = 1 << N;
        double gamma = 1.0 / Q;
        var profile = Enumerable.Repeat(gamma, N).ToArray();
        var chain = new ChainSystem(N, 2.0 * Qh, 1.0, HamiltonianType.XY, TopologyKind.Chain);
        var H = chain.BuildHamiltonian();

        // 1. the survivor sector: the (p,p) diagonal block carrying the global-slowest soft mode.
        int survivorP = 1;
        double bestRe = double.NegativeInfinity;
        for (int p = 1; p < N; p++)
        {
            var lamP = PerBlockLiouvillianBuilder.BuildBlockZ(H, profile, SectorFlat(N, p, p)).Evd().EigenValues;
            double slowP = lamP.Where(z => z.Real < -KernelTol).Select(z => z.Real).DefaultIfEmpty(double.NegativeInfinity).Max();
            if (slowP > bestRe) { bestRe = slowP; survivorP = p; }
        }
        _survivorP = survivorP;
        var flat = SectorFlat(N, survivorP, survivorP);
        var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, profile, flat);

        // 2. the slowest right-eigenmode + its biorthogonal dual (R^-1 row): <l|r> = 1 by construction.
        var evd = block.Evd();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();
        var lam = evd.EigenValues;
        int idx = -1; double best = double.NegativeInfinity;
        for (int i = 0; i < lam.Count; i++)
            if (lam[i].Real < -KernelTol && lam[i].Real > best) { best = lam[i].Real; idx = i; }
        var rcol = R.Column(idx);
        var lrow = Rinv.Row(idx);

        // 3. the density profile n(j) of the slow mode (its diagonal occupation), phase-fixed to real.
        var rUnit = rcol / rcol.L2Norm();
        var full = new Complex[d * d];
        for (int k = 0; k < flat.Length; k++) full[flat[k]] = rUnit[k];
        Complex anchor = Complex.Zero;
        for (int s = 0; s < d; s++) if (full[s * d + s].Magnitude > anchor.Magnitude) anchor = full[s * d + s];
        Complex phase = anchor.Magnitude > 1e-15 ? Complex.Conjugate(anchor) / anchor.Magnitude : Complex.One;
        var n = new double[N];
        for (int s = 0; s < d; s++)
        {
            double diag = (phase * full[s * d + s]).Real;
            if (diag != 0.0)
                for (int j = 0; j < N; j++) if (((s >> j) & 1) != 0) n[j] += diag;
        }
        _density = n;

        // density-mode character: ‖M − diag(M)‖ / ‖M‖ of the embedded mode. full is unit-norm (rUnit), so
        // ‖M‖² = Σ|full|²; the off-diagonal weight is sqrt(1 − diagonal-weight). ~0 = a pure density mode.
        double totSq = 0.0, diagSq = 0.0;
        for (int i = 0; i < full.Length; i++) totSq += full[i].Real * full[i].Real + full[i].Imaginary * full[i].Imaginary;
        for (int s = 0; s < d; s++) { var z = full[s * d + s]; diagSq += z.Real * z.Real + z.Imaginary * z.Imaginary; }
        _offDiag = totSq > 1e-300 ? Math.Sqrt(Math.Max(0.0, 1.0 - diagSq / totSq)) : 0.0;

        // 4. per bond: |dRe(b)| = |Re (R^-1 row . V_b . R col)| and grad^2 = (n(j)-n(j+1))^2.
        var bonds = new List<GradientBond>(chain.Bonds.Count);
        foreach (var bond in chain.Bonds)
        {
            var Vfull = BondPerturbation.Build(N, bond.Site1, bond.Site2, BondPerturbation.Kind.XY);
            var Vblock = ComplexMatrix.Build.Dense(flat.Length, flat.Length, (a, b) => Vfull[flat[a], flat[b]]);
            double dRe = (lrow.DotProduct(Vblock * rcol)).Real;     // biorthogonal first-order shift
            double grad = n[bond.Site1] - n[bond.Site2];
            double gradSq = grad * grad;
            double ratio = gradSq > 1e-12 ? Math.Abs(dRe) / gradSq : double.NaN;
            bonds.Add(new GradientBond(bond.Site1, bond.Site2, Math.Abs(dRe), gradSq, ratio));
        }
        return bonds;
    }

    /// <summary>The flat (row*d+col) indices of the (pCol,pRow) joint-popcount sector (the
    /// SectorReductionWitness recipe, as in <see cref="StoneSurvivorClosureWitness"/>).</summary>
    private static int[] SectorFlat(int n, int pCol, int pRow)
    {
        var decomp = JointPopcountSectorBuilder.Build(n);
        var s = decomp.SectorRanges.First(r => r.PCol == pCol && r.PRow == pRow);
        var flat = new int[s.Size];
        for (int k = 0; k < s.Size; k++) flat[k] = decomp.Permutation[s.Offset + k];
        return flat;
    }

    public string DisplayName => $"SurvivorDiffusionGradientWitness (N={N}, Q={Q.ToString("0.##", Inv)})";

    public string Summary =>
        "(D) the felt_time closure functional: the survivor's first-order bond rate shift dRe(b) ~ " +
        "(density-mode gradient at bond b)^2 -- the diffusion Rayleigh quotient (amplitude^2). The slow " +
        $"survivor (sector ({SurvivorP},{SurvivorP})) is a DENSITY/diffusion mode; a delta-J defect perturbs the " +
        "local diffusion coefficient, so dRe = dlambda/dD_b ~ (n(j)-n(j+1))^2, ~0 at the no-flux chain ends. " +
        $"dRe/grad^2 bond-independent (CV={RatioCv.ToString("0.000", Inv)}), log-log slope dRe vs |grad| = " +
        $"{PowerSlope.ToString("0.00", Inv)} ({(LawHolds ? "amplitude^2 CONFIRMED" : "law not clean")}). The " +
        "eigenvalue-level functional underlying the PTF closure (inspect --root stone); right power, wrong " +
        "wave for the old single-particle phi*phi guess.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("the diffusion-Rayleigh law",
                summary: $"dRe(b)/grad^2 over the {Bonds.Count} bonds has CV={RatioCv.ToString("0.000", Inv)} " +
                         $"(bond-independent => LINEAR in grad^2); log-log slope dRe vs |grad| = {PowerSlope.ToString("0.00", Inv)} " +
                         $"(~2 => amplitude^2). {(LawHolds ? "LAW HOLDS." : "law not clean.")}");
            yield return new InspectableNode("the density standing wave",
                summary: $"n(j) = [{string.Join(", ", Density.Select(x => x.ToString("+0.00;-0.00", Inv)))}] " +
                         $"(off-diag weight {OffDiagonalWeight.ToString("0.000", Inv)}: " +
                         $"{(OffDiagonalWeight < 0.1 ? "a pure density mode (strong-dephasing limit)" : "dressed with coherence at this Q")}; " +
                         "antisymmetric, quiet at the reflecting ends).");
            foreach (var b in Bonds)
                yield return new InspectableNode($"bond ({b.SiteA},{b.SiteB})",
                    summary: $"|dRe| = {b.RateShift.ToString("0.0000", Inv)}, grad^2 = {b.GradSq.ToString("0.00000", Inv)}, " +
                             $"ratio = {(double.IsNaN(b.Ratio) ? "n/a" : b.Ratio.ToString("0.000", Inv))}");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
