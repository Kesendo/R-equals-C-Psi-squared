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
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>One mode's closure readout: the PTF painter closure Sum_i ln(alpha_i) over the reliable
/// sites, its sign-coherence coh = |Sum f|/Sum|f| (1 = a global RATE shift, ~0 = a redistribution),
/// and whether it sits in the +-0.05 window.</summary>
public readonly record struct ModeClosure(
    string Label, int PCol, int PRow, double Closure, double SignCoherence, int Reliable, int Sites, bool InWindow)
{
    /// <summary>A genuine K_decay (rate) shift: OUT of the window AND sign-coherent, so it is a global
    /// rate change, not an asymmetric per-site redistribution (which would also move Sum ln alpha).</summary>
    public bool IsRateShift => !InWindow && SignCoherence > 0.8;
}

/// <summary>THE STONE (felt_time_dimensions arc, step B), graduated to a live witness
/// (<c>inspect --root stone</c>): the C# landing of <c>simulations/stone_survivor_alpha_closure.py</c>.
///
/// <para>For the near-stationary MODE-ISOLATING probe rho_0 = I/d + eps*Herm(mode) (I/d is the stationary
/// steady state L(I/d)=0, so the single-site purity is driven almost only by the chosen eigenmode), the
/// PTF painter closure Sum_i ln(alpha_i) -- computed through the CANONICAL Symphony FitAlpha / PuritySites
/// -- reads the mode's first-order RATE shift under a delta-J bond defect: it BREAKS (out of the +-0.05
/// window) AND is sign-coherent for the soft survivor interior (2,2) mode (Re lambda moves), and HOLDS
/// (in window) for the rigid (0,1) band edge (Re = -2gamma, frozen). This is the TRAJECTORY-level dual of
/// the eigenvalue-level value/vector split (A) -- a CONSTRUCTIVE confirmation for this probe, NOT a
/// universal trajectory law.</para>
///
/// <para>Review-pinned scope (two-lens panel 2026-06-19): the break is probe-state-dependent (a survivor-
/// dominated but POLARIZED state holds); Sum ln alpha != 0 certifies a rate shift only via the
/// SIGN-COHERENCE of the reliable per-site f (asserted: <see cref="ModeClosure.IsRateShift"/>); the
/// magnitude is scaling + sign, not quantitative. The standing-wave NODE site (odd N) drops out of the
/// reliable set (|f|>10), exactly as Symphony's guard intends.</para></summary>
public sealed class StoneSurvivorClosureWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Qh = 0.5;             // H = Qh*(XX+YY) (the carbon J=1 convention; ChainSystem J=2*Qh)
    private const double ClosureWindow = 0.05; // EQ-014 (= PaintersMovement.ClosureWindow)
    private const double FMax = 10.0, ConsistencyTol = 0.5, FFloor = 0.5;  // = PaintersMovement reliability guard
    private const double KernelTol = 1e-7;

    public int N { get; }
    public double Q { get; }
    public double DeltaJ { get; }

    public StoneSurvivorClosureWitness(int n = 4, double q = 1.5, double deltaJ = 0.02)
    {
        if (n < 4 || n > 5)
            throw new ArgumentOutOfRangeException(nameof(n), n, "N in 4..5 (interior (2,2) needs N>=4; the canonical alpha-fit path is N<=5)");
        if (q <= 0) throw new ArgumentOutOfRangeException(nameof(q), q, "Q must be > 0");
        if (deltaJ <= 0 || deltaJ > 0.1) throw new ArgumentOutOfRangeException(nameof(deltaJ), deltaJ, "deltaJ in (0, 0.1] (perturbative)");
        N = n; Q = q; DeltaJ = deltaJ;
    }

    private ModeClosure? _survivor, _bandEdge;
    /// <summary>The soft interior incompleteness survivor (2,2): OUT + sign-coherent = a rate shift.</summary>
    public ModeClosure Survivor { get { _survivor ??= ComputeMode("survivor interior", 2, 2); return _survivor.Value; } }
    /// <summary>The rigid (0,1) band edge (Re = -2gamma): IN window = frozen.</summary>
    public ModeClosure BandEdge { get { _bandEdge ??= ComputeMode("band edge", 0, 1); return _bandEdge.Value; } }

    private ModeClosure ComputeMode(string label, int pc, int pr)
    {
        int d = 1 << N;
        double gamma = 1.0 / Q;
        var profile = Enumerable.Repeat(gamma, N).ToArray();
        var chain = new ChainSystem(N, 2.0 * Qh, 1.0, HamiltonianType.XY, TopologyKind.Chain);
        var H = chain.BuildHamiltonian();

        // 1. the (pc,pr) block slowest right-eigenmode, embedded as a d x d operator
        //    (the SectorReductionWitness.DensityMode recipe; flat = row*d + col).
        var flat = SectorFlat(N, pc, pr);
        var block = PerBlockLiouvillianBuilder.BuildBlockZ(H, profile, flat);
        var slow = PhaseRigidity.Compute(block)
            .Where(m => m.Lambda.Magnitude > KernelTol)
            .OrderByDescending(m => m.Lambda.Real).First();
        double reMode = slow.Lambda.Real;
        var full = ComplexVector.Build.Dense(d * d);
        for (int k = 0; k < flat.Length; k++) full[flat[k]] = slow.Right[k];
        var M = ComplexMatrix.Build.Dense(d, d, (a, b) => full[a * d + b]);

        // 2. the mode-isolating probe rho_0 = I/d + eps*Herm(mode): Hermitian, traceless, unit-Frobenius
        //    Herm(mode); eps keeps rho_0 PSD. (I/d is stationary, so only this mode drives the dynamics.)
        var I = ComplexMatrix.Build.DenseIdentity(d);
        var Mh = 0.5 * (M + M.ConjugateTranspose());
        Mh -= (Mh.Trace() / (double)d) * I;
        Mh /= Mh.FrobeniusNorm();
        double minEig = Mh.Evd().EigenValues.Select(z => z.Real).Min();
        double eps = 0.5 / d / Math.Max(Math.Abs(minEig), 1e-9);
        var rho0 = (Complex)(1.0 / d) * I + (Complex)eps * Mh;

        // 3. clean + defected (deltaJ) + guard (deltaJ/2) full Liouvillians; defect on a non-central bond.
        var L_A = PauliDephasingDissipator.BuildZ(H, profile);
        var bond = chain.Bonds[chain.Bonds.Count / 3];
        var V = BondPerturbation.Build(N, bond.Site1, bond.Site2, BondPerturbation.Kind.XY);
        var L_B = L_A + (Complex)DeltaJ * V;
        var L_guard = L_A + (Complex)(DeltaJ / 2.0) * V;

        // 4. site Paulis + a time grid that captures the mode's decay (5 lifetimes, 60 points).
        var sitePaulis = new (ComplexMatrix X, ComplexMatrix Y, ComplexMatrix Z)[N];
        for (int i = 0; i < N; i++)
            sitePaulis[i] = (PauliString.SiteOp(N, i, PauliLetter.X),
                             PauliString.SiteOp(N, i, PauliLetter.Y),
                             PauliString.SiteOp(N, i, PauliLetter.Z));
        double tmax = 5.0 / Math.Max(Math.Abs(reMode), 1e-3);
        var tGrid = new double[60];
        for (int s = 0; s < 60; s++) tGrid[s] = tmax * s / 59.0;

        // 5. the CANONICAL Symphony propagation + alpha-fit (PaintersMovement.PuritySites / FitAlpha).
        var pA = PaintersMovement.PuritySites(N, L_A, rho0, tGrid, sitePaulis);
        var pB = PaintersMovement.PuritySites(N, L_B, rho0, tGrid, sitePaulis);
        var pGuard = PaintersMovement.PuritySites(N, L_guard, rho0, tGrid, sitePaulis);
        var alpha = new double[N];
        var alphaGuard = new double[N];
        for (int i = 0; i < N; i++)
        {
            alpha[i] = PaintersMovement.FitAlpha(tGrid, pA[i], pB[i]);
            alphaGuard[i] = PaintersMovement.FitAlpha(tGrid, pA[i], pGuard[i]);
        }

        // 6. the reliability guard + closure + sign-coherence (PaintersMovement's formula; the standing-wave
        //    NODE site is ill-conditioned -> |f| large -> drops out; coh = |sum f|/sum|f| over reliable).
        double sumLn = 0.0, sumF = 0.0, sumAbsF = 0.0;
        int nrel = 0;
        for (int i = 0; i < N; i++)
        {
            double f = (alpha[i] - 1.0) / DeltaJ;
            double fGuard = (alphaGuard[i] - 1.0) / (DeltaJ / 2.0);
            bool reliable = Math.Abs(f) <= FMax && Math.Abs(f - fGuard) <= ConsistencyTol * (Math.Abs(fGuard) + FFloor);
            if (reliable) { sumLn += Math.Log(Math.Max(alpha[i], 1e-30)); sumF += f; sumAbsF += Math.Abs(f); nrel++; }
        }
        double closure = nrel > 0 ? sumLn : double.NaN;
        double coh = sumAbsF > 0 ? Math.Abs(sumF) / sumAbsF : double.NaN;
        return new ModeClosure(label, pc, pr, closure, coh, nrel, N, Math.Abs(closure) <= ClosureWindow);
    }

    /// <summary>The flat (row*d+col) indices of the (pCol,pRow) joint-popcount sector (the private
    /// SectorReductionWitness recipe, replicated).</summary>
    private static int[] SectorFlat(int n, int pCol, int pRow)
    {
        var decomp = JointPopcountSectorBuilder.Build(n);
        var s = decomp.SectorRanges.First(r => r.PCol == pCol && r.PRow == pRow);
        var flat = new int[s.Size];
        for (int k = 0; k < s.Size; k++) flat[k] = decomp.Permutation[s.Offset + k];
        return flat;
    }

    public string DisplayName => $"StoneSurvivorClosureWitness (N={N}, Q={Q.ToString("0.##", Inv)}, dJ={DeltaJ.ToString("0.###", Inv)})";

    public string Summary =>
        "the stone (felt_time arc B): the PTF painter closure Sum_i ln(alpha_i), through the CANONICAL Symphony " +
        "FitAlpha on the mode-isolating probe rho_0=I/d+eps*Herm(mode), reads the mode's first-order RATE shift. " +
        $"survivor interior (2,2) -> {(Survivor.IsRateShift ? "OUT + sign-coherent = RATE-SHIFT" : "in window")} " +
        $"(Sum ln a={Survivor.Closure.ToString("+0.000;-0.000", Inv)}, coh={Survivor.SignCoherence.ToString("0.00", Inv)}); " +
        $"band edge (0,1) -> {(BandEdge.InWindow ? "IN = FROZEN" : "out of window")} " +
        $"(Sum ln a={BandEdge.Closure.ToString("+0.000;-0.000", Inv)}). The TRAJECTORY-level dual of the eigenvalue " +
        "value/vector split (inspect --root survivor / horizon). Review-pinned: probe-state-specific (a polarized " +
        "survivor-dominated state holds), the rate shift certified by sign-coherence, magnitude = scaling+sign.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return ModeNode(Survivor);
            yield return ModeNode(BandEdge);
        }
    }

    private InspectableNode ModeNode(ModeClosure m) =>
        new InspectableNode($"{m.Label} ({m.PCol},{m.PRow})",
            summary: $"Sum ln alpha = {m.Closure.ToString("+0.0000;-0.0000", Inv)} " +
                     $"({(m.InWindow ? "IN" : "OUT")} the +-{ClosureWindow.ToString("0.##", Inv)} window), " +
                     $"sign-coherence coh={m.SignCoherence.ToString("0.00", Inv)} ({m.Reliable}/{m.Sites} reliable sites) -> " +
                     $"{(m.IsRateShift ? "RATE-SHIFT (K_decay moves; soft darkness)" : m.InWindow ? "FROZEN (K_decay defect-invariant; rigid darkness)" : "out-of-window but sign-mixed (redistribution, NOT a rate shift)")}");

    public InspectablePayload Payload => InspectablePayload.Empty;
}
