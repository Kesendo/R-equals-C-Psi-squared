using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live lab for the coherence horizon Q*(N) (typed home intent:
/// <c>docs/ANALYTICAL_FORMULAS.md</c> F2b corollary "Coherence horizon Q*(N)"). Two thresholds live
/// here under one letter. The EP Q*(N) is where the single-excitation pair coalesces, a square-root
/// exceptional point (1, √2, 1.87874, 2.37367 at N=2..5; <see cref="EpQ"/>). The handover this witness
/// bisects on the FULL Liouvillian (<see cref="Horizon"/>) is where the slowest non-zero mode stops
/// oscillating: the handover Q_h, the clock's takeover (1, √2, 1.87854, 2.37217 on the SE block; the full-L
/// bisection reads them to its resolution). They coincide at N=2,3 and part from N=4: the full L holds
/// the (0,1) band-edge survivor at exactly −2γ, and the coalescer sits on that floor only where its
/// coherence share w2 is exactly ½ (the Absorption Theorem, Re λ = −2γ⟨n_diff⟩, so Re = −4γ·w2 on
/// this block). From N=4 the coalescer sits below the floor at every Q ≥ Q*; below Q* the pair
/// splits into two real branches, and the darker one reaches the floor, becomes the gap and freezes
/// takes the clock over at the HANDOVER Q_h = Q* − ((2w2−1)/c)² (<see cref="HandoverQ"/>, the excess-light node).
/// Why exactly at N=2,3: at an EP the pair's light is its Jordan plane's mean light (2λ_EP = tr L|_V, the
/// trace-midpoint 2×2 of <c>docs/proofs/PROOF_F86B_UNIVERSAL_SHAPE.md</c>), and a plane of one population
/// direction (light 0) and one coherence direction (light 2) has mean 1, the midpoint of two rungs. Why the
/// handover is a crossing and never an avoided one: (−1)^{n_XY} commutes with L (the Absorption proof's
/// parity-preserving mixing), so the (1,1) pair (light 0, 2) and the (0,1) survivor (light 1) never hybridize.
///
/// <para>What the Hückel side supplies, and where it stops. The single-excitation block of H in the site
/// basis IS the Hückel matrix at α=0, β=J: (J/2)(X_lX_{l+1} + Y_lY_{l+1}) = J(σ⁺_lσ⁻_{l+1} + h.c.), so
/// restricting this number-conserving chain to the one-excitation sector leaves the tridiagonal hopping matrix
/// with off-diagonal J and nothing on the diagonal, and <see cref="HuckelResidual"/> is exactly 0.0, not small.
/// The sign β = +J is a choice and not a claim: on a bipartite chain the sublattice operator diag((−1)^l)
/// conjugates A into −A, so the two signs give the same spectrum. That is the
/// whole of the inheritance, and it is a statement about the HAMILTONIAN. Both thresholds live on the
/// dephasing side, and a Hückel or Frost construction is a STATIC spectrum with no γ, no bath and no time in
/// it, so it carries no threshold to compare against; a threshold read off the Hückel side's own dephasing
/// extension is untested here rather than non-existent, and what this rests on is that no independent chemistry
/// computation was ever on the other side, which had this witness's own numbers on both.
/// N=2 (Q*=1) is the exceptional point itself (γ=J), where the ±J band mode ceases to be the gap mode.</para>
///
/// <para>The handover reuses <see cref="Symphony"/> as the spectrum engine: bisect γ at J=1 on whether
/// <see cref="Symphony.Clock"/>.Omega (max|Im λ| among the modes at the slowest non-zero decay rate)
/// is non-zero; computed once and cached (N=5 is a 1024×1024 eigendecomposition; ~16 bisection builds,
/// lazy and build-once). The EP, the handover and the coalescer's weights are read on the N²-dim
/// single-excitation block (<see cref="EpCharacterWitness.Lse"/>), where the same criterion lands on
/// the EP because that spectrum holds no survivor. The closed form of the EP: 1 and √2 exactly
/// at N=2,3 (a clean 2×2), transcendental from N=4 (the SE-block discriminant), asymptotic slope 2/π
/// derived (<c>docs/proofs/PROOF_COHERENCE_HORIZON_SLOPE.md</c>).</para></summary>
public sealed class CoherenceHorizonWitness : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The coupling the horizon is read at (J=1, so Q* = 1/γ*).</summary>
    public const double J = 1.0;

    /// <summary>The Omega threshold: |Im| above this counts as "the slowest mode still oscillates".</summary>
    private const double OmegaEps = 1e-3;

    /// <summary>γ bracket for the bisection: oscillating at γ_lo=0.3 (Q=3.33), frozen at γ_hi=1.2
    /// (Q=0.83), covering every Q* ∈ [0.83, 3.3].</summary>
    private const double GammaLo = 0.3, GammaHi = 1.2;
    private const int BisectionSteps = 16;

    /// <summary>The spectrum-reading grid: the clock reads the eigenvalues, so the time grid is
    /// resolution-independent; a small grid keeps the per-build cost down.</summary>
    private const int SpectrumPoints = 8;

    private readonly Dictionary<int, double> _horizonCache = new();

    /// <summary>The live handover on the full Liouvillian (the clock's takeover): the largest Q (= smallest
    /// γ at J=1) below which the slowest non-zero mode no longer oscillates. This is the handover Q_h(N), equal to the EP Q*(N)
    /// at N=2,3 only and below it from N=4 (<see cref="HandoverQ"/> reads the same event on the SE block;
    /// <see cref="EpQ"/> the EP): the full L holds the (0,1) survivor at −2γ, and the split pair's darker
    /// branch is the gap once it rises above that floor. Computed by bisecting γ on whether
    /// <see cref="Symphony.Clock"/>.Omega &gt; <see cref="OmegaEps"/>, then Q = J/γ*; resolution
    /// <see cref="HorizonResolution"/>. Cached per N (one bisection, ~16 dense eigendecompositions).</summary>
    public double Horizon(int n)
    {
        if (n < 2 || n > Symphony.MaxN)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"the coherence horizon needs N in 2..{Symphony.MaxN} (dense Liouvillian spectrum); got {n}");
        if (_horizonCache.TryGetValue(n, out var q)) return q;

        // γ_lo oscillates (Omega > ε), γ_hi is frozen (Omega ≈ 0); bisect to the crossover γ*.
        double glo = GammaLo, ghi = GammaHi;
        for (int step = 0; step < BisectionSteps; step++)
        {
            double mid = 0.5 * (glo + ghi);
            if (Omega(n, mid) > OmegaEps) glo = mid;   // still oscillating → need more γ (lower Q)
            else ghi = mid;                            // frozen → less γ (higher Q)
        }
        double gStar = 0.5 * (glo + ghi);
        q = J / gStar;
        _horizonCache[n] = q;
        return q;
    }

    /// <summary>Omega(N, γ) = the coherence hand at J=1: max|Im λ| among the modes sharing the slowest
    /// non-zero decay rate (the gap), read off the <see cref="Symphony"/> clock. &gt; ε ⟹ the slowest
    /// mode oscillates; ≈ 0 ⟹ it is overdamped (real).</summary>
    public double Omega(int n, double gamma) =>
        new Symphony(n: n, j: J, gamma: gamma, initialState: InitialStateKind.BellPair,
            tPoints: SpectrumPoints).Clock.Omega;

    /// <summary>The single-excitation block of H in the site basis: the N×N matrix ⟨i|H|j⟩ over the
    /// computational states carrying exactly one excitation, site l being bit N−1−l of the index.</summary>
    public static ComplexMatrix SingleExcitationBlock(int n, double j = J)
    {
        var H = new ChainSystem(n, j, 0.0).BuildHamiltonian();
        var block = ComplexMatrix.Build.Dense(n, n);
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++)
                block[r, c] = H[1 << (n - 1 - r), 1 << (n - 1 - c)];
        return block;
    }

    /// <summary>max entry-wise |H_SE − J·A|, with A the open chain's adjacency matrix (unit off-diagonal
    /// on the bonds, zero elsewhere): the Hückel matrix at α=0, β=J. This is exactly 0.0, not small.
    /// Jordan-Wigner sends (J/2)Σ(X_lX_{l+1} + Y_lY_{l+1}) onto the hopping matrix itself, so the block
    /// is that matrix rather than an approximation of it. The gate has something to catch: a ZZ term or a
    /// longitudinal field puts weight on the diagonal, where A carries none, and a non-uniform J leaves the
    /// bond entries no longer equal to one β. Measured on the Heisenberg chain, whose builder differs from
    /// the XY one in BOTH the hop (J/2 against J) and the added ZZ, the residual is 0.5; the ZZ part alone is
    /// isolated by <see cref="ZzDiagonalWeight"/>.</summary>
    public static double HuckelResidual(int n, double j = J)
    {
        var block = SingleExcitationBlock(n, j);
        double worst = 0.0;
        for (int r = 0; r < n; r++)
            for (int c = 0; c < n; c++)
            {
                double huckel = Math.Abs(r - c) == 1 ? j : 0.0;
                worst = Math.Max(worst, (block[r, c] - huckel).Magnitude);
            }
        return worst;
    }

    /// <summary>The ZZ term's own footprint on the single-excitation block, isolated from the couplings the
    /// two <see cref="HamiltonianType"/> builders happen to use: the largest |diagonal| entry of the
    /// Heisenberg chain's block. The Hückel matrix has a zero diagonal and so does the XY chain's block, at
    /// every coupling, so this weight is what the ZZ adds and nothing else. It is non-zero at every N≥2
    /// (0.25 / 0.5 / 0.25 / 0.5 at N=2..5, the bond-parity of a single flipped site), which is why the
    /// Hückel identity is a statement about the XY chain and not about spin chains in general.</summary>
    public static double ZzDiagonalWeight(int n)
    {
        var H = new ChainSystem(n, J, 0.0, HamiltonianType.Heisenberg).BuildHamiltonian();
        double worst = 0.0;
        for (int i = 0; i < n; i++)
            worst = Math.Max(worst, H[1 << (n - 1 - i), 1 << (n - 1 - i)].Magnitude);
        return worst;
    }

    /// <summary>The largest |diagonal| entry of the XY chain's single-excitation block, the quantity
    /// <see cref="ZzDiagonalWeight"/> is measured against: exactly 0.0, at every coupling.</summary>
    public static double XyDiagonalWeight(int n, double j = J)
    {
        var block = SingleExcitationBlock(n, j);
        double worst = 0.0;
        for (int i = 0; i < n; i++)
            worst = Math.Max(worst, block[i, i].Magnitude);
        return worst;
    }

    /// <summary>The low-N band-edge coincidence 2cos(π/(N+1)): equal to Q*(N) at N=2,3 only
    /// (1 = 2cos60°, √2 = 2cos45°), then departing (Q*(4)=1.8785 ≠ φ=1.618).</summary>
    public static double BandEdgeCoincidence(int n) => 2.0 * Math.Cos(Math.PI / (n + 1));

    /// <summary>The rigidity below which a gap mode counts as coalescing (an EP).</summary>
    private const double REpThreshold = 0.05;

    /// <summary>The non-zero modes within a small band of the slowest decay rate (the gap), with their
    /// phase rigidity, built at Q = J/γ on the live Liouvillian.</summary>
    private static List<PhaseRigidity.Mode> GapModes(int n, double gamma)
    {
        var L = new ChainSystem(n, J, gamma).BuildLiouvillian();
        var nz = PhaseRigidity.Compute(L).Where(m => m.Lambda.Real < -1e-6).ToList();
        double gap = nz.Max(m => m.Lambda.Real);
        return nz.Where(m => m.Lambda.Real > gap - 0.15).ToList();
    }

    /// <summary>At Q = Q*(1+δ), δ = <see cref="EpReadDelta"/>, just above the EP on the full Liouvillian:
    /// the coalescing near-gap mode (minimum phase rigidity, the {0,2}-coherence whose r → 0) with its
    /// n_diff histogram, and the co-located band-edge survivor (Im ≈ 2cos(π/(N+1)), r ≈ 1). The instrument
    /// that distinguishes the EP (erasure) from the crossing (survival). Read above the EP rather than at
    /// the handover, because at N≥4 the handover (<see cref="Horizon"/>) lies below the EP, where the pair has
    /// already split into two real branches. ReadAtQ is the Q used, so γ = J/ReadAtQ.</summary>
    public (PhaseRigidity.Mode Coalescer,
            IReadOnlyDictionary<int, double> CoalescerHist,
            double CoalescerMeanNDiff,
            PhaseRigidity.Mode BandEdge,
            double BandEdgeR,
            double ReadAtQ) EpModes(int n)
    {
        double q = EpQ(n) * (1.0 + EpReadDelta);
        var gapModes = GapModes(n, J / q);
        var coalescer = gapModes.OrderBy(m => m.Rigidity).First();
        var (mean, hist) = LiouvilleOperatorContent.NDiffHistogram(coalescer.Right, n);
        double bandIm = 2.0 * Math.Cos(Math.PI / (n + 1));
        var bandEdge = gapModes.OrderBy(m => Math.Abs(Math.Abs(m.Lambda.Imaginary) - bandIm)).First();
        return (coalescer, hist, mean, bandEdge, bandEdge.Rigidity, q);
    }

    /// <summary>√-scaling certificate of a 2nd-order EP: Im²/(Q−Q*) for the small-Im coalescer branch
    /// at Q = Q*(n)·(1+delta). Constant across deltas ⟹ Im ∝ √(Q−Q*) ⟹ a clean 2-dim EP (not a
    /// cluster). Returns NaN if no coalescer branch is resolved above Q*.</summary>
    public double SqrtScalingRatio(int n, double delta)
    {
        double qStar = EpQ(n);
        double q = qStar * (1.0 + delta);
        double bandIm = 2.0 * Math.Cos(Math.PI / (n + 1));
        var branches = GapModes(n, J / q)
            .Where(m => Math.Abs(m.Lambda.Imaginary) > 1e-6 && Math.Abs(m.Lambda.Imaginary) < bandIm - 0.2)
            .OrderBy(m => Math.Abs(m.Lambda.Imaginary))
            .ToList();
        if (branches.Count == 0) return double.NaN;
        double im = Math.Abs(branches[0].Lambda.Imaginary);
        return im * im / Math.Abs(q - qStar);
    }

    // ---- the single-excitation block: the EP, the handover, and the excess light ----

    /// <summary>The largest N the SE-block readings are offered for: at δ = <see cref="EpReadDelta"/> the
    /// coalescer branch is resolved by the picker up to here (an N²-dim eigenproblem, instant).</summary>
    public const int SeBlockMaxN = 8;

    /// <summary>The offset above the EP at which the coalescer is read, Q = Q*·(1+δ): small enough that the
    /// pair is still the near-degenerate one, large enough that its two branches are resolved.</summary>
    public const double EpReadDelta = 0.001;

    private readonly Dictionary<int, double> _epCache = new();
    private readonly Dictionary<int, double> _handoverCache = new();

    private static void CheckSeBlockN(int n)
    {
        if (n < 2 || n > SeBlockMaxN)
            throw new ArgumentOutOfRangeException(nameof(n), $"the SE-block readings need N in 2..{SeBlockMaxN}; got {n}");
    }

    /// <summary>The EP Q*(N) on the single-excitation block (<see cref="EpCharacterWitness.BisectEpQ"/>): the
    /// Q below which the slowest SE mode no longer oscillates. The same criterion as <see cref="Horizon"/> on
    /// a spectrum without the (0,1) survivor, so it lands on the EP itself; equal to the Python
    /// <c>qstar_se</c> to bisection precision. Cached per N.</summary>
    public double EpQ(int n)
    {
        CheckSeBlockN(n);
        if (_epCache.TryGetValue(n, out var q)) return q;
        q = EpCharacterWitness.BisectEpQ(n);
        _epCache[n] = q;
        return q;
    }

    /// <summary>The handover Q_h(N) on the single-excitation block: the Q at which the slowest real mode reaches
    /// the floor Re = −2γ (the D06 gap-threshold definition; up to N=8 that mode is the darker of the two real
    /// branches the pair splits into below Q*). Below Q_h that branch is slower than
    /// the floor and is the gap (the clock frozen); above it the floor's own survivor is the gap (the clock
    /// turning). Equal to the EP exactly where the EP is on the floor (N=2,3, w2 = ½); below it from N=4 by
    /// ((2w2−1)/c)² to leading order, c the square-root splitting coefficient. Bisected in [0.8·Q*, Q*], 60
    /// steps. Cached per N.</summary>
    public double HandoverQ(int n)
    {
        CheckSeBlockN(n);
        if (_handoverCache.TryGetValue(n, out var q)) return q;
        double qs = EpQ(n);
        double lo = 0.8 * qs, hi = qs;
        for (int it = 0; it < 60; it++)
        {
            double mid = 0.5 * (lo + hi);
            if (SlowestRealModeInFloorUnits(n, mid) < 1.0) lo = mid;   // a real mode above the floor: below Q_h
            else hi = mid;
        }
        q = 0.5 * (lo + hi);
        _handoverCache[n] = q;
        return q;
    }

    /// <summary>The slowest real SE mode's rate in floor units, Re/(−2γ), at Q (NaN if none is real):
    /// &lt; 1 means a real mode sits above the floor and is the gap.</summary>
    private static double SlowestRealModeInFloorUnits(int n, double q)
    {
        double g = J / q;
        var real = EpCharacterWitness.Lse(n, J, g).Evd().EigenValues
            .Where(e => e.Real < -1e-7 && Math.Abs(e.Imaginary) < 1e-9).Select(e => e.Real).ToArray();
        return real.Length == 0 ? double.NaN : real.Max() / (-2.0 * g);
    }

    /// <summary>The coalescer read at Q = Q*(1+δ) on the SE block: its rate in floor units Re/(−2γ), its
    /// coherence share w2 (the n_diff = 2 weight; populations carry n_diff = 0), the Absorption-Theorem
    /// residual Re/(−2γ) − 2·w2 (exactly 0 for every eigenvector: Re λ = −2γ⟨n_diff⟩), and |Im|. On the
    /// floor ⟺ w2 = ½. EigenResidual = ‖Lv − λv‖/‖v‖, the eigensolver's own error on the pair, which bounds
    /// the identity's residual: |Re λ/(−2γ) − ⟨n_diff⟩| ≤ EigenResidual/(2γ), since v†Lv = λ‖v‖² + v†(Lv − λv)
    /// and Re v†L_Hv = 0. Rigidity = |l†v|/(‖l‖‖v‖) with l the biorthogonal left vector (the phase rigidity,
    /// the Petermann form): on the 2×2 at N=2,3 it is √(2δ)·(1 − O(δ)) exactly, and the dressed pair from N=4
    /// sits below that, the deficit growing with N like the excess light. The picker is the Python census's:
    /// among the upper-half-plane modes within 0.25·J of the slowest oscillating rate, the one with the
    /// smallest |Im|.</summary>
    public (double ReOver2G, double W2, double AbsorptionResidual, double AbsIm, double EigenResidual, double Rigidity) CoalescerAtEp(int n, double delta = EpReadDelta)
    {
        CheckSeBlockN(n);
        double q = EpQ(n) * (1.0 + delta), g = J / q;
        var lse = EpCharacterWitness.Lse(n, J, g);
        var evd = lse.Evd();
        var lam = evd.EigenValues;
        var vecs = evd.EigenVectors;
        int dim = n * n;
        var osc = Enumerable.Range(0, dim).Where(k => lam[k].Real < -1e-7 && lam[k].Imaginary > 1e-9).ToList();
        if (osc.Count == 0)
            throw new InvalidOperationException($"N={n}: no oscillating SE mode at Q = Q*(1+{delta.ToString("0.###", Inv)})");
        double reMax = osc.Max(k => lam[k].Real);
        int pick = osc.Where(k => lam[k].Real > reMax - 0.25 * J).OrderBy(k => Math.Abs(lam[k].Imaginary)).First();
        double total = 0.0, w2 = 0.0;
        for (int i = 0; i < n; i++)
            for (int jj = 0; jj < n; jj++)
            {
                double m = vecs[i * n + jj, pick].Magnitude;
                total += m * m;
                if (i != jj) w2 += m * m;
            }
        w2 /= total;
        double re = lam[pick].Real / (-2.0 * g);
        var v = vecs.Column(pick);
        double eigenResidual = (lse * v - v * lam[pick]).L2Norm() / v.L2Norm();
        var left = vecs.Inverse().Row(pick).Conjugate();
        double rigidity = left.ConjugateDotProduct(v).Magnitude / (left.L2Norm() * v.L2Norm());
        return (re, w2, re - 2.0 * w2, Math.Abs(lam[pick].Imaginary), eigenResidual, rigidity);
    }

    /// <summary>The 2×2 closure residual, entry-wise and exact (no eigensolver): with h the single-particle
    /// hopping, P = |0⟩⟨0| − |N−1⟩⟨N−1| (a population direction, n_diff = 0) and C = [h, P]/J (a coherence
    /// direction, n_diff = 2), the plane span{P, C} is L-invariant iff [h, C] = c₂·J·P with c₂ = 4 (N=2) or
    /// 2 (N≥3). Returns max|[h, C] − c₂·J·P|: exactly 0.0 at N=2,3, where the two ends' leaks cancel on the
    /// shared middle site and the pair is the 2×2 λ² + 4γλ + c₂J² = 0 (so w2 = ½ identically and the EP sits
    /// on the floor); 2.0 from N=4 on, where the leak lands on interior sites and no γ-independent plane of a
    /// population and a coherence direction exists.</summary>
    public static double ClosureResidual(int n)
    {
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(n));
        var h = EpCharacterWitness.HSingle(n, J);
        var pMat = new double[n, n];
        pMat[0, 0] = 1.0; pMat[n - 1, n - 1] = -1.0;
        double[,] Comm(double[,] a, double[,] b)
        {
            var r = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int k = 0; k < n; k++)
                {
                    double acc = 0.0;
                    for (int m = 0; m < n; m++) acc += a[i, m] * b[m, k] - b[i, m] * a[m, k];
                    r[i, k] = acc;
                }
            return r;
        }
        var c = Comm(h, pMat);                       // [h, P] = J·C
        for (int i = 0; i < n; i++) for (int k = 0; k < n; k++) c[i, k] /= J;
        var hc = Comm(h, c);                         // [h, C]
        double c2 = n == 2 ? 4.0 : 2.0;
        double worst = 0.0;
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++)
                worst = Math.Max(worst, Math.Abs(hc[i, k] - c2 * J * pMat[i, k]));
        return worst;
    }

    /// <summary>The resolution of <see cref="Horizon"/> in Q: the γ bracket over 2^steps, mapped through
    /// Q = J/γ (dQ = J·dγ/γ*²). The instrument's own law, used to gate its agreement with the SE-block
    /// <see cref="HandoverQ"/>.</summary>
    public double HorizonResolution(int n)
    {
        double gStar = J / Horizon(n);
        return J * (GammaHi - GammaLo) / Math.Pow(2.0, BisectionSteps) / (gStar * gStar);
    }

    /// <summary>the excess light: the coalescer is on the floor Re = −2γ only where its coherence share is
    /// exactly ½ (N=2,3, the clean 2×2); from N=4 its light exceeds 1 by 2w2−1, it sits below the floor by
    /// 2γ(2w2−1) at every Q ≥ Q*, and the handover Q_h sits below the EP by the trace dressing ((2w2−1)/c)².
    /// Read on the SE block; the Absorption identity is exact in every row; the full-L handover cross-docks
    /// to Q_h within its resolution at N=2..5.</summary>
    private InspectableNode TheExcessLight()
    {
        var rows = new List<IInspectable>();
        for (int n = 2; n <= SeBlockMaxN; n++)
        {
            double qs = EpQ(n), qh = HandoverQ(n);
            var c = CoalescerAtEp(n);
            string cross = n <= 5
                ? $"; full-L handover Horizon({n}) = {Horizon(n).ToString("0.00000", Inv)}, = Q_h within its resolution " +
                  $"{HorizonResolution(n).ToString("0.0e0", Inv)}"
                : "";
            rows.Add(new InspectableNode($"N={n}: Q* = {qs.ToString("0.0000000", Inv)}, Q_h = {qh.ToString("0.0000000", Inv)}",
                summary: $"Q* − Q_h = {(qs - qh).ToString("0.000e0", Inv)}; the coalescer at Q*(1+{EpReadDelta.ToString("0.###", Inv)}): " +
                         $"Re/(−2γ) = {c.ReOver2G.ToString("0.000000", Inv)}, w2 = {c.W2.ToString("0.000000", Inv)}, " +
                         $"2w2 − 1 = {(2.0 * c.W2 - 1.0).ToString("0.000000", Inv)}, Absorption residual Re/(−2γ) − 2w2 = " +
                         $"{c.AbsorptionResidual.ToString("0.0e0", Inv)} (bound ‖Lv−λv‖/(2γ‖v‖) = {(c.EigenResidual / (2.0 * J / (qs * (1.0 + EpReadDelta)))).ToString("0.0e0", Inv)}), " +
                         $"|Im| = {c.AbsIm.ToString("0.000000", Inv)}, r/√(2δ) = {(c.Rigidity / Math.Sqrt(2.0 * EpReadDelta)).ToString("0.0000", Inv)}; " +
                         $"2×2 closure residual max|[h,C] − c₂JP| = {ClosureResidual(n).ToString("0.0", Inv)}" +
                         (n <= 3 ? " (exact: the pair IS the 2×2, w2 = ½, r = √(2δ), the EP on the floor, Q_h = Q*)" : " (no invariant 2×2; below the floor)") + cross));
        }
        return new InspectableNode("the excess light (the EP is on the floor at N=2,3 only)",
            summary: "Absorption Theorem: Re λ = −2γ⟨n_diff⟩ for every eigenvector, so on the single-excitation block " +
                     "Re = −4γ·w2 with w2 the coalescer's coherence share, and it sits on the floor Re = −2γ ⟺ w2 = ½. " +
                     "At N=2,3 the pair is a clean 2×2 (P = |0⟩⟨0| − |N−1⟩⟨N−1| and C = [h,P]/J close under h entry-wise at " +
                     "0.0, so span{P,C} is L-invariant with λ² + 4γλ + c₂J² = 0), w2 = ½ identically, the EP is on the floor and " +
                     "the handover IS the EP (Q_h = Q*, c = √(2/Q*) = √2 and 2^¼ exactly, r = √(2δ)). Why exactly there: at an EP the " +
                     "pair's light is its Jordan plane's mean light (2λ_EP = tr L|_V, the trace-midpoint 2×2 of PROOF_F86B_UNIVERSAL_SHAPE), " +
                     "and a plane of one population direction (light 0) and one coherence direction (light 2) has mean 1, the midpoint of " +
                     "two rungs, not a coincidence. From N=4 the coalescer is the {0,2}-coherence " +
                     "by support but its weights are not ½/½: it sits below the floor by 2γ(2w2−1) at every Q ≥ Q*, closest at " +
                     "the EP, its rate walking from −2γ toward the −4γ centre of the resummed dispersion λ²+8γλ+4J²q². Below Q* " +
                     "the pair is two real branches, Re/(−2γ) = 2w2 ∓ c√(Q*−Q) + O(Q*−Q); the darker one reaches the floor at " +
                     "Q_h = Q* − ((2w2−1)/c)²·(1 + O(2w2−1)), where it becomes the gap and the band edge takes the clock over (with " +
                     "the excess light and c both read in the limit δ → 0 the leading order reproduces the bisected gap to " +
                     "0.04% / 0.23% / 0.66% / 1.3% / 2.2% at N=4..8; read at the δ = 0.001 these rows print it is 0.4% .. 2.6%; the " +
                     "deficit is the next order, a negative linear term in the pair's centre). The crossing is exact, never avoided: " +
                     "(−1)^{n_XY} commutes with L (PROOF_ABSORPTION_THEOREM), so the (1,1) pair and the (0,1) survivor never hybridize. So the handover this witness bisects " +
                     "on the full L (Horizon), IncompletenessSurvivorWitness.HandoverQ and the glossary's XY Q*_gap are one event, " +
                     "the handover, and the EP is a second event from N=4 on (1.87874 vs 1.87854 at N=4, 2.37367 vs 2.37217 at N=5).",
            children: rows);
    }

    public string DisplayName =>
        $"CoherenceHorizonWitness (Q*(N) live, J={J.ToString("0.#", Inv)}, ε={OmegaEps.ToString("0.###", Inv)})";

    public string Summary =>
        "the coherence horizon Q*(N) live (typed-home intent: docs/ANALYTICAL_FORMULAS.md F2b corollary " +
        "\"Coherence horizon Q*(N)\"). Two thresholds under one letter: the EP Q*(N), where the single-excitation " +
        $"pair coalesces (SE block: Q*(2)={EpQ(2).ToString("0.#####", Inv)}, Q*(3)={EpQ(3).ToString("0.#####", Inv)}, " +
        $"Q*(4)={EpQ(4).ToString("0.#####", Inv)}, Q*(5)={EpQ(5).ToString("0.#####", Inv)}), and the handover Q_h(N), the " +
        "clock's takeover: the Q below which the slowest non-zero mode of the FULL Liouvillian stops oscillating (bisected on " +
        $"Symphony.Clock.Omega at J=1: {Horizon(2).ToString("0.#####", Inv)}, {Horizon(3).ToString("0.#####", Inv)}, " +
        $"{Horizon(4).ToString("0.#####", Inv)}, {Horizon(5).ToString("0.#####", Inv)}). Equal at N=2,3, where the EP sits on " +
        "the floor Re = −2γ; apart from N=4 by ((2w2−1)/c)², the square of the EP's excess light over the split " +
        "coefficient (see the excess-light node). Both are readings of the dephasing side; what the Hückel " +
        $"matrix supplies is the Hamiltonian, exactly (HuckelResidual(5) = {HuckelResidual(5).ToString("0.0#####", Inv)}, " +
        "see the inheritance node). " +
        "Sector overview: inspect --root blockspectrum (this zooms the (1,1) single-excitation {0,2} sector — " +
        "the low-Q EP regime; --root ceiling reads its high-Q regime).";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return TheLadder();
            yield return TheHuckelInheritance();
            yield return TheEpBase();
            yield return TheBandEdgeCoincidence();
            yield return TheEpVerdict();
            yield return TheExcessLight();
        }
    }

    /// <summary>the ladder: the computed handover Q_h(N) for N=2..5, live.</summary>
    private InspectableNode TheLadder()
    {
        var rungs = new List<IInspectable>();
        foreach (int n in new[] { 2, 3, 4, 5 })
            rungs.Add(InspectableNode.RealScalar($"Q_h({n})", Horizon(n), "0.#####"));
        return new InspectableNode("the ladder",
            summary: $"the computed handover Q_h(N) (live, bisected γ on Symphony.Clock.Omega at J=1): " +
                     $"N=2 → {Horizon(2).ToString("0.#####", Inv)}, N=3 → {Horizon(3).ToString("0.#####", Inv)} (√2), " +
                     $"N=4 → {Horizon(4).ToString("0.#####", Inv)}, N=5 → {Horizon(5).ToString("0.#####", Inv)}, each to " +
                     "this bisection's own resolution (the resolution node reads it per N). The SE-block reading of the " +
                     "same event is 1 / √2 / 1.87854 / 2.37217, and the EP a second event from N=4 (1.87874 / 2.37367).",
            children: rungs);
    }

    /// <summary>what the Hückel matrix supplies here: the Hamiltonian, exactly. The residual is 0.0.</summary>
    private InspectableNode TheHuckelInheritance()
    {
        var rows = new List<IInspectable>();
        foreach (int n in new[] { 2, 3, 4, 5, 6 })
            rows.Add(InspectableNode.RealScalar($"max |H_SE − J·A| at N={n}", HuckelResidual(n), "0.0#####"));
        return new InspectableNode("what the Hückel matrix supplies (exact)",
            summary: "the single-excitation block of H in the site basis IS the Hückel matrix at α=0, β=J, entry for " +
                     "entry: Jordan-Wigner sends (J/2)Σ(X_lX_{l+1} + Y_lY_{l+1}) onto the tridiagonal hopping matrix " +
                     $"itself, so the residual is exactly {HuckelResidual(5).ToString("0.0#####", Inv)} and not merely " +
                     "small. That is the whole of the inheritance and it is a statement about the HAMILTONIAN. It stops " +
                     "there: both Q_h and the EP are readings of the dephasing side, and a Hückel or Frost construction " +
                     "is a static spectrum with no γ, no bath and no time in it, so there is no threshold on that side to " +
                     "compare a ladder against. Anything that reads as such a comparison has this witness's own numbers " +
                     "on both sides of it. The gate can fail: a ZZ term, a non-uniform J or a longitudinal field each " +
                     "put a non-zero entry where the adjacency matrix has none.",
            children: rows);
    }

    /// <summary>the EP base (N=2): Q*(2)=1 is the exceptional point itself.</summary>
    private InspectableNode TheEpBase()
    {
        double q2 = Horizon(2);
        return new InspectableNode("the EP base (N=2)",
            summary: $"Q*(2) = {q2.ToString("0.#####", Inv)} is the exceptional point itself (γ=J, the two clocks " +
                     "merge into one). Mechanism note: at N=2 the crossover the bisection lands is NOT the pulled " +
                     "coherence hand 2√(J²−γ²) freezing. Symphony.Clock.Omega(2) reads the ±J band mode (Omega=1) right " +
                     "down to γ=J, then drops to 0, because that band mode ceases to be the slowest-decay (gap) mode " +
                     "exactly at γ=J, which coincides with the exceptional point. So Q*(2)=1 marks the ±J band mode " +
                     "ceasing to be the gap mode, coinciding with the EP, distinct from the N≥3 band-edge-vs-overdamped " +
                     "crossover (see ClockHandLadderWitness for the N=2 raw-clock subtlety).");
    }

    /// <summary>the band-edge coincidence: Q*(N) = 2cos(π/(N+1)) at N=2,3 ONLY, a low-N accident.</summary>
    private InspectableNode TheBandEdgeCoincidence()
    {
        double be3 = BandEdgeCoincidence(3);
        double phi = BandEdgeCoincidence(4); // 2cos(π/5) = φ, the golden ratio
        return new InspectableNode("the band-edge coincidence",
            summary: $"Q*(N) = 2cos(π/(N+1)) at N=2,3 ONLY: 1 = 2cos60° and √2 = 2cos45° = " +
                     $"{be3.ToString("0.#####", Inv)} match Q*(2)/Q*(3). It is a low-N accident, departing at N≥4: " +
                     $"2cos(π/5) = φ = {phi.ToString("0.#####", Inv)} but Q*(4) = {EpQ(4).ToString("0.#####", Inv)} " +
                     $"≠ φ. This is why √2 looked exact at " +
                     "N=3 while the rest awaited a clean form: the band edge is the right answer only where it happens to " +
                     "coincide with the horizon.");
    }

    /// <summary>the EP verdict, recomputed live at Q = Q*(1+δ): per N=2..5 the coalescing {0,2}-coherence
    /// (r → 0, weight on n_diff ∈ {0,2}, ½/½ at N=2,3 only) and the co-located band-edge survivor (r ≈ 1).
    /// No bifurcation at N=4: the {0,2}-coherence is the freezer at every N, the band edge the γ-protected
    /// survivor; they share the floor Re = −2γ at N=2,3 only (Absorption Theorem: the survivor has
    /// ⟨n_diff⟩ = 1 exactly, the coalescer 2w2 &gt; 1 from N=4; at N=2,3 the two share the floor).</summary>
    private InspectableNode TheEpVerdict()
    {
        var rungs = new List<IInspectable>();
        foreach (int n in new[] { 2, 3, 4, 5 })
        {
            var ep = EpModes(n);
            string h0 = ep.CoalescerHist.GetValueOrDefault(0).ToString("0.0000", Inv);
            string h2 = ep.CoalescerHist.GetValueOrDefault(2).ToString("0.0000", Inv);
            double reOver2G = ep.Coalescer.Lambda.Real / (-2.0 * J / ep.ReadAtQ);
            bool isEp = ep.Coalescer.Rigidity < REpThreshold;
            string label = isEp ? "EP" : "crossing";
            string verdict = isEp ? "genuine EP, the {0,2}-coherence coalesces" : "no EP (a crossing)";
            rungs.Add(new InspectableNode($"N={n}: {label}",
                summary: $"coalescer r = {ep.Coalescer.Rigidity.ToString("0.000", Inv)} " +
                         $"(Im = {ep.Coalescer.Lambda.Imaginary.ToString("0.000", Inv)}, hist {{0:{h0}, 2:{h2}}}, " +
                         $"mean n_diff = {ep.CoalescerMeanNDiff.ToString("0.0000", Inv)} = Re/(−2γ) = {reOver2G.ToString("0.0000", Inv)}" +
                         $"{(n <= 3 ? ", on the floor" : ", below the floor")}) → {verdict}; " +
                         $"band edge Im = {Math.Abs(ep.BandEdge.Lambda.Imaginary).ToString("0.000", Inv)} " +
                         $"(2cos(π/(N+1))), r = {ep.BandEdgeR.ToString("0.000", Inv)} → the γ-protected survivor."));
        }
        double ratio = SqrtScalingRatio(4, 0.03);
        return new InspectableNode("the EP verdict (live phase rigidity)",
            summary: "the mode that coalesces at Q*(N) is the {0,2}-coherence (population/antisymmetric block, all " +
                     "of its weight on n_diff ∈ {0,2}) at ALL N=2..5, a genuine square-root EP (phase rigidity r → 0). " +
                     "NO sector bifurcation at N=4: the band edge 2cos(π/(N+1)) is the co-located γ-protected SURVIVOR " +
                     "(r ≈ 1). The two share the floor Re = −2γ at N=2,3 ONLY, where the clean 2×2 forces the " +
                     "coalescer's coherence share to exactly ½ (Absorption Theorem: Re λ = −2γ⟨n_diff⟩, so on the floor " +
                     "⟺ w2 = ½); from N=4 the coalescer sits BELOW the floor by 2γ(2w2−1) at every Q ≥ Q*, the gap there is " +
                     "the survivor's alone, and the handover this witness bisects (the clock's takeover, Q_h) sits below the EP by " +
                     "((2w2−1)/c)² to leading order (the excess-light node). So Q*(N) is at once a {0,2}-coherence EP (the erasure point, " +
                     "which climbs the ladder) and a band-edge crossing (the clock survives). √-scaling Im²/(Q−Q*) at " +
                     $"N=4 = {ratio.ToString("0.00", Inv)} (constant ⟹ a clean 2nd-order EP). The closed form of the EP: " +
                     "1 and √2 exactly at N=2,3, transcendental from N=4 (the SE-block discriminant), slope 2/π derived. " +
                     "Recomputed live via PhaseRigidity (the eig instrument, the F86a-misfire-prone family); supersedes " +
                     "the earlier narrated 'bifurcation at N=4'. The ARTIFACT-FREE confirmation (Riesz ‖P‖ / departure-" +
                     "from-normality / geo-vs-alg, no eig eigenvector) lives at inspect --root epcharacter " +
                     "(EpCharacterWitness): DEFECTIVE at every N=2..5, dep≈4, geo 1 < alg 2.",
            children: rungs);
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
