using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.States;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The zoom-out, the first movement: ONE open quantum system, ONE time evolution, every
/// applicable lens reading the SAME trajectory on a SHARED timeline, plus a cross-lens event axis
/// so the interplay (das Zusammenspiel) becomes visible instead of 122 separate maximum-zoom windows.
///
/// <para>Where the Object Manager's other roots each catalog one zoom, the Symphony catalogs the
/// <i>zoom-out</i>: it builds the chain's Liouvillian L = −i[H, ·] + Σ_l γ_l(Z_l ρ Z_l − ρ) once,
/// eigendecomposes it once, and unrolls the chosen initial state's density matrix ρ(t) over the
/// time grid exactly once. The lenses do not each re-evolve; they all read the one
/// <see cref="States"/> list. The mathematical license is reflections/ON_THE_ONE_DIAGONAL.md: the
/// rates, the mirror, and the verdict are one diagonal read three ways, so it is honest to hold one
/// trajectory and let every reading take its slice of it.</para>
///
/// <para><b>The lenses wired in this first movement</b> (each a child node, all over the one ρ(t)):
/// <list type="bullet">
///   <item><b>palindrome</b> — the spectrum's paired structure (F1), read off the same
///         <see cref="CarrierVectorPortfolio"/> decomposition <see cref="MirrorSystem"/> uses:
///         every decay rate r pairs with 2σ − r.</item>
///   <item><b>quarter (CΨ)</b> — CΨ(t) = purity × normalized coherence along the trajectory
///         (C = Tr(ρ²) = Σ|ρ_ij|², Ψ = ℓ₁-coherence/(d−1), the repo's canonical convention where
///         Bell+ has CΨ(0) = 1/3); curve payload, ¼-crossing times.</item>
///   <item><b>dose (K)</b> — the dimensionless dose K = γ·t along the same grid; the dose of the
///         fold is K at the first ¼ crossing.</item>
///   <item><b>light</b> — the state's light content over time, the purity-weighted mean
///         popcount(i⊕j) over the coherences |i⟩⟨j| (popcount(i⊕j) = the number of sites where ket
///         and bra differ = the Absorption-Theorem light channel); curve payload.</item>
/// </list></para>
///
/// <para><b>events</b> is the point of the movement: the cross-lens axis, every detected event from
/// every lens merged and sorted on the one timeline (a ¼ crossing, K reaching a milestone, a light
/// extremum), so coincidences between lenses show up as adjacency on a single axis.</para>
///
/// <para>Diagnostic scale: the full dense d²×d² Liouvillian is materialised and eigendecomposed, so
/// N is capped at <see cref="MaxN"/> (= 5, d²=1024). More instruments (more lenses) can join later as
/// further child nodes reading the same <see cref="States"/>; the engine is built to hand the one
/// trajectory to all of them.</para></summary>
public sealed class Symphony : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>The largest N the dense d²×d² Liouvillian eigendecomposition will run (d²=1024 at N=5).</summary>
    public const int MaxN = 5;

    public int N { get; }
    public double J { get; }
    public double Gamma { get; }
    public HamiltonianType HType { get; }
    public TopologyKind Topology { get; }
    public InitialStateKind InitialState { get; }

    /// <summary>The two sites the local-CΨ lens reduces ρ(t) onto (0-indexed, site 0 = MSB).
    /// Default (0,1) is Bell+'s carrier pair, where the lens reproduces F25 at N=2.</summary>
    public (int Site1, int Site2) CarrierPair { get; }

    public double TMax { get; }
    public int TPoints { get; }

    /// <summary>The painters' movement: if set, Symphony plays the piece a second (defected) and a
    /// guard third (δJ/2) time on a single-bond J-defect at this bond index, and wires the PTF lenses.
    /// Null = no movement, Symphony unchanged.</summary>
    public int? DefectBond { get; }

    /// <summary>The defect strength δJ (added to J on the defect bond) for the painters' movement.</summary>
    public double DeltaJ { get; }

    /// <summary>Hilbert dimension 2^N.</summary>
    public int Dim => 1 << N;

    /// <summary>σ = N·γ, the total dephasing and the palindrome centre.</summary>
    public double Sigma => N * Gamma;

    public Symphony(
        int n = 3,
        double j = 1.0,
        double gamma = 0.1,
        HamiltonianType hType = HamiltonianType.XY,
        TopologyKind topology = TopologyKind.Chain,
        InitialStateKind initialState = InitialStateKind.BellPair,
        double tMax = double.NaN,
        int tPoints = 60,
        int? defectBond = null,
        double deltaJ = 0.02,
        (int Site1, int Site2)? carrierPair = null)
    {
        if (n < 2 || n > MaxN)
            throw new ArgumentOutOfRangeException(nameof(n),
                $"Symphony needs N in 2..{MaxN} (dense d²×d² Liouvillian eigendecomposition); got {n}");
        if (gamma <= 0) throw new ArgumentOutOfRangeException(nameof(gamma), $"gamma must be positive; got {gamma}");
        if (tPoints < 2) throw new ArgumentOutOfRangeException(nameof(tPoints), $"t-points must be >= 2; got {tPoints}");

        N = n;
        J = j;
        Gamma = gamma;
        HType = hType;
        Topology = topology;
        InitialState = initialState;
        TPoints = tPoints;
        // Default window: long enough that the slow (carrier 4γ) coherence folds well past ¼ —
        // 1/γ is several carrier times; runs in milliseconds at these N.
        TMax = double.IsNaN(tMax) ? 1.0 / gamma : tMax;
        if (TMax <= 0) throw new ArgumentOutOfRangeException(nameof(tMax), $"t-max must be positive; got {TMax}");

        DefectBond = defectBond;
        DeltaJ = deltaJ;

        var pair = carrierPair ?? (0, 1);
        if (pair.Site1 < 0 || pair.Site1 >= N || pair.Site2 < 0 || pair.Site2 >= N)
            throw new ArgumentOutOfRangeException(nameof(carrierPair),
                $"carrier pair sites must be in [0, {N - 1}]; got ({pair.Site1}, {pair.Site2})");
        if (pair.Site1 == pair.Site2)
            throw new ArgumentException(
                $"carrier pair sites must be distinct; got ({pair.Site1}, {pair.Site2})", nameof(carrierPair));
        CarrierPair = pair;
    }

    // ---- The one evolution, computed once, shared by all lenses ----

    private double[]? _tGrid;
    private ComplexMatrix[]? _states;     // ρ(t) at each grid point, vec convention ρ[a,b] = vec[a*d+b]
    private Complex[]? _lambdaA;          // L_A eigenvalues, kept for the global clock (Takt + ω_mem)
    private int _evolveCount;             // how many times the trajectory was actually built (must stay 1)

    /// <summary>The L_A (clean) Liouvillian eigenvalues, available after the one evolution; the global
    /// clock (Takt gap, ω_mem) is read off these by the painters' movement.</summary>
    public IReadOnlyList<Complex> LiouvillianEigenvalues { get { EnsureEvolved(); return _lambdaA!; } }

    /// <summary>Grid-resolution fitness for the envelope lenses, from the shared L spectrum:
    /// ω = max |Im λ| (fastest coherent oscillation), samples per oscillation 2π/(ωΔt), and the
    /// conservative peak-clip floor ½·(ωΔt)²·peakScale (the raw peak-height clip; with parabolic
    /// apex the true residual is far smaller, so this over-estimates — a safe floor to surface).
    /// ω ≈ 0 (pure dephasing) ⟹ SamplesPerOscillation = +∞ and PeakClipFloor = 0.</summary>
    public (double Omega, double SamplesPerOscillation, double PeakClipFloor) GridFitness(double peakScale)
    {
        EnsureEvolved();
        double omega = 0.0;
        foreach (var e in _lambdaA!) omega = Math.Max(omega, Math.Abs(e.Imaginary));
        double dt = TPoints > 1 ? TMax / (TPoints - 1) : TMax;
        double wdt = omega * dt;
        double samples = wdt > 1e-12 ? 2.0 * Math.PI / wdt : double.PositiveInfinity;
        double floor = 0.5 * wdt * wdt * peakScale;
        return (omega, samples, floor);
    }

    /// <summary>The shared timeline: t = 0, …, TMax in <see cref="TPoints"/> equal steps.</summary>
    public IReadOnlyList<double> TimeGrid { get { EnsureEvolved(); return _tGrid!; } }

    /// <summary>The one trajectory: ρ(t) at every grid point, built once and handed to every lens.</summary>
    public IReadOnlyList<ComplexMatrix> States { get { EnsureEvolved(); return _states!; } }

    /// <summary>How many times the trajectory was built. A test guards this stays exactly 1, no matter
    /// how many lenses or how often their summaries are read: ONE evolution, many lenses.</summary>
    public int EvolveCount => _evolveCount;

    private void EnsureEvolved()
    {
        if (_states is not null) return;

        var chain = new ChainSystem(N, J, Gamma, HType, Topology);
        var H = chain.BuildHamiltonian();
        var gammas = Enumerable.Repeat(Gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammas);   // d²×d² superoperator, vec[a*d+b] = ρ[a,b]

        int d = Dim;
        int dim2 = d * d;

        // vec(ρ(t)) = R · exp(Λ t) · R⁻¹ · vec(ρ(0)). One eigendecomposition, then a mat-vec per step.
        var evd = L.Evd();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();
        var lambda = evd.EigenValues;

        var rho0 = BuildInitialState(d);
        var vec0 = ComplexVector.Build.Dense(dim2);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                vec0[a * d + b] = rho0[a, b];
        var c0 = Rinv * vec0;

        var tGrid = new double[TPoints];
        for (int i = 0; i < TPoints; i++) tGrid[i] = TMax * i / (TPoints - 1);

        var states = new ComplexMatrix[TPoints];
        var expDiag = ComplexVector.Build.Dense(dim2);
        for (int s = 0; s < TPoints; s++)
        {
            double t = tGrid[s];
            for (int i = 0; i < dim2; i++) expDiag[i] = Complex.Exp(lambda[i] * t) * c0[i];
            var vecT = R * expDiag;
            var rhoT = Matrix<Complex>.Build.Dense(d, d);
            for (int a = 0; a < d; a++)
                for (int b = 0; b < d; b++)
                    rhoT[a, b] = vecT[a * d + b];
            states[s] = rhoT;
        }

        _tGrid = tGrid;
        _states = states;
        _lambdaA = lambda.ToArray();
        _evolveCount++;
    }

    /// <summary>The chosen initial density matrix ρ(0). Bell+ is embedded on sites (0,1) with the
    /// remaining sites in |0⟩: it is the repo's hardware-confirmed ¼-crossing carrier (F25, CΨ(0)=1/3)
    /// and at N=2 reproduces the F25 closed form CΨ(t)=f(1+f²)/6, f=e^{−4γt} exactly. SingleExcitation
    /// is |10…0⟩, one excitation on site 0.</summary>
    public ComplexMatrix BuildInitialState(int d)
    {
        var psi = ComplexVector.Build.Dense(d);
        switch (InitialState)
        {
            case InitialStateKind.SingleExcitation:
                // |10…0⟩: site 0 excited (most-significant bit), rest |0⟩.
                psi[1 << (N - 1)] = Complex.One;
                break;
            case InitialStateKind.BondingMode:
            {
                // ψ₁ = Σ_l sin(π(l+1)/(N+1)) |1_l⟩: the delocalized single-excitation sine mode
                // (the F67 receiver state). The canonical carrier of the PTF painter protocol;
                // localized or multi-sector states break the rescaling picture and the guard
                // rightly flags every site unreliable.
                double norm = 0.0;
                for (int l = 0; l < N; l++)
                {
                    double a = Math.Sin(Math.PI * (l + 1) / (N + 1));
                    psi[1 << (N - 1 - l)] = a;
                    norm += a * a;
                }
                norm = Math.Sqrt(norm);
                for (int l = 0; l < N; l++)
                    psi[1 << (N - 1 - l)] /= norm;
                break;
            }
            case InitialStateKind.BellPair:
            default:
                // (|00…0⟩ + |11 0…0⟩)/√2: Bell+ on sites 0,1 (the two most-significant bits), rest |0⟩.
                double amp = 1.0 / Math.Sqrt(2.0);
                int both = (1 << (N - 1)) | (1 << (N - 2));
                psi[0] = amp;
                psi[both] = amp;
                break;
        }
        var rho = Matrix<Complex>.Build.Dense(d, d);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                rho[a, b] = psi[a] * Complex.Conjugate(psi[b]);
        return rho;
    }

    // ---- Lens primitives over the one trajectory ----

    /// <summary>Purity C = Tr(ρ²) = Σ_ij |ρ_ij|² of a density matrix.</summary>
    public static double Purity(ComplexMatrix rho)
    {
        double sum = 0.0;
        for (int i = 0; i < rho.RowCount; i++)
            for (int j = 0; j < rho.ColumnCount; j++)
            {
                double m = rho[i, j].Magnitude;
                sum += m * m;
            }
        return sum;
    }

    /// <summary>Normalized ℓ₁ coherence Ψ = (Σ_{i≠j} |ρ_ij|) / (d−1). Bell+ gives Ψ = 1/3 at d=4.</summary>
    public static double NormalizedCoherence(ComplexMatrix rho)
    {
        int d = rho.RowCount;
        if (d <= 1) return 0.0;
        double l1 = 0.0;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                if (i != j) l1 += rho[i, j].Magnitude;
        return l1 / (d - 1);
    }

    /// <summary>CΨ = purity × normalized coherence (the repo convention: Bell+ → 1/3).</summary>
    public static double Cpsi(ComplexMatrix rho) => Purity(rho) * NormalizedCoherence(rho);

    /// <summary>The local CΨ: the canonical <see cref="Cpsi"/> of the reduced 2-site density matrix on
    /// the <see cref="CarrierPair"/>. At N=2 with the default pair this equals the global CΨ exactly
    /// (the partial trace keeps both qubits), reproducing F25; for N≥3 it reads the coherence where it
    /// sits — in the carrier pair — so the fold stays audible where the global /(d−1) normalization
    /// has pushed CΨ(0) below ¼.</summary>
    public double LocalCpsi(ComplexMatrix rho)
        => Cpsi(PartialTrace.Of(rho, N, new[] { CarrierPair.Site1, CarrierPair.Site2 }));

    /// <summary>The light content of ρ: the purity-weighted mean light over the coherences, where the
    /// light of |i⟩⟨j| is popcount(i⊕j), the number of sites where ket and bra differ (the
    /// Absorption-Theorem light channel; the diagonal i=j carries light 0). Returns
    /// Σ_ij popcount(i⊕j)·|ρ_ij|² / Σ_ij |ρ_ij|² — a number in [0, N], the average number of lit
    /// sites the state's weight sits on. We picked this purity-weighted coherence light (not a
    /// per-site ⟨X/Y⟩) because it reads straight off the same ρ(t) the other lenses use, with only a
    /// popcount.</summary>
    public static double LightContent(ComplexMatrix rho)
    {
        int d = rho.RowCount;
        double weighted = 0.0, total = 0.0;
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
            {
                double m = rho[i, j].Magnitude;
                double w = m * m;
                total += w;
                weighted += w * System.Numerics.BitOperations.PopCount((uint)(i ^ j));
            }
        return total > 0 ? weighted / total : 0.0;
    }

    private double[]? _cpsi, _light, _dose, _localCpsi;

    private void EnsureLenses()
    {
        if (_cpsi is not null) return;
        EnsureEvolved();
        int n = _states!.Length;
        _cpsi = new double[n];
        _light = new double[n];
        _dose = new double[n];
        _localCpsi = new double[n];
        for (int s = 0; s < n; s++)
        {
            _cpsi[s] = Cpsi(_states[s]);
            _light[s] = LightContent(_states[s]);
            _dose[s] = Gamma * _tGrid![s];
            _localCpsi[s] = LocalCpsi(_states[s]);
        }
    }

    /// <summary>The first time CΨ(t) actually crosses ¼ (a downward/upward sign change, interpolated),
    /// and the dose K = γ·t there; null if it never crosses within the window. A state that starts
    /// already below ¼ (e.g. Bell+ at N≥3, where the global Ψ-normalization /(d−1) puts CΨ(0)=1/(d−1)<¼)
    /// has NO crossing, and this honestly returns null rather than the t=0 grid point.</summary>
    public (double Time, double Dose)? FirstQuarterCrossing(double threshold = 0.25)
    {
        var crossings = QuarterCrossingTimes(threshold);
        if (crossings.Count == 0) return null;
        double t = crossings[0];
        return (t, Gamma * t);
    }

    // ---- IInspectable: the score, the lenses, the events ----

    public string DisplayName =>
        $"Symphony (N={N}, {HType}, {Topology}, J={J.ToString("0.###", Inv)}, γ={Gamma.ToString("0.###", Inv)}, " +
        $"{InitialState}, t∈[0,{TMax.ToString("0.##", Inv)}], {TPoints} pts)";

    public string Summary
    {
        get
        {
            EnsureLenses();
            int lenses = 5;
            var events = BuildEvents();
            var cross = FirstQuarterCrossing();
            string crossClause = cross is { } c
                ? $"; first ¼ crossing at t={c.Time.ToString("0.###", Inv)} (K={c.Dose.ToString("0.####", Inv)})"
                : "";
            string movementClause = DefectBond is { } b
                ? $"; + the painters' movement (bond {b}, δJ={DeltaJ.ToString("0.###", Inv)})"
                : "";
            return $"one evolution (N={N}, J={J.ToString("0.###", Inv)}, γ={Gamma.ToString("0.###", Inv)}, " +
                   $"t∈[0,{TMax.ToString("0.##", Inv)}]); {lenses} lenses on one timeline; {events.Count} events" +
                   crossClause + movementClause + ".";
        }
    }

    public IEnumerable<IInspectable> Children
    {
        get
        {
            EnsureLenses();
            yield return ScoreNode();
            yield return PalindromeLens();
            yield return QuarterLens();
            yield return LocalQuarterLens();
            yield return DoseLens();
            yield return LightLens();
            if (DefectBond is not null)
                yield return Painters();
            yield return EventsNode();
        }
    }

    /// <summary>The painters' movement (only when --defect-bond is set): the second performance and
    /// the guard third, the PTF lenses (painters α_i, closure Σ ln α, the live K₁ chiral mirror) and
    /// the global clock. Built lazily and cached so its trajectories, like the clean one, are each
    /// built exactly once.</summary>
    private PaintersMovement? _painters;
    private PaintersMovement Painters()
    {
        EnsureEvolved();
        return _painters ??= new PaintersMovement(this, DefectBond!.Value, DeltaJ);
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    /// <summary>The shared timeline parameters: the t grid and the K = γ·t axis the lenses read.</summary>
    private InspectableNode ScoreNode()
    {
        var tGrid = _tGrid!;
        var k = tGrid.Select(t => Gamma * t).ToArray();
        return new InspectableNode("score",
            summary: $"shared timeline: {TPoints} points, t∈[0,{TMax.ToString("0.###", Inv)}], " +
                     $"K=γ·t∈[0,{(Gamma * TMax).ToString("0.####", Inv)}]; all lenses read the one ρ(t). " +
                     $"(EvolveCount={_evolveCount}: built once.)",
            children: new IInspectable[]
            {
                InspectableNode.RealScalar("t-max", TMax, "0.####"),
                InspectableNode.RealScalar("t-points", TPoints),
                InspectableNode.RealScalar("γ (Takt)", Gamma, "0.####"),
                new InspectableNode("K = γ·t axis",
                    summary: $"dose at end of window K={(Gamma * TMax).ToString("0.####", Inv)}",
                    payload: new InspectablePayload.Curve("K = γ·t", tGrid, k, "t", "K")),
            });
    }

    /// <summary>lens: palindrome — the spectrum's paired structure of THIS system (F1), read off the
    /// same Liouvillian decomposition <see cref="MirrorSystem"/> uses. Honest about whether every
    /// rate finds its 2σ − r partner.</summary>
    private InspectableNode PalindromeLens()
    {
        var chain = new ChainSystem(N, J, Gamma, HType, Topology);
        var channels = Enumerable.Range(0, N).Select(l => new ChannelRate($"q{l}", Gamma)).ToList();
        var mirror = new MirrorSystem(N, chain.BuildHamiltonian(), channels);
        var partners = mirror.PalindromePartners;
        int paired = partners.Count(p => p.PartnerPresent);
        bool holds = mirror.PalindromeHolds;
        return new InspectableNode("lens: palindrome",
            summary: $"F1 pairing verdict: holds = {holds}; {paired}/{partners.Count} rates find their " +
                     $"mirror partner 2σ − r (2σ = {(2.0 * Sigma).ToString("0.####", Inv)}).",
            children: new IInspectable[]
            {
                InspectableNode.RealScalar("slowest rate", mirror.Spectrum.SlowestRate, "0.####"),
                new InspectableNode("clock θ (deg)",
                    summary: (mirror.Rotation.Angle * 180.0 / Math.PI).ToString("0.#", Inv)),
            });
    }

    /// <summary>lens: quarter (CΨ) — the global CΨ(t), now the Envelope-Theorem witness. Reports the
    /// direction-split ¼-crossing count, the live envelope verdict (the peaks form a non-increasing
    /// sequence — proven N=2, verified N≥3, PROOF_MONOTONICITY_CPSI), and the envelope fold (the
    /// absorbing ¼ crossing; "the fold" now means THIS, never an upward oscillation).</summary>
    private InspectableNode QuarterLens()
    {
        var cpsi = _cpsi!;
        var tGrid = _tGrid!;
        double start = cpsi[0], min = cpsi.Min();

        var dirs = QuarterCrossingDirections(cpsi);
        int down = dirs.Count(d => d < 0), up = dirs.Count(d => d > 0);
        var (_, samples, floor) = GridFitness(cpsi.Max());
        string gridClause = double.IsInfinity(samples) ? "" : $" [≈{samples.ToString("0.#", Inv)} samples/oscillation]";
        string crossClause = dirs.Length == 0
            ? "no ¼ crossing in window"
            : $"{dirs.Length} ¼ crossing(s): {down}↓ + {up}↑{gridClause}";

        var env = QuarterEnvelope.Of(cpsi, tGrid);
        string envClause = env.IsNonIncreasing
            ? "envelope non-increasing ✓ (the Envelope Theorem holds live — proven N=2, verified N≥3, PROOF_MONOTONICITY_CPSI)"
            : $"envelope shows {env.RiseCount} predecessor-rise(s), max Δ={env.MaxRiseMagnitude.ToString("0.#####", Inv)} " +
              $"(peak-clip floor on this grid ≈ {floor.ToString("0.#####", Inv)}) — grid-sensitive, verify with ≥4× t-points; " +
              "a rise that SURVIVES refinement would falsify the Tier-2 verification";
        string foldClause = env.EnvelopeFoldTime is { } ft
            ? $"the fold (envelope, absorbing) at t={ft.ToString("0.###", Inv)} (K={(Gamma * ft).ToString("0.####", Inv)})"
            : "no envelope fold in window";

        return new InspectableNode("lens: quarter (CΨ)",
            summary: $"CΨ(0) = {start.ToString("0.####", Inv)}, min = {min.ToString("0.####", Inv)}; {crossClause}. " +
                     $"{envClause}. {foldClause}.",
            payload: new InspectablePayload.Curve("CΨ(t)", tGrid, cpsi, "t", "CΨ"));
    }

    /// <summary>lens: quarter (local CΨ) — CΨ of the reduced 2-site carrier-pair state along the one
    /// trajectory. Audible at N≥3 where the global lens is silent. Direction-split ¼-crossing count,
    /// plus the envelope verdict: unlike the global CΨ (theorem-bound non-increasing), the reduced open
    /// subsystem has NO such theorem, so its beat envelope can genuinely RISE — the freedom. A detected
    /// rise is reported grid-sensitive (parabolic-apex + predecessor semantics; verify under refinement).
    /// See QuarterEnvelope and PROOF_MONOTONICITY_CPSI.</summary>
    private InspectableNode LocalQuarterLens()
    {
        var local = _localCpsi!;
        var tGrid = _tGrid!;
        double start = local[0], min = local.Min(), max = local.Max();
        string pair = $"({CarrierPair.Site1},{CarrierPair.Site2})";

        var dirs = QuarterCrossingDirections(local);

        // The carrier pair begins with no shared coherence (local CΨ(0) ≈ 0) and never folds: a single
        // excitation, say, places no coherence on the pair, and although the Hamiltonian pumps some in by
        // hopping, it stays below ¼. There is no quarter-fold to report.
        if (start <= 1e-12 && dirs.Length == 0)
            return new InspectableNode("lens: quarter (local CΨ)",
                summary: $"2-site reduced ρ on carrier pair {pair}: no pair coherence in this initial " +
                         $"state (local CΨ(0) ≈ 0); the Hamiltonian pumps it to at most " +
                         $"{max.ToString("0.####", Inv)} (< ¼), so no fold.",
                payload: new InspectablePayload.Curve("local CΨ(t)", tGrid, local, "t", "local CΨ"));

        int down = dirs.Count(d => d < 0), up = dirs.Count(d => d > 0);
        string crossClause = dirs.Length == 0
            ? "no ¼ crossing in window"
            : $"{dirs.Length} ¼ crossing(s): {down}↓ + {up}↑";

        var env = QuarterEnvelope.Of(local, tGrid);
        var (_, _, floor) = GridFitness(max);
        string envClause = env.RiseCount > 0 && env.FirstRiseTime is { } rt
            ? $"envelope RISES: {env.RiseCount} predecessor-rise(s), max Δ={env.MaxRiseMagnitude.ToString("0.#####", Inv)} " +
              $"at t={rt.ToString("0.###", Inv)} — the freedom (beating; no theorem binds the reduced open subsystem; " +
              $"peak-clip floor ≈ {floor.ToString("0.#####", Inv)}). The rise is grid-sensitive: verify with ≥4× t-points"
            : "envelope non-increasing in this window (no rise resolved on this grid)";

        return new InspectableNode("lens: quarter (local CΨ)",
            summary: $"2-site reduced ρ on carrier pair {pair}: local CΨ(0) = {start.ToString("0.####", Inv)}, " +
                     $"min = {min.ToString("0.####", Inv)}, max = {max.ToString("0.####", Inv)}; {crossClause}. " +
                     $"{envClause}.",
            payload: new InspectablePayload.Curve("local CΨ(t)", tGrid, local, "t", "local CΨ"));
    }

    /// <summary>lens: dose (K) — K = γ·t marks; the dose of the fold is K at the first ¼ crossing.</summary>
    private InspectableNode DoseLens()
    {
        var cross = FirstQuarterCrossing();
        string doseClause = cross is { } c
            ? $"the dose of the fold: K = {c.Dose.ToString("0.####", Inv)} at the first ¼ crossing (t={c.Time.ToString("0.###", Inv)})"
            : "no global ¼ crossing in window, so no global fold dose";

        var localTimes = QuarterCrossingTimes(_localCpsi!, _tGrid!);
        string localClause = localTimes.Count == 0
            ? "; no local fold in window"
            : $"; the local fold: K = {(Gamma * localTimes[0]).ToString("0.####", Inv)} at the first local " +
              $"¼ crossing (carrier pair {CarrierPair.Site1},{CarrierPair.Site2}, t={localTimes[0].ToString("0.###", Inv)})";

        return new InspectableNode("lens: dose (K)",
            summary: $"K = γ·t reaches {(Gamma * TMax).ToString("0.####", Inv)} at the window end; {doseClause}{localClause}.",
            payload: new InspectablePayload.Curve("K = γ·t", _tGrid!, _dose!, "t", "K"));
    }

    /// <summary>lens: light — the state's light content over time, the purity-weighted mean
    /// popcount(i⊕j) over the coherences; curve payload.</summary>
    private InspectableNode LightLens()
    {
        var light = _light!;
        double start = light[0], end = light[^1];
        return new InspectableNode("lens: light",
            summary: $"light content (purity-weighted mean popcount(i⊕j), the lit-site count): " +
                     $"{start.ToString("0.####", Inv)} → {end.ToString("0.####", Inv)}, " +
                     $"min {light.Min().ToString("0.####", Inv)}, max {light.Max().ToString("0.####", Inv)}.",
            payload: new InspectablePayload.Curve("light(t)", _tGrid!, light, "t", "⟨popcount(i⊕j)⟩"));
    }

    /// <summary>events — the cross-lens axis: every detected event from every lens merged and sorted
    /// by time. THIS NODE IS THE POINT: interplay is coincidence made visible on one axis.</summary>
    private InspectableNode EventsNode()
    {
        var events = BuildEvents();
        var children = events.Select(e => (IInspectable)new InspectableNode(
            displayName: $"t={e.Time.ToString("0.###", Inv)}, K={e.Dose.ToString("0.####", Inv)}",
            summary: $"[{e.Lens}] {e.What}")).ToList();
        return new InspectableNode("events",
            summary: $"{events.Count} cross-lens event(s) on the one timeline, sorted by time; " +
                     "coincidences = adjacency here.",
            children: children);
    }

    // ---- Event detection across the lenses ----

    private sealed record SymphonyEvent(double Time, double Dose, string Lens, string What);

    /// <summary>All detected events merged and sorted by time: every CΨ=¼ crossing (lens quarter),
    /// K reaching {0.25, 0.5, 1.0} (lens dose), and light-content extrema / sign-changes-of-slope
    /// (lens light). The sorted list is the cross-lens axis.</summary>
    private List<SymphonyEvent> BuildEvents()
    {
        EnsureLenses();
        var tGrid = _tGrid!;
        var cpsi = _cpsi!;
        var light = _light!;
        var events = new List<SymphonyEvent>();

        // lens quarter: every ¼ crossing (linearly interpolated crossing time).
        foreach (double tc in QuarterCrossingTimes())
            events.Add(new SymphonyEvent(tc, Gamma * tc, "quarter",
                $"CΨ crosses ¼ (the fold; quantum→classical boundary)"));

        // lens local quarter: the carrier-pair fold (audible at N≥3 where the global one is silent),
        // counted in both directions (the open-subsystem heartbeat).
        var localTimes = QuarterCrossingTimes(_localCpsi!, tGrid);
        var localDirs = QuarterCrossingDirections(_localCpsi!);
        System.Diagnostics.Debug.Assert(localTimes.Count == localDirs.Length,
            "QuarterCrossingTimes and QuarterCrossingDirections must return order-aligned, equal-length results");
        for (int k = 0; k < localTimes.Count; k++)
        {
            double tc = localTimes[k];
            string dir = localDirs[k] < 0 ? "down" : "up";
            events.Add(new SymphonyEvent(tc, Gamma * tc, "local quarter",
                $"local CΨ (carrier pair {CarrierPair.Site1},{CarrierPair.Site2}) crosses ¼ ({dir})"));
        }

        // lens dose: K reaching the milestones, within the window.
        foreach (double kMark in new[] { 0.25, 0.5, 1.0 })
        {
            double tMark = kMark / Gamma;
            if (tMark <= TMax + 1e-12)
                events.Add(new SymphonyEvent(tMark, kMark, "dose",
                    $"dose K reaches {kMark.ToString("0.##", Inv)}"));
        }

        // lens light: interior extrema of the light curve (slope sign change).
        for (int s = 1; s < light.Length - 1; s++)
        {
            double dPrev = light[s] - light[s - 1];
            double dNext = light[s + 1] - light[s];
            if (dPrev == 0.0 && dNext == 0.0) continue;
            bool maxHere = dPrev > 0 && dNext <= 0;
            bool minHere = dPrev < 0 && dNext >= 0;
            if (maxHere || minHere)
                events.Add(new SymphonyEvent(tGrid[s], Gamma * tGrid[s], "light",
                    $"light content {(maxHere ? "peaks" : "bottoms")} at {light[s].ToString("0.###", Inv)}"));
        }

        // painters' movement (only when --defect-bond is set): the closure verdict as an event at the
        // window end, and the chiral-mirror deviation only if the K₁ identity is BROKEN.
        if (DefectBond is not null)
        {
            var pm = Painters();
            if (pm.HasLenses)
            {
                events.Add(new SymphonyEvent(TMax, Gamma * TMax, "painters",
                    $"closure {(pm.ClosureInWindow ? "holds" : "out of window")}: " +
                    $"Σ ln α (reliable) = {pm.ClosureSum.ToString("0.####", Inv)}"));
                if (!pm.ChiralMirrorExact)
                    events.Add(new SymphonyEvent(TMax, Gamma * TMax, "painters",
                        $"chiral mirror BROKEN: max site-wise |ΔP| = {pm.ChiralMirrorDeviation.ToString("E2", Inv)}"));
            }
        }

        events.Sort((a, b) => a.Time.CompareTo(b.Time));
        return events;
    }

    /// <summary>The ¼-crossing times of a CΨ array on its t grid, linearly interpolated between the
    /// bracketing grid points (both downward and upward crossings, in case a trajectory dips and
    /// recovers). Static so both the global and the local curves can share it.</summary>
    public static List<double> QuarterCrossingTimes(double[] cpsi, double[] tGrid, double threshold = 0.25)
    {
        var times = new List<double>();
        for (int s = 1; s < cpsi.Length; s++)
        {
            double a = cpsi[s - 1] - threshold;
            double b = cpsi[s] - threshold;
            if (a == 0.0) { times.Add(tGrid[s - 1]); continue; }
            if (a * b < 0.0)
            {
                double frac = a / (a - b);
                times.Add(tGrid[s - 1] + frac * (tGrid[s] - tGrid[s - 1]));
            }
        }
        return times;
    }

    /// <summary>The direction of each ¼ crossing of a CΨ array: −1 for a downward crossing (CΨ falls
    /// through ¼), +1 for an upward crossing (it rises back through ¼). Same order as
    /// <see cref="QuarterCrossingTimes(double[], double[], double)"/>.</summary>
    public static int[] QuarterCrossingDirections(double[] cpsi, double threshold = 0.25)
    {
        var dirs = new List<int>();
        for (int s = 1; s < cpsi.Length; s++)
        {
            double a = cpsi[s - 1] - threshold;
            double b = cpsi[s] - threshold;
            if (a == 0.0) { dirs.Add(b < 0 ? -1 : +1); continue; }
            if (a * b < 0.0) dirs.Add(a > 0 ? -1 : +1);
        }
        return dirs.ToArray();
    }

    /// <summary>The ¼-crossing times of the GLOBAL CΨ curve (the one trajectory's CΨ(t)).</summary>
    private List<double> QuarterCrossingTimes(double threshold = 0.25)
    {
        EnsureLenses();
        return QuarterCrossingTimes(_cpsi!, _tGrid!, threshold);
    }
}

/// <summary>The painters' movement: the Perspectival Time Field read off ONE Symphony, the piece
/// played a second (defected) and a guard third (δJ/2) time on a single-bond J-defect.
///
/// <para>The physics (ground truth: <c>docs/proofs/PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md</c> and
/// <c>simulations/framework/workflows/ptf.py</c>, the reference this ports):</para>
/// <list type="bullet">
///   <item><b>painters</b> — per-site purity P_i(t) of the defected run is a time-rescaling of the
///         clean run, P_B(i, t) ≈ P_A(i, α_i·t); α_i is fitted per site (bounded [0.1, 10], L2 over
///         the grid, P_A cubic-interpolated). f_i = (α_i − 1)/δJ at δJ and at the guard δJ/2; sites
///         where the two disagree (or whose α is not of sane magnitude) are unreliable and excluded.</item>
///   <item><b>closure</b> — Σ_i ln α_i over the reliable sites; ≈ 0 within ±0.05 for the canonical
///         XY/Z-dephasing/|δJ|≤0.1/real-ψ protocol (EQ-014, a perturbative empirical regularity, NOT
///         a per-site intrinsic clock).</item>
///   <item><b>chiral mirror</b> — K₁ = Π_{l odd} Z_l; P_i(t; K₁ψ) = P_i(t; ψ) exactly site-wise for
///         the defected XY chain and real ψ (<see cref="Ptf.ChiralMirrorTrajectoryClaim"/>). Reported
///         live as the max site-wise purity deviation between the ψ and K₁ψ defected runs.</item>
///   <item><b>clock</b> — the Takt gap (slowest nonzero decay rate) and ω_mem (max |Im λ| among the
///         modes at the gap) read off the shared L_A eigendecomposition.</item>
/// </list>
///
/// <para><b>Honest declines</b>: the canonical PTF reference is the XY chain; if H is not XY, or the
/// initial state is not real, the movement appears with a "silent: …" summary and no lenses.</para></summary>
public sealed class PaintersMovement : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Closure-law window: |Σ ln α_i| ≤ this is "in window" (EQ-014, ±0.05).</summary>
    public const double ClosureWindow = 0.05;
    /// <summary>The K₁ chiral mirror is EXACT iff the max site-wise |ΔP| ≤ this; else BROKEN.</summary>
    public const double ChiralExactTol = 1e-10;

    // Guard / reliability thresholds, ported from perspectives_panel.
    private const double AlphaLo = 0.1, AlphaHi = 10.0;
    private const double FMax = 10.0, ConsistencyTol = 0.5, FFloor = 0.5;

    public Symphony Parent { get; }
    public int Bond { get; }
    public double DeltaJ { get; }

    /// <summary>True when the canonical PTF protocol applies (XY chain, real initial state). When
    /// false the movement is silent (no lenses), with <see cref="DeclineReason"/> stating why.</summary>
    public bool HasLenses { get { Ensure(); return _hasLenses; } }
    public string DeclineReason { get { Ensure(); return _declineReason; } }
    private bool _hasLenses;
    private string _declineReason = "";

    // Per-trajectory build accounting: each of the four runs (clean P_A, defected P_B, guard P_guard,
    // K₁ defected) built exactly once. A test guards this.
    private int _buildCount;
    public int BuildCount => _buildCount;

    private int _n;
    private double[]? _tGrid;
    private double[][]? _pA;       // [site][t] clean per-site purity
    private double[][]? _pB;       // [site][t] defected per-site purity
    private double[]? _alpha;      // per-site α_i at δJ
    private double[]? _f;          // per-site f_i = (α−1)/δJ at δJ
    private double[]? _fGuard;     // per-site f_i at the guard δJ/2
    private bool[]? _reliable;     // per-site reliability flag (sane & linear)
    private double _closureSum;    // Σ ln α over reliable sites (NaN if none)
    private double _closureAll;    // Σ ln α over all sites
    private double _chiralDev;     // max site-wise |P_i(t;K₁ψ) − P_i(t;ψ)| on the defected run
    private double _taktGap, _omegaMem;

    public double ClosureSum { get { Ensure(); return _closureSum; } }
    public bool ClosureInWindow => HasLenses && Math.Abs(_closureSum) <= ClosureWindow;
    public double ChiralMirrorDeviation { get { Ensure(); return _chiralDev; } }
    public bool ChiralMirrorExact => HasLenses && _chiralDev <= ChiralExactTol;
    public IReadOnlyList<double> Alphas { get { Ensure(); return _alpha ?? Array.Empty<double>(); } }
    public IReadOnlyList<bool> Reliable { get { Ensure(); return _reliable ?? Array.Empty<bool>(); } }

    public PaintersMovement(Symphony parent, int bond, double deltaJ)
    {
        Parent = parent;
        Bond = bond;
        DeltaJ = deltaJ;
    }

    private bool _built;
    private void Ensure()
    {
        if (_built) return;
        _built = true;
        Build();
    }

    private void Build()
    {
        var p = Parent;
        _n = p.N;
        int d = p.Dim;

        // Decline 1: the canonical PTF reference is the XY chain.
        if (p.HType != HamiltonianType.XY)
        {
            _hasLenses = false;
            _declineReason = $"the canonical PTF protocol needs the XY chain (got {p.HType})";
            return;
        }

        // Decline 2: the K₁-conjugation + reality argument needs a real initial state.
        var rho0 = p.BuildInitialState(d);
        if (!IsReal(rho0))
        {
            _hasLenses = false;
            _declineReason = "the canonical PTF protocol needs a real initial state (got a complex ρ(0))";
            return;
        }

        // The defect bond must index into this chain's bonds.
        var chain = new ChainSystem(_n, p.J, p.Gamma, p.HType, p.Topology);
        if (Bond < 0 || Bond >= chain.Bonds.Count)
        {
            _hasLenses = false;
            _declineReason = $"defect bond {Bond} outside [0, {chain.Bonds.Count}) for this topology";
            return;
        }

        _hasLenses = true;

        var H = chain.BuildHamiltonian();
        var gammas = Enumerable.Repeat(p.Gamma, _n).ToArray();
        var L_A = PauliDephasingDissipator.BuildZ(H, gammas);

        var bondPair = chain.Bonds[Bond];
        var V = BondPerturbation.Build(_n, bondPair.Site1, bondPair.Site2, BondPerturbation.Kind.XY);
        var L_B = L_A + (Complex)DeltaJ * V;
        var L_guard = L_A + (Complex)(DeltaJ / 2.0) * V;

        // The shared timeline = the Symphony's grid.
        var tGrid = p.TimeGrid.ToArray();
        _tGrid = tGrid;

        // Per-site Pauli operators X_i, Y_i, Z_i (built once).
        var sitePaulis = new (ComplexMatrix X, ComplexMatrix Y, ComplexMatrix Z)[_n];
        for (int i = 0; i < _n; i++)
            sitePaulis[i] = (PauliString.SiteOp(_n, i, PauliLetter.X),
                             PauliString.SiteOp(_n, i, PauliLetter.Y),
                             PauliString.SiteOp(_n, i, PauliLetter.Z));

        // The clock from L_A (read once, off the SHARED clean eigendecomposition the Symphony kept).
        (_taktGap, _omegaMem) = GlobalClock(p.LiouvillianEigenvalues);

        // The four trajectories, each built exactly once.
        _pA = PuritySites(L_A, rho0, tGrid, sitePaulis);            // clean
        _pB = PuritySites(L_B, rho0, tGrid, sitePaulis);            // defected (the second performance)
        var pGuardSites = PuritySites(L_guard, rho0, tGrid, sitePaulis);   // guard third (δJ/2)

        // The K₁ chiral-mirror run: same defected system, initial state K₁ψ.
        var rhoK1 = ConjugateByK1(rho0);
        var pBk1 = PuritySites(L_B, rhoK1, tGrid, sitePaulis);
        double dev = 0.0;
        for (int i = 0; i < _n; i++)
            for (int s = 0; s < tGrid.Length; s++)
                dev = Math.Max(dev, Math.Abs(_pB[i][s] - pBk1[i][s]));
        _chiralDev = dev;

        // Per-site α fits at δJ and the guard δJ/2, then the reliability guard + closure.
        _alpha = new double[_n];
        var alphaGuard = new double[_n];
        for (int i = 0; i < _n; i++)
        {
            _alpha[i] = FitAlpha(tGrid, _pA[i], _pB[i]);
            alphaGuard[i] = FitAlpha(tGrid, _pA[i], pGuardSites[i]);
        }

        _f = new double[_n];
        _fGuard = new double[_n];
        _reliable = new bool[_n];
        for (int i = 0; i < _n; i++)
        {
            _f[i] = (_alpha[i] - 1.0) / DeltaJ;
            _fGuard[i] = (alphaGuard[i] - 1.0) / (DeltaJ / 2.0);
            bool sane = Math.Abs(_f[i]) <= FMax;
            bool linear = Math.Abs(_f[i] - _fGuard[i]) <= ConsistencyTol * (Math.Abs(_fGuard[i]) + FFloor);
            _reliable[i] = sane && linear;
        }

        double sumReliable = 0.0; int nReliable = 0;
        double sumAll = 0.0;
        for (int i = 0; i < _n; i++)
        {
            double la = Math.Log(Math.Max(_alpha[i], 1e-30));
            sumAll += la;
            if (_reliable[i]) { sumReliable += la; nReliable++; }
        }
        _closureAll = sumAll;
        _closureSum = nReliable > 0 ? sumReliable : double.NaN;
    }

    /// <summary>Per-site purity P_i(t) = ½(1 + ⟨X_i⟩² + ⟨Y_i⟩² + ⟨Z_i⟩²) along the spectral evolution
    /// of L from ρ0 on tGrid. One eigendecomposition of L, one mat-vec per step (the same propagation
    /// the Symphony uses). Counts one trajectory build.</summary>
    private double[][] PuritySites(ComplexMatrix L, ComplexMatrix rho0, double[] tGrid,
        (ComplexMatrix X, ComplexMatrix Y, ComplexMatrix Z)[] sitePaulis)
    {
        int d = rho0.RowCount;
        int dim2 = d * d;

        var evd = L.Evd();
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();
        var lambda = evd.EigenValues;

        var vec0 = ComplexVector.Build.Dense(dim2);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                vec0[a * d + b] = rho0[a, b];
        var c0 = Rinv * vec0;

        var purity = new double[_n][];
        for (int i = 0; i < _n; i++) purity[i] = new double[tGrid.Length];

        var expDiag = ComplexVector.Build.Dense(dim2);
        var rhoT = ComplexMatrix.Build.Dense(d, d);
        for (int s = 0; s < tGrid.Length; s++)
        {
            double t = tGrid[s];
            for (int k = 0; k < dim2; k++) expDiag[k] = Complex.Exp(lambda[k] * t) * c0[k];
            var vecT = R * expDiag;
            for (int a = 0; a < d; a++)
                for (int b = 0; b < d; b++)
                    rhoT[a, b] = vecT[a * d + b];
            // Hermitize for numerical safety (matches the Python reference).
            for (int a = 0; a < d; a++)
                for (int b = a + 1; b < d; b++)
                {
                    var avg = 0.5 * (rhoT[a, b] + Complex.Conjugate(rhoT[b, a]));
                    rhoT[a, b] = avg;
                    rhoT[b, a] = Complex.Conjugate(avg);
                }
            for (int i = 0; i < _n; i++)
            {
                double x = Trace(sitePaulis[i].X, rhoT);
                double y = Trace(sitePaulis[i].Y, rhoT);
                double z = Trace(sitePaulis[i].Z, rhoT);
                purity[i][s] = 0.5 * (1.0 + x * x + y * y + z * z);
            }
        }

        _buildCount++;
        return purity;
    }

    /// <summary>Re(Tr(Op·ρ)) — the expectation of a Hermitian operator in state ρ.</summary>
    private static double Trace(ComplexMatrix op, ComplexMatrix rho)
    {
        int d = op.RowCount;
        Complex tr = Complex.Zero;
        for (int a = 0; a < d; a++)
            for (int c = 0; c < d; c++)
                tr += op[a, c] * rho[c, a];
        return tr.Real;
    }

    /// <summary>Fit α ∈ [0.1, 10] minimizing Σ_t (P_A(α·t) − P_B(t))² (golden-section search; P_A is
    /// natural-cubic-spline-interpolated to match scipy's interp1d(kind='cubic'), clamped to the
    /// endpoint values off-grid, the scipy fill_value=(y0, y_end) the reference uses).</summary>
    private static double FitAlpha(double[] tGrid, double[] pA, double[] pB)
    {
        var spline = new NaturalCubicSpline(tGrid, pA);

        double Mse(double alpha)
        {
            double sum = 0.0;
            for (int s = 0; s < tGrid.Length; s++)
            {
                double dv = spline.Eval(alpha * tGrid[s]) - pB[s];
                sum += dv * dv;
            }
            return sum / tGrid.Length;
        }

        // Golden-section minimization on [AlphaLo, AlphaHi], to scipy's xatol ~1e-7.
        const double gr = 0.6180339887498949;
        double a = AlphaLo, b = AlphaHi;
        double c = b - gr * (b - a), dd = a + gr * (b - a);
        double fc = Mse(c), fd = Mse(dd);
        for (int it = 0; it < 200 && (b - a) > 1e-7; it++)
        {
            if (fc < fd) { b = dd; dd = c; fd = fc; c = b - gr * (b - a); fc = Mse(c); }
            else { a = c; c = dd; fc = fd; dd = a + gr * (b - a); fd = Mse(dd); }
        }
        return 0.5 * (a + b);
    }

    /// <summary>A natural cubic spline through (x, y), the same construction scipy's
    /// interp1d(kind='cubic') uses (second derivative zero at both ends), with the endpoint values
    /// clamped off-grid (scipy fill_value=(y[0], y[-1])). x must be strictly ascending.</summary>
    private sealed class NaturalCubicSpline
    {
        private readonly double[] _x, _y, _m;   // _m = second derivatives at the knots

        public NaturalCubicSpline(double[] x, double[] y)
        {
            _x = x; _y = y;
            int n = x.Length;
            _m = new double[n];
            if (n < 3) return;   // n<3: Eval falls back to linear

            // Solve the tridiagonal system for the second derivatives (natural BCs m[0]=m[n-1]=0).
            var sub = new double[n];   // sub-diagonal
            var diag = new double[n];
            var sup = new double[n];   // super-diagonal
            var rhs = new double[n];
            diag[0] = 1.0; rhs[0] = 0.0;
            diag[n - 1] = 1.0; rhs[n - 1] = 0.0;
            for (int i = 1; i < n - 1; i++)
            {
                double hi = x[i] - x[i - 1];
                double hip = x[i + 1] - x[i];
                sub[i] = hi;
                diag[i] = 2.0 * (hi + hip);
                sup[i] = hip;
                rhs[i] = 6.0 * ((y[i + 1] - y[i]) / hip - (y[i] - y[i - 1]) / hi);
            }
            // Thomas algorithm.
            for (int i = 1; i < n; i++)
            {
                double w = sub[i] / diag[i - 1];
                diag[i] -= w * sup[i - 1];
                rhs[i] -= w * rhs[i - 1];
            }
            _m[n - 1] = rhs[n - 1] / diag[n - 1];
            for (int i = n - 2; i >= 0; i--)
                _m[i] = (rhs[i] - sup[i] * _m[i + 1]) / diag[i];
        }

        public double Eval(double xq)
        {
            var x = _x; var y = _y;
            int n = x.Length;
            if (xq <= x[0]) return y[0];
            if (xq >= x[n - 1]) return y[n - 1];
            int lo = 0, hi = n - 1;
            while (hi - lo > 1)
            {
                int mid = (lo + hi) >> 1;
                if (x[mid] <= xq) lo = mid; else hi = mid;
            }
            int j = lo;
            double h = x[j + 1] - x[j];
            if (n < 3)   // linear fallback
                return y[j] + (y[j + 1] - y[j]) * (xq - x[j]) / h;
            double a = (x[j + 1] - xq) / h;
            double b = (xq - x[j]) / h;
            return a * y[j] + b * y[j + 1]
                 + ((a * a * a - a) * _m[j] + (b * b * b - b) * _m[j + 1]) * (h * h) / 6.0;
        }
    }

    /// <summary>K₁ = Π_{l odd} Z_l conjugation of ρ0 (the odd-sublattice Z product, 0-indexed). K₁ is
    /// diagonal ±1, so K₁ρK₁ flips the sign of ρ[a,b] when a and b differ in odd-site parity.</summary>
    private ComplexMatrix ConjugateByK1(ComplexMatrix rho0)
    {
        int d = rho0.RowCount;
        var sign = new int[d];
        for (int a = 0; a < d; a++)
        {
            int s = 1;
            // bit l (0-indexed site) is the (N-1-l)-th bit (site 0 = most-significant), matching
            // Symphony.BuildInitialState. K₁ applies Z on odd sites l = 1, 3, 5, …
            for (int l = 1; l < _n; l += 2)
            {
                int bit = (a >> (_n - 1 - l)) & 1;
                if (bit == 1) s = -s;
            }
            sign[a] = s;
        }
        var outM = ComplexMatrix.Build.Dense(d, d);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                outM[a, b] = sign[a] * sign[b] * rho0[a, b];
        return outM;
    }

    private static bool IsReal(ComplexMatrix m, double tol = 1e-12)
    {
        for (int a = 0; a < m.RowCount; a++)
            for (int b = 0; b < m.ColumnCount; b++)
                if (Math.Abs(m[a, b].Imaginary) > tol) return false;
        return true;
    }

    /// <summary>The global clock from an L spectrum: Takt gap = slowest nonzero decay rate
    /// (= min over rate = −Re λ of the rates above tol), ω_mem = max |Im λ| among the modes at the gap.</summary>
    private static (double gap, double omega) GlobalClock(IReadOnlyList<Complex> evals,
        double tol = 1e-9, double gapTol = 1e-6)
    {
        double gap = double.PositiveInfinity;
        foreach (var e in evals)
        {
            double rate = -e.Real;
            if (rate > tol && rate < gap) gap = rate;
        }
        if (double.IsInfinity(gap)) return (0.0, 0.0);   // γ = 0: the clock is stopped
        double omega = 0.0;
        foreach (var e in evals)
            if (Math.Abs(-e.Real - gap) <= gapTol) omega = Math.Max(omega, Math.Abs(e.Imaginary));
        return (gap, omega);
    }

    // ---- IInspectable ----

    public string DisplayName => "movement: painters";

    public string Summary
    {
        get
        {
            Ensure();
            if (!HasLenses)
                return $"silent: {DeclineReason}.";
            int nReliable = _reliable!.Count(r => r);
            string mirror = ChiralMirrorExact
                ? $"chiral mirror EXACT (max |ΔP| = {_chiralDev.ToString("E2", Inv)})"
                : $"chiral mirror BROKEN (max |ΔP| = {_chiralDev.ToString("E2", Inv)})";
            return $"the piece played twice (clean + defected on bond {Bond}, δJ={DeltaJ.ToString("0.###", Inv)}) " +
                   $"plus the guard third (δJ/2); painters α over {_n} sites ({nReliable} reliable); " +
                   $"closure Σ ln α = {_closureSum.ToString("0.####", Inv)} " +
                   $"({(ClosureInWindow ? "in" : "out of")} the ±{ClosureWindow.ToString("0.##", Inv)} window); {mirror}.";
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            Ensure();
            if (!HasLenses) yield break;
            yield return SecondPerformanceNode();
            yield return PaintersLens();
            yield return ClosureLens();
            yield return ChiralMirrorLens();
            yield return ClockNode();
        }
    }

    private InspectableNode SecondPerformanceNode() =>
        new InspectableNode("the second performance",
            summary: $"the piece is played twice — clean (L_A) and defected (L_B = L_A + δJ·V on bond {Bond}) " +
                     $"— plus a guard third at δJ/2 for the reliability guard. Each trajectory built exactly " +
                     $"once (BuildCount={_buildCount}: P_A, P_B, P_guard, and the K₁ run).");

    private InspectableNode PaintersLens()
    {
        var children = new List<IInspectable>();
        for (int i = 0; i < _n; i++)
            children.Add(new InspectableNode($"site {i}",
                summary: $"α = {_alpha![i].ToString("0.#####", Inv)}, f = (α−1)/δJ = {_f![i].ToString("0.###", Inv)} " +
                         $"(guard f = {_fGuard![i].ToString("0.###", Inv)}); {(_reliable![i] ? "reliable" : "UNRELIABLE")}",
                payload: new InspectablePayload.Real($"α (site {i})", _alpha[i], "0.#####")));
        int nReliable = _reliable!.Count(r => r);
        return new InspectableNode("lens: painters",
            summary: $"per-site time-rescaling α_i: P_B(i, t) ≈ P_A(i, α_i·t); {nReliable}/{_n} sites reliable " +
                     $"(sane |f|≤{FMax.ToString("0", Inv)} and linear across δJ vs δJ/2).",
            children: children);
    }

    private InspectableNode ClosureLens() =>
        new InspectableNode("lens: closure",
            summary: $"Σ_i ln α_i over the reliable sites = {_closureSum.ToString("0.######", Inv)}; " +
                     $"{(ClosureInWindow ? "IN" : "OUT OF")} the ±{ClosureWindow.ToString("0.##", Inv)} window " +
                     $"(Σ over all sites = {_closureAll.ToString("0.####", Inv)}). " +
                     "HONEST SCOPE: α_i is per-(site, purity, state, event); the closure is the EQ-014 " +
                     "perturbative empirical regularity for the canonical XY / Z-dephasing / |δJ|≤0.1 / real-ψ " +
                     "protocol, NOT an intrinsic per-site clock.");

    private InspectableNode ChiralMirrorLens() =>
        new InspectableNode("lens: chiral mirror",
            summary: $"K₁ = Π_{{l odd}} Z_l: max site-wise |P_i(t; K₁ψ) − P_i(t; ψ)| on the defected run = " +
                     $"{_chiralDev.ToString("E3", Inv)} ⟹ {(ChiralMirrorExact ? "EXACT (≤1e-10)" : "BROKEN")}. " +
                     "The PTF surviving Tier-1 law (EQ-014); see ChiralMirrorTrajectoryClaim " +
                     "(docs/proofs/PROOF_PTF_CHIRAL_MIRROR_RATE_LAW.md).");

    private InspectableNode ClockNode() =>
        new InspectableNode("clock",
            summary: $"global clock from the shared L_A spectrum: Takt gap (slowest nonzero decay rate) = " +
                     $"{_taktGap.ToString("0.#####", Inv)} (= 2γ floor for a dephasing chain), " +
                     $"ω_mem (max |Im λ| at the gap) = {_omegaMem.ToString("0.#####", Inv)}.",
            children: new IInspectable[]
            {
                InspectableNode.RealScalar("Takt gap", _taktGap, "0.######"),
                InspectableNode.RealScalar("ω_mem", _omegaMem, "0.######"),
            });
}

/// <summary>Which initial density matrix the Symphony evolves. <see cref="BellPair"/> (default) is
/// Bell+ on sites (0,1) with the rest in |0⟩: the repo's hardware-confirmed ¼-crossing carrier
/// (F25, CΨ(0)=1/3), reproducing the F25 closed form at N=2. <see cref="SingleExcitation"/> is
/// |10…0⟩, one excitation on site 0.</summary>
public enum InitialStateKind
{
    BellPair,
    SingleExcitation,
    /// <summary>The k=1 sine mode Σ_l sin(π(l+1)/(N+1))|1_l⟩ (F67 bonding receiver):
    /// the canonical carrier of the PTF painter protocol.</summary>
    BondingMode,
}
