using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
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
    public double TMax { get; }
    public int TPoints { get; }

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
        int tPoints = 60)
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
    }

    // ---- The one evolution, computed once, shared by all lenses ----

    private double[]? _tGrid;
    private ComplexMatrix[]? _states;     // ρ(t) at each grid point, vec convention ρ[a,b] = vec[a*d+b]
    private int _evolveCount;             // how many times the trajectory was actually built (must stay 1)

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

    private double[]? _cpsi, _light, _dose;

    private void EnsureLenses()
    {
        if (_cpsi is not null) return;
        EnsureEvolved();
        int n = _states!.Length;
        _cpsi = new double[n];
        _light = new double[n];
        _dose = new double[n];
        for (int s = 0; s < n; s++)
        {
            _cpsi[s] = Cpsi(_states[s]);
            _light[s] = LightContent(_states[s]);
            _dose[s] = Gamma * _tGrid![s];
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
            int lenses = 4;
            var events = BuildEvents();
            var cross = FirstQuarterCrossing();
            string crossClause = cross is { } c
                ? $"; first ¼ crossing at t={c.Time.ToString("0.###", Inv)} (K={c.Dose.ToString("0.####", Inv)})"
                : "";
            return $"one evolution (N={N}, J={J.ToString("0.###", Inv)}, γ={Gamma.ToString("0.###", Inv)}, " +
                   $"t∈[0,{TMax.ToString("0.##", Inv)}]); {lenses} lenses on one timeline; {events.Count} events" +
                   crossClause + ".";
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
            yield return DoseLens();
            yield return LightLens();
            yield return EventsNode();
        }
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

    /// <summary>lens: quarter (CΨ) — CΨ(t) along the one trajectory, curve payload, start/min and the
    /// ¼-crossing times.</summary>
    private InspectableNode QuarterLens()
    {
        var cpsi = _cpsi!;
        var tGrid = _tGrid!;
        double start = cpsi[0];
        double min = cpsi.Min();
        var crossings = QuarterCrossingTimes();
        string crossClause = crossings.Count == 0
            ? "no ¼ crossing in window"
            : $"{crossings.Count} ¼ crossing(s) at t = " +
              string.Join(", ", crossings.Select(t => t.ToString("0.###", Inv)));
        return new InspectableNode("lens: quarter (CΨ)",
            summary: $"CΨ(0) = {start.ToString("0.####", Inv)}, min = {min.ToString("0.####", Inv)}; {crossClause}.",
            payload: new InspectablePayload.Curve("CΨ(t)", tGrid, cpsi, "t", "CΨ"));
    }

    /// <summary>lens: dose (K) — K = γ·t marks; the dose of the fold is K at the first ¼ crossing.</summary>
    private InspectableNode DoseLens()
    {
        var cross = FirstQuarterCrossing();
        string doseClause = cross is { } c
            ? $"the dose of the fold: K = {c.Dose.ToString("0.####", Inv)} at the first ¼ crossing (t={c.Time.ToString("0.###", Inv)})"
            : "no ¼ crossing in window, so no fold dose";
        return new InspectableNode("lens: dose (K)",
            summary: $"K = γ·t reaches {(Gamma * TMax).ToString("0.####", Inv)} at the window end; {doseClause}.",
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

        events.Sort((a, b) => a.Time.CompareTo(b.Time));
        return events;
    }

    /// <summary>The ¼-crossing times of CΨ(t), linearly interpolated between the bracketing grid
    /// points (both downward and upward crossings, in case a trajectory dips and recovers).</summary>
    private List<double> QuarterCrossingTimes(double threshold = 0.25)
    {
        EnsureLenses();
        var tGrid = _tGrid!;
        var cpsi = _cpsi!;
        var times = new List<double>();
        for (int s = 1; s < cpsi.Length; s++)
        {
            double a = cpsi[s - 1] - threshold;
            double b = cpsi[s] - threshold;
            if (a == 0.0) { times.Add(tGrid[s - 1]); continue; }
            if (a * b < 0.0)
            {
                double frac = a / (a - b);   // linear interpolation to the crossing
                times.Add(tGrid[s - 1] + frac * (tGrid[s] - tGrid[s - 1]));
            }
        }
        return times;
    }
}

/// <summary>Which initial density matrix the Symphony evolves. <see cref="BellPair"/> (default) is
/// Bell+ on sites (0,1) with the rest in |0⟩: the repo's hardware-confirmed ¼-crossing carrier
/// (F25, CΨ(0)=1/3), reproducing the F25 closed form at N=2. <see cref="SingleExcitation"/> is
/// |10…0⟩, one excitation on site 0.</summary>
public enum InitialStateKind
{
    BellPair,
    SingleExcitation,
}
