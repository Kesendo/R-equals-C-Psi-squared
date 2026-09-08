using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.F86;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Propagation;
using RCPsiSquared.Diagnostics.Ptf;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The single-excitation population flow as a live Object-Manager GameObject: an excitation on site 0
/// of an XY chain under Z-dephasing, evolved across a Q-grid, with per-site occupation ⟨n_site⟩(τ)
/// relaxing toward the equipartitioned target 1/N for the connected Q &gt; 0 flow. The C# home of <c>simulations/post_ep_dynamics_4d.py</c>
/// and the loop of <c>experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md</c>.
///
/// <para>Two distinct objects live on this node, the two parity sectors of the drain-depth / n_XY axis:
/// the EVEN-rung single-excitation FLOW (number-conserving, the occupation trajectories ⟨n⟩(τ)
/// relaxing toward 1/N) and the ODD-rung BIRTH CHANNEL (number-changing, the longest-lived coherence
/// the N=5 open-chain-only birth/sterile surface reads through <see cref="IsInBirthCanal"/>, with maximal saturation
/// <see cref="MaxSaturationCeiling"/> = ¼). The flow never enters the birth channel
/// (<see cref="FlowEntersBirthChannel"/> = false); they are separate readings, not one timescale.</para>
///
/// <para>Built on demand from (N, Q-grid, τ-grid); the trajectories compute lazily on first access.
/// A plain <see cref="IInspectable"/> (not a Claim): a live reading, not a typed-knowledge assertion.
/// The dimensionless Liouvillian L'(N,Q) = −i(Q/2)[H_unit,·] + Σ_l(Z_lρZ_l − ρ) reuses
/// <see cref="PauliDephasingDissipator.BuildZ"/> with H = (Q/2)·H_unit and γ = 1; Q is canonical carrier Q,
/// τ = γ·t the dimensionless time. Site 0 = leftmost factor; vec is row-major (C-order), matching the
/// dissipator.</para>
///
/// <para>Numerical note: the trajectory is propagated in the Liouvillian's eigenbasis (eig once,
/// the same method as the Python prototype and <c>BlockCpsiTrajectory</c>). The sampled population
/// trajectories make no spectral-transition claim; diagonalizability failures are numerical errors,
/// not silently interpreted as transition evidence.</para></summary>
public sealed class PostEpFlowField : IInspectable
{
    public int N { get; }
    public IReadOnlyList<double> QGrid { get; }
    public IReadOnlyList<double> TauGrid { get; }

    /// <summary>The relative per-site dephasing weights (length N, all > 0). Uniform [1,…,1] is
    /// the unshaped baseline. Use <see cref="NormalizeToTotal"/> to compare shapes at a fixed Σγ.</summary>
    public IReadOnlyList<double> GammaProfile { get; }

    /// <summary>The equipartitioned fixed point of the connected Q &gt; 0 single-excitation flow:
    /// ⟨n_site⟩(∞) = 1/N. Q=0 is rejected because every diagonal population is then stationary.</summary>
    public double Target => 1.0 / N;

    /// <summary>Chain (open line, the default) or Ring (closed, the extra wrap bond N−1 ↔ 0). The
    /// ring is the aromatic substrate: benzene C₆, cyclobutadiene C₄. N ≥ 3 for a meaningful ring.</summary>
    public FlowTopology Topology { get; }

    public PostEpFlowField(int N, IReadOnlyList<double> qGrid, IReadOnlyList<double> tauGrid,
        IReadOnlyList<double>? gammaProfile = null, FlowTopology topology = FlowTopology.Chain)
    {
        if (N < 1 || N > 6) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be 1..6 (dense 4^N Liouvillian)");
        QGrid = qGrid ?? throw new ArgumentNullException(nameof(qGrid));
        TauGrid = tauGrid ?? throw new ArgumentNullException(nameof(tauGrid));
        if (qGrid.Count == 0) throw new ArgumentException("qGrid must be non-empty", nameof(qGrid));
        if (qGrid.Any(q => !double.IsFinite(q) || q <= 0.0))
            throw new ArgumentException("qGrid entries must be finite and strictly > 0; Q=0 has no connected flow to 1/N", nameof(qGrid));
        if (tauGrid.Count < 2) throw new ArgumentException("tauGrid needs >= 2 points", nameof(tauGrid));

        if (gammaProfile is null)
        {
            GammaProfile = Enumerable.Repeat(1.0, N).ToArray();
        }
        else
        {
            if (gammaProfile.Count != N)
                throw new ArgumentException($"gammaProfile must have length N={N}, got {gammaProfile.Count}", nameof(gammaProfile));
            if (gammaProfile.Any(w => w <= 0.0))
                throw new ArgumentException("gammaProfile entries must be strictly > 0 (a zero-γ site is Fall 2, out of scope)", nameof(gammaProfile));
            GammaProfile = gammaProfile.ToArray();
        }
        this.N = N;
        Topology = topology;
    }

    /// <summary>Rescale a per-site weight profile so its entries sum to <paramref name="total"/>
    /// (use N for the fixed-Σγ, mean-1 comparison). Relative ratios are preserved.</summary>
    public static double[] NormalizeToTotal(IReadOnlyList<double> profile, double total)
    {
        if (profile is null) throw new ArgumentNullException(nameof(profile));
        if (profile.Count == 0) throw new ArgumentException("profile must be non-empty", nameof(profile));
        double sum = profile.Sum();
        if (sum <= 0.0) throw new ArgumentException("profile sum must be > 0", nameof(profile));
        double scale = total / sum;
        return profile.Select(w => w * scale).ToArray();
    }

    private ComplexMatrix? _hUnit;
    /// <summary>H_unit = Σ_b (X_b X_{b+1} + Y_b Y_{b+1}), coefficient 1 per term (NOT J/2).</summary>
    private ComplexMatrix HUnit()
    {
        if (_hUnit is not null) return _hUnit;
        int d = 1 << N;
        var H = ComplexMatrix.Build.Dense(d, d);
        for (int b = 0; b < N - 1; b++)
            H = H + Bond(b, b + 1);
        if (Topology == FlowTopology.Ring && N >= 3)
            H = H + Bond(N - 1, 0);          // the wrap bond closes the ring (aromatic substrate)
        return _hUnit = H;
    }

    /// <summary>The XY hopping on one bond: X_i X_j + Y_i Y_j, coefficient 1.</summary>
    private ComplexMatrix Bond(int i, int j)
    {
        var xi = PauliString.SiteOp(N, i, PauliLetter.X);
        var xj = PauliString.SiteOp(N, j, PauliLetter.X);
        var yi = PauliString.SiteOp(N, i, PauliLetter.Y);
        var yj = PauliString.SiteOp(N, j, PauliLetter.Y);
        return xi * xj + yi * yj;
    }

    /// <summary>The dimensionless Liouvillian L'(N,Q) at canonical carrier Q (γ=1,
    /// H=(Q/2)Σ(XX+YY), hop element Q).</summary>
    public ComplexMatrix DimensionlessLiouvillian(double q)
    {
        return PauliDephasingDissipator.BuildZ(HUnit().Multiply(new Complex(q / 2.0, 0.0)), GammaProfile);
    }

    /// <summary>vec(ρ₀) for the single excitation on site 0, row-major: only entry idx·d+idx = 1,
    /// idx = 2^(N−1).</summary>
    public ComplexVector InitialStateVec()
    {
        int d = 1 << N;
        int idx = 1 << (N - 1);
        var v = ComplexVector.Build.Dense(d * d);
        v[idx * d + idx] = Complex.One;
        return v;
    }

    /// <summary>Row-major observable covector for operator O: w[b·d+a] = O[a,b], so Tr(Oρ) = w·vec(ρ).
    /// The site observables n_l are diagonal (hence symmetric), so for them this orientation is
    /// immaterial; the convention is written generally for any future non-diagonal observable.</summary>
    private static ComplexVector Covector(ComplexMatrix o)
    {
        int d = o.RowCount;
        var w = ComplexVector.Build.Dense(d * d);
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                w[b * d + a] = o[a, b];
        return w;
    }

    private IReadOnlyList<ComplexVector>? _siteObservables;
    private IReadOnlyList<ComplexVector> SiteObservables()
    {
        if (_siteObservables is not null) return _siteObservables;
        int d = 1 << N;
        var ident = ComplexMatrix.Build.DenseIdentity(d);
        var list = new List<ComplexVector>(N);
        for (int l = 0; l < N; l++)
        {
            var zl = PauliString.SiteOp(N, l, PauliLetter.Z);
            var nL = (ident - zl).Multiply(new Complex(0.5, 0.0));   // n_l = (I − Z_l)/2
            list.Add(Covector(nL));
        }
        return _siteObservables = list;
    }

    private IReadOnlyList<PostEpQFlow>? _flows;
    /// <summary>Per-Q per-site occupation trajectories ⟨n_site⟩(τ), computed lazily.</summary>
    public IReadOnlyList<PostEpQFlow> Flows => _flows ??= ComputeFlows();

    private IReadOnlyList<PostEpQFlow> ComputeFlows()
    {
        var rho0 = InitialStateVec();
        var observables = SiteObservables();
        var qFlows = new List<PostEpQFlow>(QGrid.Count);
        foreach (double q in QGrid)
        {
            var L = DimensionlessLiouvillian(q);
            var ev = SpectralPropagator.EvolveWithSpectrum(L, rho0, observables, TauGrid);
            var traj = ev.Observables;   // traj[site][τ]
            var sites = new List<PostEpSiteFlow>(N);
            for (int s = 0; s < N; s++)
            {
                bool isEdge = s == 0 || s == N - 1;
                sites.Add(new PostEpSiteFlow(s, isEdge, traj[s], NTurns(traj[s])));
            }
            qFlows.Add(new PostEpQFlow(q, HasResolvedTurns: sites.Any(s => s.Turns > 1), sites,
                GlobalSlowestNonKernelRate: SlowestNonKernelRate(ev.Eigenvalues)));
        }
        return qFlows;
    }

    /// <summary>The slowest non-kernel relaxation rate: −max{ Re λ : |λ| > tol }. Positive in every
    /// physical case (returns 0 only if L has no non-kernel mode at all). The kernel modes
    /// (|λ| ≈ 0) are excluded. This is a GLOBAL spectral reading: the winning mode may be the
    /// odd number-changing channel and have zero overlap with this even single-excitation flow,
    /// so it is not called the population trajectory's approach rate.</summary>
    private static double SlowestNonKernelRate(IReadOnlyList<Complex> eigenvalues)
    {
        const double tol = 1e-7;
        double maxRe = double.NegativeInfinity;
        foreach (var z in eigenvalues)
            if (z.Magnitude > tol && z.Real > maxRe) maxRe = z.Real;
        return double.IsNegativeInfinity(maxRe) ? 0.0 : -maxRe;
    }

    // ---- The sterile zone vs the birth canal: where the closed form holds, and where the new is born ----
    // internal (not private) so BirthCanalSurfaceWitness shares these probe points and tolerance and
    // cannot drift from them; same assembly, so no InternalsVisibleTo is needed.
    internal const double BirthCanalProbeQLow = 3.0;
    internal const double BirthCanalProbeQHigh = 2000.0;
    internal const double BirthCanalTolerance = 1e-4;

    private double SlowestRateAt(double q)
    {
        var evd = DimensionlessLiouvillian(q).Evd();
        return SlowestNonKernelRate(evd.EigenValues.ToArray());
    }

    private double? _globalSlowestRateDeviation;

    /// <summary>The Q-drift of the global slowest non-kernel rate: rate(high Q) − rate(low Q).
    /// This numerical reading is available throughout the dense N≤6 domain, without assigning a
    /// birth/sterile mechanism to the winning sector.</summary>
    public double GlobalSlowestRateDeviation =>
        _globalSlowestRateDeviation ??=
            SlowestRateAt(BirthCanalProbeQHigh) - SlowestRateAt(BirthCanalProbeQLow);

    /// <summary>Whether the verified N=5 birth/sterile surface semantics apply. At N≥6 the global
    /// slowest mode can switch to an even density sector; use the sector-resolved junction witness.</summary>
    public bool HasBirthCanalClassification => N == 5 && Topology == FlowTopology.Chain;

    private void RequireBirthCanalClassification()
    {
        if (!HasBirthCanalClassification)
            throw new InvalidOperationException(
                $"Birth/sterile classification is verified only on the N=5 open-chain surface; " +
                $"N={N}, topology={Topology} has only the " +
                "unlabelled GlobalSlowestRateDeviation and may require a min-over-sectors classification.");
    }

    /// <summary>The N=5 surface's Q-drift of the global slowest rate. On that verified surface the
    /// winner is the longest-lived NUMBER-CHANGING coherence (one light quantum, pure odd n_XY
    /// parity for the pinned uniform/peaked-V/flat-bulk-edge profiles), the birth-capable channel
    /// whose content is bounded by the maximal saturation
    /// C_block ≤ ¼ (<see cref="MaxSaturationCeiling"/>, <see cref="BlockCoherenceContent"/>). It is
    /// NOT the even-rung single-excitation flow the rest of this class propagates: that flow is
    /// number-conserving and never enters this channel (<see cref="FlowEntersBirthChannel"/> is
    /// false). Zero in the sterile zone (the birth channel's lifetime is Q-independent, the closed
    /// form holds, nothing new couples); positive in the birth canal (its lifetime is Q-modulated,
    /// the loop where the coupling drives the birth-capable coherence and new structure can form).
    /// The exploration ruler: how deep in the N=5 birth canal. Throws outside N=5 rather than
    /// projecting those semantics onto a different winning sector.</summary>
    public double BirthCanalDeviation
    {
        get
        {
            RequireBirthCanalClassification();
            return GlobalSlowestRateDeviation;
        }
    }

    /// <summary>True in the birth canal: the Q-dependent loop where the longest-lived
    /// number-changing (odd-rung, birth-capable) coherence has its lifetime modulated by the
    /// coupling (the closed form fails and the new can be born). False in the sterile zone. Scoped
    /// to N=5.</summary>
    public bool IsInBirthCanal
    {
        get
        {
            RequireBirthCanalClassification();
            return Math.Abs(GlobalSlowestRateDeviation) > BirthCanalTolerance;
        }
    }

    /// <summary>The complement of <see cref="IsInBirthCanal"/>: the Q-independent corridor (the
    /// sterile zone) where the slowest rate is a closed form, the slow mode isolated, decay and
    /// rotation decoupled. Frozen, determined, no creation.</summary>
    public bool IsInSterileZone
    {
        get
        {
            RequireBirthCanalClassification();
            return !IsInBirthCanal;
        }
    }

    /// <summary>The closed-form slowest rate, valid ONLY in the sterile zone (the Q-independent
    /// corridor). In the birth canal it throws: there the rate is Q-dependent, so there is no single
    /// closed value, only a Q-trajectory and its high-Q limit. The closed form refuses to lie where
    /// the new is born; use the per-Q rate or <see cref="BirthCanalDeviation"/> to explore the canal.</summary>
    public double ClosedFormRate => IsInBirthCanal
        ? throw new InvalidOperationException(
            $"In the birth canal (Q-dependent, deviation {BirthCanalDeviation.ToString("E2", Inv)}): no closed-form " +
            "rate exists here, only a Q-trajectory and its high-Q limit. The closed form holds only in the sterile " +
            "zone (IsInSterileZone). Use the per-Q rate or BirthCanalDeviation to explore the birth canal.")
        : SlowestRateAt(BirthCanalProbeQHigh);

    // ---- The birth channel made measurable: the between-block (number-changing) coherence ----

    /// <summary>The universal maximal-saturation ceiling on between-block coherence: any density
    /// matrix obeys C_block ≤ 1/4 (<see cref="BlockCoherenceContent"/>, Theorem 2 of
    /// PROOF_BLOCK_CPSI_QUARTER, the bilinear apex ¼ of p·(1−p)). This is the most NUMBER-CHANGING
    /// (odd-rung, birth-capable) coherence a state can hold; the Dicke superposition
    /// (|D_n⟩+|D_{n+1}⟩)/√2 saturates it. It is the ceiling of the very channel the
    /// <see cref="IsInBirthCanal"/> members track.</summary>
    public static double MaxSaturationCeiling => BlockCoherenceContent.Quarter;

    private double? _peakBetweenBlockSaturation;
    /// <summary>The peak between-block coherence the flow's state reaches over the whole Q × τ grid:
    /// max over τ, over all blocks n ∈ [0, N−1], of C_block(ρ(τ), n) (the (popcount-n, popcount-(n+1))
    /// block). This population flow is a single excitation, a NUMBER-CONSERVING (even-rung) object, so it
    /// carries no between-block coherence: this is ≈ 0 (machine zero), the LIVE measurement that the
    /// flow never enters the birth channel. Between-block coherence (the odd rung, bounded by
    /// <see cref="MaxSaturationCeiling"/> = ¼) is the birth-capable channel; reaching it needs a
    /// number-uncertain initial state (a Dicke superposition), not this flow. Lazy; reconstructs
    /// ρ(τ) via <see cref="SpectralPropagator.EvolveStateVectors"/>.</summary>
    public double PeakBetweenBlockSaturation => _peakBetweenBlockSaturation ??= ComputePeakBetweenBlockSaturation();

    /// <summary>True if the flow's state ever develops between-block (number-changing) coherence
    /// above 1e-9, i.e. enters the birth channel. False for the single-excitation flow (even-rung,
    /// number-conserving): it stays out of the birth channel by construction.</summary>
    public bool FlowEntersBirthChannel => PeakBetweenBlockSaturation > 1e-9;

    private double ComputePeakBetweenBlockSaturation()
    {
        int d = 1 << N;
        var rho0 = InitialStateVec();
        double peak = 0.0;
        foreach (double q in QGrid)
        {
            var states = SpectralPropagator.EvolveStateVectors(DimensionlessLiouvillian(q), rho0, TauGrid);
            foreach (var vt in states)
            {
                var rho = ComplexMatrix.Build.Dense(d, d);
                for (int i = 0; i < d; i++)
                    for (int j = 0; j < d; j++)
                        rho[i, j] = vt[i * d + j];        // row-major vec: ρ[i,j] = vec[i·d+j]
                for (int n = 0; n < N; n++)
                {
                    double cBlock = BlockCoherenceContent.Compute(rho, n);
                    if (cBlock > peak) peak = cBlock;
                }
            }
        }
        return peak;
    }

    /// <summary>The slowest non-kernel mode read through the WHOLE assembly at one Q: the scattered
    /// pieces brought together where the birth channel lives. Reads its spectral-edge rate and
    /// the tolerance-cluster mean separately; its drain depth
    /// ⟨popcount(i⊕j)⟩ = the light content n_XY (one axis, two names); the per-site light (the
    /// carrier vector); the Absorption-Theorem cross-check rate = 2·Σ_k γ_k·light_k (equal to the
    /// cluster mean, not necessarily the spectral edge); the parity-resolved projector support;
    /// and the maximal-saturation ceiling
    /// ¼ of the odd rail. The parity can be Mixed when the slow-rate cluster spans both exact
    /// sectors; it is not inferred from the generally noninteger mean depth.</summary>
    public FlowAssemblyReading ReadAssembly(double q)
    {
        // Basis-free since 2026-06-10: the per-site carrier is read through the orthogonal
        // projector onto the slow invariant subspace (range of the biorthogonal spectral
        // projector over the slowest-rate cluster), see SlowLightDistribution. The former
        // tolerance-cluster-averaged carrier and its gauge caveat (the |v|² average over
        // non-orthogonal eigenvectors stayed mildly basis-dependent) are retired:
        // the projector reading is frame-independent for the selected invariant subspace, keeps the absorption
        // identity machine-exact for the cluster mean, and coincides with the old average wherever
        // ClusterDimension = 1.
        var slow = SlowLightDistribution.Compute(DimensionlessLiouvillian(q), N, GammaProfile);
        return new FlowAssemblyReading(slow.SpectralEdgeRate, slow.ClusterMeanRate, slow.TotalLight,
            slow.PerSiteLight, slow.AbsorptionRate, slow.Parity, MaxSaturationCeiling,
            slow.ClusterDimension);
    }

    /// <summary>Oscillation count: strict sign changes of the first difference (local extrema).
    /// A display heuristic for the over/underdamped tag, not a bit-for-bit match of the Python
    /// prototype's n_turns (which also counts transitions through flat segments).</summary>
    private static int NTurns(IReadOnlyList<double> ys)
    {
        int turns = 0;
        for (int i = 1; i + 1 < ys.Count; i++)
        {
            int d1 = Math.Sign(ys[i] - ys[i - 1]);
            int d2 = Math.Sign(ys[i + 1] - ys[i]);
            if (d1 != 0 && d2 != 0 && d1 != d2) turns++;
        }
        return turns;
    }

    // ---- Object Manager: the single-excitation population flow as a live IInspectable node ----
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    public string DisplayName =>
        $"PostEpFlowField (N={N}, target 1/N={Target.ToString("0.0000", Inv)})";
    public string Summary =>
        $"connected-Q single-excitation flow toward 1/N; {QGrid.Count} Q × {N} sites, {TauGrid.Count} τ-points";
    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var qf in Flows)
            {
                string regime = qf.HasResolvedTurns ? "sampled population trace has resolved turns" : "no resolved turns on sampled grid";
                var siteLeaves = new List<IInspectable>(N);
                foreach (var s in qf.Sites)
                {
                    string cls = s.IsEdge ? "edge" : "bulk";
                    string tag = s.Turns <= 1 ? "monotone" : $"oscillates ({s.Turns} turns)";
                    string range = $"{s.Occupation[0].ToString("0.00", Inv)}→{s.Occupation[^1].ToString("0.00", Inv)}";
                    siteLeaves.Add(new InspectableNode(
                        displayName: $"site {s.Site} ({cls})",
                        summary: $"{range}  {tag}",
                        payload: new InspectablePayload.Curve("⟨n⟩", TauGrid, s.Occupation, "τ", "⟨n⟩")));
                }
                yield return new InspectableNode(
                    displayName: $"Q={qf.Q.ToString("0.00", Inv)}",
                    summary: $"{regime}; global slowest non-kernel rate {qf.GlobalSlowestNonKernelRate.ToString("0.0000", Inv)} (not a flow-overlap rate)",
                    children: siteLeaves);
            }
        }
    }
    public InspectablePayload Payload => InspectablePayload.Empty;
}

/// <summary>The bond graph of the flow: an open Chain (line) or a closed Ring (the aromatic
/// substrate, with the wrap bond N−1 ↔ 0). Star, complete, etc. are not flow substrates here.</summary>
public enum FlowTopology { Chain, Ring }

/// <summary>One site's occupation trajectory at one Q: the site index, edge/bulk class,
/// ⟨n_site⟩(τ) over the τ-grid, and the oscillation (turn) count.</summary>
public sealed record PostEpSiteFlow(int Site, bool IsEdge, IReadOnlyList<double> Occupation, int Turns);

/// <summary>One Q slice of the flow: canonical Q, whether the finite sampled population traces
/// contain resolved turns, the per-site trajectories, and the global slowest non-kernel rate.
/// The global winning mode need not overlap this even single-excitation preparation/readout.
/// The turn flag is a display reading, not a critical-damping or spectral-transition verdict.</summary>
public sealed record PostEpQFlow(double Q, bool HasResolvedTurns, IReadOnlyList<PostEpSiteFlow> Sites, double GlobalSlowestNonKernelRate);

/// <summary>The slowest non-kernel rate read through the whole assembly, the scattered pieces in
/// one place: <paramref name="SlowestRate"/> (−max Re λ); <paramref name="SlowClusterMeanRate"/>
/// = −(1/g)Σ Re λ over the selected tolerance cluster, the rate reproduced by the projector-trace
/// Absorption identity; <paramref name="SlowestDepth"/> = the drain
/// depth ⟨popcount(i⊕j)⟩ = the light n_XY; <paramref name="PerSiteLight"/> = the carrier vector
/// ⟨X/Y at k⟩, BASIS-FREE via the orthogonal projector onto the slow invariant subspace
/// (<see cref="RCPsiSquared.Diagnostics.Ptf.SlowLightDistribution"/>, 2026-06-10; it equals the
/// per-eigenvector reading wherever <paramref name="SlowClusterDimension"/> = 1);
/// <paramref name="AbsorptionRate"/> = 2·Σ_k γ_k·light_k (the Absorption
/// Theorem, equal to the cluster-mean rate); <paramref name="Parity"/> is derived from the slow projector's
/// exact even/odd support, not from its generally noninteger mean depth; <see cref="Rung"/> and
/// <see cref="OnBirthRail"/> are undefined when that support is mixed;
/// <paramref name="MaxSaturationCeiling"/> = ¼, the
/// rail's bilinear ceiling; <paramref name="SlowClusterDimension"/> = the number of modes selected
/// by the real-rate tolerance, not an assertion of exact degeneracy.</summary>
public sealed record FlowAssemblyReading(
    double SlowestRate,
    double SlowClusterMeanRate,
    double SlowestDepth,
    IReadOnlyList<double> PerSiteLight,
    double AbsorptionRate,
    SlowLightParity Parity,
    double MaxSaturationCeiling,
    int SlowClusterDimension)
{
    public double AbsorptionResidual => Math.Abs(SlowClusterMeanRate - AbsorptionRate);

    public int? Rung => Parity switch
    {
        SlowLightParity.Even => 0,
        SlowLightParity.Odd => 1,
        _ => null,
    };

    public bool? OnBirthRail => Parity switch
    {
        SlowLightParity.Even => false,
        SlowLightParity.Odd => true,
        _ => null,
    };
}
