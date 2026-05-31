using System.Globalization;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Propagation;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The post-EP flow as a live Object-Manager GameObject: a single excitation on site 0
/// of an XY chain under Z-dephasing, evolved across a Q-grid, with per-site occupation ⟨n_site⟩(τ)
/// relaxing to the equipartitioned target 1/N. The C# home of <c>simulations/post_ep_dynamics_4d.py</c>
/// and the loop of <c>experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md</c>.
///
/// <para>Built on demand from (N, Q-grid, τ-grid); the trajectories compute lazily on first access.
/// A plain <see cref="IInspectable"/> (not a Claim): a live reading, not a typed-knowledge assertion.
/// The dimensionless Liouvillian L'(N,Q) = −iQ[H_unit,·] + Σ_l(Z_lρZ_l − ρ) reuses
/// <see cref="PauliDephasingDissipator.BuildZ"/> with H = Q·H_unit and γ = 1; Q = J/γ is the only knob,
/// τ = γ·t the dimensionless time. Site 0 = leftmost factor; vec is row-major (C-order), matching the
/// dissipator.</para>
///
/// <para>Numerical note: the trajectory is propagated in the Liouvillian's eigenbasis (eig once,
/// the same method as the Python prototype and <c>BlockCpsiTrajectory</c>). Exactly at an
/// exceptional point the Liouvillian is defective and that eigenbasis is singular; at the
/// generic grid points away from the exact EP it is well-conditioned, which is why the
/// validated Python reference scans the same Q-grid without issue.</para></summary>
public sealed class PostEpFlowField : IInspectable
{
    public int N { get; }
    public IReadOnlyList<double> QGrid { get; }
    public IReadOnlyList<double> TauGrid { get; }

    /// <summary>The relative per-site dephasing weights (length N, all > 0). Uniform [1,…,1] is
    /// the unshaped baseline. Use <see cref="NormalizeToTotal"/> to compare shapes at a fixed Σγ.</summary>
    public IReadOnlyList<double> GammaProfile { get; }

    /// <summary>The equipartitioned fixed point every trajectory relaxes to: ⟨n_site⟩(∞) = 1/N.</summary>
    public double Target => 1.0 / N;

    public PostEpFlowField(int N, IReadOnlyList<double> qGrid, IReadOnlyList<double> tauGrid,
        IReadOnlyList<double>? gammaProfile = null)
    {
        if (N < 1 || N > 6) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be 1..6 (dense 4^N Liouvillian)");
        QGrid = qGrid ?? throw new ArgumentNullException(nameof(qGrid));
        TauGrid = tauGrid ?? throw new ArgumentNullException(nameof(tauGrid));
        if (qGrid.Count == 0) throw new ArgumentException("qGrid must be non-empty", nameof(qGrid));
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
        {
            var xb = PauliString.SiteOp(N, b, PauliLetter.X);
            var xb1 = PauliString.SiteOp(N, b + 1, PauliLetter.X);
            var yb = PauliString.SiteOp(N, b, PauliLetter.Y);
            var yb1 = PauliString.SiteOp(N, b + 1, PauliLetter.Y);
            H = H + xb * xb1 + yb * yb1;
        }
        return _hUnit = H;
    }

    /// <summary>The dimensionless Liouvillian L'(N,Q) at the given Q (γ = 1, J = Q).</summary>
    public ComplexMatrix DimensionlessLiouvillian(double q)
    {
        return PauliDephasingDissipator.BuildZ(HUnit().Multiply(new Complex(q, 0.0)), GammaProfile);
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
            qFlows.Add(new PostEpQFlow(q, Underdamped: q >= 1.0, sites, SlowestNonKernelRate(ev.Eigenvalues)));
        }
        return qFlows;
    }

    /// <summary>The slowest non-kernel relaxation rate: −max{ Re λ : |λ| > tol }. Positive in every
    /// physical case (returns 0 only if L has no non-kernel mode at all); larger means faster
    /// forgetting. The kernel modes (|λ| ≈ 0, the 1/N fixed point) are excluded. This is the
    /// profile-sensitive timescale at fixed Σγ.</summary>
    private static double SlowestNonKernelRate(IReadOnlyList<Complex> eigenvalues)
    {
        const double tol = 1e-7;
        double maxRe = double.NegativeInfinity;
        foreach (var z in eigenvalues)
            if (z.Magnitude > tol && z.Real > maxRe) maxRe = z.Real;
        return double.IsNegativeInfinity(maxRe) ? 0.0 : -maxRe;
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

    // ---- Object Manager: the post-EP flow as a live IInspectable node ----
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    public string DisplayName =>
        $"PostEpFlowField (N={N}, target 1/N={Target.ToString("0.0000", Inv)})";
    public string Summary =>
        $"single excitation → 1/N; {QGrid.Count} Q × {N} sites, {TauGrid.Count} τ-points";
    public IEnumerable<IInspectable> Children
    {
        get
        {
            foreach (var qf in Flows)
            {
                string regime = qf.Underdamped ? "underdamped (hops, remembers)" : "overdamped (diffuses, forgets)";
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
                    summary: $"{regime}; slowest rate {qf.SlowestRate.ToString("0.0000", Inv)}",
                    children: siteLeaves);
            }
        }
    }
    public InspectablePayload Payload => InspectablePayload.Empty;
}

/// <summary>One site's occupation trajectory at one Q: the site index, edge/bulk class,
/// ⟨n_site⟩(τ) over the τ-grid, and the oscillation (turn) count.</summary>
public sealed record PostEpSiteFlow(int Site, bool IsEdge, IReadOnlyList<double> Occupation, int Turns);

/// <summary>One Q slice of the flow: the Q value, whether it is above the rotation onset
/// (underdamped, Q ≥ 1), and the per-site trajectories.</summary>
public sealed record PostEpQFlow(double Q, bool Underdamped, IReadOnlyList<PostEpSiteFlow> Sites, double SlowestRate);
