using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Factorization;
using RCPsiSquared.Core.BlockSpectrum;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The carrier-as-vector reading of decoherence, held as typed objects so the
/// whole structure stands at once: a per-channel rate vector γ, every Liouvillian mode as
/// a per-channel difference-portfolio, and the Absorption-Theorem rate law that pairs them.
///
/// <para>The law (see <see cref="RCPsiSquared.Core.Symmetry.AbsorptionTheoremClaim"/>,
/// per-channel reading): for any Liouvillian mode k under site/channel-dependent
/// Z-dephasing with rate vector γ = (γ_x),</para>
///
/// <code>
///   −Re(λ_k) = 2 · Σ_x γ_x · ⟨Δ_x⟩_k
/// </code>
///
/// <para>where ⟨Δ_x⟩_k ∈ [0,1] is how much of its difference the mode stores in channel x
/// (on a computational coherence |A⟩⟨B|, simply whether bra and ket disagree at x). The
/// uniform-γ case collapses Σ_x γ_x·⟨Δ_x⟩ to γ·popcount; a non-uniform γ-vector lifts that
/// degeneracy and the portfolio sets the lifetime hierarchy directly. Exact for any
/// Hermitian H, because Herm(L) is the dephasing dissipator alone (the Hamiltonian is
/// anti-Hermitian and drops out of Re(λ)); coupling only turns the sharp bit Δ_x ∈ {0,1}
/// into the expectation ⟨Δ_x⟩ ∈ [0,1].</para>
///
/// <para>A "channel" is any independent dephasing axis: a qubit site, or a physical
/// channel of a substrate (electron spin, nuclear spin, charge, valley for a ³¹P-in-silicon
/// donor, whose clocks span decades). Companion reading to <see cref="MemoryAxisRho"/>:
/// MemoryAxisRho splits one state along the static / memory axis; this splits every mode
/// along the per-channel difference axis. Verified bit-exact in
/// <c>simulations/absorption_gamma_vector.py</c> and <c>simulations/sip_carrier_channels.py</c>.</para>
/// </summary>
public sealed record ChannelRate(string Channel, double Gamma);

/// <summary>How much of its difference a mode stores in one channel: ⟨Δ_x⟩ ∈ [0,1].</summary>
public sealed record ChannelActivity(string Channel, double Delta);

/// <summary>The dephasing carrier as a per-channel rate vector γ = (γ_x). The scalar
/// "carrier" γ₀ of the uniform reading is the special case where every entry is equal.</summary>
public sealed record CarrierVector(IReadOnlyList<ChannelRate> Channels)
{
    /// <summary>The dephasing rate of one channel; throws if the channel is unknown.</summary>
    public double Gamma(string channel) =>
        Channels.FirstOrDefault(c => c.Channel == channel)?.Gamma
        ?? throw new KeyNotFoundException($"no channel '{channel}' in carrier vector");

    /// <summary>The channel names, in declared order.</summary>
    public IReadOnlyList<string> ChannelNames => Channels.Select(c => c.Channel).ToList();

    /// <summary>True iff all channels carry the same clock (the uniform-γ / popcount case).</summary>
    public bool IsUniform(double tol = 1e-12) =>
        Channels.Count == 0 || Channels.All(c => Math.Abs(c.Gamma - Channels[0].Gamma) < tol);

    /// <summary>Σ_x γ_x: the source of the maximum rate 2·Σγ (the full-difference mode).</summary>
    public double TotalGamma => Channels.Sum(c => c.Gamma);
}

/// <summary>A Liouvillian mode's per-channel difference-portfolio: the percent of its
/// difference stored in each channel, ⟨Δ_x⟩ ∈ [0,1].</summary>
public sealed record ChannelDifferencePortfolio(IReadOnlyList<ChannelActivity> Activity)
{
    /// <summary>The dissipator coefficient (the "2" in 2γ·⟨n_XY⟩, the d=2
    /// root of the Pi2 ladder). It is not an eigenvalue-grid quantum.</summary>
    public const double DissipatorCoefficient = 2.0;

    /// <summary>The carrier-vector rate law: −Re(λ) = 2·Σ_x γ_x·⟨Δ_x⟩. Validates each
    /// activity lies in [0,1].</summary>
    public double DecayRate(CarrierVector carrier)
    {
        if (carrier is null) throw new ArgumentNullException(nameof(carrier));
        double weighted = 0.0;
        foreach (var a in Activity)
        {
            if (a.Delta < -1e-9 || a.Delta > 1.0 + 1e-9)
                throw new ArgumentOutOfRangeException(
                    nameof(Activity), a.Delta, $"⟨Δ⟩ for channel '{a.Channel}' must lie in [0,1]");
            weighted += carrier.Gamma(a.Channel) * a.Delta;
        }
        return DissipatorCoefficient * weighted;
    }

    /// <summary>Σ_x ⟨Δ_x⟩: the unweighted total activity (the uniform-γ popcount, when each
    /// activity is a sharp 0/1 bit).</summary>
    public double TotalActivity => Activity.Sum(a => a.Delta);

    /// <summary>⟨Δ⟩ stored in one channel; 0 if the channel is absent from the portfolio.</summary>
    public double ActivityIn(string channel) =>
        Activity.FirstOrDefault(a => a.Channel == channel)?.Delta ?? 0.0;

    /// <summary>The channel that contributes most to the decay rate (γ_x·⟨Δ_x⟩ largest):
    /// where this mode mostly loses its coherence.</summary>
    public string DominantChannel(CarrierVector carrier) =>
        Activity.OrderByDescending(a => carrier.Gamma(a.Channel) * a.Delta).First().Channel;

    /// <summary>True iff the mode stores (almost) no difference in the given fast channels:
    /// ⟨Δ_x⟩ ≈ 0 for each. This is the sweet-spot / clock-transition condition, hold the
    /// bra-ket agreement in the fast channels and only the slow clocks remain.</summary>
    public bool IsProtectedFrom(IEnumerable<string> fastChannels, double tol = 1e-9) =>
        fastChannels.All(c => ActivityIn(c) < tol);
}

/// <summary>One mode: its measured decay rate −Re(λ) and its difference-portfolio.
///
/// <para><see cref="Portfolio"/> is read off the one eigenvector the eigensolver returned. Where λ
/// is degenerate that vector is one arbitrary member of the eigenspace: it obeys the rate law
/// exactly, but its split across the channels is a basis choice (on a translation-invariant ring
/// it can show 26/4/4/26/40 % where symmetry demands 20 % each). The basis-free reading is
/// <see cref="EigenspacePortfolio"/>: the trace of each channel's disagreement projector over an
/// orthonormal basis of the eigenspace, divided by its dimension. Read that one whenever the
/// question is where a rate lives.</para></summary>
public sealed record CarrierMode(double ActualDecayRate, ChannelDifferencePortfolio Portfolio)
{
    /// <summary>Im(λ), the mode's oscillation frequency ω (signed; modes come in ±ω
    /// conjugate pairs). 0 for non-oscillating modes. The angular companion to
    /// ActualDecayRate = −Re(λ).</summary>
    public double OscillationFrequency { get; init; }

    private readonly ChannelDifferencePortfolio? _eigenspacePortfolio;

    /// <summary>The portfolio averaged over this mode's whole λ-eigenspace, basis-free; equal to
    /// <see cref="Portfolio"/> when λ is simple, and <see cref="Portfolio"/> itself on a mode built
    /// by hand rather than by a decomposition. Obeys the same rate law, since every vector in the
    /// eigenspace does.</summary>
    public ChannelDifferencePortfolio EigenspacePortfolio
    {
        get => _eigenspacePortfolio ?? Portfolio;
        init => _eigenspacePortfolio = value;
    }

    /// <summary>The dimension of this mode's λ-eigenspace (the geometric multiplicity, summed
    /// over the joint-popcount blocks that carry λ); 1 when λ is simple.</summary>
    public int EigenspaceDimension { get; init; } = 1;

    /// <summary>How many modes the eigensolver returned at this λ (the algebraic multiplicity).
    /// Larger than <see cref="EigenspaceDimension"/> exactly when λ is defective, an exceptional
    /// point, where the eigenvectors fail to span.</summary>
    public int EigenvalueMultiplicity { get; init; } = 1;
}

/// <summary>The whole Liouvillian spectrum read as portfolios at once: the carrier vector
/// plus every mode's portfolio, with the summary the reading raises, the lifetime
/// hierarchy, the coupling-independent decoherence budget, and the law residual.</summary>
public sealed record CarrierPortfolioSpectrum(CarrierVector Carrier, IReadOnlyList<CarrierMode> Modes)
{
    /// <summary>The slowest non-zero rate: the longest-lived coherence (the memory floor).</summary>
    public double SlowestRate =>
        Modes.Where(m => m.ActualDecayRate > 1e-12).Select(m => m.ActualDecayRate)
             .DefaultIfEmpty(0.0).Min();

    /// <summary>The fastest rate: the first coherence flushed (up to 2·Σγ).</summary>
    public double FastestRate => Modes.Select(m => m.ActualDecayRate).DefaultIfEmpty(0.0).Max();

    /// <summary>Σ_k |Re(λ_k)|, the total decoherence budget. Coupling-independent: Herm(L)
    /// is the dephasing dissipator alone, so the Hamiltonian redistributes rates among modes
    /// but never changes their sum.</summary>
    public double Budget => Modes.Sum(m => m.ActualDecayRate);

    /// <summary>Largest mismatch between a mode's measured rate and the carrier-vector law
    /// applied to its portfolio. Machine-zero when the portfolios are the true ⟨Δ_x⟩
    /// activities (the bit-exact check the Python verifiers run).</summary>
    public double MaxLawResidual =>
        Modes.Count == 0 ? 0.0 : Modes.Max(m => Math.Abs(m.ActualDecayRate - m.Portfolio.DecayRate(Carrier)));
}

/// <summary>Reads a Liouvillian into a <see cref="CarrierPortfolioSpectrum"/>: the numbers,
/// not a guess. Builds L = −i[H, ·] + Σ_l γ_l·D[Z_l] densely for the given Hamiltonian and
/// per-site dephasing rates, eigendecomposes it, and reads each mode's per-site difference
/// activity ⟨Δ_l⟩ = ⟨v|N_l|v⟩ / ‖v‖² with N_l = (I − Z_l⊗Z_l)/2 the "site-l off-diagonal"
/// projector. N_l is diagonal (Z_l is), and its diagonal , the bra-ket disagreement bit Δ_l ,
/// is read off the same siteZ that builds L, so the reading is convention-safe by construction
/// and ⟨v|N_l|v⟩ collapses to the weighted sum Σ_i Δ_l[i]·|v_i|² (no per-mode matrix-vector
/// product): Re(λ_k) = ⟨v_k|Herm(L)|v_k⟩/‖v_k‖² = −2 Σ_l γ_l·⟨Δ_l⟩_k holds for the right
/// eigenvectors regardless of the vec ordering. Diagnostic-scale (dense 4^N): for the small N
/// where we look, not the large-N engine.</summary>
public static class CarrierVectorPortfolio
{
    private static readonly ComplexMatrix PauliZ =
        ComplexMatrix.Build.DenseOfArray(new Complex[,] { { 1, 0 }, { 0, -1 } });

    /// <summary>Decompose the Liouvillian for an N-site system with Hamiltonian
    /// <paramref name="hamiltonian"/> (a 2^N × 2^N matrix) under per-site Z-dephasing given
    /// by <paramref name="siteChannels"/> (one <see cref="ChannelRate"/> per site, in site
    /// order). Returns every mode as a per-site difference-portfolio, all at once.</summary>
    public static CarrierPortfolioSpectrum Decompose(
        int N, ComplexMatrix hamiltonian, IReadOnlyList<ChannelRate> siteChannels)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be >= 1");
        if (hamiltonian is null) throw new ArgumentNullException(nameof(hamiltonian));
        if (siteChannels is null) throw new ArgumentNullException(nameof(siteChannels));
        if (siteChannels.Count != N)
            throw new ArgumentException($"need one channel per site (N={N}), got {siteChannels.Count}", nameof(siteChannels));
        int d = 1 << N;
        if (hamiltonian.RowCount != d || hamiltonian.ColumnCount != d)
            throw new ArgumentException($"Hamiltonian must be {d}×{d} for N={N}", nameof(hamiltonian));

        var id1 = ComplexMatrix.Build.DenseIdentity(2);
        var idD = ComplexMatrix.Build.DenseIdentity(d);
        var idD2 = ComplexMatrix.Build.DenseIdentity(d * d);

        ComplexMatrix SiteZ(int site)
        {
            ComplexMatrix op = site == 0 ? PauliZ : id1;
            for (int s = 1; s < N; s++)
                op = op.KroneckerProduct(s == site ? PauliZ : id1);
            return op;
        }

        var siteZ = Enumerable.Range(0, N).Select(SiteZ).ToArray();

        // N_l = (I − Z_l⊗Z_l)/2 is DIAGONAL: Z_l is diagonal, so Z_l⊗Z_l is, so N_l is. Its
        // diagonal is exactly the bra-ket disagreement bit Δ_l[a·d+b] = [Z_l(a) ≠ Z_l(b)] ∈ {0,1}.
        // Holding the diagonal (length 4^N) instead of the 4^N×4^N matrix turns ⟨v|N_l|v⟩ into a
        // weighted sum over |v_i|² , no per-mode matrix-vector product, no N dense 4^N operators
        // (the difference between minutes and milliseconds at N=6).
        int d2 = d * d;
        var zDiag = siteZ.Select(z => Enumerable.Range(0, d).Select(a => z[a, a].Real).ToArray()).ToArray();
        var nDiag = new double[N][];
        for (int l = 0; l < N; l++)
        {
            var dz = zDiag[l];
            var diag = new double[d2];
            for (int a = 0; a < d; a++)
                for (int b = 0; b < d; b++)
                    diag[a * d + b] = dz[a] * dz[b] < 0 ? 1.0 : 0.0;  // 1 iff bra and ket disagree at site l
            nDiag[l] = diag;
        }

        // L = −i (H⊗I − I⊗Hᵀ) + Σ_l γ_l (Z_l⊗Z_l − I).
        var L = (hamiltonian.KroneckerProduct(idD) - idD.KroneckerProduct(hamiltonian.Transpose()))
                .Multiply(-Complex.ImaginaryOne);
        for (int l = 0; l < N; l++)
            L += (siteZ[l].KroneckerProduct(siteZ[l]) - idD2).Multiply((Complex)siteChannels[l].Gamma);

        var evd = L.Evd();
        var vals = evd.EigenValues;
        var vecs = evd.EigenVectors;

        var ledger = new EigenspaceLedger(N);
        ledger.AddBlock(vals, vecs, (i, l) => nDiag[l][i] != 0.0);

        var modes = new List<CarrierMode>(d2);
        var absSq = new double[d2];
        for (int k = 0; k < d2; k++)
        {
            var v = vecs.Column(k);
            double normSq = 0.0;
            for (int i = 0; i < d2; i++)
            {
                var c = v[i];
                double m = c.Real * c.Real + c.Imaginary * c.Imaginary;
                absSq[i] = m;
                normSq += m;
            }
            var activity = new List<ChannelActivity>(N);
            for (int l = 0; l < N; l++)
            {
                var diag = nDiag[l];
                double w = 0.0;
                for (int i = 0; i < d2; i++) w += diag[i] * absSq[i];
                activity.Add(new ChannelActivity(siteChannels[l].Channel, w / normSq));
            }
            modes.Add(new CarrierMode(-vals[k].Real, new ChannelDifferencePortfolio(activity)) { OscillationFrequency = vals[k].Imaginary });
        }

        return new CarrierPortfolioSpectrum(new CarrierVector(siteChannels), ledger.Finish(modes, siteChannels));
    }

    /// <summary>True iff H preserves total popcount (Sz): H[a,b] = 0 whenever popcount(a) ≠
    /// popcount(b). For such H (XY, Heisenberg) the Z-dephasing Liouvillian is block-diagonal in
    /// the joint-popcount sectors, so <see cref="DecomposeBlocked"/> applies; otherwise (e.g. a
    /// non-truly Σ X_lY_{l+1}) the sectors are not closed and the full <see cref="Decompose"/> is
    /// required.</summary>
    public static bool IsNumberConserving(ComplexMatrix hamiltonian, double tol = 1e-12)
    {
        if (hamiltonian is null) throw new ArgumentNullException(nameof(hamiltonian));
        int d = hamiltonian.RowCount;
        for (int a = 0; a < d; a++)
            for (int b = 0; b < d; b++)
                if (BitOperations.PopCount((uint)a) != BitOperations.PopCount((uint)b)
                    && hamiltonian[a, b].Magnitude > tol)
                    return false;
        return true;
    }

    /// <summary>Same reading as <see cref="Decompose"/>, but block-by-block: the Z-dephasing
    /// Liouvillian of a number-conserving H is exactly block-diagonal in the joint-popcount sectors
    /// (popcount(row), popcount(col)), so one 4^N non-Hermitian eigendecomposition (an inherently
    /// sequential geev) becomes (N+1)² small ones. Bit-identical spectrum and law-consistent
    /// portfolios to <see cref="Decompose"/>; the win is wall-clock at large N.
    ///
    /// <para>Convention: <see cref="PerBlockLiouvillianBuilder.BuildBlockZ"/> takes a per-SITE
    /// carrier and does the site-to-bit re-indexing itself, and this reading (like
    /// <see cref="Decompose"/>) labels channel s at bit (N−1−s), the same book. So the carrier is
    /// handed over as it stands, and the law −Re λ = 2·Σ_s γ_s·⟨Δ_s⟩ is exact for any profile,
    /// palindromic or not. It passed the carrier REVERSED until 2026-08-20; see the inline note at
    /// the handover below.</para>
    ///
    /// <para>Requires <see cref="IsNumberConserving"/>; throws otherwise (the sectors would not be
    /// closed and the per-block spectra would be wrong).</para></summary>
    public static CarrierPortfolioSpectrum DecomposeBlocked(
        int N, ComplexMatrix hamiltonian, IReadOnlyList<ChannelRate> siteChannels)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be >= 1");
        if (hamiltonian is null) throw new ArgumentNullException(nameof(hamiltonian));
        if (siteChannels is null) throw new ArgumentNullException(nameof(siteChannels));
        if (siteChannels.Count != N)
            throw new ArgumentException($"need one channel per site (N={N}), got {siteChannels.Count}", nameof(siteChannels));
        int d = 1 << N;
        if (hamiltonian.RowCount != d || hamiltonian.ColumnCount != d)
            throw new ArgumentException($"Hamiltonian must be {d}×{d} for N={N}", nameof(hamiltonian));
        if (!IsNumberConserving(hamiltonian))
            throw new ArgumentException("DecomposeBlocked needs a number-conserving H (joint-popcount blocks); use Decompose otherwise", nameof(hamiltonian));

        // BuildBlockZ takes a per-SITE array and does the site-to-bit re-indexing itself, so the
        // carrier is handed over as it stands. The read below then takes channel s at bit N-1-s,
        // which is the same site convention, and the two agree for any profile.
        //
        // This reversed the carrier from 2026-05-29 to 2026-08-20, and it was CORRECT WHEN WRITTEN:
        // BuildBlockZ then read its argument as bit-indexed, and reversing here was what made
        // channel s carry gamma_s. Commit f3a58a0 (2026-07-25) moved the re-indexing inside
        // BuildBlockZ and left this bridge compensating for a convention that no longer existed,
        // so channel s carried gamma_{N-1-s} instead. Nothing in production saw it: every
        // MirrorSystem built by InspectCommand or Symphony passes a uniform profile, which is
        // invariant under the reversal. Six cases in MirrorSystemTests DO pass non-palindromic
        // carriers and stayed green anyway, on spectrum-level assertions that a site reversal
        // cannot disturb. The one gate that could see it, the non-palindromic {0.05, 0.13, 0.20} case of
        // DecomposeBlocked_MatchesFullDecompose_ForNumberConservingH, was red at a law residual of
        // 1.353e-3 for those 26 days, in a test project CLAUDE.md forbids running whole.
        //
        // It is the second orphan of that one commit; the first, repaired the same day, was
        // JwSlaterPairSparseLBuilder. Both are recorded in docs/CAUGHT_ERRORS.md under 2026-08-20.
        var gammasForBlock = new double[N];
        for (int l = 0; l < N; l++) gammasForBlock[l] = siteChannels[l].Gamma;

        // Partition the 4^N flat indices into closed sectors by (popcount(row), popcount(col)).
        int d2 = d * d;
        var sectors = new Dictionary<(int, int), List<int>>();
        for (int f = 0; f < d2; f++)
        {
            var key = (BitOperations.PopCount((uint)(f / d)), BitOperations.PopCount((uint)(f % d)));
            if (!sectors.TryGetValue(key, out var list)) { list = new List<int>(); sectors[key] = list; }
            list.Add(f);
        }

        var modes = new List<CarrierMode>(d2);
        var ledger = new EigenspaceLedger(N);
        foreach (var flat in sectors.Values)
        {
            var B = PerBlockLiouvillianBuilder.BuildBlockZ(hamiltonian, gammasForBlock, flat);
            var evd = B.Evd();
            var vals = evd.EigenValues;
            var vecs = evd.EigenVectors;
            int bs = flat.Count;
            ledger.AddBlock(vals, vecs, (i, s) =>
            {
                int bit = N - 1 - s;
                return (((flat[i] / d) >> bit) & 1) != (((flat[i] % d) >> bit) & 1);
            });
            var absSq = new double[bs];
            for (int k = 0; k < bs; k++)
            {
                var v = vecs.Column(k);
                double normSq = 0.0;
                for (int i = 0; i < bs; i++)
                {
                    var c = v[i];
                    double m = c.Real * c.Real + c.Imaginary * c.Imaginary;
                    absSq[i] = m;
                    normSq += m;
                }
                var activity = new List<ChannelActivity>(N);
                for (int s = 0; s < N; s++)
                {
                    int bit = N - 1 - s;
                    double w = 0.0;
                    for (int i = 0; i < bs; i++)
                    {
                        int f = flat[i];
                        if ((((f / d) >> bit) & 1) != (((f % d) >> bit) & 1)) w += absSq[i];
                    }
                    activity.Add(new ChannelActivity(siteChannels[s].Channel, w / normSq));
                }
                modes.Add(new CarrierMode(-vals[k].Real, new ChannelDifferencePortfolio(activity)) { OscillationFrequency = vals[k].Imaginary });
            }
        }

        return new CarrierPortfolioSpectrum(new CarrierVector(siteChannels), ledger.Finish(modes, siteChannels));
    }

    /// <summary>Collects, block by block, each eigenvalue cluster's eigenspace, orthonormalized,
    /// and the trace of every channel's disagreement projector over it; then merges equal λ across
    /// blocks (their supports are disjoint, so traces and dimensions add) and stamps each mode with
    /// the basis-free <see cref="CarrierMode.EigenspacePortfolio"/>.
    ///
    /// <para>Error model of the two thresholds. Eigenvalues are clustered when
    /// |λ − λ'| ≤ 1e-6·max(1, |λ|, |λ'|): a diagonalizable degenerate λ comes back from geev split
    /// by about eps·cond, near 1e-14, a defective EP2 by about √eps, near 1e-8, both far inside.
    /// An EPm splits by about eps^(1/m)·‖L‖, which for m ≥ 3 can exceed the threshold: such a
    /// cluster may fall apart into separate pieces, each read as simple, and the defective tag
    /// then does not fire. Two distinct eigenvalues closer than the threshold are merged, and the
    /// reading is then the trace over the sum of their eigenspaces, still basis-free. The span is
    /// taken from a thin QR of the unit eigenvectors and an SVD of its small R factor (the same
    /// singular values as the eigenvector matrix, with Q orthonormal to machine precision):
    /// directions with singular value below 1e-6 of the largest are dropped, which is how a
    /// defective cluster shows a dimension smaller than its multiplicity.</para></summary>
    private sealed class EigenspaceLedger
    {
        private const double ClusterTol = 1e-6;
        private const double SpanTol = 1e-6;

        private readonly int _channels;
        private readonly List<Piece> _pieces = new();
        private readonly List<int> _pieceOfMode = new();

        private sealed record Piece(Complex Lambda, int Rank, int Count, double[] Trace);

        public EigenspaceLedger(int channels) => _channels = channels;

        private static bool Close(Complex a, Complex b) =>
            (a - b).Magnitude <= ClusterTol * Math.Max(1.0, Math.Max(a.Magnitude, b.Magnitude));

        private static List<List<int>> Cluster(IReadOnlyList<Complex> lambdas)
        {
            int n = lambdas.Count;
            var parent = Enumerable.Range(0, n).ToArray();
            int Find(int x) { while (parent[x] != x) x = parent[x] = parent[parent[x]]; return x; }
            var order = Enumerable.Range(0, n).OrderBy(i => lambdas[i].Real).ToArray();
            for (int a = 0; a < n; a++)
            {
                var la = lambdas[order[a]];
                // Close() never accepts a pair whose real parts differ by more than this.
                double reach = 2.0 * ClusterTol * Math.Max(1.0, la.Magnitude);
                for (int b = a + 1; b < n && lambdas[order[b]].Real - la.Real <= reach; b++)
                    if (Close(la, lambdas[order[b]])) parent[Find(order[a])] = Find(order[b]);
            }
            return Enumerable.Range(0, n).GroupBy(Find).Select(g => g.ToList()).ToList();
        }

        public void AddBlock(MathNet.Numerics.LinearAlgebra.Vector<Complex> vals, Matrix<Complex> vecs, Func<int, int, bool> differs)
        {
            int bs = vals.Count;
            int first = _pieceOfMode.Count;
            for (int k = 0; k < bs; k++) _pieceOfMode.Add(-1);
            foreach (var cluster in Cluster(vals.ToArray()))
            {
                int k = cluster.Count;
                var V = Matrix<Complex>.Build.Dense(bs, k);
                for (int c = 0; c < k; c++)
                {
                    var col = vecs.Column(cluster[c]);
                    V.SetColumn(c, col.Divide(col.L2Norm()));
                }
                var qr = V.QR(QRMethod.Thin);
                var svd = qr.R.Svd(true);
                var sv = svd.S.Select(z => z.Real).ToArray();
                double sMax = sv.Max();
                var keep = Enumerable.Range(0, k).Where(r => sv[r] > SpanTol * sMax).ToList();
                var trace = new double[_channels];
                foreach (int r in keep)
                {
                    var u = qr.Q * svd.U.Column(r);
                    for (int i = 0; i < bs; i++)
                    {
                        var c = u[i];
                        double m = c.Real * c.Real + c.Imaginary * c.Imaginary;
                        if (m == 0.0) continue;
                        for (int l = 0; l < _channels; l++)
                            if (differs(i, l)) trace[l] += m;
                    }
                }
                Complex mean = cluster.Aggregate(Complex.Zero, (acc, i) => acc + vals[i]) / k;
                _pieces.Add(new Piece(mean, keep.Count, k, trace));
                foreach (int i in cluster) _pieceOfMode[first + i] = _pieces.Count - 1;
            }
        }

        public IReadOnlyList<CarrierMode> Finish(IReadOnlyList<CarrierMode> modes, IReadOnlyList<ChannelRate> channels)
        {
            var portfolioOfPiece = new ChannelDifferencePortfolio[_pieces.Count];
            var rankOfPiece = new int[_pieces.Count];
            var countOfPiece = new int[_pieces.Count];
            foreach (var group in Cluster(_pieces.Select(p => p.Lambda).ToList()))
            {
                int rank = group.Sum(p => _pieces[p].Rank);
                int count = group.Sum(p => _pieces[p].Count);
                var activity = Enumerable.Range(0, _channels)
                    .Select(l => new ChannelActivity(channels[l].Channel, group.Sum(p => _pieces[p].Trace[l]) / rank))
                    .ToList();
                var portfolio = new ChannelDifferencePortfolio(activity);
                foreach (int p in group) { portfolioOfPiece[p] = portfolio; rankOfPiece[p] = rank; countOfPiece[p] = count; }
            }
            return modes.Select((m, i) =>
            {
                int p = _pieceOfMode[i];
                return m with
                {
                    EigenspacePortfolio = portfolioOfPiece[p],
                    EigenspaceDimension = rankOfPiece[p],
                    EigenvalueMultiplicity = countOfPiece[p],
                };
            }).ToList();
        }
    }
}
