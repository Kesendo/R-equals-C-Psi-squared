using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
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
/// <c>simulations/_absorption_gamma_vector.py</c> and <c>simulations/_sip_carrier_channels.py</c>.</para>
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
    /// <summary>The absorption quantum (Absorption Theorem; the "2" in 2γ·⟨n_XY⟩, the d=2
    /// root of the Pi2 ladder). Kept as a named constant so the law reads literally.</summary>
    public const double AbsorptionQuantum = 2.0;

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
        return AbsorptionQuantum * weighted;
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

/// <summary>One mode: its measured decay rate −Re(λ) and its difference-portfolio.</summary>
public sealed record CarrierMode(double ActualDecayRate, ChannelDifferencePortfolio Portfolio);

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
            modes.Add(new CarrierMode(-vals[k].Real, new ChannelDifferencePortfolio(activity)));
        }

        return new CarrierPortfolioSpectrum(new CarrierVector(siteChannels), modes);
    }
}
