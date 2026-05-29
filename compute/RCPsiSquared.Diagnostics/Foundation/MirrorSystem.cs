using System;
using System.Collections.Generic;
using System.Linq;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The conductor's stand: one open quantum system held live, with the F-formula
/// readings glued on as computed properties. Give it (N, H, per-site γ); each reading is a
/// property derived from that one input, so moving the input moves every reading together,
/// the whole box at once, instead of one stop-and-restart script per mechanism.
///
/// <para>This is the first object of a growing object-manager: today one system with all its
/// mirror-readings side by side; later, many systems held live at once. We are the conductor;
/// this is the stand the score sits on.</para>
///
/// <para><b>Voices wired so far</b> (each a property on this same object, all from the one input):
/// <list type="bullet">
///   <item><see cref="Spectrum"/> , the inner law (<see cref="CarrierVectorPortfolio.Decompose"/>):
///         every mode's per-channel difference-portfolio and decay rate.</item>
///   <item><see cref="PalindromePartners"/> , F1: every decay rate r pairs with 2σ − r
///         (Π·L·Π⁻¹ = −L − 2σ·I), read live off the spectrum.</item>
///   <item><see cref="Evolve"/> , time: unroll the static spectrum over K carrier-ticks into the
///         decay each mode shows. This is the connection from what is (the spectrum) to what is
///         measured (the FID); see its own note on naming time and the connection.</item>
/// </list></para>
///
/// <para><b>Voices still to add</b> (each a future property on this object): the 90° memory
/// rotation H ↔ M (F80, energy ↔ time), the bit_a/bit_b sectors (F61/F63, cavity/transport),
/// the inside-observable Q = J/γ, and the measurement wrapper (spectrum → fitted observable).
/// Each grows the stand without changing what is already on it.</para>
///
/// <para>Diagnostic-scale (dense 4^N via <see cref="CarrierVectorPortfolio.Decompose"/>): the
/// small N where the whole box fits in one view, not the large-N engine.</para></summary>
public sealed class MirrorSystem
{
    /// <summary>Number of sites.</summary>
    public int N { get; }

    /// <summary>The Hamiltonian (2^N × 2^N). Hermitian in the physical case, but the readings
    /// hold for any H the Absorption Theorem covers.</summary>
    public ComplexMatrix Hamiltonian { get; }

    /// <summary>The per-site dephasing carrier vector, one <see cref="ChannelRate"/> per site.</summary>
    public IReadOnlyList<ChannelRate> Channels { get; }

    public MirrorSystem(int N, ComplexMatrix hamiltonian, IReadOnlyList<ChannelRate> channels)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be >= 1");
        Hamiltonian = hamiltonian ?? throw new ArgumentNullException(nameof(hamiltonian));
        Channels = channels ?? throw new ArgumentNullException(nameof(channels));
        if (channels.Count != N)
            throw new ArgumentException($"need one channel per site (N={N}), got {channels.Count}", nameof(channels));
        int d = 1 << N;
        if (hamiltonian.RowCount != d || hamiltonian.ColumnCount != d)
            throw new ArgumentException($"Hamiltonian must be {d}×{d} for N={N}", nameof(hamiltonian));
        this.N = N;
    }

    /// <summary>σ = Σ_l γ_l, the total dephasing and the palindrome centre.</summary>
    public double TotalDephasing => Channels.Sum(c => c.Gamma);

    private CarrierPortfolioSpectrum? _spectrum;

    /// <summary>The inner law: every Liouvillian mode as a per-channel difference-portfolio and
    /// decay rate. Computed once, lazily, from this system's one input.</summary>
    public CarrierPortfolioSpectrum Spectrum =>
        _spectrum ??= CarrierVectorPortfolio.Decompose(N, Hamiltonian, Channels);

    /// <summary>F1 palindrome read as a live property: each mode's decay rate r is paired with
    /// its mirror partner 2σ − r. <see cref="PalindromePairing.PartnerPresent"/> is true when a
    /// mode at that rate exists in the spectrum, which it always should by Π·L·Π⁻¹ = −L − 2σ·I.</summary>
    public IReadOnlyList<PalindromePairing> PalindromePartners
    {
        get
        {
            double twoSigma = 2.0 * TotalDephasing;
            var rates = Spectrum.Modes.Select(m => m.ActualDecayRate).ToList();
            return rates.Select(r =>
            {
                double partner = twoSigma - r;
                bool present = rates.Any(s => Math.Abs(s - partner) < 1e-7);
                return new PalindromePairing(r, partner, present);
            }).ToList();
        }
    }

    /// <summary>True iff every mode has its F1 mirror partner present in the spectrum, i.e. the
    /// palindrome Π·L·Π⁻¹ = −L − 2σ·I holds for this system. A live check, not a stored fact.</summary>
    public bool PalindromeHolds => PalindromePartners.All(p => p.PartnerPresent);

    /// <summary>Unroll the static spectrum into time: the connection by which what *is* (the
    /// spectrum, the inner observation) becomes what is *measured* (the decay over time, the
    /// outer observation). Returns each mode's survival e^(−(rate/σ)·K) after K elapsed
    /// carrier-ticks. At K=0 every mode is fully present; as K grows the fast modes vanish and
    /// only the slow (memory) modes remain.
    ///
    /// <para><b>On time and the connection, named on purpose.</b> K = γ₀·t is the dimensionless
    /// time, the count of carrier-ticks (here in units of σ = Σγ, the carrier scale and the
    /// palindrome centre). The absolute time t and the carrier γ₀ cancel into K. That
    /// cancellation is clean, and it robs the view: it hides that a real <i>time</i> flows here,
    /// and that this method is a real <i>connection</i>, the channel where the inner spectrum
    /// turns into the outer measurement. We name them so the wegkürzen does not blind us:
    /// K is the time, in ticks; Evolve is the connection, spectrum → measurement. From inside,
    /// only K is readable (γ₀ is the silent unit), so raw local t would carry no statement; K does.</para>
    ///
    /// <para>A measured FID is a wrapper-weighted sum of these survivals (the weights come from the
    /// preparation and the observable, the outer side); Evolve gives the inner envelope each weight
    /// rides on. This is the spectrum's own decay, before any prep or observable is chosen.</para></summary>
    public IReadOnlyList<ModeSurvival> Evolve(double K)
    {
        if (K < 0) throw new ArgumentOutOfRangeException(nameof(K), K, "K (carrier-ticks elapsed) must be >= 0");
        double scale = TotalDephasing > 0 ? TotalDephasing : 1.0;
        return Spectrum.Modes
            .Select(m => new ModeSurvival(m.ActualDecayRate, Math.Exp(-(m.ActualDecayRate / scale) * K), m.Portfolio))
            .ToList();
    }

    /// <summary>Move the carrier and get a fresh reading of the whole system (the box moves):
    /// a new <see cref="MirrorSystem"/> with the given per-site dephasing and the same H. Every
    /// property recomputes from the moved input.</summary>
    public MirrorSystem WithChannels(IReadOnlyList<ChannelRate> channels) =>
        new(N, Hamiltonian, channels);
}

/// <summary>One mode's F1 mirror pairing: its decay rate, the partner rate 2σ − r the palindrome
/// requires, and whether that partner is present in the spectrum.</summary>
public sealed record PalindromePairing(double Rate, double PartnerRate, bool PartnerPresent);

/// <summary>One mode at K elapsed carrier-ticks: its decay rate, its survival factor
/// e^(−(rate/σ)·K), and its channel-difference portfolio. The portfolio says <i>where</i> the mode
/// lives (which channels carry its bra-ket difference); the survival says <i>how much</i> of it is
/// still there at time K. Together they are the spectrum seen as it decays, not as a static list.</summary>
public sealed record ModeSurvival(double Rate, double Survival, ChannelDifferencePortfolio Portfolio);
