using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

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
///   <item><see cref="MemoryRotation"/> , F80: the 90° turn H ↔ M. The mirror-defect
///         M = Π·L·Π⁻¹ + L + 2σ·I is 0 for truly H (perfect mirror) and carries H's spectrum
///         rotated a quarter turn (Spec(M) = ±2i·Spec(H)) for non-truly chain Π²-odd H , energy
///         on the real axis, memory on the imaginary one. The wave remembers across the turn.</item>
/// </list></para>
///
/// <para><b>Voices still to add</b> (each a future property on this object): the bit_a/bit_b sectors
/// (F61/F63, cavity/transport), the inside-observable Q = J/γ, and the measurement wrapper
/// (spectrum → fitted observable). Each grows the stand without changing what is already on it.</para>
///
/// <para>Diagnostic-scale (dense 4^N via <see cref="CarrierVectorPortfolio.Decompose"/>): the
/// small N where the whole box fits in one view, not the large-N engine.</para></summary>
public sealed class MirrorSystem : IInspectable
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
    /// decay rate. Computed once, lazily, from this system's one input. For a number-conserving H
    /// (XY, Heisenberg) the Z-dephasing Liouvillian is block-diagonal in the joint-popcount sectors,
    /// so this uses <see cref="CarrierVectorPortfolio.DecomposeBlocked"/> ((N+1)² small Evds instead
    /// of one 4^N one); a non-conserving H (e.g. the xybond) falls back to the full decomposition.
    /// Either way the spectrum and portfolios are the same.</summary>
    public CarrierPortfolioSpectrum Spectrum =>
        _spectrum ??= CarrierVectorPortfolio.IsNumberConserving(Hamiltonian)
            ? CarrierVectorPortfolio.DecomposeBlocked(N, Hamiltonian, Channels)
            : CarrierVectorPortfolio.Decompose(N, Hamiltonian, Channels);

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

    private MemoryRotationReading? _memoryRotation;

    /// <summary>The 90° memory rotation (F80): build this system's mirror-defect
    /// M = Π·L·Π⁻¹ + L + 2σ·I and read how it relates to the Hamiltonian.
    ///
    /// <para>For a truly Hamiltonian (Heisenberg, the XY model) under Z-dephasing the F1 palindrome
    /// holds exactly and M = 0: a perfect mirror, no defect (<see cref="MemoryRotationReading.PerfectMirror"/>).
    /// For a non-truly chain Π²-odd Hamiltonian the palindrome breaks, and the defect is not noise:
    /// Spec(M) = ±2i·Spec(H). H's real energies (the measured, outer axis) reappear on M's imaginary
    /// axis (the time / decay / memory axis), turned a quarter and scaled by 2. The wave remembers by
    /// sharing its spectrum across the 90° turn; H is the distance between the two mirror sectors.</para>
    ///
    /// <para>The quarter turn is the third step of a short lineage: Π² = I (the mirror exists) → F1
    /// (Π·L·Π⁻¹ = −L − 2σ·I, the palindrome) → F80 (when the palindrome breaks, the break carries H,
    /// rotated 90°). <see cref="MemoryRotationReading.MemoryCarriesEnergy"/> is the live F80 check via
    /// the Frobenius signature ‖M‖²_F = 4·‖H‖²_F·2^N. Computed once, lazily, from this system's input.</para></summary>
    public MemoryRotationReading MemoryRotation => _memoryRotation ??= ComputeMemoryRotation();

    private MemoryRotationReading ComputeMemoryRotation()
    {
        var gammas = Channels.Select(c => c.Gamma).ToList();
        var L = PauliDephasingDissipator.BuildZ(Hamiltonian, gammas);
        var M = PalindromeResidual.Build(L, N, TotalDephasing);
        double defect = M.FrobeniusNorm();

        // H's energies (real, since H is Hermitian) and their 90° images 2λ on the memory axis.
        var energies = Hamiltonian.Evd().EigenValues
            .Select(z => Math.Round(z.Real, 10)).Distinct().OrderBy(x => x).ToList();
        var rotation = energies.Select(e => new EnergyMemoryImage(e, 2.0 * e)).ToList();

        // F80 Frobenius signature: M carries all of H's spectrum (rotated 90°) iff ‖M‖²_F = 4‖H‖²_F·2^N.
        double hNormSq = Hamiltonian.FrobeniusNorm();
        hNormSq *= hNormSq;
        double predicted = 4.0 * hNormSq * (1 << N);
        bool carries = hNormSq > 1e-12 && Math.Abs(defect * defect - predicted) < 1e-6 * Math.Max(1.0, predicted);
        bool perfectMirror = defect < 1e-9;

        return new MemoryRotationReading(defect, perfectMirror, carries, rotation);
    }

    private TaktReading? _takt;

    /// <summary>The Takt: the beat <see cref="Evolve"/> silently counts in. K = γ₀·t is a count of
    /// carrier-ticks, and the carrier γ₀ that sets the tick cancels out of the view; this voice names
    /// the unit so the wegkürzen does not leave the clock invisible.
    ///
    /// <para>γ is the timekeeper: it is what makes time tick at all. The floor of felt motion
    /// (<see cref="TaktReading.Gap"/>, the 2γ floor at small N) is the spectrum's own
    /// <see cref="CarrierPortfolioSpectrum.SlowestRate"/> , the framework's shared memory floor, so the
    /// Takt and the spectrum read the one number instead of two; its reciprocal is the longest breath
    /// before the slow mode reaches the centre (<see cref="TaktReading.Tau"/>). Read live off that
    /// floor, so it shows what the modes actually do, not an assumed 2γ. At γ=0 the floor is 0 and the
    /// clock stands still (<see cref="TaktReading.Stopped"/>): no felt time, only frozen Hamiltonian
    /// oscillation.</para>
    ///
    /// <para>The clock is a circle: a single mode winds as e^(λt) = e^(−αt)·e^(iωt), a radius e^(−αt)
    /// (decay, set by γ) turning an angle e^(iωt) (rotation, set by J) into a logarithmic spiral. This
    /// is the <i>radial</i> hand, the inward winding. The <i>angular</i> hand, the rotation ω=2J and the
    /// angle θ=arctan(Q), is a later voice; the circle closes then. Computed once, lazily, from this
    /// system's input.</para></summary>
    public TaktReading Takt => _takt ??= ComputeTakt();

    private TaktReading ComputeTakt()
    {
        double gap = Spectrum.SlowestRate;
        return gap <= 0.0
            ? new TaktReading(true, 0.0, double.PositiveInfinity)
            : new TaktReading(false, gap, 1.0 / gap);
    }

    private RotationReading? _rotation;

    /// <summary>The Rotation: the angular hand of the circular quantum clock. A single mode winds as
    /// e^(λt) = e^(−αt)·e^(iωt); the <see cref="Takt"/> reads the radius e^(−αt) (the inward decay set
    /// by γ), this voice reads the angle e^(iωt) (the turning set by J), on the same memory mode the
    /// Takt's Tau tracks. Together the two hands trace the memory mode's logarithmic spiral.
    ///
    /// <para>The reading is the memory mode's rotation ω = Im(λ) (<see cref="RotationReading.Frequency"/>,
    /// the slowest mortal mode's oscillation) and its angle θ = arctan(ω/Gap)
    /// (<see cref="RotationReading.Angle"/>), the F95 angle (= arctan(Q) for the 2-level case), the seam
    /// between the radial bank (past, {I,Z}) and the angular bank (future, {X,Y}). θ is stored in
    /// radians, canonical and composing with F95; the CLI renders it in degrees.</para>
    ///
    /// <para>Pairs with the Takt: when γ=0 the radial hand stops (<see cref="TaktReading.Stopped"/>) and
    /// the angle reaches its maximum θ = π/2, the pure circle turning forever with no inward pull. When
    /// J=0 nothing rotates (all ω=0), so the hand does not turn (<see cref="RotationReading.Turning"/> is
    /// false) and θ = 0, pure radial decay. Computed once, lazily, from this system's input.</para></summary>
    public RotationReading Rotation => _rotation ??= ComputeRotation();

    private RotationReading ComputeRotation()
    {
        const double tol = 1e-9;
        double gap = Spectrum.SlowestRate;
        if (gap > tol)
        {
            // The mortal memory mode's rotation: max |ω| over the modes sitting at the slowest rate
            // (deterministic across a degenerate slowest rate / ±ω conjugate pair).
            double frequency = Spectrum.Modes
                .Where(m => Math.Abs(m.ActualDecayRate - gap) <= tol)
                .Select(m => Math.Abs(m.OscillationFrequency))
                .DefaultIfEmpty(0.0).Max();
            return new RotationReading(frequency > tol, frequency, Math.Atan2(frequency, gap));
        }

        // γ=0, the pure-circle limit (Takt.Stopped): no decay floor, read the turning over all modes.
        double freq = Spectrum.Modes
            .Select(m => Math.Abs(m.OscillationFrequency)).DefaultIfEmpty(0.0).Max();
        bool turning = freq > tol;
        return new RotationReading(turning, freq, turning ? Math.PI / 2.0 : 0.0);
    }

    /// <summary>Move the carrier and get a fresh reading of the whole system (the box moves):
    /// a new <see cref="MirrorSystem"/> with the given per-site dephasing and the same H. Every
    /// property recomputes from the moved input.</summary>
    public MirrorSystem WithChannels(IReadOnlyList<ChannelRate> channels) =>
        new(N, Hamiltonian, channels);

    // ---- Object Manager: the conductor's stand as a live IInspectable node ----
    // MirrorSystem is the first live-data object of the growing object-manager (see the class
    // summary), not a registry Claim: built on demand from (N, H, γ), readings computed lazily. The
    // voices that fall out of the one spectrum decomposition , the slow modes, the F1 palindrome, the
    // clock's two hands , surface as a tree. Numbers format InvariantCulture so the reading is the
    // same whoever renders it. The heavy F80 MemoryRotation (it builds the full 4^N defect) is not in
    // the default tree; read it through the MemoryRotation property.

    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    /// <summary>Short label: N, the palindrome centre σ, and the Liouvillian dimension 4^N.</summary>
    public string DisplayName =>
        $"MirrorSystem (N={N}, σ={TotalDephasing.ToString("0.####", Inv)}, dim {1L << (2 * N)})";

    /// <summary>One-line fingerprint: the slowest surviving rate, whether F1 holds, the clock angle.</summary>
    public string Summary =>
        $"slowest rate {Spectrum.SlowestRate.ToString("0.####", Inv)}, " +
        $"F1 palindrome holds = {PalindromeHolds}, " +
        $"clock θ = {(Rotation.Angle * 180.0 / Math.PI).ToString("0.#", Inv)}°";

    /// <summary>The voices that fall out of the one spectrum decomposition: the slow modes, the F1
    /// palindrome check, and the clock's two hands. A pure container, no leaf payload of its own.</summary>
    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return SpectrumNode();
            yield return PalindromeNode();
            yield return ClockNode();
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    /// <summary>The slowest <paramref name="top"/> distinct-rate modes (the most-rotating
    /// representative at each rate, matching the Rotation voice), each carrying its per-channel
    /// Δ-portfolio as a Vector payload , the bars the renderer draws under --draw.</summary>
    private InspectableNode SpectrumNode(int top = 8)
    {
        var modes = Spectrum.Modes
            .Where(m => m.ActualDecayRate > 1e-9)
            .GroupBy(m => m.ActualDecayRate.ToString("0.000000", Inv))
            .OrderBy(g => g.First().ActualDecayRate)
            .Take(top)
            .Select(g => g.OrderByDescending(m => Math.Abs(m.OscillationFrequency)).First())
            .Select(ModeNode)
            .ToList();
        return InspectableNode.Group("Spectrum", modes, modes.Count);
    }

    private static IInspectable ModeNode(CarrierMode m)
    {
        double theta = Math.Atan2(Math.Abs(m.OscillationFrequency), m.ActualDecayRate) * 180.0 / Math.PI;
        var labels = m.Portfolio.Activity.Select(a => a.Channel).ToList();
        var deltas = ComplexVector.Build.DenseOfEnumerable(
            m.Portfolio.Activity.Select(a => new Complex(a.Delta, 0.0)));
        string portfolio = string.Join("  ",
            m.Portfolio.Activity.Select(a => $"{a.Channel} {a.Delta.ToString("0%", Inv)}"));
        return new InspectableNode(
            displayName: $"rate {m.ActualDecayRate.ToString("0.0000", Inv)}  θ {theta.ToString("0.0", Inv)}°",
            summary: portfolio,
            payload: new InspectablePayload.Vector("Δ-portfolio", deltas, labels));
    }

    private InspectableNode PalindromeNode()
    {
        var sample = PalindromePartners.FirstOrDefault(p => p.Rate > 1e-9 && p.PartnerPresent);
        IInspectable[] children = sample is null
            ? Array.Empty<IInspectable>()
            : new IInspectable[]
            {
                new InspectableNode(
                    $"e.g. rate {sample.Rate.ToString("0.0000", Inv)} pairs with {sample.PartnerRate.ToString("0.0000", Inv)}",
                    summary: $"partner present: {sample.PartnerPresent}"),
            };
        return new InspectableNode("F1 palindrome",
            summary: $"holds = {PalindromeHolds}; rate r pairs with 2σ − r (2σ = {(2.0 * TotalDephasing).ToString("0.####", Inv)})",
            children: children);
    }

    private InspectableNode ClockNode()
    {
        var takt = Takt;
        var rot = Rotation;
        var children = new IInspectable[]
        {
            InspectableNode.RealScalar("Takt gap (2γ floor)", takt.Gap, "0.####"),
            InspectableNode.RealScalar("Takt τ (longest breath)", takt.Tau, "0.####"),
            InspectableNode.RealScalar("Rotation ω", rot.Frequency, "0.####"),
            InspectableNode.RealScalar("Rotation angle θ (deg)", rot.Angle * 180.0 / Math.PI, "0.#"),
        };
        return InspectableNode.Group("Clock", children);
    }
}

/// <summary>One mode's F1 mirror pairing: its decay rate, the partner rate 2σ − r the palindrome
/// requires, and whether that partner is present in the spectrum.</summary>
public sealed record PalindromePairing(double Rate, double PartnerRate, bool PartnerPresent);

/// <summary>One mode at K elapsed carrier-ticks: its decay rate, its survival factor
/// e^(−(rate/σ)·K), and its channel-difference portfolio. The portfolio says <i>where</i> the mode
/// lives (which channels carry its bra-ket difference); the survival says <i>how much</i> of it is
/// still there at time K. Together they are the spectrum seen as it decays, not as a static list.</summary>
public sealed record ModeSurvival(double Rate, double Survival, ChannelDifferencePortfolio Portfolio);

/// <summary>The F80 reading of this system's mirror-defect M = Π·L·Π⁻¹ + L + 2σ·I.
/// <see cref="DefectNorm"/> is ‖M‖_F, the size of the break from the F1 palindrome (0 = perfect
/// mirror). <see cref="PerfectMirror"/> is true for truly Hamiltonians (M = 0).
/// <see cref="MemoryCarriesEnergy"/> is the live F80 identity check (Spec(M) = ±2i·Spec(H), via the
/// Frobenius signature): true when the defect carries the full Hamiltonian spectrum, rotated 90°.
/// <see cref="Rotation"/> lists each energy and its quarter-turn image on the memory axis.</summary>
public sealed record MemoryRotationReading(
    double DefectNorm,
    bool PerfectMirror,
    bool MemoryCarriesEnergy,
    IReadOnlyList<EnergyMemoryImage> Rotation);

/// <summary>One rung shared between energy and memory across the 90° turn: a real energy eigenvalue
/// λ of H (the measured, outer axis) and its image 2λ on the imaginary memory axis (the
/// imaginary part of M's eigenvalue 2iλ). When the F80 identity holds, these are the same spectrum
/// seen from the two sides of the quarter turn.</summary>
public sealed record EnergyMemoryImage(double Energy, double MemoryAxisValue);

/// <summary>The Takt reading: the unit of felt time <see cref="MirrorSystem.Evolve"/> counts in.
/// <see cref="Stopped"/> is true at γ=0, when there is no decay clock, only frozen Hamiltonian
/// oscillation. <see cref="Gap"/> is the spectrum's <see cref="CarrierPortfolioSpectrum.SlowestRate"/>
/// (the framework's shared memory floor, the 2γ floor of felt motion at small N), <c>0.0</c> when
/// stopped. <see cref="Tau"/> = 1/Gap is the longest felt duration, the slowest breath before the
/// slow mode reaches the centre; <see cref="double.PositiveInfinity"/> when stopped. This is the
/// radial hand of the circular quantum clock e^(λt) = e^(−αt)·e^(iωt).</summary>
public sealed record TaktReading(bool Stopped, double Gap, double Tau);

/// <summary>The Rotation reading: the angular hand of the circular quantum clock
/// e^(λt) = e^(−αt)·e^(iωt), the companion to <see cref="TaktReading"/>'s radial hand.
/// <see cref="Turning"/> is true when the memory mode rotates (ω ≠ 0, i.e. J ≠ 0). <see cref="Frequency"/>
/// is ω = |Im(λ)| of the memory mode, its turning rate. <see cref="Angle"/> is θ = arctan(ω/Gap) in
/// radians, the F95 angle (= arctan(Q) for the 2-level case): θ → π/2 at the pure circle (γ=0, the
/// radial hand stopped, turning forever) and θ = 0 for pure radial decay (J=0, no turning).</summary>
public sealed record RotationReading(bool Turning, double Frequency, double Angle);
