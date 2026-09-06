using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The centre-line fixed object in the Object Manager: not the system, an object INSIDE it.
/// The SYSTEM (a <see cref="MirrorSystem"/>) bears the x/y/z frame , z = the watched axis (Re λ),
/// x-y = the coherent-motion plane (Im λ), and σ = Σγ , while this object inherits that frame.
///
/// <para>The map must be named precisely. Linear F1 sends λ ↦ −2σ − λ and fixes only the single
/// complex point λ = −σ. Hermiticity preservation gives the Lindbladian spectrum conjugate closure
/// λ ↦ conj(λ). Their composite sends λ ↦ −2σ − conj(λ) and fixes the full line Re λ = −σ.
/// This live object counts that <b>conjugate-composite fixed line</b>, and separately reports the
/// linear-F1 fixed-point subset. It never calls every centre-line mode its own linear F1 partner.</para>
///
/// <para>Not a Claim, a live reading. Both counts are recomputed from the inherited spectrum.</para></summary>
public sealed class SelfMirrorObject : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double RelativeTolerance = 1e-7;
    private const double ClosedSystemRelativeTolerance = 1e-12;
    private const double MinimumResolvedDecayToFrequencyRatio = 1e-8;
    private const double MachineEpsilon = 2.2204460492503131e-16;
    private const double EigenvalueResolutionFactor = 256.0;

    private readonly MirrorSystem _system;   // INHERITED FROM: the x/y/z frame lives here, not on the object.

    public SelfMirrorObject(MirrorSystem system) =>
        _system = system ?? throw new ArgumentNullException(nameof(system));

    /// <summary>Inherited from the system, not owned.</summary>
    public int N => _system.N;

    /// <summary>Inherited from the system: σ = Σγ, the total dephasing (and the frame's center scale).</summary>
    public double Sigma => _system.TotalDephasing;

    /// <summary>The object's place in the inherited frame: Re λ = −σ, the palindrome center line.</summary>
    public double Center => -_system.TotalDephasing;

    private double MaxDecayScale => _system.Spectrum.Modes
        .Select(m => Math.Abs(m.ActualDecayRate)).Append(Math.Abs(Sigma)).Max();

    private double MaxFrequencyScale => _system.Spectrum.Modes
        .Select(m => Math.Abs(m.OscillationFrequency)).DefaultIfEmpty(0.0).Max();

    /// <summary>Whether a nonzero dissipative centre is resolved against the coherent scale.
    /// Sigma=0 is a separately known closed-system case; otherwise an unresolved scale separation
    /// is surfaced instead of returning a tolerance-dependent integer.</summary>
    public bool IsFixedSetResolved => Sigma == 0.0 || MaxFrequencyScale == 0.0 ||
        MaxDecayScale / MaxFrequencyScale >= MinimumResolvedDecayToFrequencyRatio;

    /// <summary>A decay-axis tolerance. Oscillation frequencies do not set this scale: letting
    /// a large Hamiltonian widen a real-part test would classify off-centre modes as fixed.
    /// A common rescaling of H and all channel rates still rescales the decay rates and this
    /// tolerance together.</summary>
    private double DecayAxisTolerance
    {
        get
        {
            if (Sigma == 0.0)
                return ClosedSystemRelativeTolerance * Math.Max(MaxFrequencyScale, 1.0);
            return RelativeTolerance * MaxDecayScale;
        }
    }

    /// <summary>The Im λ = 0 gate uses an eigensolver-resolution estimate, not the much wider
    /// decay-axis membership window. In a dissipative problem it is tied to the decay scale, not
    /// the largest unrelated Hamiltonian branch; otherwise one large frequency could erase another
    /// exact, small frequency. A common change of energy units still rescales the estimate. In the
    /// separate closed-system case, where there is no decay reference, the frequency scale is used.</summary>
    private double FrequencyZeroTolerance =>
        EigenvalueResolutionFactor * MachineEpsilon *
        (Sigma == 0.0 ? MaxFrequencyScale : MaxDecayScale);

    private static string F(double value) =>
        value != 0.0 && Math.Abs(value) < 1e-4
            ? value.ToString("0.####E+0", Inv)
            : value.ToString("0.####", Inv);

    /// <summary>Multiplicity on the fixed line of the conjugate-composite map
    /// λ ↦ −2σ − conj(λ), equivalently Re λ = −σ. Read live from the inherited spectrum.</summary>
    public int CompositeFixedLineCount
    {
        get
        {
            EnsureResolved();
            double tolerance = DecayAxisTolerance;
            return _system.Spectrum.Modes.Count(m =>
                Math.Abs(m.ActualDecayRate - Sigma) <= tolerance);
        }
    }

    /// <summary>Multiplicity at the fixed point of the linear F1 map λ ↦ −2σ − λ.
    /// Linear F1 fixes λ only when λ = −σ, so both the centre-rate and zero-frequency conditions
    /// are required. This is a subset of <see cref="CompositeFixedLineCount"/>.</summary>
    public int LinearF1FixedPointCount
    {
        get
        {
            EnsureResolved();
            double decayTolerance = DecayAxisTolerance;
            double frequencyTolerance = FrequencyZeroTolerance;
            return _system.Spectrum.Modes.Count(m =>
                Math.Abs(m.ActualDecayRate - Sigma) <= decayTolerance &&
                Math.Abs(m.OscillationFrequency) <= frequencyTolerance);
        }
    }

    private void EnsureResolved()
    {
        if (!IsFixedSetResolved)
            throw new InvalidOperationException("fixed-set count unresolved: dissipative scale is below the eigensolver resolution relative to the coherent scale");
    }

    public string DisplayName =>
        $"SelfMirrorObject (conjugate-composite fixed line Re λ = −σ = {F(Center)}; ⊂ MirrorSystem N={N})";

    public string Summary =>
        IsFixedSetResolved
        ? $"an object INSIDE the system: {CompositeFixedLineCount} modes fixed by the composite " +
        $"λ ↦ −2σ − conj(λ); its {LinearF1FixedPointCount}-mode subset at λ = −σ is fixed by linear F1 " +
        $"λ ↦ −2σ − λ. Everything else (x/y/z and σ = {F(Sigma)}) is INHERITED."
        : $"fixed-set count UNRESOLVED: σ = {F(Sigma)} is below numerical resolution relative to the coherent scale; no integer count is reported.";

    public IEnumerable<IInspectable> Children
    {
        get
        {
            // The jump, made visible: x/y/z is the system's inheritance. Drill into the parent itself.
            yield return new InspectableNode(
                displayName: "↑ inherited frame (x/y/z) , the SYSTEM, not the object",
                summary: "z = the watched axis (Re λ), x-y = the coherent-motion plane (Im λ), σ = the center: " +
                         "the object does NOT own these. They are the system's structure; the object only inherits a place in them.",
                children: new IInspectable[] { _system });

            // The object's own delta: the composite fixed line, with the smaller linear-F1 subset fenced.
            yield return IsFixedSetResolved
                ? new InspectableNode(
                    displayName: "the object itself: conjugate-composite fixed line",
                    summary: $"{CompositeFixedLineCount} modes at Re λ = −σ are fixed by λ ↦ −2σ − conj(λ). " +
                             $"Only {LinearF1FixedPointCount} of them also have Im λ = 0 and are fixed by the linear F1 map " +
                             $"λ ↦ −2σ − λ, whose fixed-point equation is λ = −σ.")
                : new InspectableNode(
                    displayName: "the object itself: fixed-set count unresolved",
                    summary: "the dissipative axis is below numerical resolution relative to H; rerun at a resolvable ratio or use an exact/block-aware calculation");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
