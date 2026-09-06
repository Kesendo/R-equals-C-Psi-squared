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

    private readonly MirrorSystem _system;   // INHERITED FROM: the x/y/z frame lives here, not on the object.

    public SelfMirrorObject(MirrorSystem system) =>
        _system = system ?? throw new ArgumentNullException(nameof(system));

    /// <summary>Inherited from the system, not owned.</summary>
    public int N => _system.N;

    /// <summary>Inherited from the system: σ = Σγ, the total dephasing (and the frame's center scale).</summary>
    public double Sigma => _system.TotalDephasing;

    /// <summary>The object's place in the inherited frame: Re λ = −σ, the palindrome center line.</summary>
    public double Center => -_system.TotalDephasing;

    /// <summary>A spectral-scale tolerance. A common rescaling of H and all channel rates
    /// rescales both this tolerance and every fixed-set residual by the same factor.</summary>
    private double FixedSetTolerance
    {
        get
        {
            double scale = Math.Abs(Sigma);
            foreach (var mode in _system.Spectrum.Modes)
                scale = Math.Max(scale,
                    Math.Max(Math.Abs(mode.ActualDecayRate), Math.Abs(mode.OscillationFrequency)));
            return RelativeTolerance * scale;
        }
    }

    /// <summary>Multiplicity on the fixed line of the conjugate-composite map
    /// λ ↦ −2σ − conj(λ), equivalently Re λ = −σ. Read live from the inherited spectrum.</summary>
    public int CompositeFixedLineCount
    {
        get
        {
            double tolerance = FixedSetTolerance;
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
            double tolerance = FixedSetTolerance;
            return _system.Spectrum.Modes.Count(m =>
                Math.Abs(m.ActualDecayRate - Sigma) <= tolerance &&
                Math.Abs(m.OscillationFrequency) <= tolerance);
        }
    }

    public string DisplayName =>
        $"SelfMirrorObject (conjugate-composite fixed line Re λ = −σ = {Center.ToString("0.####", Inv)}; ⊂ MirrorSystem N={N})";

    public string Summary =>
        $"an object INSIDE the system: {CompositeFixedLineCount} modes fixed by the composite " +
        $"λ ↦ −2σ − conj(λ); its {LinearF1FixedPointCount}-mode subset at λ = −σ is fixed by linear F1 " +
        $"λ ↦ −2σ − λ. Everything else (x/y/z and σ = {Sigma.ToString("0.####", Inv)}) is INHERITED.";

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
            yield return new InspectableNode(
                displayName: "the object itself: conjugate-composite fixed line",
                summary: $"{CompositeFixedLineCount} modes at Re λ = −σ are fixed by λ ↦ −2σ − conj(λ). " +
                         $"Only {LinearF1FixedPointCount} of them also have Im λ = 0 and are fixed by the linear F1 map " +
                         $"λ ↦ −2σ − λ, whose fixed-point equation is λ = −σ.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
