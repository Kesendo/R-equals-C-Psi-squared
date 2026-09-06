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
/// <para>Not a Claim, a live reading. Counts are emitted only on the explicitly certified
/// zero-Hamiltonian, uniform-rate branch; floating-spectrum membership remains unresolved.</para></summary>
public sealed class SelfMirrorObject : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private readonly MirrorSystem _system;   // INHERITED FROM: the x/y/z frame lives here, not on the object.

    public SelfMirrorObject(MirrorSystem system) =>
        _system = system ?? throw new ArgumentNullException(nameof(system));

    /// <summary>Inherited from the system, not owned.</summary>
    public int N => _system.N;

    /// <summary>Inherited from the system: σ = Σγ, the total dephasing (and the frame's center scale).</summary>
    public double Sigma => _system.TotalDephasing;

    /// <summary>The object's place in the inherited frame: Re λ = −σ, the palindrome center line.</summary>
    public double Center => -_system.TotalDephasing;

    /// <summary>Exact fixed-point multiplicity cannot be certified from floating eigenvalues alone:
    /// a genuinely small spectral gap is indistinguishable from eigensolver error. This object
    /// therefore reports counts only for the algebraically closed zero-Hamiltonian, uniform-rate
    /// case. All other inputs remain explicitly unresolved until a certified rank/gap route exists.</summary>
    public bool IsFixedSetResolved => IsCertifiedZeroHamiltonianUniformRateCase;

    private bool IsCertifiedZeroHamiltonianUniformRateCase
    {
        get
        {
            if (_system.Hamiltonian.Enumerate().Any(z => z.Real != 0.0 || z.Imaginary != 0.0))
                return false;
            double gamma = _system.Channels[0].Gamma;
            return gamma >= 0.0 && _system.Channels.All(channel => channel.Gamma == gamma);
        }
    }

    private static int Binomial(int n, int k)
    {
        k = Math.Min(k, n - k);
        long value = 1;
        for (int i = 1; i <= k; i++)
            value = checked(value * (n - k + i) / i);
        return checked((int)value);
    }

    private int CertifiedFixedCount
    {
        get
        {
            if (Sigma == 0.0)
                return Enumerable.Repeat(4, N).Aggregate(1, (value, factor) => checked(value * factor));
            if ((N & 1) != 0)
                return 0;
            return checked(Enumerable.Repeat(2, N).Aggregate(1, (value, factor) => checked(value * factor)) *
                           Binomial(N, N / 2));
        }
    }

    private static string F(double value) =>
        value != 0.0 && Math.Abs(value) < 1e-4
            ? value.ToString("0.####E+0", Inv)
            : value.ToString("0.####", Inv);

    /// <summary>Multiplicity on the fixed line of the conjugate-composite map
    /// λ ↦ −2σ − conj(λ), equivalently Re λ = −σ. Reported only from the certified
    /// algebraic branch described by <see cref="IsFixedSetResolved"/>.</summary>
    public int CompositeFixedLineCount
    {
        get
        {
            EnsureResolved();
            return CertifiedFixedCount;
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
            return CertifiedFixedCount;
        }
    }

    private void EnsureResolved()
    {
        if (!IsFixedSetResolved)
            throw new InvalidOperationException("fixed-set count unresolved: floating eigenvalues do not certify exact line or point membership");
    }

    public string DisplayName =>
        $"SelfMirrorObject (conjugate-composite fixed line Re λ = −σ = {F(Center)}; ⊂ MirrorSystem N={N})";

    public string Summary =>
        IsFixedSetResolved
        ? $"an object INSIDE the system: {CompositeFixedLineCount} modes fixed by the composite " +
        $"λ ↦ −2σ − conj(λ); its {LinearF1FixedPointCount}-mode subset at λ = −σ is fixed by linear F1 " +
        $"λ ↦ −2σ − λ. Everything else (x/y/z and σ = {F(Sigma)}) is INHERITED."
        : $"fixed-set count UNRESOLVED: floating eigenvalues do not certify exact membership in the composite line " +
          $"λ ↦ −2σ − conj(λ) or the linear F1 point λ = −σ; no integer count is reported (σ = {F(Sigma)}).";

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
                    summary: "floating eigenvalues alone cannot certify exact line or point membership; use an exact/block-aware calculation with a certified gap");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
