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
/// <para>Not a Claim, a live reading, and an exact one: both counts are algebraic
/// multiplicities of the characteristic polynomial of L for the inputs as given (every double is a
/// dyadic rational), computed without an eigensolver by <see cref="CentreLineExactCount"/>. N = 2
/// under γ = 0.1 on every site: the Heisenberg chain (J/4)Σσ·σ gives 10 on the line and 4 at the
/// point at J = 1, 10³, 10⁴ and 10⁸, the XY chain (J/2)(XX+YY) gives 10 and 0; at γ = 0 and J = 1
/// the line holds all 16 and the point 10 (Heisenberg) or 6 (XY). H = 0 at a uniform rate γ &gt; 0
/// gives 2^N·C(N, N/2) on both at even N and 0 at odd N.
/// A floating window cannot tell a mode on the line from one 5·10⁻⁸ beside it, and the exact
/// route never needs to.</para></summary>
public sealed class SelfMirrorObject : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private readonly MirrorSystem _system;   // INHERITED FROM: the x/y/z frame lives here, not on the object.
    private CentreLineExactCount.Result? _exact;
    private bool _exactTried;

    public SelfMirrorObject(MirrorSystem system) =>
        _system = system ?? throw new ArgumentNullException(nameof(system));

    /// <summary>Inherited from the system, not owned.</summary>
    public int N => _system.N;

    /// <summary>Inherited from the system: σ = Σγ, the total dephasing (and the frame's center scale).</summary>
    public double Sigma => _system.TotalDephasing;

    /// <summary>The object's place in the inherited frame: Re λ = −σ, the palindrome center line.</summary>
    public double Center => -_system.TotalDephasing;

    private CentreLineExactCount.Result? Exact
    {
        get
        {
            if (!_exactTried)
            {
                _exact = CentreLineExactCount.TryCount(
                    N, _system.Hamiltonian, _system.Channels.Select(c => c.Gamma).ToArray());
                _exactTried = true;
            }
            return _exact;
        }
    }

    /// <summary>True when the exact route ran: every connected block of L fits
    /// <see cref="CentreLineExactCount.MaxBlockDimension"/>. False is a cost bound, not a
    /// precision verdict; the counts then throw rather than guess.</summary>
    public bool IsFixedSetResolved => Exact is not null;

    private static string F(double value) =>
        value != 0.0 && Math.Abs(value) < 1e-4
            ? value.ToString("0.####E+0", Inv)
            : value.ToString("0.####", Inv);

    /// <summary>Algebraic multiplicity on the fixed line of the conjugate-composite map
    /// λ ↦ −2σ − conj(λ), equivalently Re λ = −σ.</summary>
    public int CompositeFixedLineCount => EnsureResolved().LineCount;

    /// <summary>Algebraic multiplicity at the fixed point of the linear F1 map λ ↦ −2σ − λ, which
    /// is λ = −σ: both the centre rate and zero frequency. Always a subset of
    /// <see cref="CompositeFixedLineCount"/>; the two coincide when every centre-line mode has
    /// Im λ = 0 (H = 0, for one) and part when H moves some of them off the real axis.</summary>
    public int LinearF1FixedPointCount => EnsureResolved().PointCount;

    private CentreLineExactCount.Result EnsureResolved() =>
        Exact ?? throw new InvalidOperationException(
            $"fixed-set count not computed: a connected block of L exceeds the exact route's " +
            $"{CentreLineExactCount.MaxBlockDimension}-dimensional cost bound");

    public string DisplayName =>
        $"SelfMirrorObject (conjugate-composite fixed line Re λ = −σ = {F(Center)}; ⊂ MirrorSystem N={N})";

    public string Summary =>
        IsFixedSetResolved
        ? $"an object INSIDE the system: {CompositeFixedLineCount} modes on the composite fixed line " +
          $"Re λ = −σ of λ ↦ −2σ − conj(λ); {LinearF1FixedPointCount} of them sit at λ = −σ, the fixed point " +
          $"of linear F1 λ ↦ −2σ − λ (exact algebraic multiplicities, no eigensolver). Everything else " +
          $"(x/y/z and σ = {F(Sigma)}) is INHERITED."
        : $"fixed-set count not computed: a connected block of L exceeds the exact route's " +
          $"{CentreLineExactCount.MaxBlockDimension}-dimensional cost bound (σ = {F(Sigma)}).";

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
                             $"{LinearF1FixedPointCount} of them also have Im λ = 0 and are fixed by the linear F1 map " +
                             $"λ ↦ −2σ − λ, whose fixed-point equation is λ = −σ; the other " +
                             $"{CompositeFixedLineCount - LinearF1FixedPointCount} sit on the line at Im λ ≠ 0 and pair " +
                             $"with their conjugates. Counted exactly: a Gaussian-integer characteristic polynomial per " +
                             $"connected block of L, roots on the line by Sturm sequences, none rounded.")
                : new InspectableNode(
                    displayName: "the object itself: fixed-set count not computed",
                    summary: $"a connected block of L exceeds {CentreLineExactCount.MaxBlockDimension} dimensions, the exact route's cost bound");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
