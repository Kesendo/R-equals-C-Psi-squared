using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The self-mirror fixed point AS AN OBJECT in the Object Manager: not the system, an object
/// INSIDE it. This is the jump that was missing. The SYSTEM (a <see cref="MirrorSystem"/>) bears the
/// x/y/z frame , z = the watched axis (Re λ), x-y = the coherent-motion plane (Im λ), the center
/// σ = Σγ , and the object INHERITS that frame; it does not own it. The object's only delta, the one
/// thing it adds, is BEING the F1 fixed point: the modes at Re λ = −σ that are their own mirror
/// partner under λ ↦ −2σ − λ. x/y/z is the system's inheritance, not the object's property; the object
/// is built FROM a system and never carries the frame itself.
///
/// <para>Not a Claim, a live reading. The fixed sector it names is the k=N/2 self-mirror rung,
/// populated for even N and empty for odd N (half-integer w_XY = N/2). It breadcrumbs the already-typed
/// facets (F1PalindromeIdentity, HalfIntegerMirrorClaim, HalfAsStructuralFixedPointClaim) rather than
/// re-deriving them.</para></summary>
public sealed class SelfMirrorObject : IInspectable
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;
    private const double Tol = 1e-7;

    private readonly MirrorSystem _system;   // INHERITED FROM: the x/y/z frame lives here, not on the object.

    public SelfMirrorObject(MirrorSystem system) =>
        _system = system ?? throw new ArgumentNullException(nameof(system));

    /// <summary>Inherited from the system, not owned.</summary>
    public int N => _system.N;

    /// <summary>Inherited from the system: σ = Σγ, the total dephasing (and the frame's center scale).</summary>
    public double Sigma => _system.TotalDephasing;

    /// <summary>The object's place in the inherited frame: Re λ = −σ, the palindrome center line.</summary>
    public double Center => -_system.TotalDephasing;

    /// <summary>The object's own delta: how many modes are their OWN F1 mirror , decay rate r = σ
    /// (Re λ = −σ), the fixed sector of r ↦ 2σ − r. Even N populates the k=N/2 rung; odd N leaves it
    /// empty (half-integer w_XY = N/2). Read live from the inherited spectrum.</summary>
    public int SelfPairedCount =>
        _system.Spectrum.Modes.Count(m => Math.Abs(m.ActualDecayRate - Sigma) < Tol);

    public string DisplayName =>
        $"SelfMirrorObject (the F1 fixed point at Re λ = −σ = {Center.ToString("0.####", Inv)}; ⊂ MirrorSystem N={N})";

    public string Summary =>
        $"an object INSIDE the system. Its only delta: it IS the self-paired sector , {SelfPairedCount} modes at " +
        $"Re λ = −σ, each its own mirror. Everything else (x/y/z: z = watched Re λ, x-y = motion Im λ, σ = " +
        $"{Sigma.ToString("0.####", Inv)}) is INHERITED from the system, not owned.";

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

            // The object's own delta: the fixed point. This, and only this, IS the object.
            yield return new InspectableNode(
                displayName: "the object itself: the F1 fixed point (its only delta)",
                summary: $"{SelfPairedCount} modes with decay rate r = σ (Re λ = −σ): each is its OWN mirror under " +
                         $"λ ↦ −2σ − λ. The fixed sector of the F1 reflection = the k=N/2 self-mirror rung (even N populates " +
                         $"it, odd N leaves it empty per the half-integer w_XY = N/2). This is what the object IS; where it sits is inherited.");
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
