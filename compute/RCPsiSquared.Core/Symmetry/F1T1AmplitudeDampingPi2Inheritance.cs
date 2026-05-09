using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F1's palindrome residual under T1 amplitude damping (per
/// <c>memory/project_palindrome_frobenius_scaling.md</c> + F1 open-question Item 2)
/// has the empirical closed form:
///
/// <code>
///   ‖M‖² = 2^(N+2) · n_YZ · ‖H‖²_F  +  4^(N−1) · [3 · Σγ_T1²  +  4 · (Σγ_T1)²]
/// </code>
///
/// <para>verified bit-exact at N = 3..6. Each scaling factor and integer multiplier
/// in this formula sits on the Pi2-Foundation. The analytic derivation is the
/// open Tier 1 promotion path (F1 open-question Item 2); this claim makes the
/// Pi2-Foundation anchoring of the existing empirical formula explicit.</para>
///
/// <list type="bullet">
///   <item><b>"2^(N+2)" prefactor of the H-part</b>: <c>a_{−(N+1)}</c> on the Pi2
///         dyadic ladder. At N=2: 2⁴ = 16; at N=3: 2⁵ = 32. The exponent shift
///         <c>−(N+1)</c> matches the F1 derivation overhead one Π squaring above
///         the F49 power-factor ladder index <c>5−2N</c>.</item>
///   <item><b>"4^(N−1)" prefactor of the T1-only part</b>: <c>a_{3−2N}</c> on the
///         Pi2 dyadic ladder, equivalently d² for (N−1) qubits via
///         <see cref="Pi2OperatorSpaceMirrorClaim"/>. At N=2: 4; at N=3: 16; at
///         N=4: 64. This matches F49's own factor offset by one qubit (F49 uses
///         <c>4^(N−2)</c>; this T1 piece uses <c>4^(N−1)</c>).</item>
///   <item><b>"4" multiplier</b> in <c>4·(Σγ_T1)²</c>: exactly <c>a_{−1}</c> on the
///         Pi2 dyadic ladder = d² for 1 qubit. Same anchor used by F86 t_peak and
///         F49 (via Pi2OperatorSpaceMirror).</item>
///   <item><b>"3" multiplier</b> in <c>3·Σγ_T1²</c>: small integer from the T1
///         dissipator algebra; not directly a Pi2 anchor but algebraically clean.
///         Documented as such; not claimed as Pi2-Foundation inheritance.</item>
/// </list>
///
/// <para>Tier outcome: <b>Tier1Candidate</b>. The closed form itself is currently
/// Tier1Candidate (verified N=3..6 empirically; analytic derivation open per
/// F1 open-question Item 2). The Pi2-Anker readings of its constants are
/// algebraic-trivial composition once the formula is granted. When the analytic
/// derivation lands, this claim's status (and the underlying formula's status)
/// promotes to Tier1Derived together.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F82 + F84 (T1-related closed
/// forms; the F82/F84 closed forms cover the 2^(N−1) F81-anti-component, this
/// claim's H-part 2^(N+2) and full-‖M‖² scaling are the open extension) +
/// <c>compute/RCPsiSquared.Core/F1/F1OpenQuestions.cs</c> (Item 2: T1 amplitude
/// damping full closed form) +
/// <c>compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs</c> (the parent F1
/// identity whose residual ‖M‖² is being scaled) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs</c>.</para></summary>
public sealed class F1T1AmplitudeDampingPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly Pi2OperatorSpaceMirrorClaim _mirror;

    /// <summary>The "4" multiplier in <c>4·(Σγ_T1)²</c>. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(−1) = <c>a_{−1}</c> = 4 = d² for N=1.</summary>
    public double FourMultiplier => _ladder.Term(-1);

    /// <summary>The small integer multiplier "3" in <c>3·Σγ_T1²</c>. From the T1
    /// dissipator algebra; documented but not a direct Pi2-Foundation anchor.</summary>
    public const int ThreeMultiplier = 3;

    /// <summary>The "2^(N+2)" prefactor of the H-part. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(<c>−(N+1)</c>) =
    /// <c>2^{1−(−(N+1))} = 2^{N+2}</c>.</summary>
    public double HPartPrefactor(int N)
    {
        if (N < 2)
            throw new ArgumentOutOfRangeException(nameof(N), N, "T1 amplitude damping scaling requires N ≥ 2.");
        return _ladder.Term(HPartLadderIndex(N));
    }

    /// <summary>The Pi2 ladder index where the H-part prefactor lands: <c>−(N+1)</c>.
    /// At N=3 this is −4 (a_{−4} = 32); at N=4 it is −5 (a_{−5} = 64).</summary>
    public int HPartLadderIndex(int N) => -(N + 1);

    /// <summary>The "4^(N−1)" prefactor of the T1-only part. Exactly equal to
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(<c>3−2N</c>) =
    /// <c>2^{1−(3−2N)} = 2^{2N−2} = 4^{N−1}</c>. Matches d² for (N−1) qubits via
    /// <see cref="Pi2OperatorSpaceMirrorClaim"/>.</summary>
    public double T1PartPrefactor(int N)
    {
        if (N < 2)
            throw new ArgumentOutOfRangeException(nameof(N), N, "T1 amplitude damping scaling requires N ≥ 2.");
        return _ladder.Term(T1PartLadderIndex(N));
    }

    /// <summary>The Pi2 ladder index where the T1-only prefactor lands: <c>3−2N</c>.
    /// At N=3 this is −3 (a_{−3} = 16 = d² for N=2 qubits); at N=4 it is −5
    /// (a_{−5} = 64 = d² for N=3 qubits).</summary>
    public int T1PartLadderIndex(int N) => 3 - 2 * N;

    /// <summary>The qubit count whose operator-space dimension equals the T1-only
    /// prefactor: <c>N − 1</c>. Connects to <see cref="Pi2OperatorSpaceMirrorClaim"/>.
    /// At N_chain=3 the T1-prefactor is d² for 2 qubits; at N_chain=8 it is d²
    /// for 7 qubits (N=7 not in pinned mirror table; only N=1..6 pinned).</summary>
    public int T1PartOperatorSpaceQubitCount(int N) => N - 1;

    public F1T1AmplitudeDampingPi2Inheritance(Pi2DyadicLadderClaim ladder, Pi2OperatorSpaceMirrorClaim mirror)
        : base("F1 T1-amplitude-damping ‖M‖² scaling factors inherit from Pi2-Foundation (analytic derivation: Tier 1 candidate)",
               Tier.Tier1Candidate,
               "docs/ANALYTICAL_FORMULAS.md F82 + " +
               "compute/RCPsiSquared.Core/F1/F1OpenQuestions.cs + " +
               "compute/RCPsiSquared.Core/F1/F1PalindromeIdentity.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2OperatorSpaceMirrorClaim.cs")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _mirror = mirror ?? throw new ArgumentNullException(nameof(mirror));
    }

    public override string DisplayName =>
        "F1 T1-amplitude-damping closed-form factors as Pi2-Foundation inheritance";

    public override string Summary =>
        $"‖M‖² = 2^(N+2)·n_YZ·‖H‖² + 4^(N−1)·[3·Σγ² + 4·(Σγ)²]: 2^(N+2) = a_{{−(N+1)}}, 4^(N−1) = a_{{3−2N}} = d² für (N−1) qubits, 4 = a_{{−1}}; analytic derivation open ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("empirical closed form",
                summary: "‖M‖² = 2^(N+2)·n_YZ·‖H‖²_F + 4^(N−1)·[3·Σγ_T1² + 4·(Σγ_T1)²]; verified bit-exact N=3..6");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "2^(N+2) = a_{−(N+1)} (ladder); 4^(N−1) = a_{3−2N} (ladder, = d² für (N−1) qubits via OperatorSpaceMirror); 4 = a_{−1}");
            yield return new InspectableNode("3 multiplier",
                summary: "small integer from T1 dissipator algebra; documented but not a direct Pi2-Foundation anchor");
            yield return InspectableNode.RealScalar("FourMultiplier (= a_{-1})", FourMultiplier);
            yield return InspectableNode.RealScalar("ThreeMultiplier", ThreeMultiplier);
            for (int N = 2; N <= 6; N++)
            {
                yield return new InspectableNode(
                    $"chain N={N}",
                    summary: $"H-part prefactor 2^(N+2) = {HPartPrefactor(N)} (a_{HPartLadderIndex(N)}); " +
                             $"T1-part prefactor 4^(N−1) = {T1PartPrefactor(N)} (a_{T1PartLadderIndex(N)} = d² for {T1PartOperatorSpaceQubitCount(N)} qubits)");
            }
            yield return new InspectableNode("Tier1 promotion path",
                summary: "F1 open-question Item 2: derive the closed form analytically (currently empirical N=3..6); when proven, this claim and parent formula promote to Tier1Derived together");
        }
    }
}
