using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F172, the clock mirror (Tier 1 derived, <c>docs/proofs/PROOF_QUDIT_CLOCK_MIRROR.md</c>): N sites of d
/// levels, the clock Z|j⟩ = ω^j|j⟩ (ω = e^{2πi/d}, the Weyl clock operator, not MirrorWorld's <c>Clock</c>), the shift
/// X|j⟩ = |j + 1⟩, U = (X^{d/2})^⊗N.
///
/// <para>What the page proves. Under the clock pair {Z_l, Z_l†} at one rate γ, one site's letters X^aZ^b decay at
/// r(a) = 2γ(1 − cos(2πa/d)), the column alone, and the half-shift a ↦ a + d/2 pairs every rate with the one the
/// palindrome needs exactly when d is even. At even d, for every Hermitian H commuting with U, ρ ↦ Uρ carries L onto
/// −L† − 4Nγ (Prosen's PT relation; an instance of the Lemma of his Brief Report of 2012), and for real symmetric H,
/// ρ ↦ ρᵀU gives −L − 4Nγ, the shape of F1's Π. Writing j = (d/2)·b + j', Z = Z_b ⊗ diag(ω^{j'}) and
/// X^{d/2} = X_b ⊗ 1: the element is the X of a qubit every even qudit carries. F158's two-end count carries over to
/// unitary jumps whose dissipator is self-adjoint, so for the pair at even d the palindrome holds exactly when some
/// invertible UΔ, Δ diagonal, commutes with H. The lone jump Z pairs about −Nγ when H commutes with U and meets one of
/// four routes ([Q, H] = 0; Q H̄ Q = H; Y H Y⁻¹ = −H with Y = (Z^{d/2})^⊗N; H̄ = −H), and with the counts agreeing its
/// palindrome is exactly spec L(H, Z) = spec L(H, Z†). At odd d no Hamiltonian makes the spectrum palindromic, with the
/// pair or the lone jump, and one-rate full-Cartan dephasing at d ≥ 3 never: the numerical range puts every eigenvalue
/// above −N·r_max, while 0's partner would sit at −2N times the mean rate.</para>
///
/// <para>Typed parent: <see cref="PalindromeTwoEndCountClaim"/> (F158), whose count the page carries from Hermitian
/// involutions to unitary jumps with a self-adjoint dissipator.</para>
///
/// <para>Live witness: <c>inspect --root clockmirror</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/QuditClockMirrorWitness.cs</c>), which reads the mirror identities
/// exactly in ℤ[ζ_m] on every matrix unit, a different road from the page's gate
/// <c>simulations/qudit_clock_mirror_gate.py</c> (superoperators modulo two primes under a norm bound).</para></summary>
public sealed class QuditClockMirrorClaim : Claim
{
    public PalindromeTwoEndCountClaim TwoEndCount { get; }

    /// <summary>The clock pair's rate of column a in units of γ: 2(1 − cos(2πa/d)).</summary>
    public static double PairRateOverGamma(int a, int d)
    {
        if (d < 2) throw new ArgumentOutOfRangeException(nameof(d), d, "d ≥ 2");
        return 2.0 * (1.0 - Math.Cos(2.0 * Math.PI * a / d));
    }

    /// <summary>The element exists exactly at even d: ℤ_d has an element of order two.</summary>
    public static bool HasHalfShift(int d) => d >= 2 && d % 2 == 0;

    /// <summary>The pair's centre, −2Nγ; the reflection is λ ↦ −λ − 4Nγ.</summary>
    public static long PairTwiceCentre(int n, long gamma) => checked(4L * n * gamma);

    /// <summary>The lone jump's centre, −Nγ; the reflection is λ ↦ −λ − 2Nγ.</summary>
    public static long LoneTwiceCentre(int n, long gamma) => checked(2L * n * gamma);

    public QuditClockMirrorClaim(PalindromeTwoEndCountClaim twoEndCount)
        : base(
            "Clock dephasing of d-level chains: at even d the half-shift U = (X^{d/2})^⊗N, the X of a qubit every even qudit " +
            "carries, reflects L onto −L† − 4Nγ for every H commuting with U, and the palindrome holds exactly when some " +
            "invertible UΔ commutes with H (F158's count for unitary jumps with a self-adjoint dissipator); the lone jump pairs " +
            "on four routes, its palindrome being spec L(H, Z) = spec L(H, Z†) where the counts agree; at odd d, and under " +
            "one-rate full-Cartan dephasing at d ≥ 3, no Hamiltonian makes the spectrum palindromic",
            Tier.Tier1Derived,
            "docs/proofs/PROOF_QUDIT_CLOCK_MIRROR.md + docs/ANALYTICAL_FORMULAS.md (F172, the index entry)")
    {
        TwoEndCount = twoEndCount ?? throw new ArgumentNullException(nameof(twoEndCount));
    }

    public override string DisplayName =>
        "F172: the clock mirror, a qubit inside every even qudit";

    public override string Summary =>
        "clock jumps {Z, Z†}: palindromic at even d exactly when an invertible UΔ commutes with H, U = (X^{d/2})^⊗N the hidden " +
        "qubit's X; the lone jump Z on four routes, reduced to kicks forward against kicks back; odd d and one-rate full-Cartan " +
        $"at d ≥ 3 never, by the numerical range ({Tier.Label()}; live witness: inspect --root clockmirror)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the rates on a ring",
                summary: "r(a) = 2γ(1 − cos(2πa/d)) on the column a of the letters X^aZ^b; a ↦ a + d/2 sends r to 4γ − r, the element " +
                         "of order two of ℤ_d, present exactly at even d");
            yield return new InspectableNode("the element and its widening",
                summary: "ρ ↦ Uρ gives −L† − 4Nγ for [H, U] = 0; the far end is U times the diagonals, so the palindrome holds exactly " +
                         "when an invertible UΔ commutes with H; for real symmetric H, ρ ↦ ρᵀU gives −L − 4Nγ");
            yield return new InspectableNode("the lone jump",
                summary: "four routes, one of each kind of symmetry: [Q, H] = 0, Q H̄ Q = H, Y H Y⁻¹ = −H, H̄ = −H, each with [H, U] = 0; " +
                         "where the counts agree the palindrome is spec L(H, Z) = spec L(H, Z†), and the count alone does not decide");
            yield return new InspectableNode("odd d",
                summary: "Re λ ≥ −N·r_max with r_max = 2γ(1 + cos(π/d)) < 4γ, while a palindrome needs 0's partner at −4Nγ: no " +
                         "Hamiltonian; the same bound excludes one-rate full-Cartan dephasing at every d ≥ 3");
        }
    }
}
