using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F167, the sodium chain (Tier 1 derived, <c>experiments/THE_SODIUM_CHAIN.md</c>): n spins
/// 3/2, each a qubit pair (|3/2⟩, |1/2⟩, |−1/2⟩, |−3/2⟩ = |00⟩, |01⟩, |10⟩, |11⟩) relaxing by its
/// own field gradient, per nucleus H_k = ω_Q,k·T₂₀ + ω_k·I_z and the jumps ZZ at 3J₀,k/2, ZX and ZY
/// at 3J₁,k/2, XI and YI at 3J₂,k/2 (<see cref="Spin32QuadrupolarMirrorsClaim"/>), coupled by a real
/// two-body Hamiltonian on any graph. Q = ⊗_k ZZ_k, S(ρ) = Q·ρᵀ; P_k the projector on the cells
/// bright for nucleus k.
///
/// <para>What the page proves, for any n and any graph. The lemma: if H is real and [H, Q] = 0,
/// S·L(J₀)·S⁻¹ = −L(−J₀) − 6·Σ_k(J₁,k + J₂,k), equivalently −L + 2·L_{J₀} − 6Σ_k(J₁,k + J₂,k) with
/// L_{J₀} = −3·Σ_k J₀,k·P_k: the mirror image is the chain with its J₀ rates reversed, a negative
/// density. The ladder: if H commutes with every ZZ_k, the part whose bright nuclei form the set B
/// is palindromic about −s_B, s_B = 3·Σ_k(J₁,k + J₂,k) + 3·Σ_{k∈B} J₀,k (<see cref="TwiceCentre"/>),
/// for identical rates s_m = 3n(J₁ + J₂) + 3m·J₀ (<see cref="TwiceRung"/>); the Ising coupling and
/// the secular dipolar coupling between nuclei of unequal, non-zero splittings qualify. Equal
/// splittings: the coupling keeps Q but moves brightness, and K_m, the largest L-invariant subspace
/// inside the cells with m bright nuclei (identical rates), is palindromic about −s_m. The full
/// secular dipolar coupling: its eight √3 strings anticommute with Q and break S; at J₀ = 0,
/// J₁, J₂ > 0 the global mirror about −σ breaks (F158: the far span is Q alone and Q does not
/// commute with H, for any non-zero coefficient of the √3 strings); U_m, the largest L-invariant
/// subspace inside V = {x : H₃x = 0 = xH₃} and the m-bright cells, is palindromic about −s_m for any
/// real H₃ anticommuting with Q; K₀ is not palindromic about −s₀ (the identity has no partner).</para>
///
/// <para>Proved at two nuclei with identical rates, one offset and an equal splitting: K₁ = U₁ and
/// K₂ = U₂, both spaces written down (E₂ = (Z₁ + Z₃)·(X₂X₄ + Y₂Y₄), E₁ the order ±5 block and four
/// exchange-antisymmetric order ±4 operators), at every parameter off a proper algebraic set, since
/// E_m lies in U_m at every parameter and one exact point has dim K_m = dim E_m; the set contains
/// d = 0. What the page reads modulo two
/// primes, not proved: the dimensions of K_m (at most 24, 128, 24 at two nuclei for the equal
/// truncation, 64, 140, 84, 16 at three) and, for the full coupling, of K₀ (10) and U₀ (5); float
/// spectra of the U_m are a further reading. The live witness reads the two-nucleus dimensions
/// exactly at its parameters and gates E₁, E₂ against every parameter piece of L.</para>
///
/// <para>Adjacency the page names and leaves open: the lemma's kept part is the dark remainder of
/// the corner-shift row of <c>docs/THE_ONE_SQUARE.md</c> §9, F137's kept jump, and MirrorWorld's
/// GammaFold gain turn read on the J₀ channel alone; the ladder and the sharp/blurred split are the
/// blockwise centres of <c>experiments/THE_MIRROR_UNDER_A_SHARED_FIELD.md</c> in another
/// dissipator; the hulls K_m and U_m are the construction of F157's blind subspace. Whether one
/// statement holds them is not decided, and no typed edge joins them here.</para>
///
/// <para>Typed parents: <see cref="Spin32QuadrupolarMirrorsClaim"/>, the n = 1 case (its dark and
/// bright centres are s_∅ and s_{1}); <see cref="PalindromeTwoEndCountClaim"/> (F158), the route of
/// the full coupling's global break.</para>
///
/// <para>Live witness: <c>inspect --root sodiumchain</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/SodiumChainWitness.cs</c>), two nuclei, exact over
/// ℚ(i), the coupling read from the spin matrices, the √3 of the full coupling carried by a diagonal
/// similarity that makes the generator rational; a different road from the page's verifier
/// <c>simulations/sodium_chain_mirrors.py</c>.</para></summary>
public sealed class SodiumChainMirrorLadderClaim : Claim
{
    public Spin32QuadrupolarMirrorsClaim Nucleus { get; }
    public PalindromeTwoEndCountClaim TwoEndCount { get; }

    /// <summary>6·Σ_k(J₁,k + J₂,k), the constant of the lemma: twice the reflected rate, each of the four
    /// reflected jumps of nucleus k paying −2·(3J/2).</summary>
    public static long LemmaConstant(IReadOnlyList<(long J0, long J1, long J2)> rates) =>
        checked(6 * rates.Sum(r => r.J1 + r.J2));

    /// <summary>2s_B = 6·Σ_k(J₁,k + J₂,k) + 6·Σ_{k∈B} J₀,k, bit k of <paramref name="brightMask"/> set
    /// when nucleus k is in B: the part with bright set B pairs about −s_B.</summary>
    public static long TwiceCentre(IReadOnlyList<(long J0, long J1, long J2)> rates, int brightMask)
    {
        long s = LemmaConstant(rates);
        for (int k = 0; k < rates.Count; k++)
            if ((brightMask >> k & 1) != 0) s = checked(s + 6 * rates[k].J0);
        return s;
    }

    /// <summary>2s_m = 6n(J₁ + J₂) + 6m·J₀ for n nuclei with identical rates, m of them bright.</summary>
    public static long TwiceRung(int n, int m, long j0, long j1, long j2) => checked(6 * n * (j1 + j2) + 6 * m * j0);

    public SodiumChainMirrorLadderClaim(Spin32QuadrupolarMirrorsClaim nucleus, PalindromeTwoEndCountClaim twoEndCount)
        : base(
            "Coupled spin-3/2 nuclei, each a qubit pair with its own quadrupolar relaxation: for a real H commuting with " +
            "Q = ⊗ZZ_k, S(ρ) = Q·ρᵀ gives S·L(J₀)·S⁻¹ = −L(−J₀) − 6Σ_k(J₁,k + J₂,k); if H commutes with every ZZ_k, the part " +
            "with bright set B pairs about −s_B, s_B = 3Σ_k(J₁,k + J₂,k) + 3Σ_{k∈B} J₀,k; with equal splittings the ladder holds " +
            "on the sharp subspaces K_m; the full secular dipolar coupling breaks S at its √3 strings and, at J₀ = 0, the " +
            "global mirror, while the invariant part of the √3 coupling's two-sided kernel pairs on every rung",
            Tier.Tier1Derived,
            "experiments/THE_SODIUM_CHAIN.md + docs/ANALYTICAL_FORMULAS.md (F167, the index entry)")
    {
        Nucleus = nucleus ?? throw new ArgumentNullException(nameof(nucleus));
        TwoEndCount = twoEndCount ?? throw new ArgumentNullException(nameof(twoEndCount));
    }

    public override string DisplayName =>
        "F167: the sodium chain, a ladder of mirrors one rung per set of bright nuclei";

    public override string Summary =>
        "S(ρ) = (⊗ZZ_k)·ρᵀ maps the chain onto its J₀-reversed image, S·L(J₀)·S⁻¹ = −L(−J₀) − 6Σ(J₁ + J₂); brightness-keeping " +
        "couplings give one centre −s_B per bright set, the equal-splitting flip-flop leaves it on K_m, the full dipolar " +
        $"coupling breaks the map at its √3 strings ({Tier.Label()}; live witness: inspect --root sodiumchain)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the lemma",
                summary: "H real, [H, Q] = 0: S·L(J₀)·S⁻¹ = −L(−J₀) − 6Σ_k(J₁,k + J₂,k); the transpose reverses −i[H, ·], Q keeps it; " +
                         "a string jump P commuting with Q (ZZ_k) is kept, one anticommuting (ZX, ZY, XI, YI) is reflected at −2γ; any n, any graph");
            yield return new InspectableNode("the ladder",
                summary: "[H, ZZ_k] = 0 for every k: P_k commutes with L, and the part with bright set B pairs about −s_B, " +
                         "s_B = 3Σ_k(J₁,k + J₂,k) + 3Σ_{k∈B} J₀,k; identical rates: n + 1 rungs 3J₀ apart, the parent's two at n = 1; " +
                         "the Ising coupling and the unequal-splitting truncation of the dipolar coupling qualify");
            yield return new InspectableNode("equal splittings and the full coupling",
                summary: "equal splittings: the satellite flip-flop keeps Q, moves m by 0 or ±2, and K_m pairs about −s_m; the full " +
                         "secular dipolar coupling: eight √3 strings anticommuting with Q break S, the global mirror about −σ breaks at " +
                         "J₀ = 0 (F158), U_m inside the √3 coupling's two-sided kernel pairs on every rung, K₀ does not pair about −s₀");
        }
    }
}
