using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F167, the spin-3/2 nucleus and its two mirrors (Tier 1 derived,
/// <c>experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md</c>): one spin 3/2, levels m = 3/2, 1/2, −1/2,
/// −3/2, under the Redfield generator of quadrupolar relaxation, secular in the Larmor frequency with the
/// static splitting kept whole, at infinite temperature,
/// dissipator Σ_m J_|m|·(T₂ₘρT₂ₘ† − ½{T₂ₘ†T₂ₘ, ρ}) (m = −2 … 2, T₂₀ = (3I_z² − I(I + 1))/√6,
/// T₂,±₁ = ∓(I_zI_± + I_±I_z)/2, T₂,±₂ = I_±²/2), Hamiltonian ω_Q·T₂₀ + ω·I_z.
///
/// <para>What the page proves: the five Hermitian components Γ₀ = T₂₀, Γ₁ = T₂,₁ + T₂,₁†,
/// Γ₂ = i(T₂,₁ − T₂,₁†), Γ₃ = T₂,₂ + T₂,₂†, Γ₄ = i(T₂,₂ − T₂,₂†) anticommute pairwise and square to
/// 3/2, 3, 3, 3, 3 times 1 (five Dirac matrices; in the labelling |3/2⟩, |1/2⟩, |−1/2⟩, |−3/2⟩ =
/// |00⟩, |01⟩, |10⟩, |11⟩ the strings ZZ, ZX, ZY, XI, YI up to scalars), and among the spins I ≥ 1
/// only at 3/2, so the dissipator is Pauli-string dephasing. With Γ₀ scaled to diag(1, −1, −1, 1)
/// the one map S(X) = Γ₀·Xᵀ gives S·L·S⁻¹ = −L + 2·L_{J₀} − 6(J₁ + J₂), L_{J₀} = −3J₀·P (P the
/// projector on the bright cells |m⟩⟨m′|, |m| ≠ |m′|, which commutes with L and with S), at every
/// ω_Q and ω: the dark modes pair about −c₁, c₁ = 3(J₁ + J₂), the bright ones about −c₂,
/// c₂ = c₁ + 3J₀. The global-centre theorem: the sixteen eigenvalues have one centre −s,
/// 2s = 3J₀ + 6J₁ + 6J₂ (fixed by the trace), exactly when J₀ = 0, or ω_Q = ω = 0 and J₁·J₂ = 0,
/// for every real J (<see cref="HasGlobalCentre"/>). With a splitting, q = ω_Q·√6/2 the coefficient of
/// diag(1, −1, −1, 1) in H, each bright pair is a closed
/// 2 × 2 block whose two rates sum to 2c₂, and at 2q = 3J₂ (3J₁ for the double quantum), that density
/// non-zero, it is one 2 × 2 Jordan block on Re λ = −c₂: a defective exceptional point where the
/// splitting, which carries each bright string into its own mirror image, makes a mode and its mirror
/// image one, an exceptional point of the eigenvalues; where the lines of a spectrum separate depends on
/// the observable (the page's rule). The dynamic frequency shift, the Lamb-shift Hamiltonian −3K₂·ZI − 3K₁·IZ
/// of the shifts K(m·ω₀) = Im Γ, is real and diagonal: both mirrors and the exceptional points survive it, and
/// the global centre follows <see cref="HasGlobalCentreWithShift"/>.</para>
///
/// <para>Shape and neighbours. At J₀ = 0 the map S is an F1 palindromizer built in F118's shape
/// (<c>docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md</c>: Π_Y(ρ) = X^⊗N·ρᵀ), with Γ₀ in place of
/// X^⊗N. F121's qudit family lives in the same d = 4 Hilbert space under another dissipator
/// (equidistant dephasing, 8 of 16 modes about one centre against 16 of 16 about two here); the
/// two do not contradict each other, and no typed edge joins them. Prior art: the rates are Souza
/// et al.'s (arXiv:1105.2539, Eqs. (15) to (27)), and so is the reading of the four levels as two
/// qubits; the Dirac algebra of the five quadrupole components is Wu, Hu and Zhang's
/// (Phys. Rev. Lett. 91, 186402 (2003)).</para>
///
/// <para>Typed parent: <see cref="PalindromeTwoEndCountClaim"/> (F158), which for non-negative
/// densities gives the global-centre theorem a second, independent route (the normalised jumps are
/// Hermitian strings squaring to 1, σ = s). The coupled-nuclei ladder of
/// <c>experiments/THE_SODIUM_CHAIN.md</c>, indexed by F167 as well, is typed as its child
/// <see cref="SodiumChainMirrorLadderClaim"/>.</para>
///
/// <para>Live witness: <c>inspect --root sodium</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/SodiumTwoMirrorsWitness.cs</c>), which builds L
/// from the spin matrices over ℚ(i), the √3 of I_± factored out exactly and the √6 of T₂₀ taken
/// into the integer splitting parameter q = ω_Q·√6/2, and reads the theorem by exact kernel ranks
/// and an exact characteristic polynomial, a different road from the page's sympy verifier
/// <c>simulations/sodium_two_mirrors.py</c>.</para></summary>
public sealed class Spin32QuadrupolarMirrorsClaim : Claim
{
    public PalindromeTwoEndCountClaim TwoEndCount { get; }

    /// <summary>The dark centre c₁ = 3(J₁ + J₂): the dark modes pair about −c₁.</summary>
    public static long DarkCentre(long j1, long j2) => checked(3 * (j1 + j2));

    /// <summary>The bright centre c₂ = 3(J₀ + J₁ + J₂) = c₁ + 3J₀: the bright modes pair about −c₂.</summary>
    public static long BrightCentre(long j0, long j1, long j2) => checked(3 * (j0 + j1 + j2));

    /// <summary>2s = 3J₀ + 6J₁ + 6J₂, the global candidate centre −s fixed by the trace
    /// (Tr L = −24J₀ − 48J₁ − 48J₂, the Hamiltonian part being traceless).</summary>
    public static long TwiceGlobalCentre(long j0, long j1, long j2) => checked(3 * j0 + 6 * j1 + 6 * j2);

    /// <summary>The global-centre rule for real densities and frequencies: one centre for the whole
    /// spectrum exactly when J₀ = 0, or ω_Q = ω = 0 and J₁·J₂ = 0. The splitting enters as
    /// q = ω_Q·√6/2, the coefficient of diag(1, −1, −1, 1) in H, so ω_Q = 0 exactly when q = 0.</summary>
    public static bool HasGlobalCentre(long j0, long j1, long j2, long q, long omega) =>
        j0 == 0 || (q == 0 && omega == 0 && (j1 == 0 || j2 == 0));

    /// <summary>The global-centre rule with the dynamic frequency shift, the Lamb-shift Hamiltonian
    /// K₁(T₂,₁†T₂,₁ − T₂,₋₁†T₂,₋₁) + K₂(T₂,₂†T₂,₂ − T₂,₋₂†T₂,₋₂) of the shifts K(m·ω₀) = Im Γ, odd in m: one centre
    /// exactly when J₀ = 0, or q = 0 and J₁ = 0 with ω = 6K₁ (the satellites on resonance), or q = 0 and J₂ = 0 with
    /// ω = 3K₂ (the double quanta on resonance). At K₁ = K₂ = 0 it is <see cref="HasGlobalCentre"/>.</summary>
    public static bool HasGlobalCentreWithShift(long j0, long j1, long j2, long q, long omega, long k1, long k2) =>
        j0 == 0 || (q == 0 && ((j1 == 0 && omega == checked(6 * k1)) || (j2 == 0 && omega == checked(3 * k2))));

    public Spin32QuadrupolarMirrorsClaim(PalindromeTwoEndCountClaim twoEndCount)
        : base(
            "Spin-3/2 quadrupolar relaxation: the five Hermitian quadrupole components are Dirac matrices (only at " +
            "I = 3/2), so the dissipator is Pauli-string dephasing; S(X) = Γ₀·Xᵀ gives S·L·S⁻¹ = −L + 2·L_{J₀} − 6(J₁ + J₂) " +
            "at every ω_Q and ω, the dark modes pairing about −3(J₁ + J₂) and the bright ones about a centre 3J₀ further " +
            "out; the whole spectrum has one centre exactly when J₀ = 0, or ω_Q = ω = 0 and J₁·J₂ = 0",
            Tier.Tier1Derived,
            "experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md + docs/ANALYTICAL_FORMULAS.md (F167, the index entry)")
    {
        TwoEndCount = twoEndCount ?? throw new ArgumentNullException(nameof(twoEndCount));
    }

    public override string DisplayName =>
        "F167: the spin-3/2 nucleus carries two exact mirrors, one map S(X) = Γ₀·Xᵀ";

    public override string Summary =>
        "five Dirac matrices at I = 3/2 make quadrupolar relaxation Pauli-string dephasing; S·L·S⁻¹ = −L + 2·L_{J₀} − 6(J₁ + J₂), " +
        "dark modes about −3(J₁ + J₂), bright ones about −3(J₀ + J₁ + J₂); one global centre −s, 2s = 3J₀ + 6J₁ + 6J₂, " +
        $"exactly when J₀ = 0, or ω_Q = ω = 0 and J₁·J₂ = 0 ({Tier.Label()}; live witness: inspect --root sodium)";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("five Dirac matrices",
                summary: "Γ₀ = T₂₀, Γ₁ = T₂,₁ + T₂,₁†, Γ₂ = i(T₂,₁ − T₂,₁†), Γ₃ = T₂,₂ + T₂,₂†, Γ₄ = i(T₂,₂ − T₂,₂†) anticommute " +
                         "pairwise and square to 3/2, 3, 3, 3, 3 times 1 at I = 3/2, the strings ZZ, ZX, ZY, XI, YI up to scalars; " +
                         "at every other spin I ≥ 1, Γ₀² is no multiple of 1 (3m² − I(I + 1) does not take two opposite values)");
            yield return new InspectableNode("one map, two mirrors",
                summary: "S(X) = Γ₀·Xᵀ, Γ₀ = diag(1, −1, −1, 1): S·L·S⁻¹ = −L + 2·L_{J₀} − 6(J₁ + J₂), L_{J₀} = −3J₀·P, P the bright " +
                         "projector commuting with L and S; the transpose reverses the real H and keeps every jump term, Γ₀ commutes " +
                         "with H, reflects the J₁ and J₂ dissipators and keeps the J₀ one; at J₀ = 0 an F1 palindromizer in F118's shape");
            yield return new InspectableNode("the bright pairs' exceptional point",
                summary: "with a splitting q = ω_Q·√6/2 each bright pair is a closed 2 × 2 block, the satellites at −c₂ − iω ± √(9J₂² − 4q²), " +
                         "the double quantum at −c₂ − 2iω ± √(9J₁² − 4q²), so the two rates of a pair sum to 2c₂ at every splitting; at " +
                         "2q = 3J₂ (resp. 3J₁), that density non-zero, the pair is one 2 × 2 Jordan block on Re λ = −c₂, a defective " +
                         "exceptional point of the eigenvalues (where lines separate depends on the observable); at " +
                         "every q the splitting q·diag(1, −1, −1, 1) carries each bright string X into its own mirror image S(X) = ±Γ₀X, and the coalesced " +
                         "vector is self-orthogonal, vᵀv = 0");
            yield return new InspectableNode("the dynamic frequency shift",
                summary: "the shifts K(m·ω₀) = Im Γ, odd in m, add the Lamb-shift Hamiltonian " +
                         "−3K₂·ZI − 3K₁·IZ = −3K₂·I_z − (3/2)(2K₁ − K₂)·IZ, real and diagonal: S·L·S⁻¹ = −L + 2·L_{J₀} − 6(J₁ + J₂) " +
                         "still holds, each bright pair moves as a whole (6K₁, 6K₂), so the exceptional points stay, the central " +
                         "transition moves 6(2K₁ − K₂) off the satellites, and one global centre exists exactly when J₀ = 0, or " +
                         "q = 0 and J₁ = 0 with ω = 6K₁, or q = 0 and J₂ = 0 with ω = 3K₂");
            yield return new InspectableNode("global-centre theorem",
                summary: "one centre −s, 2s = 3J₀ + 6J₁ + 6J₂, exactly when J₀ = 0, or ω_Q = ω = 0 and J₁·J₂ = 0, every real J; the " +
                         "coefficients of χ(−x − 2s) − χ(x) generate the ideal (J₀·ω_Q², J₀·(9J₁² + ω²)·(9J₂² + 4ω²)); for J ≥ 0 " +
                         "F158 gives it again as dim ker L = dim ker(L + 2s)");
        }
    }
}
