using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F170, the XY ring's frozen band at odd N: on the R₉₀ locus with γ̄ ≠ 0 every diagonal block
/// (p, p), 1 ≤ p ≤ N − 1, carries the frozen root −4γ̄ with geometric and algebraic multiplicity exactly
/// ⌊N/2⌋ for (J, δ) off a proper algebraic subset, and N − 1 at the uniform point for all but finitely many
/// J, with the wrap bond priced by occupation (the spin ring) or by order (the free periodic ring).
/// Recomputed live by <see cref="RingFrozenBandWitness"/> (<c>inspect --root ringband</c>).
///
/// <para><b>The chain of reasons, in four steps.</b> Each rung of the spin ring is a rung of a ring with
/// one fixed flux, s = (−1)^(p−1) (<see cref="SectorWrapSign"/>), because the wrap hop passes the other
/// p − 1 particles. Inside a fixed-flux ring one single-particle matrix serves every sector, so the η ladder
/// Φ and its adjoint Ψ commute with L, and the frozen space is an sl(2) module along the diagonal: the rung
/// p carries Σ_{ℓ ≤ min(p, N−p)} n_ℓ, n_ℓ the frozen modes of lowest weight on the rung ℓ. At odd N no two
/// allowed momenta add up to π (<see cref="NoResonance"/>), which turns the cell counts into the identity
/// N·D̂ = ℓ² − Dmom on LW_ℓ ∩ V₀ and the floor 𝒦 ≥ ℓ(N − ℓ)/N above 1 for ℓ ≥ 2 once N ≥ 5
/// (<see cref="CellFloorExceedsOne"/>), so n_ℓ = 0 above the bottom rung for generic parameters. What is left
/// is the corner, whose frozen modes at large J tend to one per pair of opposite momenta.</para>
///
/// <para><b>Typed parents, and what each gives.</b> <see cref="FrozenDivisorClaim"/> (F140) supplies the
/// FLOOR only. It is typed as a statement on the open chains; its proof uses only that h is real, symmetric
/// and invariant under the reversal l ↦ N − 1 − l, which both fixed-flux rings with uniform bonds are, so on
/// the ring it gives at least ⌊N/2⌋ at every J. Its tightness does not carry (it needs a nondegenerate
/// spectrum, and the ring's is doubly degenerate); the exact count on the ring is the proof's
/// Proposition 6, and the number ⌊N/2⌋ is borrowed from <see cref="FrozenDivisorClaim.FrozenMultiplicity"/>
/// as a number, not as that claim's tightness. <see cref="SpectatorIntertwinerClaim"/> (F125) supplies the
/// LADDER: its W is Φ, an exact intertwiner for any quadratic particle-conserving H, hence inside each
/// fixed-flux ring, and its reverse W† is the Ψ whose kernel is the lowest-weight space. The typed claim
/// says nothing about rings; the ring parity of the spectator, that W carries a rung of the spin ring to a
/// rung of the ring with the opposite flux, is PROOF_CODIM1_BY_ADDITIVITY §5 and the untyped gate
/// <c>SpectatorRingParityGateTests</c>, and the witness checks its adjoint on the diagonal. Neither parent
/// alone gives the ceiling: the first has the floor without the climb, the second the climb without the
/// count.</para>
///
/// <para><b>What is not typed here.</b> The proof also rests on PROOF_FROZEN_BAND_SO4's sum formula
/// (Corollary 4.2) and pencil count (Propositions 5.1 and 5.3), which carry no typed claim. The identity and
/// its floor need the commutant of H in the plane waves and stay in the proof and in the gate's R3/R4; the
/// even-N reads, the single-pair profiles at composite N and the defective corner pair at N = 5 stay reads
/// (PROOF_RING_FROZEN_BAND §8). The ring corner's commutant rates 2(N − 1)/N and 2(N − 2)/N, which
/// <c>StructuralCeilingClaim</c> carries, are derived in the proof's §5.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F170 + <c>docs/proofs/PROOF_RING_FROZEN_BAND.md</c> +
/// <c>simulations/ring_frozen_band.py</c> (the Python gate, 40 checks).</para></summary>
public sealed class RingFrozenBandClaim : Claim
{
    /// <summary>The corner floor, F140: typed on the open chains, and valid on both fixed-flux rings by its
    /// proof (h real, symmetric, reversal-invariant), an argument this claim carries.</summary>
    public FrozenDivisorClaim FrozenDivisor { get; }

    /// <summary>The η ladder W = Φ and its reverse W† = Ψ, F125.</summary>
    public SpectatorIntertwinerClaim Spectator { get; }

    /// <summary>The fixed-flux sign of the spin ring's rung p: (−1)^(p−1) (Lemma 1).</summary>
    public static int SectorWrapSign(int p) =>
        p >= 1 ? (p % 2 == 1 ? 1 : -1) : throw new ArgumentOutOfRangeException(nameof(p), p, "rungs start at p = 1");

    /// <summary>Lemma 2: no two allowed momenta of either fixed-flux ring add up to π (mod 2π) exactly when N
    /// is odd. In units of π/N two momenta sum to an even integer and π is N.</summary>
    public static bool NoResonance(int n) =>
        n >= 2 ? n % 2 == 1 : throw new ArgumentOutOfRangeException(nameof(n), n, "a ring needs N ≥ 2");

    /// <summary>Corollary 4's floor on LW_ℓ ∩ V₀ as an exact fraction: 𝒦 ≥ ℓ(N − ℓ)/N, at odd N, where the
    /// identity behind it holds.</summary>
    public static (int Numerator, int Denominator) CellFloor(int n, int l)
    {
        if (n < 3 || n % 2 == 0)
            throw new ArgumentOutOfRangeException(nameof(n), n, "the floor is an odd-N statement, N ≥ 3");
        if (l < 1 || l > n / 2)
            throw new ArgumentOutOfRangeException(nameof(l), l, $"lowest weights sit on rungs 1..⌊N/2⌋ = {n / 2}");
        return (l * (n - l), n);
    }

    /// <summary>Whether the floor ℓ(N − ℓ)/N exceeds 1 on every rung 2 ≤ ℓ ≤ ⌊N/2⌋, compared in integers:
    /// ℓ(N − ℓ) &gt; N. True at every odd N ≥ 5; at N = 3 there is no such rung, which is vacuously true.
    /// Integer arithmetic that cannot fail at odd N ≥ 5; the odd-N content is in the identity, not here.
    /// </summary>
    public static bool CellFloorExceedsOne(int n)
    {
        for (int l = 2; l <= n / 2; l++)
        {
            var (num, den) = CellFloor(n, l);
            if (num <= den) return false;
        }
        return true;
    }

    /// <summary>The count on every rung 1 ≤ p ≤ N − 1 at odd N: ⌊N/2⌋ generically on the locus (the number
    /// read off <see cref="FrozenDivisorClaim.FrozenMultiplicity"/>, the tightness from Proposition 6), and
    /// N − 1 at the uniform point for all but finitely many J.</summary>
    public static int RungCount(int n, bool uniform)
    {
        if (n < 3 || n % 2 == 0)
            throw new ArgumentOutOfRangeException(nameof(n), n, "F170 is an odd-N statement, N ≥ 3");
        return uniform ? n - 1 : FrozenDivisorClaim.FrozenMultiplicity(n);
    }

    public RingFrozenBandClaim(FrozenDivisorClaim frozenDivisor, SpectatorIntertwinerClaim spectator)
        : base(
            "F170: on the XY ring at odd N, uniform bonds, on the R₉₀ locus with γ̄ ≠ 0, every diagonal block " +
            "(p, p), 1 ≤ p ≤ N − 1, carries the frozen root −4γ̄ with geometric and algebraic multiplicity " +
            "exactly ⌊N/2⌋ for (J, δ) generic and N − 1 at the uniform point for all but finitely many J, the " +
            "wrap bond priced by occupation or by order. Each rung of the spin ring is the fixed-flux ring with " +
            "s = (−1)^(p−1); inside it the η ladder commutes with L; at odd N no two plane-wave momenta add up " +
            "to π, so N·D̂ = ℓ² − Dmom on LW_ℓ ∩ V₀ and 𝒦 ≥ ℓ(N − ℓ)/N > 1 above the seed: nothing new " +
            "is seeded above the bottom rung, and every rung carries what the corner carries, whose frozen modes " +
            "tend at large J to one per pair of opposite momenta. SCOPE: odd N; at γ̄ = 0 (F140's zero-mean " +
            "stratum) false; even N outside this claim (the proof's §7 derives that no rung above the corner carries frozen modes of its own at N ≥ 8, generically, so each carries its ring's corner, untyped; N = 4 overfull on its middle rung); the exceptional set holds " +
            "single-pair profiles at composite N, couplings at which the corner's root turns defective, and at " +
            "the uniform point couplings at which a higher rung gains lowest-weight frozen modes (at N = 5 the " +
            "rung 2 at J = (2/√3)·γ̄).",
            Tier.Tier1Derived,
            "docs/ANALYTICAL_FORMULAS.md F170 + docs/proofs/PROOF_RING_FROZEN_BAND.md + " +
            "simulations/ring_frozen_band.py + " +
            "compute/RCPsiSquared.Diagnostics/Foundation/RingFrozenBandWitness.cs (inspect --root ringband)")
    {
        FrozenDivisor = frozenDivisor ?? throw new ArgumentNullException(nameof(frozenDivisor));
        Spectator = spectator ?? throw new ArgumentNullException(nameof(spectator));
    }

    public override string DisplayName =>
        "F170: the XY ring's frozen band at odd N, every rung carrying what the corner carries";

    public override string Summary =>
        "at odd N no two plane-wave momenta add up to π, so nothing new is seeded above the bottom rung: every " +
        "diagonal rung carries ⌊N/2⌋ generically on the R₉₀ locus and N − 1 at the uniform point for all but " +
        "finitely many J, for both pricings of the wrap bond; live as kernel dimensions modulo primes meeting " +
        $"the floor, and as the nullity of the square where the block is small ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the four steps",
                summary: "sector identity s = (−1)^(p−1); the ladder commutes inside a fixed-flux ring, so the " +
                         "rung p carries Σ_{ℓ ≤ min(p, N−p)} n_ℓ; no resonance at odd N, hence " +
                         "N·D̂ = ℓ² − Dmom on LW_ℓ ∩ V₀ and n_ℓ = 0 for ℓ ≥ 2 generically; the corner, " +
                         "⌊N/2⌋ by F140's floor and a second-order lift, N − 1 at the uniform point");
            yield return new InspectableNode("the floor at N = 7",
                summary: $"ℓ(N − ℓ)/N = {string.Join(", ", Enumerable.Range(1, 3).Select(l => $"{CellFloor(7, l).Numerator}/7"))} " +
                         $"on the rungs 1, 2, 3; above 1 from ℓ = 2 on: {CellFloorExceedsOne(7)}");
            yield return new InspectableNode("what stays outside",
                summary: "even N (the proof's §7: from N = 8 on no rung above the corner carries its own, the corner exact at the uniform point and read N/2 on the locus to N = 20; " +
                         "N = 4 reads [2, 4, 2]; the periodic ring's uniform point at N = 6 carries lowest weight on rung 2), " +
                         "the single-pair profiles ⌊N/2⌋ + (gcd(2l + 1, N) − 1)/2 at composite N (the number of missed pairs " +
                         "is F157's blind Heisenberg count, the proof's §8, and the excess read matches it), the defective corner pair at N = 5, and the large-J corner " +
                         "span{P_k − P_{−k}}: reads and float checks in the gate, not typed here");
            yield return FrozenDivisor;
            yield return Spectator;
        }
    }
}
