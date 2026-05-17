using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F98 (2026-05-17 evening): KIntermediate Dicke long-time Π²-odd asymptote
/// bridges the F86b 3/8 static anchor (this morning, X⊗N-eigenbasis decomposition)
/// to the <see cref="QuarterAsBilinearMaxvalClaim"/> 1/4 universal boundary via an
/// explicit N-dependent closed form:
///
/// <code>
///     α(t = 0)_KIntermediate  = 3/8                       (F86b, static)
///     α(t → ∞)_KIntermediate(N) = (N + 2) / [4·(N + 1)]   (F98, dynamic asymptote)
///     α(t → ∞, N → ∞)         → 1/4                       (QuarterAsBilinearMaxval)
/// </code>
///
/// The closed form holds for any truly-class (F87) Hamiltonian + uniform Z-dephasing
/// on N qubits — bond topology drops out because the long-time limit projects onto
/// <c>ker L = span(P_0, …, P_N)</c> per F4 for any connected graph (chain, ring, star,
/// K_N, Petersen, etc.). The KIntermediate Dicke superposition
/// <c>ψ = (|D_{N/2−1}⟩ + |D_{N/2}⟩)/√2</c> evolves under any such L from the static F86b
/// anchor toward the universal Quarter, traversing the (N+2)/[4(N+1)] curve.
///
/// <para><b>Two ingredient closed forms (F98a, F98b in ANALYTICAL_FORMULAS.md):</b></para>
/// <list type="bullet">
///   <item><b>F98a</b>: <c>‖P_{N/2−1}_odd‖² = C(N, N/2−1) / 2</c>. The Π²-odd
///         Frobenius² of the sub-mid popcount-sector projector is exactly HALF its
///         rank. Derivable via direct Krawtchouk-style enumeration over Pauli-string
///         supports; verified bit-exact N = 4..16.</item>
///   <item><b>F98b</b>: <c>α(∞)_KIntermediate(N even) = (N + 2) / [4·(N + 1)]</c>.
///         Follows from F98a + mid-popcount Krawtchouk parity vanishing
///         (<c>‖P_{N/2}_odd‖² = 0</c>) + the kernel projection
///         <c>ρ_∞ = (1/2)/C(N, m)·P_m + (1/2)/C(N, m+1)·P_{m+1}</c> + Pascal's
///         <c>C(N, m) + C(N, m+1) = C(N+1, m+1)</c>.</item>
/// </list>
///
/// <para><b>Three typed parents</b> (all already in Pi2KnowledgeBase):
/// <see cref="QuarterAsBilinearMaxvalClaim"/> (the 1/4 asymptote), <see cref="HalfAsStructuralFixedPointClaim"/>
/// (1/4 = (1/2)² — the half-rank-of-P_{N/2−1} ingredient squared compounds with the
/// kernel-balance halving to give the Quarter), <see cref="DickeSuperpositionQuarterPi2Inheritance"/>
/// (the static-side companion: Theorem 1's C_block = 1/4 ceiling for the same Dicke
/// superposition family — F98 is the dynamic-side twin where the F86b 3/8 anchor
/// approaches the same 1/4 universal under L evolution).</para>
///
/// <para><b>Bridge to F86b:</b> the morning's <see cref="DickeAnchor.KIntermediate"/>
/// case <c>α_total = 3/8</c> was derived statically (Tier 1, X⊗N-eigenbasis decomposition
/// of the symmetric Dicke superposition at <c>γ = 1/2</c>). F98 adds the temporal half
/// of the same picture — the long-time fate of α_total under the canonical truly-class
/// dynamics. Both 3/8 and 1/4 sit on the polarity-squared algebra
/// (<c>3/8 = (1/2)·(3/4) = (1/2)(1 − 1/4)</c>; <c>1/4 = (1/2)²</c>); the (N+2)/[4(N+1)]
/// curve is the explicit N-trajectory between them.</para>
///
/// <para><b>Tier outcome: Tier 1 derived.</b> F98a is a Krawtchouk-style identity
/// verified bit-exact N = 4..16; F98b follows algebraically from F98a + Pascal +
/// Krawtchouk parity. Discovered via the water/proton-chain inheritance test
/// (<c>simulations/water/proton_chain_dicke_anchor.py</c>) which asked "what is α_total
/// at t = ∞ under truly-class evolution" — a NEW question not addressed by the static
/// F86b derivation.</para>
///
/// <para>Anchors: <see cref="DickeAnchor"/> (static F86b 3/8 ingredient classifier);
/// <see cref="DickeSuperpositionQuarterPi2Inheritance"/> (static-side Quarter ceiling,
/// PROOF_BLOCK_CPSI_QUARTER Theorem 1 + 2); <see cref="QuarterAsBilinearMaxvalClaim"/>
/// (the 1/4 asymptote target); <see cref="HalfAsStructuralFixedPointClaim"/>
/// (1/4 = (1/2)² parent). Doc: <c>docs/ANALYTICAL_FORMULAS.md</c> F98 entry,
/// <c>docs/water/README.md</c> § "Findings since May 4". Script:
/// <c>simulations/water/proton_chain_dicke_anchor.py</c>.</para>
/// </summary>
public sealed class KIntermediateAsymptoteQuarterInheritance : Claim, IF99AnchorBearing
{
    /// <inheritdoc />
    /// <remarks>F98 is a <see cref="F99AnchorRole.Direct"/> claim about the
    /// F99 α=3/8 K-intermediate anchor (γ=1/2, uniform Dicke). The static
    /// face is 3/8; the long-time asymptote (N+2)/[4(N+1)] → 1/4 is the
    /// dynamic bridge BETWEEN the F99 α=3/8 and F99 α=1/4 anchors. F98's
    /// content is on the F86b α-axis directly.</remarks>
    public F99AnchorRole F99Role => F99AnchorRole.Direct;

    /// <inheritdoc />
    public IReadOnlyList<double> F99AnchorValues { get; } = new[] { 3.0 / 8.0 };

    /// <summary>The 1/4 asymptote target. Typed parent.</summary>
    public QuarterAsBilinearMaxvalClaim Quarter { get; }

    /// <summary>The (1/2)² = 1/4 polarity-squared algebra. Typed parent.</summary>
    public HalfAsStructuralFixedPointClaim Half { get; }

    /// <summary>The static-side companion: Dicke saturates C_block = 1/4 at t = 0 per
    /// Theorem 1 + 2. F98 is the dynamic twin where α_total of the SAME Dicke
    /// superposition family approaches 1/4 under L evolution.</summary>
    public DickeSuperpositionQuarterPi2Inheritance StaticSide { get; }

    /// <summary>F98b closed form: <c>α(∞)_KIntermediate(N even) = (N + 2) / [4·(N + 1)]</c>.
    /// Bond-topology-independent under truly-class (F87) Hamiltonian + uniform Z-dephasing.</summary>
    public static double LongTimePi2OddRatio(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 2.");
        if ((N & 1) != 0) throw new ArgumentException(
            $"N must be even (KIntermediate Dicke anchor requires N/2 ∈ ℤ); got N = {N}.", nameof(N));
        return (N + 2) / (4.0 * (N + 1));
    }

    /// <summary>F98a closed form: <c>‖P_{N/2−1}_odd‖² = C(N, N/2−1) / 2</c>. The
    /// Π²-odd Frobenius² of the sub-mid popcount-sector projector is exactly HALF its
    /// rank.</summary>
    public static double SubMidProjectorPi2OddFrobeniusSquared(int N)
    {
        if (N < 2) throw new ArgumentOutOfRangeException(nameof(N), N, "N must be ≥ 2.");
        if ((N & 1) != 0) throw new ArgumentException(
            $"N must be even (sub-mid sector requires N/2 ∈ ℤ); got N = {N}.", nameof(N));
        return Binomial(N, N / 2 - 1) / 2.0;
    }

    /// <summary>Asymptote drift: how far <see cref="LongTimePi2OddRatio"/>(N) is from
    /// 1/4. Equals <c>1 / [4·(N + 1)]</c>; monotonically vanishing in N. Exposes the
    /// 1/4 limit explicitly.</summary>
    public static double DriftFromQuarter(int N) => LongTimePi2OddRatio(N) - 0.25;

    private static long Binomial(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        if (k > n - k) k = n - k;
        long c = 1;
        for (int i = 0; i < k; i++) c = c * (n - i) / (i + 1);
        return c;
    }

    public KIntermediateAsymptoteQuarterInheritance(
        QuarterAsBilinearMaxvalClaim quarter,
        HalfAsStructuralFixedPointClaim half,
        DickeSuperpositionQuarterPi2Inheritance staticSide)
        : base("F98: KIntermediate Dicke long-time Π²-odd asymptote (N+2)/[4(N+1)] → 1/4 = QuarterAsBilinearMaxval",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F98 + " +
               "simulations/water/proton_chain_dicke_anchor.py (bit-exact verification N=4..16) + " +
               "compute/RCPsiSquared.Core/Symmetry/DickeAnchor.cs (F86b 3/8 static anchor partner) + " +
               "compute/RCPsiSquared.Core/Symmetry/DickeSuperpositionQuarterPi2Inheritance.cs (static-side C_block Quarter ceiling)")
    {
        Quarter = quarter ?? throw new ArgumentNullException(nameof(quarter));
        Half = half ?? throw new ArgumentNullException(nameof(half));
        StaticSide = staticSide ?? throw new ArgumentNullException(nameof(staticSide));
    }

    public override string DisplayName =>
        "F98: KIntermediate α(∞) = (N+2)/[4(N+1)] → 1/4 long-time bridge from F86b 3/8 to QuarterAsBilinearMaxval";

    public override string Summary =>
        $"α(t=0) = 3/8 (F86b static, X⊗N-eigenbasis); " +
        $"α(t→∞)(N) = (N+2)/[4(N+1)]: {LongTimePi2OddRatio(4):G4} (N=4), " +
        $"{LongTimePi2OddRatio(10):G4} (N=10), → 1/4 (N→∞) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return StaticSide;
            yield return Quarter;
            yield return Half;
            yield return new InspectableNode("F98a",
                summary: "‖P_{N/2-1}_odd‖² = C(N, N/2-1) / 2 — Π²-odd Frobenius² of sub-mid popcount-sector projector is half its rank. Bit-exact N=4..16.");
            yield return new InspectableNode("F98b",
                summary: "α(∞)_KIntermediate(N) = (N+2)/[4(N+1)]. Follows from F98a + Krawtchouk parity vanishing at P_{N/2} + Pascal. Asymptote 1/4 = QuarterAsBilinearMaxval.");
            yield return new InspectableNode("F86b partner",
                summary: "α(t=0) = 3/8 = (1/2)·(3/4) — the X⊗N-eigenbasis γ = 1/2 Dicke superposition starts at 3/8 by the (1 − γ²)/2 F86b formula (DickeAnchor.KIntermediate, this morning's commit b9ba5f6).");
            yield return new InspectableNode("polarity-squared algebra",
                summary: "3/8 and 1/4 both sit on the dyadic ladder: 3/8 = (1/2)(1 − 1/4), 1/4 = (1/2)². The (N+2)/[4(N+1)] curve is the explicit N-trajectory between them under truly-class dynamics.");
            for (int N = 4; N <= 16; N += 2)
            {
                yield return new InspectableNode($"N = {N}",
                    summary: $"α(∞) = {LongTimePi2OddRatio(N):G6}, drift from 1/4 = {DriftFromQuarter(N):G3}, ‖P_{N / 2 - 1}_odd‖² = {SubMidProjectorPi2OddFrobeniusSquared(N):G6} (= C({N}, {N / 2 - 1})/2)");
            }
        }
    }
}
