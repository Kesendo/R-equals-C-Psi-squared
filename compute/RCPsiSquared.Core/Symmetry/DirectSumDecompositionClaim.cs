using RCPsiSquared.Core.F1;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The direct-sum decomposition of the Liouvillian
/// (<c>docs/proofs/DIRECT_SUM_DECOMPOSITION.md</c>, corollary of F1 + F61):
///
/// <code>
///   L = L_even ⊕ L_odd            by n_XY parity (bit_a)
///
///   1. dim(V_even) = dim(V_odd) = 2^(2N−1)         (half of 4^N each)
///   2. N odd:  Π exchanges sectors, and
///              L_odd = −Π · L_even · Π⁻¹ − 2Σγ · I  (mirror-image dynamics)
///   3. N even: Π preserves sectors; each is independently self-palindromic
///   4. Superselection charge P_XY = (−1)^{n_XY} = Π²_X with [P_XY, L] = 0
/// </code>
///
/// <para>The N-parity mechanism: Π maps XY-weight w to N − w, and
/// (N − w) − w = N − 2w has the parity of N. Odd N flips the sector bit,
/// even N preserves it. The two-sector split is the coarsest layer of the
/// sector structure: the mod-2 shadow of the U(1) joint-popcount grading,
/// refined further by bit_b (F63) into the Klein-V₄ four-sector split.</para>
///
/// <para>Selective breaking (the two standard perturbations break OPPOSITE
/// halves; <c>simulations/direct_sum_scope_probe.py</c>): a transverse field
/// (odd-n_XY Hamiltonian term, one-sided in the commutator) destroys the
/// direct sum and the superselection charge while the palindrome relation
/// survives; amplitude damping (bilinear sandwich, σ∓ bit_a-homogeneous)
/// preserves the direct sum EXACTLY while breaking the palindrome content
/// of statement 2/3 (the sectors persist but stop being mirror images).
/// The combinatorial parts (dimension count, sector-exchange arithmetic)
/// are properties of the Pauli basis and never break.</para>
///
/// <para>Structural parallel (independent development): the direct-sum
/// quantum theory of Gaztañaga, Kumar, Marto (Class. Quantum Grav. 43,
/// 015023, 2026): two sectors related by a discrete transformation carrying
/// opposite arrows of time, with a superselection rule.</para>
///
/// <para>Tier1Derived: corollary of the two Tier1Derived parents; no new
/// numerical input (the dimension count and parity arithmetic are exact).
/// Live witness: <c>inspect --root directsum</c>.</para></summary>
public sealed class DirectSumDecompositionClaim : Claim
{
    /// <summary>The palindrome Π·L·Π⁻¹ = −L − 2Σγ·I this claim restricts to sectors.</summary>
    public F1PalindromeIdentity F1 { get; }

    /// <summary>The bit_a (n_XY parity) grading that defines V_even/V_odd and the charge P_XY = Π²_X.</summary>
    public F61BitAParityPi2Inheritance F61 { get; }

    public DirectSumDecompositionClaim(F1PalindromeIdentity f1, F61BitAParityPi2Inheritance f61)
        : base("Direct-sum decomposition L = L_even ⊕ L_odd with odd-N Π sector exchange",
               Tier.Tier1Derived,
               "docs/proofs/DIRECT_SUM_DECOMPOSITION.md + " +
               "docs/ANALYTICAL_FORMULAS.md F61 + " +
               "simulations/direct_sum_scope_probe.py")
    {
        F1 = f1 ?? throw new ArgumentNullException(nameof(f1));
        F61 = f61 ?? throw new ArgumentNullException(nameof(f61));
    }

    /// <summary>Sector dimension: dim(V_even) = dim(V_odd) = 2^(2N−1), exactly half of 4^N.
    /// From dim(V_k) = C(N,k)·2^N and Σ_{k even} C(N,k) = 2^(N−1).</summary>
    public long SectorDimension(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "The direct sum requires N ≥ 1.");
        if (N > 31) throw new ArgumentOutOfRangeException(nameof(N), N, "2^(2N−1) overflows long for N > 31.");
        return 1L << (2 * N - 1);
    }

    /// <summary>True iff Π exchanges the sectors (V_even ↔ V_odd): exactly the odd-N case,
    /// because Π maps w → N − w and N − 2w has the parity of N.</summary>
    public bool PiExchangesSectors(int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "The direct sum requires N ≥ 1.");
        return N % 2 == 1;
    }

    /// <summary>The superselection charge: P_XY = (−1)^{n_XY} = Π²_X (conjugation by Z⊗N),
    /// F61's bit_a parity operator. [P_XY, L] = 0 exactly.</summary>
    public string SuperselectionCharge => "P_XY = (−1)^{n_XY} = Π²_X (Z⊗N conjugation)";

    public override string DisplayName =>
        "Direct sum L = L_even ⊕ L_odd: equal halves 2^(2N−1); odd N mirror-exchanged, even N self-palindromic";

    public override string Summary =>
        $"L = L_even ⊕ L_odd by n_XY parity; dim 2^(2N−1) each; odd N: L_odd = −Π·L_even·Π⁻¹ − 2Σγ·I (sector exchange); " +
        $"even N: each sector self-palindromic; charge P_XY = Π²_X; T1 preserves the sum and breaks the mirror, " +
        $"a transverse field the reverse ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("statement 1: equal dimensions",
                summary: "dim(V_even) = dim(V_odd) = 2^(2N−1); from dim(V_k) = C(N,k)·2^N and the half-binomial identity; never breaks (pure Pauli combinatorics)");
            yield return new InspectableNode("statements 2/3: N-parity dichotomy",
                summary: "Π maps w → N−w; N−2w has N's parity, so odd N exchanges the sectors (two time-reversed copies joined by Π), even N preserves them (each self-palindromic)");
            yield return new InspectableNode("statement 4: superselection",
                summary: $"[P_XY, L] = 0 with {SuperselectionCharge}; the two sectors are dynamically disconnected");
            yield return new InspectableNode("selective breaking",
                summary: "transverse field: kills the sum, keeps the palindrome; amplitude damping: keeps the sum EXACTLY, kills the palindrome (direct_sum_scope_probe.py); the wall and the mirror are separate prisoners");
            yield return new InspectableNode("finer structure",
                summary: "V_even/V_odd = mod-2 coarsening of the U(1) joint-popcount grading; bit_b (F63) refines to the Klein-V₄ four-sector split 4^(N−1) each");
            for (int n = 2; n <= 6; n++)
            {
                yield return new InspectableNode(
                    $"N={n}",
                    summary: $"per-sector dim = {SectorDimension(n)} (of {1L << (2 * n)} total); Π {(PiExchangesSectors(n) ? "EXCHANGES V_even ↔ V_odd" : "preserves each sector")}");
            }
        }
    }
}
