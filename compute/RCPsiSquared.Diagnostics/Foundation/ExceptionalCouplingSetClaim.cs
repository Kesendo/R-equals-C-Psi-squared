using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Symmetry;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The exceptional couplings of F50 (<c>docs/proofs/PROOF_WEIGHT1_DEGENERACY.md</c>, "The count at
/// exceptional couplings"): the count d_real(Re = −2γ) = 2N (+ δ_G, the weight-1 commutant's dimension) holds at
/// every ratio γ/J outside a FINITE set E(N, G), and is exceeded on it. Live witness:
/// <see cref="ExceptionalCouplingWitness"/> (<c>inspect --root exceptional</c>, N ≤ 4), which recovers every block
/// polynomial exactly over ℤ[i], counts and isolates its positive roots by Sturm's theorem, certifies the pinned
/// minimal polynomials by exact division, measures the count at −2γ generically and at every point of E with the
/// nullity beside it, and reads the chain's gap on both sides of min E.
///
/// <para><b>The theorem.</b> L is block diagonal in the joint popcount (p, q) and its dissipator is −2γ times
/// the Hamming distance of the cell. By the Absorption Theorem Re λ = −2γ⟨n_XY⟩ on every eigenvector, so a
/// real eigenvalue −2γ needs ⟨n_XY⟩ = 1. Odd p − q: every distance is ≥ 1, so the mode is PURE weight 1, a
/// commutant mode by F50's Step 1, and a left eigenvector too (anti-Hermitian L_H, Hermitian dissipator), so
/// no Jordan chain grows on it. Even p ≠ q: every distance is ≥ 2, impossible. Diagonal (p, p): the distances
/// are even, never 1, so P_p(γ) = det(B_pp(γ) + 2γ·I) has leading coefficient Π(2 − 2·Hamming) ≠ 0 and finitely
/// many roots. E(N, G) is the set of positive roots over all p. The extra modes at a point of E mix XY-weight 0
/// with weight ≥ 2.</para>
///
/// <para><b>Scope and fences.</b> Of H only Hermiticity and U(1) are used, so the localisation and the
/// finiteness hold verbatim for XXZ and XY, where off E the count is that model's weight-1 commutant. E fences
/// REAL eigenvalues only; totals that include oscillating modes (F33's 14 over the −2γ and −4γ rungs together at
/// N = 3) are measured off E, not covered. The argument sorts the blocks by parity at every rung, and at k = 1 both kinds close;
/// at a rung k ≥ 2 the blocks of p − q ≡ k (mod 2) that hold cells at distance k are not fenced (their real modes
/// are a question about H's commutant, the weight-k commutant and DEGENERACY_PALINDROME's weight-mixing modes on the chain), while the
/// blocks of p − q ≢ k whose distances straddle k reach the rung only on a finite set E_k by the same leading coefficient
/// (THE_ONE_SQUARE §9, gate cube_moves_gate.py M7). The F1 palindrome carries the exceptions to the mirror rung N − 1,
/// E_{N−k} = E_k.</para>
///
/// <para><b>What is proved beside the theorem, and what is measured.</b> E is non-empty on every graph computed (N = 2: γ/J = 2, a
/// defective double root; N = 3 chain: √((√17−1)/2) and √3; K₃: √3; N = 4: 8 chain, 3 ring, 4 star, 1
/// complete; N = 5 chain: 13). Block (p, p) carries exactly C(N, p) − c_p − m_p points counted with multiplicity (Theorem D of the plane-crossing
/// section, the inertia identity n(γ) = #{r_j(γ) &lt; 1} on the Schur complement onto the populations): c_p the components of
/// the exclusion graph, m_p the Krein debt, the block's non-stationary modes inside the half-plane Re λ &gt; −2γ at the
/// Hamiltonian end, so that no mode returns below the plane on any graph and the 2γ regime exists iff every debt is zero.
/// The witness reads the debt EXACTLY from the Sturm count (<see cref="ExceptionalCouplingWitness.Debts"/>) with the
/// eigensolver's count beside it. On the chain the debt is zero and the count C(N, p) − 1: exactly through N = 5, N − 1
/// at p = 1 for every N (proved), in float through N = 6 by the crossings and through N = 16 at the Hamiltonian end from the
/// compression W_p = Σ n_l ∘ n_l of the eigenstate populations (second eigenvalue below p − ½, open for every N; W_p is one
/// Clebsch–Gordan-dressed reduced matrix per N whose spin-shell block is exact, heights 0, 2, …, 2·min(p, N − p), so the orbit
/// sum T_{2p} over the cells of Hamming 2p is in its kernel for p ≤ N/2 at every N with a simple sector spectrum, and at half
/// filling W_p is mirrored μ ↔ N/2 − μ, the open statement there the second-smallest eigenvalue above ½; in the full space
/// h ≥ 4/3 for every traceless SU(2)-invariant operator by isotropy, the sector projector paying the difference, and W_p is the
/// Δ-axis handover's rate matrix R = 4W_p − 4p·I, gap(R) = 2 being the same threshold along Δ), and proved
/// on the uniform XY chain at every N by F141's ladder, F143's rung and F144's floor (PROOF_WEIGHT1_DEGENERACY § The
/// count as a plane crossing); the debts (m₁, m₂, m₃) are (1, 3, 1) on the N = 4 ring, (1, 0, 1) on the star and (3, 2, 3) on K₄. The chain's smallest point is 1/Q*_gap(N) at N = 2 to 5 by two routes (the exact root
/// against the spectrum on both sides of it) and by the theorem: every mode slower than 2γ is a real, semisimple
/// mode of a diagonal block (the Krein bound, no complex pair sets the gap), so wherever the 2γ regime exists its
/// threshold is 1/min E: closed forms Q*_gap = 1/2, √((1+√17)/8), 1/x₀ with x₀ the root 0.7450215 of
/// 9x¹² + 132x¹⁰ + 68x⁸ − 1696x⁶ − 2240x⁴ + 1280x² + 256, and at N = 5 an algebraic number of degree 48. The
/// event is the handover <see cref="HandoverFloorClaim"/> types for the XY chain, and at N = 2 the point is the
/// coherence horizon's exceptional point of <see cref="CoherenceHorizonClaim"/> (carrier Q = 1). Units: the PAULI
/// book H = J·Σ(XX+YY+ZZ); the spin book's J is four times this.</para></summary>
public sealed class ExceptionalCouplingSetClaim : Claim
{
    /// <summary>Typed parent: the floor Re λ = −2γ⟨n_XY⟩ that makes a real −2γ mode a ⟨n_XY⟩ = 1 mode.</summary>
    public AbsorptionTheoremClaim Absorption { get; }

    /// <summary>Typed parent: the count 2N whose converse this claim fences.</summary>
    public F50WeightOneDegeneracyPi2Inheritance WeightOne { get; }

    /// <summary>Typed parent: the (N+1)² joint-popcount grading the localisation is stated in.</summary>
    public JointPopcountSectors Sectors { get; }

    public ExceptionalCouplingSetClaim(AbsorptionTheoremClaim absorption, F50WeightOneDegeneracyPi2Inheritance weightOne, JointPopcountSectors sectors)
        : base("F50 at exceptional couplings: for the isotropic Heisenberg H = J·Σ(XX+YY+ZZ) on a graph G under " +
               "uniform Z-dephasing, let E(N, G) be the set of positive γ/J at which some diagonal joint-popcount block " +
               "(p, p), 1 ≤ p ≤ N−1, has −2γ as an eigenvalue. E(N, G) is FINITE: P_p(γ) = det(B_pp(γ) + 2γ·I) has " +
               "leading coefficient Π(2 − 2·Hamming) ≠ 0 because a diagonal block holds only even distances. Off E, " +
               "every eigenvector and generalized eigenvector at a real −2γ is a weight-1 commutant mode (Absorption: " +
               "⟨n_XY⟩ = 1; odd p−q forces pure weight 1 and a left eigenvector; even p≠q is impossible), so " +
               "d_real(Re = −2γ) is the weight-1 commutant's dimension, 2N on every tested connected graph and 2N+2 on " +
               "K₃. On E the count is exceeded by modes mixing XY-weight 0 with weight ≥ 2. Only Hermiticity and U(1) " +
               "of H are used (XXZ, XY verbatim). FENCES: real eigenvalues only; rung k = 1 only (at k = 2 the leading " +
               "coefficient vanishes); the F1 mirror carries E to the rung N−1. PROVED beside it (the count as a plane " +
               "crossing): every block (p,p) of a connected graph holds at most C(N,p) − 1 points (the Schur complement onto " +
               "the populations, each eigencurve crossing 1 at most once); exactly N − 1 at p = 1 on the chain for every N; " +
               "the Krein bound keeps every complex and defective mode on or below the line, so a gap below 2γ is a real " +
               "mode of a diagonal block, and min E(chain) = 1/Q*_gap(N) is a theorem consequence wherever the 2γ regime " +
               "exists; on the uniform XY chain #E_p = C(N,p) − 1 at every p and every N (F141's ladder, F143's rung and F144's floor put every " +
               "non-stationary mode above the plane at γ → 0⁺). MEASURED: E non-empty on every graph computed; C(N,p) − 1 points " +
               "at p ≥ 2 on the chain through N = 6, and through N = 16 by m_p = 0 at the Hamiltonian end (no non-stationary mode " +
               "below the plane at γ → 0⁺, read from the compression W_p = Σ n_l ∘ n_l of the eigenstate populations); " +
               "the identification at N = 2..5 by two routes, the Heisenberg chain's handover (closed forms at " +
               "N = 2, 3, 4; degree 48 at N = 5). Pauli book: the spin book's J is four times this one.",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_WEIGHT1_DEGENERACY.md (The count at exceptional couplings) + " +
               "docs/proofs/PROOF_ABSORPTION_THEOREM.md + " +
               "docs/proofs/derivations/D06_SPECTRAL_GAP.md (Q*_gap, the chain's smallest point) + " +
               "simulations/f50_exceptional_couplings.py (gates G1 to G6, --n5 for the N = 5 chain) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/ExceptionalCouplingWitness.cs (inspect --root exceptional)")
    {
        Absorption = absorption ?? throw new ArgumentNullException(nameof(absorption));
        WeightOne = weightOne ?? throw new ArgumentNullException(nameof(weightOne));
        Sectors = sectors ?? throw new ArgumentNullException(nameof(sectors));
    }

    public override string DisplayName =>
        "F50 at exceptional couplings: the count 2N holds off a finite set E(N, G) of γ/J";

    public override string Summary =>
        "a real eigenvalue −2γ needs ⟨n_XY⟩ = 1 (Absorption); odd blocks force pure weight 1, even off-diagonal " +
        "blocks are impossible, and the diagonal blocks (p,p) hold −2γ only at the finitely many positive roots of " +
        "det(B_pp + 2γ), leading coefficient Π(2 − 2·Hamming) ≠ 0; off E the count is the weight-1 commutant's " +
        "(2N; K₃: 2N+2), on E it is exceeded; min E(chain) = 1/Q*_gap(N) measured at N = 2..5 and a consequence wherever the 2γ regime exists; live at " +
        $"inspect --root exceptional ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("where the exceptions can be",
                summary: "only the diagonal joint-popcount blocks (p,p), 1 ≤ p ≤ N−1: Absorption puts a real −2γ mode at " +
                         "⟨n_XY⟩ = 1, odd p−q forces pure weight 1 (F50's commutant, no Jordan chain), even p≠q holds " +
                         "distance ≥ 2 only");
            yield return new InspectableNode("why finitely many",
                summary: "P_p(γ) = det(B_pp(γ) + 2γ·I) is a polynomial whose leading coefficient Π(2 − 2·Hamming) is " +
                         "nonzero, since no cell of a diagonal block has distance 1; E(N, G) = its positive roots over p");
            yield return new InspectableNode("the live witness at N = 3 (chain)",
                summary: LiveN3Chain(),
                provenance: NodeProvenance.Live);
            yield return new InspectableNode("proved beside the theorem (the count as a plane crossing)",
                summary: "every block (p,p) carries exactly C(N,p) − c_p − m_p points counted with multiplicity, c_p the components of the exclusion graph and m_p the " +
                         "Krein debt, the block's modes inside the half-plane at the Hamiltonian end (Theorem D, the inertia identity " +
                         "n(γ) = #{r_j(γ) < 1} on the Schur complement onto the populations: Haynsworth, the constant inertia of the coherence " +
                         "corner, Theorem A); so no mode returns below the plane on any graph, and the 2γ regime exists iff every debt is zero, " +
                         "a count read exactly by the witness; nothing slower than 2γ is complex or defective (the transpose " +
                         "as a Krein form), so a gap below 2γ is a real mode of a diagonal block and wherever the 2γ regime exists its " +
                         "threshold is 1/min E; on the chain exactly N − 1 points at p = 1 for every N; the (1,1) ring sectors carry " +
                         "the coherence-horizon dispersion");
            yield return new InspectableNode("measured",
                summary: "E non-empty on every graph computed (N = 2: γ/J = 2 defective; N = 3 chain: 1.249621, 1.732051; " +
                         "K₃: √3 with count 12; N = 4: 8 chain, 3 ring, 4 star, 1 complete; N = 5 chain: 13); chain block " +
                         "(p,p) carries C(N,p) − 1 points at p ≥ 2: exactly through N = 5, in float through N = 6, and m_p = 0 at the Hamiltonian end through N = 16 (the second eigenvalue of W_p = Σ n_l ∘ n_l below p − ½); " +
                         "the debts (m₁, m₂, m₃) are (1, 3, 1) on the N = 4 ring, (1, 0, 1) on the star and (3, 2, 3) on K₄, so there the 2γ regime does not exist; min E(chain) = 1/Q*_gap(N) at N = 2..5 by two " +
                         "routes (the Heisenberg chain's handover, typed for XY by HandoverFloorClaim, and at " +
                         "N = 2 the coherence horizon's point of CoherenceHorizonClaim); totals with oscillating modes only " +
                         "measured off E; Pauli J throughout");
            yield return new InspectableNode("fences",
                summary: "real eigenvalues only; rung k = 1 only (at a rung k ≥ 2 the blocks of p − q ≡ k mod 2 with cells at distance k are not fenced, their real modes being a question about H's commutant, " +
                         "while the blocks of p − q ≢ k whose distances straddle k reach the rung only on a finite E_k, THE_ONE_SQUARE §9); uniform γ; the F1 " +
                         "mirror carries E to rung N−1; XXZ and XY verbatim with their own commutant off E");
            yield return Absorption;   // typed parent edge
            yield return WeightOne;    // typed parent edge
            yield return Sectors;      // typed parent edge
        }
    }

    private static string LiveN3Chain()
    {
        var w = new ExceptionalCouplingWitness(3, "chain");
        return $"E(3, chain) = {{{string.Join(", ", w.ExceptionalSet.Select(r => r.ToString("0.000000", System.Globalization.CultureInfo.InvariantCulture)))}}} " +
               $"(pinned {(w.PinnedCoversAllRoots ? "exactly" : "NOT certified")}); count at −2γ {w.GenericCount} generic, " +
               $"[{string.Join(", ", w.CountsAtPoints)}] at E (nullities [{string.Join(", ", w.NullitiesAtPoints)}]); J/γ at min E = {w.HandoverQ!.Value.ToString("0.000000", System.Globalization.CultureInfo.InvariantCulture)} " +
               $"against Q*_gap(3) = √((1+√17)/8) = {Math.Sqrt((1 + Math.Sqrt(17)) / 8).ToString("0.000000", System.Globalization.CultureInfo.InvariantCulture)}, " +
               $"the gap read on both sides: {(w.GapSidesHold ? "2γ below, a real mode below 2γ above (measured)" : "THE GAP READING FAILS")}";
    }

    /// <summary>Shared instance with its parents built from the Pi2-Foundation root, mirroring
    /// <see cref="PinnedBlockFloorClaim.Shared"/>.</summary>
    public static ExceptionalCouplingSetClaim Shared { get; } =
        new ExceptionalCouplingSetClaim(
            new AbsorptionTheoremClaim(new Pi2DyadicLadderClaim()),
            new F50WeightOneDegeneracyPi2Inheritance(new Pi2DyadicLadderClaim()),
            new JointPopcountSectors());
}
