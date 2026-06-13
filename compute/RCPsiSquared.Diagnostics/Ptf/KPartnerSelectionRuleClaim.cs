using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.States;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>The K-partner selection rule (the reading-grammar arc's first DERIVED result):
/// the carrier never leaks into its K-partner under any bond defect. The single-excitation
/// band-edge mode (the carrier) ψ_1 and the top mode ψ_N are K-partners, ψ_N = K₁ψ_1 with
/// K₁ = Π_{l odd} Z_l (equivalently ψ_N(i) = (−1)^i ψ_1(i)); and for EVERY bond
/// V_b = ½(X_bX_{b+1} + Y_bY_{b+1}) the matrix element ⟨ψ_N|V_b|ψ_1⟩ = 0. The carrier never
/// couples to its K-partner through any bond defect.
///
/// <para>Two-line derivation, both ingredients borrowed from <see cref="ChiralMirrorTrajectoryClaim"/>:
/// ⟨ψ_N|V_b|ψ_1⟩ = ⟨K₁ψ_1|V_b|ψ_1⟩ = ⟨ψ_1|K₁V_b|ψ_1⟩ = −⟨ψ_1|V_bK₁|ψ_1⟩ = −⟨ψ_1|V_b|ψ_N⟩
/// = −⟨ψ_N|V_b|ψ_1⟩ (the modes are real, V_b symmetric), so the real number equals its own
/// negative ⟹ 0. The two facts used (K₁ψ_1 = ψ_N from the parent's Step 4; K₁V_bK₁ = −V_b,
/// V_b K-odd, from the parent's Step 1) ARE the parent claim.</para>
///
/// <para>Consequence: rank(location dictionary) = N−2. The carrier couples to N−1 other modes
/// (k = 2..N) under bond defects, but the K-partner channel k = N is selection-rule-forbidden,
/// so the location channels span only N−2 dimensions (machine-verified rank = N−2, N = 3..8).
/// The <see cref="Foundation.DefectDecoder"/>'s sign-location ambiguity (an edge bond reads
/// like the complementary interior bond, residual ratio ≈ 1.5 at N=5) IS this: the K-partner
/// pairs are the dictionary's null direction.</para>
///
/// <para>Anchor: <c>hypotheses/HANDSHAKE_GEOMETRY.md</c> + <c>simulations/_k_partner_selection_rule.py</c>
/// + <c>compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs</c> +
/// <c>compute/RCPsiSquared.Diagnostics/Ptf/ChiralMirrorTrajectoryClaim.cs</c> (the parent).</para></summary>
public sealed class KPartnerSelectionRuleClaim : Claim
{
    /// <summary>One pure-algebra self-check on small dense single-excitation matrices.</summary>
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    private readonly ChiralMirrorTrajectoryClaim _chiralMirror;

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public KPartnerSelectionRuleClaim(ChiralMirrorTrajectoryClaim chiralMirror)
        : base("K-partner selection rule: ⟨ψ_N|V_b|ψ_1⟩ = 0 for every bond defect V_b (the carrier never " +
               "couples to its K-partner ψ_N = K₁ψ_1), a two-line corollary of ChiralMirrorTrajectoryClaim; " +
               "the location dictionary has rank N−2 and the DefectDecoder's sign-location ambiguity IS this null direction",
               Tier.Tier1Derived,
               "hypotheses/HANDSHAKE_GEOMETRY.md + " +
               "simulations/_k_partner_selection_rule.py (the selection-rule + rank N−2 probe, machine-exact N = 3..8) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs (the decoder whose rank-(N−2) ambiguity this explains) + " +
               "compute/RCPsiSquared.Diagnostics/Ptf/ChiralMirrorTrajectoryClaim.cs (the parent: both ingredients)")
    {
        _chiralMirror = chiralMirror ?? throw new ArgumentNullException(nameof(chiralMirror));
        Cases = BuildBattery();
    }

    public string Partnership =>
        "The carrier ψ_1 (the single-excitation band-edge mode) and the top mode ψ_N are K-partners: " +
        "ψ_N = K₁ψ_1 with K₁ = Π_{l odd} Z_l (the odd-sublattice Z product), equivalently ψ_N(i) = (−1)^i ψ_1(i). " +
        "The sine identity makes the map exact and sign-free.";

    public string SelectionRule =>
        "For every bond b the hopping/defect V_b = ½(X_bX_{b+1} + Y_bY_{b+1}) has ⟨ψ_N|V_b|ψ_1⟩ = 0: " +
        "the carrier never couples to its K-partner through any bond defect, exactly and for all N.";

    public string Derivation =>
        "Two lines, both ingredients from ChiralMirrorTrajectoryClaim: ⟨ψ_N|V_b|ψ_1⟩ = ⟨K₁ψ_1|V_b|ψ_1⟩ = " +
        "⟨ψ_1|K₁V_b|ψ_1⟩ = −⟨ψ_1|V_bK₁|ψ_1⟩ (V_b is K₁-odd: K₁V_bK₁ = −V_b, parent Step 1) = −⟨ψ_1|V_b|ψ_N⟩ " +
        "= −⟨ψ_N|V_b|ψ_1⟩ (ψ_N = K₁ψ_1, parent Step 4; modes real, V_b symmetric). A real number equal to its " +
        "own negative is 0.";

    public string RankConsequence =>
        "Build the location dictionary M[b, k] = ⟨ψ_k|V_b|ψ_1⟩ over bonds b = 0..N−2 and modes k = 2..N (k = 1 " +
        "is the strength channel). The selection rule kills the k = N column entirely, so M has rank N−2 " +
        "(machine-verified N = 3..8): the carrier couples to N−1 modes under bond defects, but the K-partner " +
        "channel is forbidden, leaving only N−2 independent location channels.";

    public string DecoderMeaning =>
        "The DefectDecoder's sign-location ambiguity (an edge bond reads almost like the complementary interior " +
        "bond, residual ratio ≈ 1.5 at N=5, AmbiguityFactor flagged) IS the rank deficit: the K-partner pairs " +
        "are the dictionary's null direction. The decoder cannot cleanly separate sign from location precisely " +
        "where the forbidden channel would have lifted the degeneracy.";

    public string Source =>
        "Sibling/source: ChiralMirrorTrajectoryClaim (Tier1Derived). This claim is its first downstream corollary " +
        "in the reading-grammar arc: K₁ψ_1 = ψ_N (parent Step 4) and K₁V_bK₁ = −V_b (parent Step 1) are exactly " +
        "the two facts the two-line derivation consumes. ChiralKClaim is the grandparent (eigenvalue side).";

    public override string DisplayName =>
        "K-partner selection rule (⟨ψ_N|V_b|ψ_1⟩ = 0, location dictionary rank N−2, Tier1Derived)";

    public override string Summary =>
        "The carrier ψ_1 never couples to its K-partner ψ_N = K₁ψ_1 through any bond defect: the selection rule " +
        "⟨ψ_N|V_b|ψ_1⟩ = 0 (a two-line corollary of ChiralMirrorTrajectoryClaim) makes the location dictionary " +
        $"rank N−2; the DefectDecoder's sign-location ambiguity is this K-partner null direction; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("K-partnership (ψ_N = K₁ψ_1, the carrier and its mirror)", summary: Partnership);
            yield return new InspectableNode("Selection rule (⟨ψ_N|V_b|ψ_1⟩ = 0 for every bond)", summary: SelectionRule);
            yield return new InspectableNode("Derivation (two lines, both ingredients from the parent)", summary: Derivation);
            yield return new InspectableNode("Rank consequence (location dictionary rank N−2)", summary: RankConsequence);
            yield return new InspectableNode("Decoder meaning (the sign-location ambiguity IS the null direction)", summary: DecoderMeaning);
            yield return new InspectableNode("Source: ChiralMirrorTrajectoryClaim (the two facts consumed)", summary: Source);
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return _chiralMirror; // typed parent edge
        }
    }

    /// <summary>Pure-algebra battery on dense single-excitation matrices at N = 4 AND N = 5 (no
    /// Liouvillian, no propagation): (a) ψ_N = K₁ψ_1 exactly; (b) the selection rule
    /// max_b |⟨ψ_N|V_b|ψ_1⟩| = 0; (c) the location dictionary M[b, k] = ⟨ψ_k|V_b|ψ_1⟩
    /// (b = 0..N−2, k = 2..N) has rank N−2.</summary>
    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        const double tol = 1e-9;
        string Zero(double v) => v < tol ? "0" : v.ToString("E3", CultureInfo.InvariantCulture);

        var cases = new List<BatteryCase>();
        foreach (int N in new[] { 4, 5 })
        {
            ComplexMatrix BondTerm(int l) =>
                0.5 * (PauliString.SiteOp(N, l, PauliLetter.X) * PauliString.SiteOp(N, l + 1, PauliLetter.X)
                     + PauliString.SiteOp(N, l, PauliLetter.Y) * PauliString.SiteOp(N, l + 1, PauliLetter.Y));

            // K₁ = Π_{l odd} Z_l (0-indexed odd sublattice).
            ComplexMatrix K1 = PauliString.SiteOp(N, 1, PauliLetter.Z);
            for (int l = 3; l < N; l += 2)
                K1 *= PauliString.SiteOp(N, l, PauliLetter.Z);

            var psi1 = BondingMode.Build(N, 1);
            var psiN = BondingMode.Build(N, N);

            // (a) ψ_N = K₁ψ_1 exactly.
            cases.Add(new BatteryCase(
                Name: $"ψ_N = K₁ψ_1 exactly (N={N})",
                Detail: $"‖K₁·ψ_1 − ψ_N‖₂ at N={N} (odd-sublattice Z product)",
                Expected: "0",
                Actual: Zero((K1 * psi1 - psiN).L2Norm())));

            // ⟨ψ_a|V_b|ψ_c⟩ = ψ_a* · (V_b · ψ_c) (modes real).
            double Element(ComplexVector bra, ComplexMatrix Vb, ComplexVector ket) =>
                bra.Conjugate().DotProduct(Vb * ket).Real;

            // (b) selection rule: max_b |⟨ψ_N|V_b|ψ_1⟩| = 0.
            double worstSel = 0.0;
            for (int b = 0; b < N - 1; b++)
                worstSel = Math.Max(worstSel, Math.Abs(Element(psiN, BondTerm(b), psi1)));
            cases.Add(new BatteryCase(
                Name: $"selection rule max_b |⟨ψ_N|V_b|ψ_1⟩| = 0 (N={N})",
                Detail: $"the carrier never couples to its K-partner through any of the {N - 1} bond defects",
                Expected: "0",
                Actual: Zero(worstSel)));

            // (c) location dictionary M[b, k] = ⟨ψ_k|V_b|ψ_1⟩, b = 0..N−2, k = 2..N; rank = N−2.
            int bonds = N - 1;          // b = 0..N−2
            int locModes = N - 1;       // k = 2..N
            var M = ComplexMatrix.Build.Dense(bonds, locModes);
            for (int b = 0; b < bonds; b++)
            {
                var Vb = BondTerm(b);
                for (int kIdx = 0; kIdx < locModes; kIdx++)
                {
                    int k = kIdx + 2;
                    M[b, kIdx] = new Complex(Element(BondingMode.Build(N, k), Vb, psi1), 0);
                }
            }
            int rank = M.Svd().S.Count(s => s.Real > tol);
            cases.Add(new BatteryCase(
                Name: $"rank(location dictionary) = N−2 (N={N})",
                Detail: $"M[b, k] = ⟨ψ_k|V_b|ψ_1⟩, b = 0..{N - 2}, k = 2..{N}; K-partner channel k={N} forbidden",
                Expected: (N - 2).ToString(CultureInfo.InvariantCulture),
                Actual: rank.ToString(CultureInfo.InvariantCulture)));
        }

        return cases;
    }
}
