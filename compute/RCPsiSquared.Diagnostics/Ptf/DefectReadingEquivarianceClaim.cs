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

/// <summary>M3 of the reading-grammar arc: the defect-reading map is equivariant under the
/// geometric spatial reflection R (i → N−1−i, the chain mirror, a ℤ₂ with R² = I). The location
/// dictionary M[b,k] = ⟨ψ_k|V_b|ψ_1⟩ (carrier ψ_1, modes k = 2..N, V_b = ½(X_bX_{b+1}+Y_bY_{b+1}))
/// satisfies the exact identity M[N−2−b, k] = (−1)^{k−1} M[b, k]: reflecting the bond reflects the
/// reading, mode by mode.
///
/// <para>Two-line derivation (this claim's own contribution): R V_b R = V_{N−2−b}; the carrier is
/// R-even, Rψ_1 = +ψ_1; mode k has reflection parity Rψ_k = (−1)^{k−1}ψ_k. Hence
/// M[N−2−b,k] = ⟨ψ_k|V_{N−2−b}|ψ_1⟩ = ⟨Rψ_k|V_b|Rψ_1⟩ = (−1)^{k−1}⟨ψ_k|V_b|ψ_1⟩. The
/// sign-location confusability follows in closed form as a parity-weighted mode sum:
/// cos(b, N−2−b) = Σ_{k=2..N} (−1)^{k−1} w_k, w_k = M[b,k]²/‖M[b,·]‖²; negative (anti-collinear)
/// because the R-odd channels (above all the seesaw k=2) carry the net location weight.</para>
///
/// <para>Two structures on one dictionary: the within-feature stabilizer is the K-partner null
/// (k=N forbidden column, <see cref="KPartnerSelectionRuleClaim"/>, the load-bearing parent that
/// defines M[b,k] and rank N−2); the cross-feature group action is this spatial-reflection
/// equivariance. The R here is NOT the coherence-space mirror group's R: the
/// <see cref="RCPsiSquared.Core.Symmetry.MirrorGroupD4Claim"/>'s R is the ket-flip I⊗X^⊗N (which
/// does not even preserve the single-excitation sector), and the geometric spatial mirror is
/// deliberately outside that D₄ (PROOF_PI_FACTORS_AS_R_TIMES_D §5). Sibling mirrors in two spaces.</para>
///
/// <para>Carrier-parity note: the (−1)^{k−1} sign is specific to the R-even carrier ψ_1; a carrier
/// of parity (−1)^{c−1} gives (−1)^{k−c}. The canonical carrier is ψ_1. Honest scope: the bare
/// cosine is −0.33 at N=5, NOT −0.97; the −0.97 is the painted (propagated α-profile) instance
/// (handshake_gram_metric.py, Q=20), where the readout concentrates weight on the R-odd seesaw.</para>
///
/// <para>Anchor: hypotheses/HANDSHAKE_GEOMETRY.md + simulations/handshake_reading_equivariance.py
/// + compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs (the decoder whose sign-location
/// ambiguity this explains structurally) + compute/RCPsiSquared.Diagnostics/Ptf/KPartnerSelectionRuleClaim.cs
/// (the parent: the dictionary, the null column, rank N−2).</para></summary>
public sealed class DefectReadingEquivarianceClaim : Claim
{
    public readonly record struct BatteryCase(string Name, string Detail, string Expected, string Actual)
    {
        public bool Passes => string.Equals(Expected, Actual, StringComparison.Ordinal);
    }

    private readonly KPartnerSelectionRuleClaim _kPartner;

    public IReadOnlyList<BatteryCase> Cases { get; }
    public int PassCount => Cases.Count(c => c.Passes);

    public DefectReadingEquivarianceClaim(KPartnerSelectionRuleClaim kPartner)
        : base("Defect-reading spatial-reflection equivariance: the location dictionary " +
               "M[b,k] = ⟨ψ_k|V_b|ψ_1⟩ satisfies M[N−2−b,k] = (−1)^{k−1}M[b,k] exactly (the reading " +
               "commutes with the geometric chain mirror R: i→N−1−i), so the sign-location confusability " +
               "is the closed-form parity-weighted mode sum cos(b,N−2−b) = Σ(−1)^{k−1}w_k; the same " +
               "dictionary's within-feature stabilizer is the K-partner null (KPartnerSelectionRuleClaim)",
               Tier.Tier1Derived,
               "hypotheses/HANDSHAKE_GEOMETRY.md + " +
               "simulations/handshake_reading_equivariance.py (R-parity + equivariance + cosine, machine-exact N = 4,5,6) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs (the sign-location ambiguity this explains) + " +
               "compute/RCPsiSquared.Diagnostics/Ptf/KPartnerSelectionRuleClaim.cs (the parent: dictionary, null column, rank N−2)")
    {
        _kPartner = kPartner ?? throw new ArgumentNullException(nameof(kPartner));
        Cases = BuildBattery();
    }

    public string Equivariance =>
        "The geometric spatial reflection R (i → N−1−i) sends bond b to bond N−2−b (R V_b R = V_{N−2−b}); " +
        "the carrier is R-even (Rψ_1 = +ψ_1) and mode k has reflection parity Rψ_k = (−1)^{k−1}ψ_k, so the " +
        "location dictionary obeys M[N−2−b,k] = (−1)^{k−1}M[b,k] exactly, for all N.";

    public string Confusability =>
        "The mirror-image bonds read with cosine cos(b,N−2−b) = Σ_{k=2..N} (−1)^{k−1} w_k (w_k = M[b,k]²/‖M[b,·]‖²), " +
        "a closed-form parity-weighted mode sum. It is negative (anti-collinear) because the R-odd channels, above " +
        "all the seesaw k=2, carry the net location weight: this IS the DefectDecoder's sign-location ambiguity.";

    public string Derivation =>
        "M[N−2−b,k] = ⟨ψ_k|V_{N−2−b}|ψ_1⟩ = ⟨ψ_k|R V_b R|ψ_1⟩ = ⟨Rψ_k|V_b|Rψ_1⟩ = (−1)^{k−1}⟨ψ_k|V_b|ψ_1⟩. " +
        "Only R (a ℤ₂, R² = I) is used; this is R-equivariance, not the full coherence-space D₄.";

    public string Source =>
        "Load-bearing parent: KPartnerSelectionRuleClaim (Tier1Derived) defines M[b,k], the K-partner null k=N " +
        "column, and rank N−2 (the within-feature stabilizer). See-also (NOT a typed parent): MirrorGroupD4Claim — " +
        "the coherence-space mirror group whose R is the ket-flip I⊗X^⊗N, which deliberately EXCLUDES this geometric " +
        "spatial mirror (PROOF_PI_FACTORS_AS_R_TIMES_D §5). Two distinct mirrors in two spaces.";

    public override string DisplayName =>
        "Defect-reading spatial-reflection equivariance (M[N−2−b,k] = (−1)^{k−1}M[b,k], Tier1Derived)";

    public override string Summary =>
        "The defect-reading map's spatial-reflection equivariance: the reading commutes with the geometric " +
        "chain mirror R, M[N−2−b,k] = (−1)^{k−1}M[b,k] exactly, so the sign-location confusability is the " +
        "closed-form parity-weighted mode sum " +
        $"cos = Σ(−1)^{{k−1}}w_k; the dictionary's K-partner null is the within-feature stabilizer; {PassCount}/{Cases.Count} PASS ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("Equivariance (M[N−2−b,k] = (−1)^{k−1}M[b,k], the reading reflects)", summary: Equivariance);
            yield return new InspectableNode("Confusability (cos = Σ(−1)^{k−1}w_k, anti-collinear via the R-odd seesaw)", summary: Confusability);
            yield return new InspectableNode("Derivation (two lines; only the ℤ₂ reflection R, not D₄)", summary: Derivation);
            yield return new InspectableNode("See-also: MirrorGroupD4Claim (sibling mirror in coherence space, NOT this geometric R)",
                summary: Source);
            foreach (var c in Cases)
                yield return new InspectableNode(c.Name,
                    summary: $"{c.Detail}; expected {c.Expected}, got {c.Actual}, " + (c.Passes ? "PASS" : "FAIL"));
            yield return _kPartner; // the single typed parent edge
        }
    }

    /// <summary>Pure single-excitation algebra at N = 4, 5, 6 (no Liouvillian, no propagation):
    /// (1) R-parity Rψ_k = (−1)^{k−1}ψ_k; (2) equivariance M[N−2−b,k] = (−1)^{k−1}M[b,k];
    /// (3) cosine = Σ(−1)^{k−1}w_k; (4) anti-collinear mirror-pair cos &lt; 0.</summary>
    private static IReadOnlyList<BatteryCase> BuildBattery()
    {
        const double tol = 1e-9;
        string Zero(double v) => v < tol ? "0" : v.ToString("E3", CultureInfo.InvariantCulture);
        static int Sign(int k) => ((k - 1) % 2 == 0) ? 1 : -1;

        var cases = new List<BatteryCase>();
        foreach (int N in new[] { 4, 5, 6 })
        {
            ComplexMatrix BondTerm(int l) =>
                0.5 * (PauliString.SiteOp(N, l, PauliLetter.X) * PauliString.SiteOp(N, l + 1, PauliLetter.X)
                     + PauliString.SiteOp(N, l, PauliLetter.Y) * PauliString.SiteOp(N, l + 1, PauliLetter.Y));
            double Element(ComplexVector bra, ComplexMatrix Vb, ComplexVector ket) =>
                bra.Conjugate().DotProduct(Vb * ket).Real;

            var psi1 = BondingMode.Build(N, 1);
            int bonds = N - 1;       // b = 0..N−2
            int locModes = N - 1;    // k = 2..N

            // location dictionary M[b][kIdx] = ⟨ψ_k|V_b|ψ_1⟩
            var M = new double[bonds][];
            for (int b = 0; b < bonds; b++)
            {
                var Vb = BondTerm(b);
                M[b] = new double[locModes];
                for (int kIdx = 0; kIdx < locModes; kIdx++)
                    M[b][kIdx] = Element(BondingMode.Build(N, kIdx + 2), Vb, psi1);
            }

            // (1) R-parity: BondingMode is big-endian (site j → index 1<<(N−1−j)); the SE amplitude
            // a_j = ψ_k[1<<(N−1−j)] must satisfy a_{N−1−j} = (−1)^{k−1} a_j.
            double worstRparity = 0.0;
            for (int k = 1; k <= N; k++)
            {
                var psik = BondingMode.Build(N, k);
                for (int j = 0; j < N; j++)
                {
                    double aj = psik[1 << (N - 1 - j)].Real;
                    double aMirror = psik[1 << j].Real;          // amplitude at the reflected site N−1−j
                    worstRparity = Math.Max(worstRparity, Math.Abs(aMirror - Sign(k) * aj));
                }
            }
            cases.Add(new BatteryCase(
                $"R-parity Rψ_k = (−1)^(k−1)ψ_k (N={N})",
                "single-excitation amplitudes reflect with sign (−1)^(k−1) under the chain mirror i→N−1−i",
                "0", Zero(worstRparity)));

            // (2) equivariance M[N−2−b,k] = (−1)^{k−1} M[b,k]
            double worstEq = 0.0;
            for (int b = 0; b < bonds; b++)
                for (int kIdx = 0; kIdx < locModes; kIdx++)
                    worstEq = Math.Max(worstEq, Math.Abs(M[bonds - 1 - b][kIdx] - Sign(kIdx + 2) * M[b][kIdx]));
            cases.Add(new BatteryCase(
                $"equivariance M[N−2−b,k] = (−1)^(k−1)M[b,k] (N={N})",
                "the location dictionary commutes with the spatial reflection, mode by mode",
                "0", Zero(worstEq)));

            // (3) cosine = Σ(−1)^{k−1} w_k, and (4) cos < 0, over each unordered mirror pair.
            double worstCos = 0.0;
            bool allNeg = true;
            int pairCount = 0;
            double sampleCos = 0.0;
            for (int b = 0; b < bonds; b++)
            {
                int bm = bonds - 1 - b;
                if (bm <= b) continue;                            // each pair once; skip the self-mirror bond
                double dot = 0, nb = 0, nbm = 0;
                for (int kIdx = 0; kIdx < locModes; kIdx++)
                {
                    dot += M[b][kIdx] * M[bm][kIdx];
                    nb += M[b][kIdx] * M[b][kIdx];
                    nbm += M[bm][kIdx] * M[bm][kIdx];
                }
                double cosMeasured = dot / Math.Sqrt(nb * nbm);
                double pred = 0;
                for (int kIdx = 0; kIdx < locModes; kIdx++)
                    pred += Sign(kIdx + 2) * (M[b][kIdx] * M[b][kIdx]) / nb;
                worstCos = Math.Max(worstCos, Math.Abs(cosMeasured - pred));
                allNeg &= cosMeasured < 0;
                sampleCos = cosMeasured;
                pairCount++;
            }
            cases.Add(new BatteryCase(
                $"cosine = Σ(−1)^(k−1)w_k (N={N})",
                $"measured mirror-pair cosine equals the parity-weighted mode sum, {pairCount} pair(s)",
                "0", Zero(worstCos)));
            cases.Add(new BatteryCase(
                $"anti-collinear mirror-pair cos < 0 (N={N})",
                $"the R-odd seesaw carries the location; sample cos = {sampleCos.ToString("F3", CultureInfo.InvariantCulture)}",
                "True", allNeg.ToString()));
        }
        return cases;
    }
}
