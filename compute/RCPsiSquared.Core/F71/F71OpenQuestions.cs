using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.F71;

/// <summary>The substantive open theoretical items for the F71 family: directions where
/// the bond-mirror identity is expected to break, partially survive, or generalise. Anchored
/// at the F71 "Breaks for" clause in <c>docs/ANALYTICAL_FORMULAS.md</c>.
///
/// <para>Closed: "non-uniform J_b" was resolved 2026-05-20 by F100
/// (<see cref="C1QPeakMirrorJParity"/>, proof
/// <c>docs/proofs/PROOF_F100_C1_QPEAK_MIRROR_J_PARITY.md</c>). "non-uniform γ_i:
/// site-dependent dephasing" was resolved 2026-05-21 by F101
/// (<see cref="C1MirrorGammaParity"/>, proof
/// <c>docs/proofs/PROOF_F101_C1_MIRROR_GAMMA_PARITY.md</c>): the c₁ bond-mirror
/// deviation is exactly odd in the F71-anti-palindromic component of the per-site γ
/// profile, zero for palindromic γ, leading-order linear. Graceful breakdown, not a
/// hard violation; the symmetric component of γ still pairs bonds.</para></summary>
public static class F71OpenQuestions
{
    private const string Anchor = "docs/ANALYTICAL_FORMULAS.md F71 \"Breaks for\" clause";

    public static IReadOnlyList<OpenQuestion> Standard { get; } = new[]
    {
        new OpenQuestion(
            "asymmetric initial states: per-site purity reflection",
            "F71 requires ρ₀ R-symmetric in per-site purities (not in ρ itself). " +
            "Open: classify initial states for which the per-site purity reflection fails, and " +
            "characterise the residual mirror structure (if any) for those states.",
            "Sweep ρ₀ by Bloch-vector parameter; look for hidden conserved quantities under R.",
            Anchor),
        new OpenQuestion(
            "per-F71-orbit substructure: derive Q_peak(orbit, c, N)",
            "F71 pairs bonds but does not predict the value. The per-F71-orbit substructure observed " +
            "in F86 (e.g. c=2 N=6: central b=2 → 1.440, flanking b=1, b=3 → 1.648) needs an analytical " +
            "model. Likely depends on the multi-particle XY structure of the (n, n+1) block.",
            "OBC sine-mode algebra applied to the per-bond M_H matrix elements; deriving g_eff(b, c, N) " +
            "would close this jointly with F86 Item 1'.",
            Anchor),
    };
}
