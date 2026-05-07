namespace RCPsiSquared.Core.F86;

/// <summary>Shared anchor strings for the F86 family. Mirrors the
/// <see cref="F71.F71Anchors"/> pattern: one place to update when a proof reference moves,
/// instead of re-typing the path in every claim's ctor.</summary>
internal static class F86Anchors
{
    public const string ProofQPeak = "docs/proofs/PROOF_F86_QPEAK.md";
    public const string Statement2Plus3 = "docs/proofs/PROOF_F86_QPEAK.md Statement 2 + Statement 3";
}
