using RCPsiSquared.Core.Knowledge;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.Foundation;
using RCPsiSquared.Diagnostics.Knowledge;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>F127 (CrossTripleOrthogonalityClaim / CrossTripleOrthogonalityWitness):
/// (a) the live witness passes (𝔉-slice zeros, core-identity T zeros, and the sheet lattice
/// combinatorics, each with an off-variety / structural control);
/// (b) registry wiring: the claim resolves from the default registry, Tier1Candidate (the
/// code-trust caveat is the named reason), with its typed parent SeedExistenceCountingClaim;
/// (c) the claim's constants match the F127 registry entry;
/// (d) the residue-collapse extension: the sheet lattice (§2) and the core identity T (§3).</summary>
public class CrossTripleOrthogonalityTests
{
    // the first 30-bit wall prime (core_grid.gen_primes(17)[0]; ≡ 1 mod 4, so i = √−1 exists)
    private const long WallPrime = 1073741833L;

    [Fact]
    public void Witness_SlicePasses()
    {
        var w = new CrossTripleOrthogonalityWitness();
        Assert.StartsWith("PASS", w.Summary);
        Assert.Contains("controls nonzero", w.Summary);
    }

    [Fact]
    public void Witness_SummaryReportsAllThreeChecks()
    {
        var w = new CrossTripleOrthogonalityWitness();
        Assert.Contains("core-T", w.Summary);
        Assert.Contains("lattice", w.Summary);
    }

    [Fact]
    public void SheetLattice_HasTheExpectedResidueCombinatorics()
    {
        // PIECE 1 (§2): exact integer combinatorics of the 72 atoms → 288 pole events → 32 sheets,
        // every coefficient in {±1}, exactly 9 events per sheet, one per (i,j) block.
        var report = SheetLattice.Analyze();
        Assert.Equal(72, report.Atoms);
        Assert.Equal(288, report.Events);
        Assert.True(report.AllCoefficientsPmOne);
        Assert.Equal(32, report.DistinctSheets);
        Assert.Equal(9, report.MinEventsPerSheet);
        Assert.Equal(9, report.MaxEventsPerSheet);
        Assert.True(report.EachSheetOneEventPerBlock);
    }

    [Fact]
    public void CoreIdentity_VanishesOnVariety_ControlNonzero()
    {
        // PIECE 2 (§3): T ≡ 0 on the three-constraint variety {Qz, Qw, product = 1}, live in GF(p);
        // the off-variety control (re-randomized w₁) must be nonzero at nearly every point.
        var (onVariety, bad, ctrlNonzero, ctrlEval) =
            CrossFormCertificate.CertifyCoreIdentitySlice(WallPrime, 20);
        Assert.True(onVariety >= 20, $"expected ≥ 20 on-variety points, got {onVariety}");
        Assert.Equal(0, bad);
        Assert.True(ctrlNonzero >= (int)(0.8 * ctrlEval),
            $"controls: {ctrlNonzero}/{ctrlEval} nonzero (expected ≥ 80%)");
    }

    [Fact]
    public void ClosedForm_HoldsAtGenericPoints_AndCorruptionControlBreaks()
    {
        // PIECE 3 (§3 closed form, 2026-07-14): T·P = ⅛[2cos s((e₁−f₁)² − 2sin²s)
        // + sin s(Σsin2a + Σsin2b)]·Va·Vb is an UNCONDITIONAL identity, so it must hold at
        // GENERIC half-angle points in GF(p) (no variety construction); the corruption control
        // (the sin s term sign-flipped) must break at nearly every point, else the check is vacuous.
        var (points, mismatches, ctrlMismatches, ctrlEvals) =
            CrossFormCertificate.CertifyClosedFormSlice(WallPrime, 24);
        Assert.True(points >= 24, $"expected ≥ 24 generic points, got {points}");
        Assert.Equal(0, mismatches);
        Assert.True(ctrlMismatches >= (int)(0.8 * ctrlEvals),
            $"corruption control: {ctrlMismatches}/{ctrlEvals} broke (expected ≥ 80%)");
    }

    [Fact]
    public void SharperLocus_EqualNonzeroCosineSums_OnSheet_TVanishes()
    {
        // PIECE 4 (Corollary 2 of the closed form, NEW surface): T = 0 already on
        // {Σcos a = Σcos b ≠ 0, sheet} — equality of the two cosine sums suffices. The Σ ≠ 0
        // requirement makes the slice discriminating against the old variety (where both are 0);
        // the control re-randomizes w₁ (breaking sheet and equality) and must be nonzero.
        var (points, bad, ctrlNonzero, ctrlEval) =
            CrossFormCertificate.CertifySharperLocusSlice(WallPrime, 20);
        Assert.True(points >= 20, $"expected ≥ 20 sharper-locus points, got {points}");
        Assert.Equal(0, bad);
        Assert.True(ctrlNonzero >= (int)(0.8 * ctrlEval),
            $"controls: {ctrlNonzero}/{ctrlEval} nonzero (expected ≥ 80%)");
    }

    [Fact]
    public void Witness_SummaryReportsClosedFormAndSharperLocus()
    {
        var w = new CrossTripleOrthogonalityWitness();
        Assert.Contains("closed form", w.Summary);
        Assert.Contains("sharper", w.Summary);
    }

    [Fact]
    public void Claim_ResolvesFromDefaultRegistry_WithTypedParent()
    {
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<CrossTripleOrthogonalityClaim>();
        Assert.NotNull(claim);
        Assert.Equal(Tier.Tier1Candidate, claim.Tier);
        Assert.NotNull(claim.SeedExistence);
        Assert.Contains("F127", claim.DisplayName);
    }

    [Fact]
    public void Claim_ConstantsMatchTheRegistryEntry()
    {
        Assert.Equal(527, CrossTripleOrthogonalityClaim.CertifiedTasks);
        Assert.Equal(17, CrossTripleOrthogonalityClaim.PrimeCount);
        Assert.Equal(510.0, CrossTripleOrthogonalityClaim.Log2PrimeProduct);
        Assert.True(CrossTripleOrthogonalityClaim.AssemblyDIsSymbolic);
    }

    [Fact]
    public void Claim_AnchorNamesTheWallAndTheAssembly()
    {
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<CrossTripleOrthogonalityClaim>();
        Assert.Contains("grid_proof_sweep.py", claim.Anchor);
        Assert.Contains("assembly_d_symbolic.py", claim.Anchor);
        Assert.Contains("docs/ANALYTICAL_FORMULAS.md", claim.Anchor);
    }
}
