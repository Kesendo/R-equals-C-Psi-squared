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
        // {Σcos a = Σcos b ≠ 0, sheet}: equality of the two cosine sums suffices. The Σ ≠ 0
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
    public void FlipLemma_OddProjectionAnnihilatesTheNumerator_ExactOverZ()
    {
        // PIECE 5 (F128 flip lemma, 2026-07-14 evening): the integer trig polynomial
        // cos s · B · V_a V_b P̃ (8640 monomials in half-angle units) is annihilated by the
        // (ℤ/2)⁶ signed character sum, exactly over ℤ; the projector self-test guards against
        // a projection that annihilates everything; the 28 Murnaghan-Nakayama-shifted exponent
        // sets of proof (B) all die structurally (0/repeat/±pair).
        var r = FlipLemma.Analyze();
        Assert.True(r.ProjectorSelfTestOk, "projector self-test failed");
        Assert.Equal(8640, r.NumeratorMonomials);
        Assert.Equal(0, r.SurvivingMonomials);
        Assert.Equal(28, r.ShiftedSets);
        Assert.True(r.AllShiftedSetsDie, "a shifted exponent set survives the death criteria");
        Assert.True(r.AllShiftedProjectionsVanish,
            "a shifted alternant's odd projection is not the zero polynomial");
    }

    [Fact]
    public void F128Factorization_HoldsAtGenericPoints_AndCorruptionControlBreaks()
    {
        // PIECE 6 (F128, the factorization): 𝔉 = −(e₁−f₁)²·𝒪[cos s·cot s·V_a V_b/P] at GENERIC
        // half-angle points in GF(p); LHS via the committed 𝔉 transcription (Evaluate), RHS via
        // the literal 64-flip sum: disjoint code paths. The corruption control replaces
        // (e₁−f₁)² by (e₁+f₁)² and must break at nearly every point.
        var (points, mismatches, ctrlMismatches, ctrlEvals) =
            CrossFormCertificate.CertifyF128FactorizationSlice(WallPrime, 12);
        Assert.True(points >= 12, $"expected ≥ 12 generic points, got {points}");
        Assert.Equal(0, mismatches);
        Assert.True(ctrlMismatches >= (int)(0.8 * ctrlEvals),
            $"corruption control: {ctrlMismatches}/{ctrlEvals} broke (expected ≥ 80%)");
    }

    [Fact]
    public void F128SharperLocus_EqualNonzeroCosineSums_NoSheet_CrossFormVanishes()
    {
        // PIECE 7 (F128 corollary, the 𝔉-scoped sharper locus): 𝔉 = 0 already on {e₁ = f₁ ≠ 0},
        // ONE constraint, no sheet: strictly wider than F127's V (both sums zero) and than the
        // T-scoped sharper slice (which needs the sheet). The control re-randomizes w₃
        // (breaking the equality) and must be nonzero.
        var (points, bad, ctrlNonzero, ctrlEval) =
            CrossFormCertificate.CertifyF128SharperLocusSlice(WallPrime, 12);
        Assert.True(points >= 12, $"expected ≥ 12 sharper-locus points, got {points}");
        Assert.Equal(0, bad);
        Assert.True(ctrlNonzero >= (int)(0.8 * ctrlEval),
            $"controls: {ctrlNonzero}/{ctrlEval} nonzero (expected ≥ 80%)");
    }

    [Fact]
    public void Witness_SummaryReportsF128()
    {
        var w = new CrossTripleOrthogonalityWitness();
        Assert.Contains("F128", w.Summary);
        Assert.Contains("flip lemma", w.Summary);
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

    [Fact]
    public void CyclotomicRing_BasicIdentitiesHoldExactly()
    {
        // PIECE 8 (the exact layer under F129/F130): φ(30) = 8; the order-3 root sum
        // 1 + ζ^10 + ζ^20 = 0 at m = 30; the zero mode ζ^10 + ζ^30 = i + (−i) = 0 at m = 40;
        // and ζ^m = 1 (exponent wrap).
        Assert.Equal(8, CyclotomicRing.Degree(30));
        var r3 = CyclotomicRing.Zero(30);
        CyclotomicRing.AddRootPower(r3, 30, 0, 1);
        CyclotomicRing.AddRootPower(r3, 30, 10, 1);
        CyclotomicRing.AddRootPower(r3, 30, 20, 1);
        Assert.True(CyclotomicRing.IsZero(r3), "1 + ζ^10 + ζ^20 must vanish at m = 30");
        var zm = CyclotomicRing.Zero(40);
        CyclotomicRing.AddRootPower(zm, 40, 10, 1);
        CyclotomicRing.AddRootPower(zm, 40, 30, 1);
        Assert.True(CyclotomicRing.IsZero(zm), "ζ^10 + ζ^30 must vanish at m = 40");
        var one = CyclotomicRing.Zero(30);
        CyclotomicRing.AddRootPower(one, 30, 30, 1);
        CyclotomicRing.AddRootPower(one, 30, 0, -1);
        Assert.True(CyclotomicRing.IsZero(one), "ζ^30 must equal 1 at m = 30");
    }

    [Fact]
    public void F129Census_NonFiringInjective_FiringCollide_LawEquality()
    {
        // PIECE 9 (F129, the exact census n ≤ 60): at every non-firing n all clean-triple
        // level VECTORS (ℤ[ζ_2n] free-basis form) are pairwise distinct — a PROOF of
        // injectivity at that n; at every firing n (3|n ≥ 9, 10|n ≥ 20) an exact collision
        // exists; the firing set equals the predicted set exactly.
        var r = LevelCollisionCensus.Analyze();
        Assert.Equal(LevelCollisionCensus.LiveMaxN, r.MaxN);
        Assert.True(r.AllNonFiringInjective, "a non-firing n has an exact collision");
        Assert.True(r.AllFiringCollide, "a firing n has no exact collision");
        Assert.True(r.LawEquality, "the firing set does not equal the predicted set");
        Assert.True(r.CleanTriplesChecked > 400_000,
            $"census unexpectedly small: {r.CleanTriplesChecked} clean triples");
    }

    [Fact]
    public void F129Census_MechanismAnchorsAreExact()
    {
        // PIECE 10 (F129 §4 anchors): n = 15 (8,12,14)~(9,11,13) collide at a NONZERO level and
        // the ± root set splits into four rotated R₃ cycles each summing to zero; n = 20
        // (1,7,9)~(3,5,10) is the R₅ conjugate pair plus the zero mode. All in ℤ[ζ_2n].
        var r = LevelCollisionCensus.Analyze(maxN: 5);
        Assert.True(r.AnchorsExact, "a mechanism anchor failed its exact decomposition");
    }

    [Fact]
    public void F129Census_CorruptedFiringPredicateBreaksTheLaw()
    {
        // Discrimination guard: a deliberately wrong firing predicate (n = 9 declared
        // non-firing) must break the law equality — proving the census actually measures
        // collisions rather than echoing the predicate.
        var r = LevelCollisionCensus.Analyze(maxN: 15,
            fireOverride: n => n != 9 && LevelCollisionCensus.Fires(n));
        Assert.False(r.LawEquality, "the corrupted predicate went unnoticed: the census is vacuous");
        Assert.False(r.AllNonFiringInjective);
    }

    [Fact]
    public void F130Decoupling_AllCollisionPairsVanishExactly()
    {
        // PIECE 11 (F130): at the named pairs covering every cell of the proof (disjoint
        // ε=−1 at nonzero level, the pentagon door, the non-clean pair, overlap-1 both signs,
        // the resonant special case, the level-free ε=+1 cell), both scaled half-Grams
        // Ê± vanish exactly in ℤ[ζ_2n]: U⁺ = U⁻ = 0, no floats.
        var r = CollisionDecoupling.Analyze();
        Assert.Equal(7, r.Collisions.Count);
        Assert.True(r.AllCollisionsDecouple,
            "a collision pair has a nonzero half-Gram: " +
            string.Join("; ", r.Collisions.Where(c => !c.PlusZero || !c.MinusZero).Select(c => c.Label)));
        Assert.True(r.LevelFlagsAsExpected, "an exact level-equality flag contradicts the pair table");
    }

    [Fact]
    public void F130Decoupling_ControlsNonzero_AndLemma3SignHolds()
    {
        // Discrimination guard + the unconditional pin: one unequal-level control per
        // level-sensitive cell must be NONZERO (unequal level does not imply B ≠ 0 in general,
        // so this is verified, not assumed), and Lemma 3's Ê⁺ = ε·Ê⁻ must hold exactly on
        // every pair including the controls.
        var r = CollisionDecoupling.Analyze();
        Assert.Equal(3, r.Controls.Count);
        Assert.True(r.AllControlsNonzero,
            "a control pair vanished: " +
            string.Join("; ", r.Controls.Where(c => c.PlusZero || c.MinusZero).Select(c => c.Label)));
        Assert.True(r.AllLemma3SignsOk, "the exact Lemma-3 sign identity failed on a pair");
    }

    [Fact]
    public void F130Decoupling_CorruptedPairIsNonzero()
    {
        // Corruption control: replace the n = 15 collision partner by a disjoint ε=−1 triple
        // at a DIFFERENT level (the bump must keep ε = −1: a bump into the disjoint ε=+1
        // cell would vanish level-free by Lemmas 3+4 and prove nothing) — the half-Grams
        // must NOT vanish, so the zero at the true pair is a measurement, not an artifact.
        var r = CollisionDecoupling.Evaluate("corrupted", 15,
            new[] { 8, 12, 14 }, new[] { 7, 9, 13 });
        Assert.False(r.LevelsEqual);
        Assert.False(r.PlusZero, "the corrupted pair's Ê⁺ vanished: the evaluator may be vacuous");
        Assert.False(r.MinusZero);
        Assert.True(r.Lemma3SignOk, "Lemma 3 must hold even on the corrupted pair");
    }

    [Fact]
    public void Witness_SummaryReportsF129AndF130()
    {
        var w = new CrossTripleOrthogonalityWitness();
        Assert.Contains("F129", w.Summary);
        Assert.Contains("F130", w.Summary);
        Assert.StartsWith("PASS", w.Summary);
    }

    [Fact]
    public void Claim_CarriesF129F130ConstantsAndAnchors()
    {
        Assert.Equal(60, CrossTripleOrthogonalityClaim.F129LiveCensusMaxN);
        Assert.Equal(7, CrossTripleOrthogonalityClaim.F130CollisionPairs);
        Assert.Equal(3, CrossTripleOrthogonalityClaim.F130ControlPairs);
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<CrossTripleOrthogonalityClaim>();
        Assert.Contains("f129_level_collision_law.py", claim.Anchor);
        Assert.Contains("f130_collision_decoupling.py", claim.Anchor);
    }

    [Fact]
    public void F129Inventory_ClosedFormsTieToExactCensus_RowForRow()
    {
        // PIECE 12 (the F129 family inventory, 2026-07-15): the thirteen derived closed
        // forms sum (total AND d-split) to the exact ℤ[ζ_2n] census at every row n ≤ 60,
        // firing and non-firing alike (the formulas vanish at the silent doors n = 6, 10
        // themselves), plus the n = 70 capstone (L's door: C 120 + L 20, all disjoint)
        // and the internal M-split sum.
        var r = CollisionFamilyInventory.Analyze();
        Assert.Equal(CollisionFamilyInventory.LiveMaxN, r.MaxN);
        Assert.True(r.AllRowsTied, "a closed-form sum disagrees with the exact census");
        Assert.True(r.Capstone70Tied, "the n = 70 capstone (L's door) broke");
        Assert.True(r.MSplitSumsToCountM, "40 + 0 + 60 does not sum to Count(M)");
        Assert.Equal((140L, 140L, 0L), CollisionFamilyInventory.CensusAt(70));
    }

    [Fact]
    public void F129Inventory_Capstone105_SevenFamiliesSum8858()
    {
        // PIECE 13 (M's door, m = 210): at n = 105 exactly seven families co-fire (doors
        // 3, 15, 21, 105: A 5457, D 1152, E 612, F 901, H 576, I 60, M 100) and the exact
        // census reproduces the sum 8858 with the d-split D + H = 1728 shared-mode pairs.
        Assert.Equal(5457, CollisionFamilyInventory.Count(CollisionFamilyInventory.Family.A, 105));
        Assert.Equal(100, CollisionFamilyInventory.Count(CollisionFamilyInventory.Family.M, 105));
        var census = CollisionFamilyInventory.CensusAt(105);
        Assert.Equal(CollisionFamilyInventory.ClosedFormsOf(105), census);
        Assert.Equal((8858L, 7130L, 1728L), census);
    }

    [Fact]
    public void F129Inventory_CorruptedClosedFormBreaksTheTie()
    {
        // Discrimination guard (two-sided): a deliberately wrong family D formula
        // (12(n−8) instead of 12(n−9)) shifts Total and Overlap1 at every 15|n ≤ 60 while
        // the census side is recomputed from the triples — the tie must break, proving the
        // tie actually measures the census rather than echoing the formulas.
        static long Corrupted(CollisionFamilyInventory.Family f, int n) =>
            f == CollisionFamilyInventory.Family.D && n % 15 == 0 && n >= 15
                ? 12L * (n - 8)
                : CollisionFamilyInventory.Count(f, n);
        var r = CollisionFamilyInventory.Analyze(countsOverride: Corrupted);
        Assert.False(r.AllRowsTied, "a corrupted D closed form failed to break the tie");
    }

    [Fact]
    public void F129Inventory_OnsetZerosAndSubLawDoors()
    {
        // The F129 thresholds are zeros of the formulas (A and B at the silent 3|6, C at
        // the silent 10|10 — A needs floor division, C# truncation would give A(6) = 1);
        // off-door every family is 0 (no door divides 25); and every d = 2 family's door
        // carries the factor 3 (the sub-law, visible piece by piece).
        Assert.Equal(0, CollisionFamilyInventory.Count(CollisionFamilyInventory.Family.A, 6));
        Assert.Equal(0, CollisionFamilyInventory.Count(CollisionFamilyInventory.Family.B, 6));
        Assert.Equal(0, CollisionFamilyInventory.Count(CollisionFamilyInventory.Family.C, 10));
        foreach (CollisionFamilyInventory.Family f in Enum.GetValues<CollisionFamilyInventory.Family>())
        {
            Assert.Equal(0, CollisionFamilyInventory.Count(f, 25));
            if (CollisionFamilyInventory.SharesAMode(f))
                Assert.Equal(0, CollisionFamilyInventory.Door(f) % 3);
        }
    }

    [Fact]
    public void CyclotomicRing_Phi210CarriesThePlusTwoCoefficient()
    {
        // The inventory capstones take the ring past its old m ≤ 120 envelope: Φ₁₄₀ still
        // has coefficients in {−1, 0, 1}, but Φ₂₁₀ (odd part 105 = 3·5·7) is the classic
        // first break — its x⁷ coefficient is +2 (cross-checked against sympy's
        // cyclotomic_poly(210)). Correctness rests on DivideExact + checked arithmetic,
        // never on the coefficient bound; this pin asserts the extended envelope knowingly.
        Assert.Equal(48, CyclotomicRing.Degree(140));
        Assert.Equal(48, CyclotomicRing.Degree(210));
        var phi140 = CyclotomicRing.Phi(140);
        Assert.All(phi140, c => Assert.InRange(c, -1, 1));
        var phi210 = CyclotomicRing.Phi(210);
        Assert.Equal(2, phi210[7]);
        Assert.Equal(2, phi210.Max());
        Assert.Equal(-1, phi210.Min());
    }

    [Fact]
    public void Claim_CarriesInventoryConstantsAndAnchors()
    {
        Assert.Equal(13, CrossTripleOrthogonalityClaim.F129InventoryFamilies);
        Assert.Equal(40, CrossTripleOrthogonalityClaim.F129MSplitFannedR5Pair);
        Assert.Equal(0, CrossTripleOrthogonalityClaim.F129MSplitMiddleImpossible);
        Assert.Equal(60, CrossTripleOrthogonalityClaim.F129MSplitFixedVertexR5);
        Assert.Equal(
            CollisionFamilyInventory.Count(CollisionFamilyInventory.Family.M, 105),
            CrossTripleOrthogonalityClaim.F129MSplitFannedR5Pair
            + CrossTripleOrthogonalityClaim.F129MSplitMiddleImpossible
            + CrossTripleOrthogonalityClaim.F129MSplitFixedVertexR5);
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<CrossTripleOrthogonalityClaim>();
        Assert.Contains("f129_family_inventory.py", claim.Anchor);
        Assert.Contains("docs/proofs/PROOF_F129_FAMILY_INVENTORY_COUNTS.md", claim.Anchor);
        Assert.Contains("experiments/F129_FAMILY_INVENTORY.md", claim.Anchor);
    }

    // ---- F133: the symplectic closed form of the F128 cofactor W ----

    [Fact]
    public void F133SinSLemma_OddProjectionAnnihilatesSinSTimesDelta_ExactOverZ()
    {
        // (F133 gate G2): the integer trig polynomial (2i sin s)·(2i)¹⁵ Δ̂ is annihilated by the
        // (ℤ/2)⁶ signed character sum, exactly over ℤ; the projector self-test guards against a
        // projection that annihilates everything (a generic monomial must expand to 64 ±1 terms).
        var r = WSymplecticClosedForm.AnalyzeSinSLemma();
        Assert.Equal(720, r.DeltaMonomials);          // (2i)¹⁵ Δ̂ = the A₅ Weyl denominator
        Assert.Equal(0, r.SurvivingMonomials);
        Assert.True(r.ProjectorSelfTestOk, "projector self-test failed");
    }

    [Fact]
    public void F133Readoff_HalvesHaveTheCommittedSizes()
    {
        // (F133 gate G5): the meet-in-the-middle split X = P₁·P₂ reproduces the committed monomial
        // counts (P₁ = Δ̂ · first 15 sheets, P₂ = the other 16); cancellation keeps both small.
        var h = WSymplecticClosedForm.BuildHalves();
        Assert.Equal(590016, h.P1Size);
        Assert.Equal(5817, h.P2Size);
    }

    [Fact]
    public void F133Readoff_FullWindowMatchesTable_WithAggregateChecksums()
    {
        // (F133 gate G5): n_raw(λ) = 2·n_λ across the FULL 557-λ even-degree window, exactly the
        // 143 committed terms nonzero, with the aggregate checksums (max|n_λ| = 8, Σ|n_λ| = 359,
        // n_() = −1). The committed table is embedded in the witness; a transcription slip fails here.
        var r = WSymplecticClosedForm.AnalyzeReadoff();
        Assert.Equal(590016, r.P1Size);
        Assert.Equal(5817, r.P2Size);
        Assert.Equal(557, r.WindowCount);
        Assert.Equal(0, r.Mismatches);
        Assert.True(r.AllMatchTable, "the read-off disagrees with the committed table");
        Assert.Equal(143, r.NonzeroCount);
        Assert.Equal(8, r.MaxAbsN);
        Assert.Equal(359, r.SumAbsN);
        Assert.Equal(-1, r.NAtEmpty);
    }

    [Fact]
    public void F133GfpCertificate_ClosedFormIdentityHolds_BothPrimes_CorruptionBreaks()
    {
        // (F133, the GF(p) certificate): C₆·SP == 2⁻⁵·(2i)⁻⁴⁶·Σ_λ n_λ·A_{λ+ρ} at random points over
        // BOTH 30-bit NTT primes, using the READ-OFF-DERIVED n_λ (a code path disjoint from the
        // embedded table and from the 64-flip LHS). The corruption control bumps one n_λ and must
        // break the identity at every point, else the check is vacuous.
        var derived = WSymplecticClosedForm.DeriveCoefficients(WSymplecticClosedForm.BuildHalves());
        Assert.Equal(143, derived.Count);
        foreach (long p in new long[] { 2013265921L, 1811939329L })
        {
            var g = WSymplecticClosedForm.CertifyGfpSlice(p, 6, derived);
            Assert.True(g.Points >= 6, $"expected >= 6 points at p = {p}, got {g.Points}");
            Assert.Equal(0, g.Mismatches);
            Assert.True(g.ControlChecks > 0, "no corruption controls were evaluated");
            Assert.Equal(g.ControlChecks, g.ControlMismatches);
        }
    }

    [Fact]
    public void Witness_SummaryReportsF133()
    {
        var w = new CrossTripleOrthogonalityWitness();
        Assert.Contains("F133", w.Summary);
        Assert.Contains("W closed form", w.Summary);
        Assert.StartsWith("PASS", w.Summary);
    }

    [Fact]
    public void Claim_CarriesF133ConstantsAndAnchors()
    {
        Assert.Equal(143, CrossTripleOrthogonalityClaim.F133CharacterTerms);
        Assert.Equal(8, CrossTripleOrthogonalityClaim.F133MaxCoefficient);
        var claim = KnowledgeRegistryFactory.BuildDefault().Get<CrossTripleOrthogonalityClaim>();
        Assert.Contains("docs/proofs/PROOF_F133_W_SYMPLECTIC_CLOSED_FORM.md", claim.Anchor);
        Assert.Contains("simulations/f133_w_closed_form.py", claim.Anchor);
    }
}
