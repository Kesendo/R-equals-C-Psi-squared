using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Numerics;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>The live F127 witness (<c>inspect --root crosstriple</c>): recomputes, at inspect time,
/// an INDEPENDENT C# slice of the cross-triple orthogonality along the residue-collapse chain of
/// <c>docs/proofs/PROOF_F127_RESIDUE_COLLAPSE.md</c> and the F128 factorization of
/// <c>docs/proofs/PROOF_F128_FLIP_SUM_FACTORIZATION.md</c>. Five recomputations, each with its own
/// discriminating control:
/// <list type="bullet">
/// <item>the full form 𝔉 vanishes on the two-constraint variety over split primes, with an
/// off-variety nonzero control (<see cref="CrossFormCertificate.CertifySlice"/>);</item>
/// <item>the §3 core function T vanishes on the three-constraint variety, with an off-variety
/// nonzero control (<see cref="CrossFormCertificate.CertifyCoreIdentitySlice"/>);</item>
/// <item>the §2 sheet lattice: 72 atoms → 288 pole events → 32 sheets, all coefficients ±1, nine
/// events per sheet, exact integer combinatorics (<see cref="SheetLattice"/>);</item>
/// <item>the §3 closed form (2026-07-14): T·P = ⅛[2cos s((e₁−f₁)² − 2sin²s) + sin s(Σsin2a +
/// Σsin2b)]·V_a·V_b holds at generic GF(p) points with a corruption control, and Corollary 2's
/// SHARPER locus (Σcos a = Σcos b ≠ 0 on the sheet) kills T live
/// (<see cref="CrossFormCertificate.CertifyClosedFormSlice"/>,
/// <see cref="CrossFormCertificate.CertifySharperLocusSlice"/>);</item>
/// <item>F128 (2026-07-14 evening, PROOF_F128_FLIP_SUM_FACTORIZATION): the flip lemma exactly
/// over ℤ (8640 monomials annihilated by the signed character sum, <see cref="FlipLemma"/>),
/// the factorization 𝔉 = −(e₁−f₁)²·𝒪[cos s·cot s·V_a V_b/P] at generic GF(p) points (LHS/RHS
/// disjoint code paths, corruption control), and the 𝔉-scoped sharper locus {Σcos a = Σcos b
/// ≠ 0}, ONE constraint, no sheet (<see cref="CrossFormCertificate.CertifyF128FactorizationSlice"/>,
/// <see cref="CrossFormCertificate.CertifyF128SharperLocusSlice"/>);</item>
/// <item>F129 (2026-07-14 night, PROOF_F129_LEVEL_COLLISION_LAW): the exact ℤ[ζ_2n] level
/// census to n ≤ 60 (non-firing n injective, firing n collide, law equality, mechanism
/// anchors term-exact; a corrupted-predicate discrimination guard lives in the tests;
/// <see cref="LevelCollisionCensus"/>);</item>
/// <item>F130 (2026-07-14 night, PROOF_F130_COLLISION_DECOUPLING): Ê± = 0 exactly in ℤ[ζ_2n]
/// at named collision pairs covering every proof cell, unequal-level nonzero controls, and
/// the exact Lemma-3 sign Ê⁺ = ε·Ê⁻ on every pair (<see cref="CollisionDecoupling"/>);</item>
/// <item>the F129 family inventory (2026-07-15, PROOF_F129_FAMILY_INVENTORY_COUNTS): the
/// thirteen derived closed-form counts tie their SUMS (total + d-split) to the exact census
/// at every n ≤ 60 (eleven families exercised live) plus the n = 70 capstone (L's door);
/// M and its 40+0+60 split stay with the committed gate's I5
/// (<see cref="CollisionFamilyInventory"/>).</item>
/// </list>
/// This is the second implementation the F127/F128 code-trust caveat asks for, at sample scale; the
/// proof objects remain the 527/527 grid+CRT wall and the committed <c>f12*_*</c> gates (F129's
/// n ≤ 210 census stays with <c>simulations/f129_level_collision_law.py</c>; the named corner
/// was closed as empty 2026-07-15, PROOF_F129 §4).
/// A couple of seconds (the exact F129 census dominates; computed once per process); fully
/// deterministic (xorshift streams seeded from the prime; the exact layers have no randomness).</summary>
public sealed class CrossTripleOrthogonalityWitness : IInspectable
{
    // the wall's own first two 30-bit split primes (core_grid.gen_primes(17)[0..1]; both = 1 mod 4)
    private static readonly long[] Primes = { 1073741833L, 1073741857L };
    private const int SamplesPerPrime = 6;
    private const int CoreSamplesPerPrime = 12;
    private const int ClosedFormSamplesPerPrime = 12;
    private const int SharperSamplesPerPrime = 10;
    private const int F128SamplesPerPrime = 6;
    private const int F128SharperSamplesPerPrime = 10;

    // F133 (the symplectic closed form of the F128 cofactor W): the two 30-bit NTT split primes
    // of the committed gate (both = 1 mod 4). The live subset is fast (sin-s lemma + a ~20-λ read-off
    // spot check + the GF(p) certificate at 2 points/prime); the full 557-λ read-off is a test.
    private static readonly long[] F133Primes = { 2013265921L, 1811939329L };
    private const int F133GfpSamplesPerPrime = 2;

    // building P₁ (~1 s) is shared across Summary and the F133 node (built once per process).
    private static readonly Lazy<WSymplecticClosedForm.Halves> F133Halves =
        new(WSymplecticClosedForm.BuildHalves);

    // ~20 window λ spot-checked live against the embedded table (incl. the three max-|n| = 8 entries,
    // the empty partition n = −1, and several in-window zeros). Expected values come from the table.
    private static readonly int[][] F133SpotLambdas =
    {
        Array.Empty<int>(), new[]{3,3,2,2,1,1}, new[]{4,2,2,2,1,1}, new[]{4,4,2,2,1,1},
        new[]{1,1}, new[]{2}, new[]{4,2}, new[]{5,3}, new[]{3,3}, new[]{4,4}, new[]{8}, new[]{10},
        new[]{2,2,2,2,2,2}, new[]{3,3,3,3}, new[]{4,4,4}, new[]{6,4},
        new[]{2,2}, new[]{6,6}, new[]{4,4,4,4}, new[]{8,8},
    };

    // F139 (the seam identity): exact ℤ, deterministic, no primes; computed once per process.
    private static readonly Lazy<WSymplecticClosedForm.SeamIdentityReport> Seam =
        new(WSymplecticClosedForm.AnalyzeSeamIdentity);

    // the two exact ℤ[ζ_2n] layers are deterministic and independent of the primes;
    // computed once per process and shared between Summary and their nodes (no double run).
    private static readonly Lazy<LevelCollisionCensus.Report> Census =
        new(() => LevelCollisionCensus.Analyze());
    private static readonly Lazy<CollisionDecoupling.Report> Decoupling =
        new(CollisionDecoupling.Analyze);
    private static readonly Lazy<CollisionFamilyInventory.Report> Inventory =
        new(() => CollisionFamilyInventory.Analyze());

    public string DisplayName =>
        "F127/F128/F129/F130/F133/F134/F139 cross-triple orthogonality (live: 𝔉-slice + core-identity T + closed form + " +
        "sheet lattice + F128 flip lemma/factorization/sharper locus + F129 exact level census + " +
        "F129 family-inventory sum tie + F130 exact collision decoupling + F133 symplectic closed form of W + " +
        "F134 two-row reflection law + F139 seam identity)";

    /// <summary>The live F133 subset (fast): the sin-s lemma over ℤ, a ~20-λ read-off spot check vs the
    /// embedded table, and the GF(p) closed-form certificate at a few points/prime (table coefficients
    /// for speed; the disjoint read-off-derived path and the full 557-λ sweep are the tests).</summary>
    private static (bool SinSOk, int Surviving, int P1, int P2, int SpotMismatch,
                    int GfpPts, int GfpMism, int GfpCtrl, int GfpCtrlBroke,
                    WSymplecticClosedForm.TwoRowReflectionReport Refl) F133Evaluate()
    {
        var lem = WSymplecticClosedForm.AnalyzeSinSLemma();
        bool sinSOk = lem.SurvivingMonomials == 0 && lem.ProjectorSelfTestOk;

        var h = F133Halves.Value;
        var table = WSymplecticClosedForm.TableCoefficients();
        var tableMap = new Dictionary<string, long>();
        foreach (var (lam, n) in table) tableMap[string.Join(",", lam)] = n;

        int spotMismatch = 0;
        foreach (var lam in F133SpotLambdas)
        {
            long want = 2 * (tableMap.TryGetValue(string.Join(",", lam), out long n) ? n : 0);
            if (WSymplecticClosedForm.NRaw(lam, h) != want) spotMismatch++;
        }

        int gfpPts = 0, gfpMism = 0, gfpCtrl = 0, gfpBroke = 0;
        foreach (long p in F133Primes)
        {
            var g = WSymplecticClosedForm.CertifyGfpSlice(p, F133GfpSamplesPerPrime, table);
            gfpPts += g.Points; gfpMism += g.Mismatches;
            gfpCtrl += g.ControlChecks; gfpBroke += g.ControlMismatches;
        }
        var refl = WSymplecticClosedForm.AnalyzeTwoRowReflection(h);
        return (sinSOk, lem.SurvivingMonomials, h.P1Size, h.P2Size, spotMismatch,
                gfpPts, gfpMism, gfpCtrl, gfpBroke, refl);
    }

    public string Summary
    {
        get
        {
            int onVariety = 0, bad = 0, controls = 0, samples = 0;
            int coreEvals = 0, coreBad = 0, coreCtrl = 0, coreCtrlEval = 0;
            int cfPts = 0, cfBad = 0, cfCtrl = 0, cfCtrlEval = 0;
            int shPts = 0, shBad = 0, shCtrl = 0, shCtrlEval = 0;
            int fxPts = 0, fxBad = 0, fxCtrl = 0, fxCtrlEval = 0;
            int flPts = 0, flBad = 0, flCtrl = 0, flCtrlEval = 0;
            foreach (long p in Primes)
            {
                var (ev, b, c, s) = CrossFormCertificate.CertifySlice(p, SamplesPerPrime);
                onVariety += ev; bad += b; controls += c; samples += s;
                var (cev, cb, cc, cce) = CrossFormCertificate.CertifyCoreIdentitySlice(p, CoreSamplesPerPrime);
                coreEvals += cev; coreBad += cb; coreCtrl += cc; coreCtrlEval += cce;
                var (fp, fm, fc, fce) = CrossFormCertificate.CertifyClosedFormSlice(p, ClosedFormSamplesPerPrime);
                cfPts += fp; cfBad += fm; cfCtrl += fc; cfCtrlEval += fce;
                var (sp, sb, sc, sce) = CrossFormCertificate.CertifySharperLocusSlice(p, SharperSamplesPerPrime);
                shPts += sp; shBad += sb; shCtrl += sc; shCtrlEval += sce;
                var (xp, xm, xc, xce) = CrossFormCertificate.CertifyF128FactorizationSlice(p, F128SamplesPerPrime);
                fxPts += xp; fxBad += xm; fxCtrl += xc; fxCtrlEval += xce;
                var (lp, lb, lc, lce) = CrossFormCertificate.CertifyF128SharperLocusSlice(p, F128SharperSamplesPerPrime);
                flPts += lp; flBad += lb; flCtrl += lc; flCtrlEval += lce;
            }
            var flip = FlipLemma.Analyze();
            bool flipOk = flip.NumeratorMonomials == 8640 && flip.SurvivingMonomials == 0 &&
                          flip.ProjectorSelfTestOk && flip.ShiftedSets == 28 &&
                          flip.AllShiftedSetsDie && flip.AllShiftedProjectionsVanish;
            var lattice = SheetLattice.Analyze();
            bool latticeOk = lattice.Atoms == SheetLattice.AtomCount &&
                             lattice.Events == SheetLattice.EventCount &&
                             lattice.AllCoefficientsPmOne &&
                             lattice.DistinctSheets == SheetLattice.SheetCount &&
                             lattice.MinEventsPerSheet == SheetLattice.EventsPerSheet &&
                             lattice.MaxEventsPerSheet == SheetLattice.EventsPerSheet &&
                             lattice.EachSheetOneEventPerBlock;
            var census = Census.Value;
            bool censusOk = census.AllNonFiringInjective && census.AllFiringCollide &&
                            census.LawEquality && census.AnchorsExact;
            var dec = Decoupling.Value;
            bool decOk = dec.AllCollisionsDecouple && dec.AllControlsNonzero &&
                         dec.AllLemma3SignsOk && dec.LevelFlagsAsExpected;
            var inv = Inventory.Value;
            bool invOk = inv.AllRowsTied && inv.Capstone70Tied && inv.MSplitSumsToCountM;
            var f133 = F133Evaluate();
            bool f133Ok = f133.SinSOk && f133.P1 == 590016 && f133.P2 == 5817 &&
                          f133.SpotMismatch == 0 && f133.GfpMism == 0 &&
                          f133.GfpCtrlBroke == f133.GfpCtrl && f133.GfpCtrl > 0;
            var refl = f133.Refl;
            bool f134Ok = refl.PairsChecked == 36 && refl.LivePairs == 14 && refl.Mismatches == 0 &&
                          refl.L2PairsChecked == 16 && refl.L2Breaks == 8 &&
                          refl.L2CenterHoldsNonzero &&
                          refl.SpotPairsChecked == refl.LivePairs && refl.SpotMismatches == 0;
            var seam = Seam.Value;
            bool f139Ok = seam.SkewOk && seam.DegreeLemmaOk && seam.TableMismatches == 0 &&
                          seam.TableCellsChecked == 27 && seam.QuotientsMatch &&
                          seam.RemainderBoundsOk && seam.K0RemainderZero &&
                          seam.FenceL1Ok && seam.FenceL2Overflow && seam.CorruptionBroke;
            bool pass = bad == 0 && controls >= (int)(0.9 * samples) &&
                        coreBad == 0 && coreCtrl >= (int)(0.8 * coreCtrlEval) &&
                        cfBad == 0 && cfCtrl >= (int)(0.8 * cfCtrlEval) &&
                        shBad == 0 && shCtrl >= (int)(0.8 * shCtrlEval) && latticeOk &&
                        flipOk && fxBad == 0 && fxCtrl >= (int)(0.8 * fxCtrlEval) &&
                        flBad == 0 && flCtrl >= (int)(0.8 * flCtrlEval) &&
                        censusOk && decOk && invOk && f133Ok && f134Ok && f139Ok;
            return $"{(pass ? "PASS" : "FAIL")}: 𝔉-slice {onVariety} zeros ({bad} bad), " +
                   $"{controls}/{samples} controls nonzero; core-T {coreEvals} zeros ({coreBad} bad), " +
                   $"{coreCtrl}/{coreCtrlEval} controls nonzero; closed form {cfPts} generic points " +
                   $"({cfBad} mismatch), {cfCtrl}/{cfCtrlEval} corruptions broke; sharper locus {shPts} " +
                   $"points ({shBad} bad), {shCtrl}/{shCtrlEval} controls nonzero; lattice " +
                   $"{lattice.DistinctSheets} sheets × {lattice.MaxEventsPerSheet} events {(latticeOk ? "OK" : "BAD")}; " +
                   $"F128 flip lemma {flip.NumeratorMonomials}→{flip.SurvivingMonomials} " +
                   $"{(flipOk ? "OK" : "BAD")}, factorization {fxPts} points ({fxBad} mismatch), " +
                   $"{fxCtrl}/{fxCtrlEval} corruptions broke, 𝔉-sharper locus {flPts} points " +
                   $"({flBad} bad), {flCtrl}/{flCtrlEval} controls nonzero; F129 exact census n ≤ " +
                   $"{census.MaxN} ({census.CleanTriplesChecked} clean triples) {(censusOk ? "OK" : "BAD")}; " +
                   $"F130 exact decoupling {dec.Collisions.Count} collision pairs / {dec.Controls.Count} " +
                   $"controls {(decOk ? "OK" : "BAD")}; F129 inventory {inv.RowsChecked} rows tied + " +
                   $"capstone 70 + M-split sum {(invOk ? "OK" : "BAD")}; F133 W closed form: sin-s lemma " +
                   $"→{f133.Surviving}, read-off P₁={f133.P1}·P₂={f133.P2} spot {f133.SpotMismatch} mismatch, " +
                   $"GF(p) {f133.GfpPts} pts ({f133.GfpMism} mismatch), {f133.GfpCtrlBroke}/{f133.GfpCtrl} " +
                   $"corruptions broke {(f133Ok ? "OK" : "BAD")}; F134 two-row reflection " +
                   $"{refl.Mismatches}/{refl.PairsChecked} mismatch ({refl.LivePairs} live, spot " +
                   $"{refl.SpotMismatches} mismatch), l=2 fence {refl.L2Breaks}/{refl.L2PairsChecked} " +
                   $"breaks {(f134Ok ? "OK" : "BAD")}; F139 seam identity: a-priori table " +
                   $"{seam.TableMismatches}/{seam.TableCellsChecked} mismatch, quotients = (−1)^k P_k " +
                   $"{(seam.QuotientsMatch ? "OK" : "BAD")}, deg R_k ≤ 4+k {(seam.RemainderBoundsOk ? "OK" : "BAD")}, " +
                   $"fence l1/l2 {(seam.FenceL1Ok && seam.FenceL2Overflow ? "OK" : "BAD")} " +
                   $"→ {(f139Ok ? "OK" : "BAD")}";
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return SheetLatticeNode();
            yield return CoreIdentityNode();
            yield return ClosedFormNode();
            yield return F128Node();
            yield return F129CensusNode();
            yield return F129InventoryNode();
            yield return F130DecouplingNode();
            yield return F133ClosedFormNode();
            yield return F134ReflectionNode();
            yield return F139SeamNode();
            yield return CrossFormNode();
            yield return GatesNode();
            yield return new InspectableNode("What this is",
                summary: "the residue-collapse chain recomputed live in C#: the §2 lattice (exact integer), " +
                         "the §3 core identity T (GF(p), independent of 𝔉), the §3 closed form + its sharper " +
                         "locus, the F128 flip lemma (exact ℤ) + factorization + 𝔉-sharper locus, the " +
                         "𝔉-slice, the F129 exact ℤ[ζ_2n] level census (n ≤ 60; the n ≤ 210 census stays " +
                         "with the Python gate, the corner closed empty 2026-07-15), the F129 family-inventory " +
                         "sum tie (eleven families live, L at the n = 70 capstone; membership and the M split " +
                         "stay with the gate's I1-I5), the F130 exact collision decoupling, and the F133 " +
                         "symplectic closed form of the F128 cofactor W (the sin-s lemma over ℤ, a read-off " +
                         "spot check vs the embedded 143-term table, and the GF(p) alternant certificate), " +
                         "and the F134 two-row reflection law n_(j,k) = n_(10−j,k) (table + read-off spot + " +
                         "the l = 2 domain fence), and the F139 seam identity deriving F134 from F133 " +
                         "(Φ_k = S₁₀·(−1)^k P_k + R_k, deg R_k ≤ 4+k: the wall as a Chebyshev divisor); " +
                         "a SLICE of F127/F128 (not the wall) and an independent exact re-derivation of " +
                         "F129/F130/F133's content on the witness range");
        }
    }

    private static InspectableNode SheetLatticeNode()
    {
        var r = SheetLattice.Analyze();
        return new InspectableNode("§2 sheet lattice (exact integer)",
            summary: $"{r.Atoms} atoms → {r.Events} pole events → {r.DistinctSheets} sheets; " +
                     $"all coefficients ±1: {r.AllCoefficientsPmOne}; " +
                     $"{r.MinEventsPerSheet}..{r.MaxEventsPerSheet} events per sheet; " +
                     $"one event per (i,j) block per sheet: {r.EachSheetOneEventPerBlock}",
            provenance: NodeProvenance.Live);
    }

    private static InspectableNode CoreIdentityNode()
    {
        int evals = 0, bad = 0, ctrl = 0, ctrlEval = 0;
        var perPrime = new List<IInspectable>();
        foreach (long p in Primes)
        {
            var (ev, b, c, ce) = CrossFormCertificate.CertifyCoreIdentitySlice(p, CoreSamplesPerPrime);
            evals += ev; bad += b; ctrl += c; ctrlEval += ce;
            perPrime.Add(new InspectableNode($"p = {p}",
                summary: $"{ev} on-variety points, T = 0 (nonzero: {b}, must be 0); " +
                         $"{c}/{ce} off-variety controls nonzero (must be ≥ 80%)",
                provenance: NodeProvenance.Live));
        }
        return new InspectableNode("§3 core identity T (live GF(p))",
            summary: $"{evals} three-constraint variety points, {bad} nonzero (must be 0); " +
                     $"{ctrl}/{ctrlEval} off-variety controls nonzero; T independent of 𝔉",
            children: perPrime, provenance: NodeProvenance.Live);
    }

    private static InspectableNode ClosedFormNode()
    {
        int pts = 0, bad = 0, ctrl = 0, ctrlEval = 0;
        int shPts = 0, shBad = 0, shCtrl = 0, shCtrlEval = 0;
        var perPrime = new List<IInspectable>();
        foreach (long p in Primes)
        {
            var (fp, fm, fc, fce) = CrossFormCertificate.CertifyClosedFormSlice(p, ClosedFormSamplesPerPrime);
            pts += fp; bad += fm; ctrl += fc; ctrlEval += fce;
            var (sp, sb, sc, sce) = CrossFormCertificate.CertifySharperLocusSlice(p, SharperSamplesPerPrime);
            shPts += sp; shBad += sb; shCtrl += sc; shCtrlEval += sce;
            perPrime.Add(new InspectableNode($"p = {p}",
                summary: $"identity at {fp} generic points ({fm} mismatch, must be 0), {fc}/{fce} corruptions " +
                         $"broke; sharper locus {sp} points ({sb} nonzero, must be 0), {sc}/{sce} controls nonzero",
                provenance: NodeProvenance.Live));
        }
        return new InspectableNode("§3 closed form (live GF(p), 2026-07-14)",
            summary: $"T·P = ⅛[2cos s((e₁−f₁)² − 2sin²s) + sin s(Σsin2a + Σsin2b)]·V_a·V_b at {pts} generic " +
                     $"points ({bad} mismatch); corruption control {ctrl}/{ctrlEval} broke; Corollary 2's sharper " +
                     $"locus (Σcos a = Σcos b ≠ 0, sheet): T = 0 at {shPts} points ({shBad} bad), " +
                     $"{shCtrl}/{shCtrlEval} controls nonzero",
            children: perPrime, provenance: NodeProvenance.Live);
    }

    private static InspectableNode F128Node()
    {
        var flip = FlipLemma.Analyze();
        bool flipOk = flip.NumeratorMonomials == 8640 && flip.SurvivingMonomials == 0 &&
                      flip.ProjectorSelfTestOk && flip.ShiftedSets == 28 &&
                      flip.AllShiftedSetsDie && flip.AllShiftedProjectionsVanish;
        var children = new List<IInspectable>
        {
            new InspectableNode("The flip lemma (exact over ℤ)",
                summary: $"𝒪[cos s·B·V_aV_bP̃]: {flip.NumeratorMonomials} monomials → " +
                         $"{flip.SurvivingMonomials} after the signed character sum (must be 0); " +
                         $"projector self-test {(flip.ProjectorSelfTestOk ? "OK" : "BAD")}; " +
                         $"{flip.ShiftedSets} shifted alternant exponent sets, all die " +
                         $"(0/repeat/±pair): {flip.AllShiftedSetsDie}, and each one's odd " +
                         $"projection is the zero polynomial: {flip.AllShiftedProjectionsVanish}",
                provenance: NodeProvenance.Live),
        };
        int fxPts = 0, fxBad = 0, fxCtrl = 0, fxCtrlEval = 0;
        int flPts = 0, flBad = 0, flCtrl = 0, flCtrlEval = 0;
        foreach (long p in Primes)
        {
            var (xp, xm, xc, xce) = CrossFormCertificate.CertifyF128FactorizationSlice(p, F128SamplesPerPrime);
            fxPts += xp; fxBad += xm; fxCtrl += xc; fxCtrlEval += xce;
            var (lp, lb, lc, lce) = CrossFormCertificate.CertifyF128SharperLocusSlice(p, F128SharperSamplesPerPrime);
            flPts += lp; flBad += lb; flCtrl += lc; flCtrlEval += lce;
            children.Add(new InspectableNode($"p = {p}",
                summary: $"factorization at {xp} generic points ({xm} mismatch, must be 0), {xc}/{xce} " +
                         $"corruptions broke; 𝔉-sharper locus {lp} points ({lb} nonzero, must be 0), " +
                         $"{lc}/{lce} controls nonzero",
                provenance: NodeProvenance.Live));
        }
        return new InspectableNode("F128 factorization + sharper locus (live, 2026-07-14 evening)",
            summary: $"𝔉 = −(e₁−f₁)²·𝒪[cos s·cot s·V_aV_b/P]: flip lemma {(flipOk ? "OK" : "BAD")} " +
                     $"(exact ℤ); factorization at {fxPts} generic GF(p) points ({fxBad} mismatch), " +
                     $"{fxCtrl}/{fxCtrlEval} corruptions broke; 𝔉 = 0 on {{Σcos a = Σcos b ≠ 0}} at " +
                     $"{flPts} points ({flBad} bad), {flCtrl}/{flCtrlEval} controls nonzero: ONE " +
                     $"constraint, no sheet; F127's V is the codim-2 special case",
            children: children, provenance: NodeProvenance.Live);
    }

    private static InspectableNode CrossFormNode()
    {
        int onVariety = 0, bad = 0, ctrl = 0, samples = 0;
        var perPrime = new List<IInspectable>();
        foreach (long p in Primes)
        {
            var (ev, b, c, s) = CrossFormCertificate.CertifySlice(p, SamplesPerPrime);
            onVariety += ev; bad += b; ctrl += c; samples += s;
            long iRoot = CrossFormCertificate.SqrtMinusOne(p);
            perPrime.Add(new InspectableNode($"p = {p}",
                summary: $"i = {iRoot}; {ev} on-variety evaluations, {b} nonzero (must be 0); " +
                         $"{c}/{s} off-variety controls nonzero (must be ≥ 90%)",
                provenance: NodeProvenance.Live));
        }
        return new InspectableNode("𝔉-slice (live GF(p))",
            summary: $"{onVariety} on-variety evaluations, {bad} nonzero (must be 0); " +
                     $"{ctrl}/{samples} off-variety controls nonzero",
            children: perPrime, provenance: NodeProvenance.Live);
    }

    private static InspectableNode F129CensusNode()
    {
        var r = Census.Value;
        bool ok = r.AllNonFiringInjective && r.AllFiringCollide && r.LawEquality && r.AnchorsExact;
        return new InspectableNode("F129 level-collision census (exact ℤ[ζ_2n], live)",
            summary: $"n ≤ {r.MaxN}, {r.CleanTriplesChecked} clean triples: {r.NonFiringChecked} non-firing n " +
                     $"all INJECTIVE (exact level vectors pairwise distinct): {r.AllNonFiringInjective}; " +
                     $"{r.FiringChecked} firing n (3|n ≥ 9, 10|n ≥ 20) all carry an exact collision: " +
                     $"{r.AllFiringCollide}; firing set = predicted set: {r.LawEquality}; the n=15 R₃-cycle and " +
                     $"n=20 R₅+zero-mode anchors term-exact: {r.AnchorsExact} → {(ok ? "OK" : "BAD")}. " +
                     "Scope: the law's n ≤ 210 census stays with " +
                     "simulations/f129_level_collision_law.py (corner closed empty 2026-07-15, PROOF_F129 §4)",
            provenance: NodeProvenance.Live);
    }

    private static InspectableNode F129InventoryNode()
    {
        var r = Inventory.Value;
        bool ok = r.AllRowsTied && r.Capstone70Tied && r.MSplitSumsToCountM;
        var (fanned, middle, fixedVertex) = CollisionFamilyInventory.MSplit;
        return new InspectableNode("F129 family inventory (exact ℤ[ζ_2n] sum tie, live)",
            summary: $"thirteen derived closed forms (PROOF_F129_FAMILY_INVENTORY_COUNTS.md): their sums " +
                     $"(total + d-split) tie to the exact census at all {r.RowsChecked} rows n ≤ {r.MaxN} " +
                     $"(eleven families A..K exercised live): {r.AllRowsTied}; the n = 70 capstone pins L, " +
                     $"the corner-closure's second mechanism (C 120 + L 20, all disjoint): {r.Capstone70Tied}; " +
                     $"M split {fanned}+{middle}+{fixedVertex} sums to Count(M) (internal consistency of " +
                     $"adopted constants, not a from-below leg): {r.MSplitSumsToCountM} " +
                     $"→ {(ok ? "OK" : "BAD")}. " +
                     "Scope: the tie certifies the SUMS; per-family membership (first-fit decomposition, " +
                     "labels) and M's split reconstruction stay with simulations/f129_family_inventory.py " +
                     "(I2/I3/I5); the n = 105 capstone (total 8858) is pinned in the tests",
            provenance: NodeProvenance.Live);
    }

    private static InspectableNode F130DecouplingNode()
    {
        var r = Decoupling.Value;
        bool ok = r.AllCollisionsDecouple && r.AllControlsNonzero && r.AllLemma3SignsOk && r.LevelFlagsAsExpected;
        var children = new List<IInspectable>();
        foreach (var c in r.Collisions)
            children.Add(new InspectableNode(c.Label,
                summary: $"ε = {c.Eps:+0;-0}; levels equal (exact): {c.LevelsEqual}; " +
                         $"Ê⁺ = 0: {c.PlusZero}, Ê⁻ = 0: {c.MinusZero}; Lemma-3 sign: {c.Lemma3SignOk}",
                provenance: NodeProvenance.Live));
        foreach (var c in r.Controls)
            children.Add(new InspectableNode(c.Label,
                summary: $"ε = {c.Eps:+0;-0}; levels equal (exact): {c.LevelsEqual}; " +
                         $"nonzero as required: {!c.PlusZero && !c.MinusZero}; Lemma-3 sign: {c.Lemma3SignOk}",
                provenance: NodeProvenance.Live));
        return new InspectableNode("F130 collision decoupling (exact ℤ[ζ_2n], live)",
            summary: $"equal level ⟹ B(τ,σ) = 0: {r.Collisions.Count} named pairs (every proof cell) vanish " +
                     $"exactly: {r.AllCollisionsDecouple}; {r.Controls.Count} unequal-level controls nonzero: " +
                     $"{r.AllControlsNonzero}; Lemma-3 sign Ê⁺ = ε·Ê⁻ on every pair: {r.AllLemma3SignsOk}; " +
                     $"level flags as designed: {r.LevelFlagsAsExpected} → {(ok ? "OK" : "BAD")}. " +
                     "No floats, no primes: ŝ(a) = ζ^a − ζ^{−a}, norms cancel from the zero claim",
            children: children, provenance: NodeProvenance.Live);
    }

    private static InspectableNode F133ClosedFormNode()
    {
        var f = F133Evaluate();
        bool ok = f.SinSOk && f.P1 == 590016 && f.P2 == 5817 && f.SpotMismatch == 0 &&
                  f.GfpMism == 0 && f.GfpCtrlBroke == f.GfpCtrl && f.GfpCtrl > 0;
        var children = new List<IInspectable>
        {
            new InspectableNode("The sin-s lemma (exact over ℤ)",
                summary: $"𝒪[(2i sin s)·(2i)¹⁵ Δ̂]: {f.Surviving} monomials survive the (ℤ/2)⁶ signed " +
                         $"character sum (must be 0); projector self-test {(f.SinSOk ? "OK" : "BAD")}",
                provenance: NodeProvenance.Live),
            new InspectableNode("The alternant read-off (meet-in-the-middle, spot check)",
                summary: $"X = P₁·P₂ with |P₁| = {f.P1} (Δ̂ · 15 sheets), |P₂| = {f.P2} (16 sheets); " +
                         $"n_raw(λ) = 2·n_λ on {F133SpotLambdas.Length} spot-checked window λ (incl. the three " +
                         $"|n| = 8 entries, λ = () → −1, and in-window zeros): {f.SpotMismatch} mismatch " +
                         $"(must be 0). The full 557-λ sweep + the aggregate checksums (143, max 8, sum 359) " +
                         "are the tests",
                provenance: NodeProvenance.Live),
        };
        foreach (long p in F133Primes)
        {
            var g = WSymplecticClosedForm.CertifyGfpSlice(p, F133GfpSamplesPerPrime,
                WSymplecticClosedForm.TableCoefficients());
            children.Add(new InspectableNode($"p = {p}",
                summary: $"C₆·SP == 2⁻⁵·(2i)⁻⁴⁶·Σ n_λ A_{{λ+ρ}} at {g.Points} points ({g.Mismatches} mismatch, " +
                         $"must be 0); {g.ControlMismatches}/{g.ControlChecks} corruptions (one n_λ bumped) broke",
                provenance: NodeProvenance.Live));
        }
        return new InspectableNode("F133 symplectic closed form of W (live, 2026-07-17)",
            summary: $"W = −2⁹·(Π sin x_u)·V_c(a)·V_c(b)·K/SP, K = 2⁻³⁰·Σ_λ n_λ·χ^{{C₆}}_λ (143 terms): " +
                     $"sin-s lemma →{f.Surviving} (exact ℤ), read-off P₁={f.P1}·P₂={f.P2} spot " +
                     $"{f.SpotMismatch} mismatch, GF(p) certificate {f.GfpPts} points ({f.GfpMism} mismatch), " +
                     $"{f.GfpCtrlBroke}/{f.GfpCtrl} corruptions broke → {(ok ? "OK" : "BAD")}. " +
                     "The GF(p) sum uses the read-off n_λ (disjoint from the flip sum); the live check uses " +
                     "the embedded table for speed, the disjoint read-off-derived path is the test " +
                     "(docs/proofs/PROOF_F133_W_SYMPLECTIC_CLOSED_FORM.md, simulations/f133_w_closed_form.py)",
            children: children, provenance: NodeProvenance.Live);
    }

    private static InspectableNode F134ReflectionNode()
    {
        var r = F133Evaluate().Refl;
        bool ok = r.PairsChecked == 36 && r.LivePairs == 14 && r.Mismatches == 0 &&
                  r.L2PairsChecked == 16 && r.L2Breaks == 8 && r.L2CenterHoldsNonzero &&
                  r.SpotPairsChecked == r.LivePairs && r.SpotMismatches == 0;
        return new InspectableNode("F134 two-row reflection law (live, 2026-07-17)",
            summary: $"n_(j,k) = n_(10−j,k) (μ₁ ↦ 22−μ₁): {r.Mismatches}/{r.PairsChecked} table mismatches " +
                     $"({r.LivePairs} live pairs), read-off spot on every live pair {r.SpotMismatches} mismatch " +
                     $"(disjoint from the table); domain fence: l = 2 breaks {r.L2Breaks}/{r.L2PairsChecked} " +
                     $"(must be exactly 8), the two nonzero l = 2 holds are the center-fixed j = 5 weights: " +
                     $"{r.L2CenterHoldsNonzero} → {(ok ? "OK" : "BAD")}. Certificate grade (the F127-wall class); " +
                     "the structural mechanism is F139's seam identity, the sibling node below " +
                     "(docs/proofs/PROOF_F134_TWO_ROW_REFLECTION_LAW.md, simulations/f134_two_row_reflection_law.py)",
            provenance: NodeProvenance.Live);
    }

    private static InspectableNode F139SeamNode()
    {
        var s = Seam.Value;
        bool ok = s.SkewOk && s.DegreeLemmaOk && s.TableMismatches == 0 && s.TableCellsChecked == 27 &&
                  s.QuotientsMatch && s.RemainderBoundsOk && s.K0RemainderZero &&
                  s.FenceL1Ok && s.FenceL2Overflow && s.CorruptionBroke;
        return new InspectableNode("F139 seam identity (live, exact ℤ, 2026-07-21)",
            summary: $"Φ_k = S₁₀·(−1)^k P_k + R_k, deg R_k ≤ 4+k (the F134 wall as a Chebyshev divisor): " +
                     $"ψ-skew {(s.SkewOk ? "OK" : "BAD")}, deg Φ_k ≤ 15 {(s.DegreeLemmaOk ? "OK" : "BAD")}, " +
                     $"a-priori table {s.TableMismatches}/{s.TableCellsChecked} mismatch vs the embedded " +
                     $"two-row table, quotients = (−1)^k·P_k {(s.QuotientsMatch ? "OK" : "BAD")}, remainder " +
                     $"bounds {(s.RemainderBoundsOk ? "OK" : "BAD")} (k = 0 exactly zero: {s.K0RemainderZero}), " +
                     $"fence l=1 small/l=2 overflow {(s.FenceL1Ok && s.FenceL2Overflow ? "OK" : "BAD")}, " +
                     $"corruption broke: {s.CorruptionBroke} → {(ok ? "OK" : "BAD")}. " +
                     "Derives F134 from F133 (docs/proofs/PROOF_F139_SEAM_IDENTITY.md, " +
                     "simulations/f139_seam_identity.py); own five-variable engine, disjoint from the " +
                     "six-variable read-off and the Python gate",
            provenance: NodeProvenance.Live);
    }

    private static InspectableNode GatesNode() =>
        new("The proof chain (committed gates)",
            summary: "docs/proofs/PROOF_F127_RESIDUE_COLLAPSE.md + PROOF_F128_FLIP_SUM_FACTORIZATION.md " +
                     "+ eight committed gate scripts",
            children: new IInspectable[]
            {
                new InspectableNode("f128_flip_sum_factorization.py",
                    summary: "F128: the flip lemma over ℤ (G1) + Weyl folding (G2) + raising shifts (G3) + " +
                             "the 28 deaths (G4) + factorization pin (G5) + sharper locus from below (G6) + " +
                             "restriction hygiene (G7)"),
                new InspectableNode("PROOF_F127_RESIDUE_COLLAPSE.md",
                    summary: "the derived chain: sheet lattice (§2), core identity by the closed form, prior " +
                             "resultant route retained (§3), transport (§4), the window + oddness (§5), the mirror anchor (§6)"),
                new InspectableNode("f127_closed_form.py",
                    summary: "§3 primary: the closed form T·P = RHS exact over ℚ(i) (864 monomials each side), " +
                             "the bordered Frobenius determinant, the three corollaries"),
                new InspectableNode("f127_core_identity.py",
                    summary: "§3 cross-check: naked-T falsification, N assembly (540 terms), Res | E exact, branch check"),
                new InspectableNode("f127_core_locus_patch.py",
                    summary: "§3.2: forcing chain, gcd(Res, S_w−m₀S_z) = 1, degenerate sublocus r₀ = r₁ = 0"),
                new InspectableNode("f127_exact_gates.py",
                    summary: "S1 (72 atoms), S2 (288 events), S3/S5/S6 (oddness bijections), S4 (mirror gate); LEM/tanLEM"),
                new InspectableNode("f127_transport_lemma.py",
                    summary: "§4: incidence table, 31 discriminants, norm identity (62), step-2 divisibility (384)"),
                new InspectableNode("f127_relabel_pin.py",
                    summary: "§5: the literal w₂-distinguished rerun (inventory, all poles simple, window ⊆ {−1,0,1})"),
                new InspectableNode("f127_sheet_lattice.py",
                    summary: "§2 numerics: the 32×9 lattice survey and the residue-bookkeeping gates"),
                new InspectableNode("f129_level_collision_law.py",
                    summary: "F129: the reduction exact (G1), injectivity via mod-p distinctness at all " +
                             "non-firing n ≤ 210 (G2), one exact collision per firing n (G3), the law " +
                             "equality (G4), both mechanism anchors term-exact (G5)"),
                new InspectableNode("f130_collision_decoupling.py",
                    summary: "F130: assembly (D) off-resonance at 400 pairs (G1), B = 0 at all 8452 " +
                             "equal-level pairs (G2), generic unequal-level controls (G3), the level-free " +
                             "cell (G4), overlap-2 vacuity (G5); float grade, the exactness lives here"),
                new InspectableNode("f129_family_inventory.py",
                    summary: "the F129 inventory: every pair decomposes (I1), the thirteen closed forms " +
                             "at every firing n ≤ 140 + capstones 150/210 (I2), the families partition " +
                             "the census (I3), family ⟺ door (I4), the M W-sets rebuilt from committed " +
                             "substitution recipes at n = 105 and 210 (I5)"),
                new InspectableNode("f139_seam_identity.py + PROOF_F139_SEAM_IDENTITY.md",
                    summary: "F139: the slice identity vs the committed read-off (G1), the ψ-skew (G2), " +
                             "the degree lemma (G3), the a-priori table 27/27 (G4), the division with " +
                             "Q_k = (−1)^k P_k and deg R_k ≤ 4+k (G5), the l = 1/l = 2 fence (G6), " +
                             "corruption control (G7)"),
                new InspectableNode("f133_w_closed_form.py + PROOF_F133_W_SYMPLECTIC_CLOSED_FORM.md",
                    summary: "F133: the pair escape (G1), the sin-s lemma 𝒪[sin s·Δ̂] ≡ 0 (G2), SP " +
                             "flip-invariance (G3), X alternates (G4), the 143 n_λ read off by " +
                             "meet-in-the-middle (G5, --full sweeps all 18564 dominant λ), the C₆ Weyl " +
                             "denominator (G6), the χ^C/m-basis table cross-link (G7), the end-to-end " +
                             "numeric pin (G8)"),
            });
}
