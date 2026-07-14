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
/// the exact Lemma-3 sign Ê⁺ = ε·Ê⁻ on every pair (<see cref="CollisionDecoupling"/>).</item>
/// </list>
/// This is the second implementation the F127/F128 code-trust caveat asks for, at sample scale; the
/// proof objects remain the 527/527 grid+CRT wall and the committed <c>f12*_*</c> gates (F129's
/// n ≤ 210 census and named corner stay with <c>simulations/f129_level_collision_law.py</c>).
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

    // the two exact ℤ[ζ_2n] layers are deterministic and independent of the primes;
    // computed once per process and shared between Summary and their nodes (no double run).
    private static readonly Lazy<LevelCollisionCensus.Report> Census =
        new(() => LevelCollisionCensus.Analyze());
    private static readonly Lazy<CollisionDecoupling.Report> Decoupling =
        new(CollisionDecoupling.Analyze);

    public string DisplayName =>
        "F127/F128/F129/F130 cross-triple orthogonality (live: 𝔉-slice + core-identity T + closed form + " +
        "sheet lattice + F128 flip lemma/factorization/sharper locus + F129 exact level census + " +
        "F130 exact collision decoupling)";

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
            bool pass = bad == 0 && controls >= (int)(0.9 * samples) &&
                        coreBad == 0 && coreCtrl >= (int)(0.8 * coreCtrlEval) &&
                        cfBad == 0 && cfCtrl >= (int)(0.8 * cfCtrlEval) &&
                        shBad == 0 && shCtrl >= (int)(0.8 * shCtrlEval) && latticeOk &&
                        flipOk && fxBad == 0 && fxCtrl >= (int)(0.8 * fxCtrlEval) &&
                        flBad == 0 && flCtrl >= (int)(0.8 * flCtrlEval) &&
                        censusOk && decOk;
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
                   $"controls {(decOk ? "OK" : "BAD")}";
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
            yield return F130DecouplingNode();
            yield return CrossFormNode();
            yield return GatesNode();
            yield return new InspectableNode("What this is",
                summary: "the residue-collapse chain recomputed live in C#: the §2 lattice (exact integer), " +
                         "the §3 core identity T (GF(p), independent of 𝔉), the §3 closed form + its sharper " +
                         "locus, the F128 flip lemma (exact ℤ) + factorization + 𝔉-sharper locus, the " +
                         "𝔉-slice, the F129 exact ℤ[ζ_2n] level census (n ≤ 60; the n ≤ 210 census + the " +
                         "named corner stay with the Python gate), and the F130 exact collision decoupling; " +
                         "a SLICE of F127/F128 (not the wall) and an independent exact re-derivation of " +
                         "F129/F130's content on the witness range");
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
                     "Scope: the law's n ≤ 210 census and its one named corner stay with " +
                     "simulations/f129_level_collision_law.py",
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
            });
}
