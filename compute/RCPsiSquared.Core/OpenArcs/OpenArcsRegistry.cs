namespace RCPsiSquared.Core.OpenArcs;

/// <summary>The open-arcs ledger: research threads that were started, reached a first
/// exemplar, then parked, indexed by name. The world's antidote to its chronic failure
/// mode (victory declared at the first exemplar, then forgotten): each entry records where
/// it stopped and the concrete next move, so a returning session looks the arc up instead
/// of re-discovering its incompleteness. Mirrors <see cref="Confirmations.ConfirmationsRegistry"/>
/// in shape; surfaced as the inspect "arcs" section via <see cref="Inspection.OpenArcsInspectableNode"/>.</summary>
public static class OpenArcsRegistry
{
    private static readonly IReadOnlyList<OpenArc> _all = new[]
    {
        new OpenArc(
            Name: "diabolic_over_higher_n",
            Opened: "2026-06-27",
            Origin: "the path-k diabolic investigation, now solid at two N (a diabolic = a SEMISIMPLE eigenvalue " +
                "coalescence, eigenvectors stay independent, silent/pass-through; as opposed to a DEFECTIVE EP, " +
                "where the eigenvectors also coalesce and the amplitude diverges; the whole arc is about counting " +
                "the diabolics and showing integrability creates them). Path-3 (N=4): the within-block " +
                "twin-pair-onto-fold diabolic at a REAL q_EP=0.659, λ_EP=−4+2iJ, traced from below (gmscan " +
                "--trace) and Δ-flipped defective by XXZ anisotropy (DIABOLIC_BY_INTEGRABILITY). Path-4 (N=5): " +
                "11 character-verified diabolics, ALL at complex q (none at physical real q), Δ-VERIFIED " +
                "integrability-protected (XxzCoherenceBlock + XxzDeltaFlipTests: they flip defective / lift under " +
                "Δ≠0, a defective control stays put). The synthesis (integrability = existence of the diabolics, " +
                "the N=4 self-fold [an antiunitary mirror the watched (SE,DE)=single-excitation/double-excitation " +
                "coherence block has ONLY at N=4] = what put ONE diabolic on the real axis): now Tier-1 for the " +
                "mechanism, full write-up + the term defined in experiments/F89_PATH_K_DIABOLIC.md. Tom's " +
                "question, the arc: does the structure " +
                "(complex-q-only, integrability-protected, and HOW MANY) generalize over N, and what would it " +
                "take to pursue it at path-5 (N=6), path-6 (N=7), and beyond?",
            ParkedAt: "[SUPERSEDED 2026-06-30 - this ParkedAt is the ORIGINAL two-N snapshot; the arc has since run " +
                "N=4 through N=9, the onset is grounded (parity-gated odd N>=7, dimension-mismatch/sector-swap), and the C# " +
                "witness is built. Read NextStep's 'CURRENT STATE' first; this field is kept as the historical why-it-looked-hard.] " +
                "TWO N DONE; the generalization not run. WHAT GENERALIZES FOR FREE (the good news): compute " +
                "is NOT the obstacle. The (SE,DE) coherence block is POLYNOMIAL in N, dim = N·C(N,2) ~ N³/2 " +
                "(24/50/90/147 for N=4..7; residual F_d = 8/18/32/53 for path-3..6), TINY beside the 4^N " +
                "Liouvillian, so an EVD per q-point is sub-millisecond at every N reachable here; the (SE,DE) " +
                "restriction is the enabler. And the Δ-test tooling is ALREADY N-general: XxzCoherenceBlock " +
                "(BuildFull/BuildSym/SeDeSymSpectrum/CharacterAtDiabolicNear/TrackDiabolicUnderDelta all take n), " +
                "plus gmscan --trace and pkmono --delta-flip. So path-5/6 need a RUN, not a rebuild, for the " +
                "Δ-test. WHAT NEEDS EXTENSION: the diabolic LOCATOR (FindDiabolics / pkmono --diabolic) rides " +
                "PathKMonodromyScout, valid for the k the F89PathKFdOracle covers (≤6); for k>6 extend the oracle " +
                "+ F89AtFactorReconstruction.ForPathK. (Caveat: the AT-removal [AT = the Absorption Theorem, the " +
                "closed-form decay-rate strands at Re lambda = -2<n_XY> in gamma units] is needed only for the monodromy/" +
                "Galois layer; the diabolic hunt itself is coalescences in the FULL symmetric block and does not " +
                "require it, so XxzCoherenceBlock can locate-and-test at any N on its own.) THE REAL OBSTACLES " +
                "(the genuinely hard part, get harder with N): (1) COMPLETENESS. The exact discriminant " +
                "disc_λ(F_d) factorization (Route B, the only way to COUNT diabolics definitively, as " +
                "(3q⁴+q²−1)² did at path-3) is already infeasible at F_18 (degree-52 disc) and hopeless at " +
                "F_32/F_53; beyond path-3 the count is answerable only by bounded-region NUMERICAL scans with a " +
                "stated coverage box, never provably complete. (2) The scan region/resolution must scale: more " +
                "strands → more inter-sector crossings → more candidates, plus the dense-spectrum hazards already " +
                "met (loop contamination needing the small intrinsic radius; the √-cusp needing GapRefine-from-" +
                "box-min for LIFT-vs-defect). Engineering, not fundamental. (3) The naive Slater additivity " +
                "(E_DE = pairwise SE sums) bookkeeping FAILS from path-5 (the F_b modes are not simple sums; the " +
                "AT factor needs the rate-confined reconstruction), BUT free-fermion integrability itself " +
                "(Jordan-Wigner) holds at every N, so the diabolics persist and stay Δ-killable; only the AT-" +
                "bookkeeping is harder. The deferred Q4 cross-fold edge (do the complex-q diabolics pair across " +
                "the cross-block fold (SE,DE)↔(SE,w_{N-2})?) also generalizes over N and is untried.",
            NextStep: "CURRENT STATE (2026-06-30, read this first; the dated layers below are the journal). " +
                "RESUMING IN ONE LINE: this arc is DONE except one minor edge - all four Moves AND the Move-4 follow-on (the " +
                "cross-fold is integrability-independent and docks onto F1 as Pi's bra leg) are closed; the one OPEN edge is the " +
                "within-odd EXACT threshold's CLOSED FORM (the threshold is empirically N=7, grounded on residual density; a closed " +
                "form is residual-density driven and likely does NOT exist), so the arc is effectively complete; the residual " +
                "grounding is already DONE (in DiabolicReflectionParityWitness), and the one bounded probe IF you still want to push " +
                "it is under 'REMAINING EDGES (a)' below. ARC TERMS (all glossed " +
                "in experiments/F89_PATH_K_DIABOLIC.md, its 'What this is about' + 'Terms used here' sections): path-k = the " +
                "(k+1)-site chain, so path-6 = N=7; (SE,DE) = the single-excitation/double-excitation coherence block; AT = the " +
                "Absorption-Theorem decay rate, a POSITIVE rate 2*gamma*n_diff, so the eigenvalue sits at Re lambda = " +
                "-2*gamma*n_diff = -2<n_XY> (gamma=1), with <n_XY> = n_diff the per-coherence XY-disagreement count (same quantity; " +
                "the sign is just rate-vs-eigenvalue); the exponent = the gap-scaling " +
                "exponent that classifies a coalescence (~1 for a diabolic / semisimple, ~0.5 for a defective EP / sqrt-branch); " +
                "EP = exceptional point (a non-Hermitian spectral degeneracy); F_d = the degree-d H_B-mixed " +
                "RESIDUAL polynomial (H_B = the XY bond Hamiltonian; F_18/F_53/F_116 are DEGREES 18/53/116, NOT F-registry " +
                "numbers like F89d); R = the site-reflection i->N-1-i (the order-2 symmetry 'S2') the block is built on; Sigma " +
                "(capital) = the antiunitary realness operator, Sigma*L*Sigma = L-dagger (the adjoint, the 'L+' below) - NOT " +
                "lowercase sigma_even/sigma_odd, which are the SPECTRA of the R-even/R-odd sectors. NOW: N=7 (path-6) " +
                "is DONE, AND its Δ-test (Move 2) is now DONE too - so the diabolic-CHARACTER question is CLOSED at N=7 " +
                "(all four real-q diabolics are integrability-protected, a PLACEMENT mechanism not a new species; details " +
                "in the 'MOVE 2 IS NOW DONE' block below). The real-q ONSET has since been pinned to a LAW by the " +
                "N=8/N=9 discriminator (2026-06-30, the 'ONSET RESOLVED' block below): it is PARITY-GATED, odd N >= 7 " +
                "(N=8 even has NO on-axis real-lambda diabolic; N=9 odd has >=3, each a clean isolated on-axis point) - " +
                "refuting both the threshold (N>=7) and the one-off readings. The CAUSE of that parity is now GROUNDED from " +
                "below (2026-06-30, the 'MECHANISM GROUNDED' block): the DIMENSION-MISMATCH / SECTOR-SWAP. R (i->N-1-i, the S2 " +
                "the block is built on) splits the (SE,DE) block into R-even/R-odd; the realness antiunitarity Sigma*L*Sigma=L+ " +
                "maps even<->odd at EVEN N (equal dims, sigma_even = conj(sigma_odd) EXACTLY) so neither sector is self-conjugate " +
                "and a real-axis collision is the generic pseudo-Hermitian DEFECTIVE EP; but the odd-N reflection-fixed central " +
                "site makes dim(even)-dim(odd)=(N-1)/2 != 0, forbidding the cross-pairing, so each sector is SELF-conjugate and two " +
                "REAL R-even eigenvalues can cross semisimply = the real-q diabolic. The C# witness for the mechanism is now BUILT " +
                "(DiabolicReflectionParityWitness, 'inspect --root diabolicparity', TDD-gated). The within-odd THRESHOLD has since been " +
                "GROUNDED ON THE RESIDUAL (2026-06-30): reading the parity directly on the AT-stripped residual (ResidualRootsExact at " +
                "real q, the F_d strands the diabolics live in) shows even N carries ZERO real residual eigenvalues anywhere in q in [0.2,3] " +
                "(no real strands to host a real-q diabolic), while odd N carries a real population that GROWS with N (up to 4/5/8 at " +
                "N=5/7/9, degrees F_18/F_53/F_116); the onset tracks this residual density (F_18 too sparse, F_53 the first to host a " +
                "real-q coalescence in the window). The recorded 'two crossing real residual strands' picture was REFINED: the onset " +
                "coalescences come in TWO geometries - a real-real CROSSING (two real strands meet, first at N=7, q~2.628) and a " +
                "CONJUGATE-PAIR TANGENCY (a (lambda,conj) pair whose Im->0 only at q*, e.g. N=9 at q~0.4755) - so counting real-real " +
                "crossings is NOT a proxy for the real-q diabolic count (that stays FindDiabolicsExact's job), and N=9 has real-q diabolics " +
                "yet NO real-real crossing. Both readings are in the witness (TDD-gated). The exact within-odd threshold N still has no " +
                "closed form (residual-density driven). MOVE 4 (the cross-fold) is now DONE (2026-06-30): the (SE,DE)=(w1,w2) <-> " +
                "(SE,w_{N-2})=(w1,N-2) fold (the F89c bra bit-flip rho[a,b]->rho[a,bbar]) is an EXACT antiunitary similarity at the " +
                "MATRIX level - L(1,N-2)(qbar) = -P conj(L(1,2)(q)) P^T - 2N*I to MACHINE ZERO (exact arithmetic), N=4..9, all q real " +
                "AND complex (P = bra-complement permutation). This is STRONGER than foldcross's spectrum match: an antiunitary " +
                "similarity preserves Jordan structure, so every (SE,DE) diabolic at (q,lambda) pairs with a (SE,w_{N-2}) diabolic at " +
                "(qbar, -lambdabar-2N) with IDENTICAL character and gap, for all N and all q at once - no enumeration. N=4 is the " +
                "degenerate partner=self within-block self-fold; for N>=5 the N=4 on-line 'zeros' become cross-block mirror partners. " +
                "Verified: the N=7 real-q diabolic (q=1.1264, lam=-4.942, gap 4.19e-5) pairs with (1,5) at -9.058, same gap. LANDED: a " +
                "shared Core builder WeightCoherenceBlock (general (wKet,wBra) block + BraComplementPermutation, promoted from " +
                "FoldCrossCommand.BuildBlock which now delegates to it), the live witness CrossFoldSimilarityWitness " +
                "('inspect --root crossfold', TDD-gated CrossFoldSimilarityWitnessTests), and REGISTERED as F89d " +
                "(docs/ANALYTICAL_FORMULAS.md, right after the F89c lemma it extends) + typed as F89CrossFoldSimilarityClaim " +
                "(parents F1PalindromeIdentity + F89BranchLocusPalindromeClaim, Tier1Derived, wiring-audited; verify via " +
                "'knowledge ancestors F89CrossFoldSimilarityClaim'). SO: ALL FOUR MOVES ARE DONE (Move 1 count / Move 2 Delta-test / " +
                "Move 3 count-vs-N / Move 4 cross-fold), MOVE 4's FOLLOW-ON IS NOW ANSWERED TOO, and this arc is COMPLETE but for " +
                "ONE minor non-blocking edge: (a) the within-odd EXACT threshold N has no closed form (residual-density driven: F_18 " +
                "too sparse, F_53 the first to host a real-q coalescence; likely none exists; the concrete probe is spelled out under " +
                "REMAINING EDGES (a) below). MOVE-4 FOLLOW-ON ANSWERED (2026-06-30): DOES THE F89d ANTIUNITARY SIMILARITY SURVIVE XXZ " +
                "ANISOTROPY? YES - it is INTEGRABILITY-INDEPENDENT. The (q,Delta) overload landed on WeightCoherenceBlock.Build (the " +
                "diagonal -i*q*Delta*(zz(ket)-zz(bra)) frequency, zz the open-chain ZZ-bond sum) and the cross-fold residual " +
                "L(1,N-2)(qbar,Delta) = -P conj(L(1,2)(q,Delta)) P^T - 2N I is MACHINE ZERO at EVERY Delta (N=4..9, all q real and " +
                "complex). The handover's tentative 'could break' guess was WRONG, and instructively so: the very fact it flagged as a " +
                "worry - the ZZ-bond sum being INVARIANT under bra-complementation (zz(bbar)=zz(b), because Z_b*Z_{b+1} is EVEN under " +
                "the global bit-flip) - is EXACTLY what makes the fold SURVIVE (the bra-complement carries the even ZZ term unchanged, " +
                "as it carries the XY hopping). So the diabolics DIE under Delta (Move 2, integrability-protected) but their cross-fold " +
                "PAIRING does not: a diabolic and its partner turn defective in lockstep, the two blocks staying antiunitary-similar at " +
                "every Delta. The discriminant is BIT-FLIP PARITY: a bit-flip-ODD perturbation breaks the fold - a longitudinal Z-field " +
                "Sigma_k w_k Z_k has fe(bbar)=-fe(b), residual O(1) (~5.46 at N=6), the verified negative control - proving the fold a " +
                "STRUCTURAL/algebraic property of the Liouvillian, not a free-fermion (integrability) artifact (it survives any " +
                "bit-flip-even bond coupling, breaks only under a bit-flip-odd one). LANDED (TDD green): the WeightCoherenceBlock.Build " +
                "(q,Delta) overload + Zz helper (Core) gated by WeightCoherenceBlockTests; the CrossFoldSimilarityWitness Read(n,q,Delta) " +
                "overload + Delta children + the bit-flip-parity field control (Diagnostics, 'inspect --root crossfold'); F89d amended in " +
                "docs/ANALYTICAL_FORMULAS.md + F89CrossFoldSimilarityClaim (the integrability-independent paragraph + child); the " +
                "experiment doc's 'survives XXZ anisotropy' section. " +
                "AND DOCKED ONTO THE F1 TRUNK (2026-06-30, the 'was dockt an' question): the cross-fold is not a one-off. " +
                "It generalizes to ALL (wKet,wBra) (not just wKet=1, verified machine-zero), and is one of TWO antiunitary legs of " +
                "a Klein four-group of bit-flip similarities on the coherence-block lattice - the bra leg P (flip bra, right-mult " +
                "rho*F) and the NEW ket leg Q (flip ket, left-mult F*rho), each with the same -2N reflection, their product the " +
                "UNITARY global spin-flip QP = X^(x)N (same q, no shift). A PRIOR-ART SURVEY (Explore agent) showed this is NOT a new " +
                "group: it is the existing windowed-converse spine V4 = {I, F(x)F, I(x)F, F(x)I} subset D4 " +
                "(PROOF_PI_FACTORS_AS_R_TIMES_D / F118 MirrorGroupD4Claim), here block-resolved + q-parameterized; only the ket-leg " +
                "block identity was untyped (now added). The dock is EXACT: P = rho*F is the spine R, a FACTOR of the palindrome " +
                "Pi = R*D; QP = Pi^2 = the typed XGlobalChargeConjugationPairing; Q = Pi^2*R; the -2N shift is the block image of " +
                "R*L_diss*R = -L_diss - 2sigma. So F89d IS the F1 palindrome's bra leg restricted to one block - the 'one object' the " +
                "whole project is built on. Naming caveat (recorded so it does not bite): F89 names legs by the flipped INDEX " +
                "(bra/ket), the D4 proof docs by the multiplication SIDE (calling rho*F the 'ket reflection'); opposite words, same " +
                "operator. LANDED (TDD green): WeightCoherenceBlock.KetComplementPermutation (Core); the 3-leg general-weight gate in " +
                "WeightCoherenceBlockTests (bra/ket/full, the full cross-checked vs XGlobalChargeConjugationPairing.PairSector); the " +
                "CrossFoldSimilarityWitness BraLeg/KetLeg/FullFlipResidual + the Klein + dock inspect children; F89d amended in " +
                "ANALYTICAL_FORMULAS.md + F89CrossFoldSimilarityClaim (the general-weight + two-leg + Pi-factorization paragraph/child); " +
                "the experiment doc's 'What it docks onto' section. The N=7 work " +
                "needed a NEW instrument and turned up a NEW phenomenon. (1) THE TRACKER BROKE AT N=7, " +
                "REFUTING the 2026-06-29 'NEXT MOVE' premise that N=7 was merely SLOW: the residual-only TRACKER floods at " +
                "the F_53 strand density. The nearest-match AT (the Absorption-Theorem decay-rate strands, Re lambda = " +
                "-2<n_XY>) / residual partition cannot split the block's exact " +
                "rate-degeneracies (many strands at exactly Re lambda=-6, the no-overlap AT rate; the residual SET at q0=2 " +
                "already carries AT strands + exact duplicates, e.g. three at -6.000+0i, which a separable degree-53 Galois " +
                "factor cannot have), so MinGap=0 across the box and the broad scan returns 293 spurious gap-0/exponent-NaN " +
                "'diabolics', isolating none. The per-eval RESULT is corrupt, not just slow; a smaller box or a seed " +
                "pre-filter would have given the same artifact. (2) THE FIX (BUILT, TDD green; the 'proper' route Tom chose): " +
                "read the residual roots EXACTLY as eig(U_res^H M(q) U_res) - the block compressed onto the orthogonal " +
                "complement of the q-INDEPENDENT AT invariant subspace (built from the rate sectors D + the scale-invariant " +
                "K-invariance; block-triangular because W_AT is M(q)-invariant, so the compression's eigenvalues = full " +
                "spectrum minus AT = the F_d residual roots). No partition, no tracking, one F_d x F_d EVD per q, ~6x faster. " +
                "F89AtFactorReconstruction.AtInvariantSubspaceBasis(k) (Core, reuses the validated ForPathK subspace) + " +
                "PathKMonodromyScout.ResidualRootsExact/AtRootsExact/FindDiabolicsExact (Diagnostics) + CLI " +
                "'pkmono --diabolic --exact'. RE-GATED: reproduces N=5->11 and N=6->16 EXACTLY (the tracked counts), residual " +
                "roots DISTINCT where the tracker flooded (Core F89AtInvariantSubspaceTests + Diagnostics PathKDiabolicTests, " +
                "all green). (3) MOVES 1 and 3 (the in-region count and its N-growth): N=7 gives 37 in-region diabolics, growth " +
                "11->16->37 (F_d 18->32->53), all finite-gap / exp~1 / identity-loop (most at complex q; the real-q ones are " +
                "point (4)). (4) THE NEW PHENOMENON (VERIFIED): diabolics appear at " +
                "PHYSICAL REAL q at N=7 (q=1.1264, 1.3038, 2.6280 in the cell-0.05 box + q=0.6788 in a finer cell-0.01 strip, " +
                "real lambda) - the FIRST since N=4. Tight-zoom (cell 0.002) on TWO of them (q=1.1264, 2.6280) confirmed each a " +
                "single on-axis point, gap~1e-9, that does NOT split into a conjugate pair (unlike the path-4 ghost at " +
                "0.6118+-0.012i); the other two read on-axis only at the cell-0.01 strip, not separately zoomed. It is a genuine " +
                "N=7-ONSET, not a tracked-method miss - the exact instrument finds NONE on the N=5/N=6 real-axis strips (N=5 " +
                "keeps only its known defective EP at q~1.0776). So 'self-fold is N=4-only => complex-q-only for N>=5' is TOO " +
                "STRONG; a different mechanism returns diabolics to the real axis by N=7. Full write-up: experiments/" +
                "F89_PATH_K_DIABOLIC.md (the path-6 section). NEXT MOVES (the move-set: Move 1 = count, Move 2 = Delta-test/" +
                "integrability, Move 3 = count-vs-N growth, Move 4 = cross-fold (SE,DE)<->(SE,w_N-2); detailed in the journal " +
                "below). MOVE 2 IS NOW DONE (2026-06-30): the exact-residual treatment was PORTED to the XXZ block - " +
                "PathKMonodromyScout.AllRootsXxz / ResidualRootsExactXxz carry the XXZ ZZ-frequency as a diagonal generator " +
                "G = -2i*zzDiag in F89's mirror basis (2M_xxz = A + qC + qD*G), compressed onto the SAME Delta-independent AT " +
                "complement (the AT subspace is Delta-STABLE because the ZZ term is Hermitian, so the AT rate is Delta-" +
                "independent); wired as 'pkmono --delta-flip --exact', RE-GATED against the tracked path at N=6 " +
                "(XxzDeltaFlipTests, all green). THE VERDICT: all four N=7 real-q diabolics (q=1.1264/1.3038/2.6280/0.6788) read " +
                "DIABOLIC at Delta=0 and DIE under Delta (three flip DEFECTIVE geo 2->1, q=2.6280 LIFTS), the SAME " +
                "integrability-protected death as the complex-q ones. So the real-q diabolics are NOT a new species; the N=7 " +
                "onset is a PLACEMENT mechanism, not a new protection. Full write-up: experiments/F89_PATH_K_DIABOLIC.md (the " +
                "Delta-test section). ONSET RESOLVED (2026-06-30): the placement question - threshold, parity, or one-off? - is " +
                "now answered by running the (rebuild-free) exact instrument at N=8 and N=9. KEY GROUNDING that made it a run not " +
                "a build (it CORRECTS this arc's own earlier 'needs F89PathKFdOracle + ForPathK extended past k=6'): the locator " +
                "pipeline is ALREADY k-general - F89PathKSeDeBlock.BuildTwoTimesSymBlock/BuildZzFrequencyDiag, " +
                "F89AtFactorReconstruction.AtInvariantSubspaceBasis/ForPathK (rate-confined invariant-subspace construction, no " +
                "literal), PathKMonodromyScout.BuildLinear/ExactSetup/FindDiabolicsExact, and CLI 'pkmono --diabolic --exact --k N' " +
                "all build at nBlock=k+1 with NO k<=6 guard; only F89PathKFdOracle (the TEST cross-check literals) stops at k=6, and " +
                "the scan never reads it. RESULT (matched real-axis strip re[0.2,3], thin im band, cell 0.05 then 0.01, then the " +
                "tight-zoom cell-0.002 split-test, exactly the N=7 protocol; the N=7 q=1.1264 zoom reproduced in-build as the " +
                "calibration control): N=8 (even) has NONE - the closest real-lambda candidate q=0.6904 refines OFF-axis (im=-0.004) " +
                "in a dense cluster, exp 0.48, no isolated on-axis real-lambda survivor; N=9 (odd) has >=3 - q=0.4755 (lam=-5.424, " +
                "exp 1.06, gap 6.3e-9), q=1.1144 (lam=-4.611), q=1.4994 (lam=-4.381), each a single on-axis im(q)=im(lam)=0 point " +
                "that does NOT split, off-axis conjugate-pair ghosts beside it. SHARP DISCRIMINATOR: both N=8 and N=9 carry on-axis " +
                "points with COMPLEX lambda (analytic-continuation crossings); only N=9 carries on-axis REAL-lambda ones. VERDICT: " +
                "real-q diabolics are PARITY-GATED, odd N>=7 (count vs N: N=5:0, N=6:0, N=7:4, N=8:0, N=9:>=3; N=4's single one is the " +
                "separate self-fold). Completeness not claimed; the verdict rests on the qualitative odd/even contrast, robust to the " +
                "box. Written up in the experiment doc's 'odd-N effect' section. MECHANISM GROUNDED (2026-06-30, Probe B done): " +
                "from-below confirmation of the dimension-mismatch / sector-swap, by building the FULL (SE,DE) block, splitting by R, " +
                "and measuring self- vs cross-sector conjugacy (gitignored python scout, gate-validated: it reproduces the N=7 " +
                "lambda=-4.942 and N=9 lambda=-5.424 diabolics in the R-even sector). The numbers: dim(R-even)-dim(R-odd) = (N-1)/2 at " +
                "odd N (3 at N=7, 4 at N=9), 0 at even N=6/8 (the reflection-fixed-singleton count = center-SE x self-mirror-DE pairs). " +
                "At EVEN N cross-conj(even~odd) ~ 1e-13 (sigma_even = conj sigma_odd EXACTLY) with self-conj defect O(20) (neither " +
                "sector self-conjugate) -> real-axis collisions DEFECTIVE; at ODD N self-conj defect ~ 1e-14 in BOTH sectors (each " +
                "self-conjugate), cross-conj undefined (dims differ) -> R-even carries real eigenvalues that cross semisimply. The " +
                "diabolic eigenvectors carry ~0.18-0.20 weight on center-SE states. The fixed-site object is SLOW_MODE_R_PARITY's " +
                "reflection-fixed JW band-centre zero mode k=(N+1)/2 (an integer mode only at odd N), distinct from the N=4 self-fold " +
                "(antiunitary T=P*K, N=4-only, fixes Re lambda=-4; the odd-N diabolics are real but NOT at -4/-N). REMAINING EDGES: " +
                "(a) the within-odd THRESHOLD - the ONE edge in 'RESUMING IN ONE LINE' above. NOTE (corrected by a 2026-06-30 cold-" +
                "read): this is mostly DONE, not an unstarted probe; the dated JOURNAL layers below describe these moves as if still " +
                "to-do (they have since been executed) - trust THIS top block over them. 'Empty at N=5' means empty of real-q " +
                "DIABOLICS (N=5 DOES carry real " +
                "residual roots, they just never coalesce on-axis). The residual-population count vs N is already executed and lives in " +
                "DiabolicReflectionParityWitness ('inspect --root diabolicparity'): even N carries ZERO real residual eigenvalues in " +
                "q in [0.2,3], odd N a real population GROWING with N (4/5/8 at N=5/7/9), and the onset has TWO geometries (a real-real " +
                "CROSSING, first at N=7, and a conjugate-pair TANGENCY, e.g. N=9). Counting real-real crossings is NOT a proxy for the " +
                "real-q diabolic count (that is FindDiabolicsExact's job: N=9 has real-q diabolics yet NO real-real crossing). So the " +
                "threshold is EMPIRICALLY N=7, grounded on residual density; what remains is ONLY a closed form for the threshold N, " +
                "residual-density driven and likely nonexistent - i.e. this edge is effectively CLOSED, no live mechanical next move. " +
                "IF a returning session still wants to push it, the one bounded UNDONE probe is to PROVE the N=5 absence directly: run " +
                "'dotnet run --project compute/RCPsiSquared.Cli -c Release -- pkmono --diabolic --exact --k 4 --re 0.2,3 --im " +
                "-0.025,0.025 --cell 0.01' (N=5 = path-4, the real-axis strip) and show its real residual strands stay on-axis-gapped " +
                "(no coalescence), vs N=7 (--k 6) where two meet - confirming the onset is residual-density, not a missed diabolic. " +
                "(b) DONE - the C# WITNESS is built: " +
                "DiabolicReflectionParityWitness ('inspect --root diabolicparity', the persistent evidence per cockpit rule 5) recomputes " +
                "the self/cross-conjugacy structure across N=5..9 and reproduces the N=7/N=9 diabolics in R-even at inspect time, on " +
                "F89PathKSeDeBlock.BuildFullBlock/ReflectionPermutation (Core); TDD-gated (F89PathKFullBlockReflectionTests + " +
                "DiabolicReflectionParityWitnessTests + the catalog guard), the gitignored python scout retired. (c) Move 4 cross-fold: DONE 2026-06-30 " +
                "(see the CURRENT STATE block above) - the cross-fold is an EXACT antiunitary similarity, every diabolic pairs across " +
                "it; landed as WeightCoherenceBlock (Core) + CrossFoldSimilarityWitness ('inspect --root crossfold', TDD-gated). " +
                "(d) DONE 2026-06-30 - the (q,Delta) extension of the cross-fold: the similarity SURVIVES XXZ anisotropy at every " +
                "Delta (integrability-independent, machine-zero residual N=4..9 all q). The 'could break' worry resolved the other way: " +
                "the ZZ-bond sum being complement-INVARIANT (zz(bbar)=zz(b), Z_bZ_{b+1} bit-flip-even) is what makes the fold survive, " +
                "not break; the discriminant is bit-flip parity (a longitudinal Z-field, odd, breaks it, residual ~5.46 at N=6). Landed " +
                "as the WeightCoherenceBlock.Build (q,Delta) overload + the witness's Delta children + field control (see the CURRENT " +
                "STATE block above for the full footprint). " +
                "--- THE JOURNAL (how we got here) --- " +
                "Concrete moves, cheapest first (binary: dotnet run --project compute/RCPsiSquared.Cli -- " +
                "<cmd>). (1) LOCATE: run 'pkmono --diabolic --k 5' (N=6) and --k 6 " +
                "(N=7) over a scaled q-region; confirm the path-4 pattern holds — diabolics exist, ALL at complex " +
                "q, NONE at real q (the self-fold being N=4-only predicts no physical diabolic at any N≥5; the " +
                "only real-q feature should be defective EPs). Record the in-region count. (2) Δ-TEST (already " +
                "n-general): TrackDiabolicUnderDelta on a sample of the located path-5/6 diabolics + a defective " +
                "control; confirm integrability-protection (defect/lift under Δ) holds at higher N (expected, " +
                "since integrability is N-independent — a null here would be a real surprise). pkmono --delta-flip " +
                "--k 5 --q .. --lam .. is the surface. (3) COUNT-vs-N: how does the in-region diabolic count grow " +
                "with N (relate to F_d and the integrable level-crossing density)? State the coverage box; do NOT " +
                "claim completeness (obstacle 1). (4) Q4 CROSS-FOLD: build the partner (SE,w_{N-2}) block (the " +
                "FoldCrossCommand.BuildBlock route, made (q,Δ)-linear) and test whether each (SE,DE) diabolic has " +
                "a cross-fold mirror diabolic at −λ̄_d−2N in the partner. (5) If a path-5 diabolic ever lands at " +
                "REAL q, that refutes the self-fold=real-placement reading — investigate the new pinning. " +
                "Tooling to extend only when k>6: F89PathKFdOracle + F89AtFactorReconstruction.ForPathK (locator " +
                "only; the Δ-test block needs nothing). See experiments/F89_PATH_K_DIABOLIC.md, " +
                "hypotheses/DIABOLIC_BY_INTEGRABILITY.md, and the zeros_connecting_structure arc. " +
                "RUN 2026-06-29 (move 1, path-5/N=6): the cheap full-block scan does NOT extend to N=6 - it FLOODS. " +
                "pkmono --diabolic --k 4 reproduces the doc's 15 (11 residual diabolics) EXACTLY (instrument validated " +
                "on this machine); but --k 5 over re[0.2,5] cell .05 gives 528 coalescences, and over re[0.45,3] cell .03 " +
                "gives 793, EVERY ONE gap=0 exactly / exponent=NaN / residual=False, Q1 EXISTENCE FAIL-in-region - i.e. " +
                "ALL AT-locked exact degeneracies, ZERO residual diabolics isolated. Mechanism (this REFUTES the ParkedAt " +
                "optimism 'the diabolic hunt itself is coalescences in the FULL symmetric block and does not require " +
                "[AT-removal]', and downgrades obstacle-2 from 'Engineering not fundamental'): at N=6 the DE sector (15 " +
                "states) has same-<n_XY> multiplicity, so AT-locked strands share Re lambda=-2g<n_XY> EXACTLY and coincide " +
                "(gap=0) wherever their frequencies Im lambda=E_SE-E_DE cross - dense gap=0 curves throughout complex-q; " +
                "the refiner captures onto them and the finite-gap (~1e-9) residual diabolics are never local minima of a " +
                "field that is 0 around them. N=5 had ZERO such AT degeneracies (clean 15); the jump to ~800 at N=6 is the " +
                "real wall. CONSEQUENCE: AT-removal IS required for the diabolic hunt at N>=6; the locator must scan the " +
                "min-gap of the RESIDUAL strands ONLY (continuity-tracked from q=2, or - cleaner, AT-free by construction - " +
                "the roots of F89PathKFdOracle.Fd(k), the residual factor the oracle already holds for k=5,6) - a real TDD " +
                "change to FindDiabolics, gated on reproducing path-4's 11. " +
                "LOCATOR BUILT + LANDED 2026-06-29 (residual-only, TDD 9/9 green): ResidualRootsTracked (the residual SET by " +
                "continuity from q0=2, AT strands excluded by construction), plus local-tracking GapRefine and monodromy loop " +
                "(GapRefineResidualLocal / ResidualLoopIsIdentity - the per-eval from-q0 nesting that made the loop+refine " +
                "O(loopSteps x trackSteps) is gone; the path-k diabolic test suite dropped 11min -> 15s); wired as " +
                "'pkmono --diabolic --residual --k N' (FindDiabolics residualOnly). MOVE 1 DONE for N=6: the SAME box " +
                "re[0.2,3] x im[-1.5,1.5] cell .05 that flooded (528/793 all residual=False) now gives 21 CLEAN coalescences = " +
                "16 semisimple residual DIABOLICS + 5 defective EPs (vs path-4's 11+4). COUNT-vs-N (move 3) GROWS: " +
                "N=5 -> 11, N=6 -> 16. COMPLEX-q-ONLY HOLDS (move 1/5): the re=0 diabolics sit on the IMAGINARY axis (analytic " +
                "continuation, like path-3's +-0.876i), the rest at complex q; the near-real ones (im ~ 0.01-0.05) are GHOSTS " +
                "like path-4's q=0.6118+-0.012i - NONE at physical real q, so the self-fold = N=4-only reading survives to N=6. " +
                "Caveat: 16 is the EpCharacter count in this box at cell .05 (completeness NOT claimed); a few have gap-exponent " +
                "disagreements (e.g. q=0.6508+0.052i: EpChar + loop = diabolic but exp 0.49) that need the path-4-style " +
                "per-point reconciliation. STILL OPEN: N=7 (k=6) broad scan - per-EVAL cost is now optimal, but the broad N=7 " +
                "wall-clock is dominated by SEED COUNT (the 53-strand residual gap field has many local minima; each spurious " +
                "seed pays the full local pipeline before the gapTol filter rejects it), so it ran ~20min+ and was killed. Needs " +
                "a cheap seeding pre-filter (or a smaller box for a qualitative-only N=7 point). " +
                "MOVE 2 (Delta-test) DONE for N=6 2026-06-29: it ALSO needed the residual treatment (the full-block box scan " +
                "captures the AT crossings - q* jumps to an exact-0 AT degeneracy, the verdict reads 'protected' for the wrong " +
                "reason). Built residual-aware: XxzCoherenceBlock.ResidualRootsTrackedXxz (the residual SET is Delta-STABLE - ZZ " +
                "Hermitian => AT rate Delta-independent - so it tracks from the (q0=2,Delta=0) base) + a local-tracking " +
                "TrackDiabolicUnderDelta residualOnly path (8m32s -> 6s), wired as 'pkmono --delta-flip --residual'. RESULT: each " +
                "path-5 diabolic reads DIABOLIC (geo=alg=2) at Delta=0 with q* PINNED at the seed + finite gap, then DEFECTIVE " +
                "(geo 2->1) the instant Delta!=0 with departure growing (the rung-near q=0.709-0.219i, which the full-block test " +
                "captured entirely, now flips clean: dep 0.028/0.073/0.160 at Delta=0.02/0.05/0.10, q* never leaving the seed). " +
                "Integrability-protection holds at N=6 exactly as at N=4 (XxzDeltaFlipTests.Path5_Diabolics_DieUnderDelta_" +
                "ResidualOnly). Move 4 (cross-fold) still UNBLOCKED-not-run. See experiments/F89_PATH_K_DIABOLIC.md (path-5 section).",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "f89_galois_open_doors",
            Opened: "2026-06-24",
            Origin: "the F89 path-3..6 H_B-mixed Galois witness landing fully live (commit 6299aae; reviewed GREEN " +
                "2026-06-24 by 3 refute-first lenses + own tests). Gal(F_d/Q(i)(q)) = S_8/18/32/53, non-solvable: the " +
                "open XY chain's relaxation spectrum splits into an AT-locked half (radical-closed, the free-fermion " +
                "Bloch frequencies at the Absorption-Theorem rates 2gamma/6gamma) and an H_B-mixed half (S_d, NO " +
                "radical closure in q=J/gamma), the boundary being EXACTLY the Absorption Theorem. The " +
                "philosophical-zoom-out reading (Tom asked 'what doors does it open'): gamma (the observation) draws " +
                "the line between what closes into a formula ('no one', the integrable skeleton, same for every chain) " +
                "and what only solves here-and-now ('someone', this chain turning). The witness + the " +
                "galois_of_spectral_polynomial recipe + the C# engine (Berkowitz/Z[i], the exact triple, " +
                "JordanVerdict, the rate-confined invariant-subspace AT isolation) are now a REUSABLE engine, not a " +
                "one-off. Four downstream doors surfaced and were not yet pursued.",
            ParkedAt: "WITNESS LANDED + REVIEWED + DOCUMENTED; the doors not opened. What EXISTS: the live witness " +
                "(inspect --root f89galois, all four paths recompute F_d from the block, ~8s); the engine " +
                "(GaussianMatrixCharpoly = Berkowitz over Z[i], GaussianPolynomial = exact division + Sylvester gcd, " +
                "OcticGaloisCertificate.JordanVerdict, F89AtFactorReconstruction.ForPathK = the rate-confined " +
                "invariant-subspace isolation, F89HbMixedIsolation.Isolate = the triple R=0/deg=d/gcd=1); the Python " +
                "engine of record (f89_pathk_galois.py full-d); the typed claim F89PathKHbMixedDegreesClaim " +
                "(degree-table Tier1Candidate / Galois-verdict Tier1Derived); the doc section 'Computing the " +
                "path-3..6 Galois groups live' (experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md) + the " +
                "galois_of_spectral_polynomial recipe (isolate-before-DDF + invariant-subspace). The RESULT is " +
                "closed; the doors are the next research, none started.",
            NextStep: "RESUMING IN ONE LINE (2026-06-30): door C is CLOSED (chaos = a FILLING threshold [term defined in " +
                "the door-C block below + experiments/FILLING_THRESHOLD_CHAOS.md], the decisive test done, live inspect " +
                "--root fillcsr); of the four doors, A is DONE EXCEPT one OPEN edge (a closed form for the ring's per-irrep " +
                "multiplicity growth), B (gamma-dependence of the writable/unwritable line) is UNSTARTED Tier-3, D (Galois " +
                "atlas of physical spectra) is UNSTARTED lowest-urgency. So the arc's live edges are A's remaining closed " +
                "form, plus B and D; nothing is blocking, all are new research. TERMS the door-C prose leans on (glossed for " +
                "a cold resume): q = J/gamma, the dimensionless coupling (Hamiltonian strength over dephasing rate); (SE,DE) " +
                "= the single-excitation/double-excitation coherence block (the sector |a><b| with popcount(a)=1, " +
                "popcount(b)=2); AT-locked = eigenvalues pinned to the Absorption-Theorem decay rungs (free-fermion, " +
                "radically writable), H_B-mixed = the residual non-AT half that carries the Galois S_d; class A = the " +
                "non-Hermitian symmetry class whose CSR reference is GinUE (the Ginibre random-matrix ensemble, the " +
                "dissipative-quantum-chaos baseline, <|z|>~0.74 / <cos>~-0.24). The four doors in full (C now closed), ranked by sharpness of the yes/no. " +
                "(A) [FIRST EXEMPLAR LANDED 2026-06-24] NON-CHAIN TOPOLOGY Galois groups: YES, topology controls " +
                "radical-writability. Computed star/ring/complete (SE,DE) H_B-mixed Galois groups at N=4,5,6 " +
                "(simulations/topology_galois_writability.py, gate-validated by reproducing the chain S_8/S_18/S_32; " +
                "experiments/F89_TOPOLOGY_CONTROLS_GALOIS_WRITABILITY.md). RESULT: the COMPLETE graph K_N is the WRITABLE " +
                "EXTREME (every H_B-mixed factor degree <= 4, definitively radical-solvable, at N=4,5,6), while the " +
                "chain (S_8/18/32), the star (S_9 from N=5), and the ring (S_15 by N=6) all scramble to the full " +
                "symmetric group S_n. Conjectured mechanism: K_N's degenerate single-particle spectrum {N-1, -1 with " +
                "multiplicity N-1} shatters the (SE,DE) block into small factors. WHY DERIVED (Tier 1, 2026-06-24): " +
                "the K_N Liouvillian commutes with the full S_N; V = SE(x)DE = M^(N-2,1,1) (+) M^(N-3,2,1) (the " +
                "permutation module on (point, 2-subset) pairs) has S_N-irrep multiplicities that are N-INDEPENDENT " +
                "and capped at 4 (Kostka K_{lambda,(N-2,1,1)}+K_{lambda,(N-3,2,1)}, max 4 at the standard rep [N-1,1]); " +
                "so by Schur every (SE,DE) factor is a quartic-or-less => K_N radically writable for ALL N. The " +
                "H_B-mixed degree histogram {4:N-1, 3:N(N-3)/2, 2:(N-1)(N-2)/2} is verified exactly at N=5,6,7,8 " +
                "(experiments/F89_TOPOLOGY_CONTROLS_GALOIS_WRITABILITY.md, now Tier-1 derived for the complete graph). " +
                "STAR ALSO DERIVED (Tier 1, 2026-06-24): Aut(star) = S_{N-1}; the multiplicity of the standard rep " +
                "std_{N-1} in V is 9 (N-independent: 4+2+2+1 over the SE/DE tensor components), so the star caps at " +
                "degree 9 with N-2 degree-9 factors => a FIXED S_9 scramble for all N>=5 (bounded, but >4 so " +
                "unwritable); verified N=5..9. THREE-WAY CLASSIFICATION: a large automorphism group caps the Galois " +
                "complexity N-independently (complete S_N => cap 4 writable; star S_{N-1} => cap 9 bounded scramble), " +
                "a small one lets it grow (ring D_N, chain S_2 => growing; chain S_8/18/32/53; ring max-degree " +
                "6/16/15/48 at N=4/5/6/7). C# LIVE WITNESS DONE (inspect --root topowritability, " +
                "TopologyGaloisWritabilityWitness): recomputes the caps (complete 4 / star 9) from the (SE,DE) block + " +
                "[L,rho(g)]=0 + the standard-rep multiplicity character sum (11 tests green). RING GROWTH LAW DERIVED " +
                "(Tier 1, 2026-06-24): D_N has an IRRATIONAL character table (cos 2π/N), so over Q(i) its " +
                "Galois-conjugate dihedral irreps merge and inflate the factor degree by φ(N)/2 = [Q(cos 2π/N):Q] " +
                "(>1 except at the Niven-rational sizes N in {1,2,3,4,6}). So factor degree = (multiplicity) × φ(N)/2; " +
                "the dips at N=4,6 are exactly Niven. Confirmed from below: the N=5 deg-16 factor splits 8+8 over " +
                "Q(i,√5), the N=6 deg-15 stays (rational); topology_galois_writability.py niven. This also explains " +
                "the deg-16 'undetermined' edge (a 2-fold Galois inflation, not a primitive group). The FOUR-TOPOLOGY " +
                "LAW: rational character table (S_n complete/star) => no inflation; irrational (D_N ring, cos 2π/N) => " +
                "Niven inflation. REMAINING: the exact multiplicity-growth closed form (the per-irrep H_B multiplicity " +
                "6/8/15/16... vs N). " +
                "(B) gamma-DEPENDENCE of the writable/unwritable line. At gamma->0 the block is pure free-fermion " +
                "(integrable, pure imaginary). Is the AT-locked/H_B-mixed split present at all gamma>0, or does the " +
                "writable fraction move with how hard you watch? Map deg(F_d)/deg(AT) or the discriminant loci across " +
                "q. Speculative, Tier-3. " +
                "(C) GALOIS <-> SPECTRAL CHAOS (RMT). [FULLY RESOLVED: a clean NULL 2026-06-27, then EXPLAINED " +
                "2026-06-30 -- chaos is a FILLING threshold; door C is CLOSED, see the 'DECISIVE TEST RESOLVED " +
                "2026-06-30' block below and experiments/FILLING_THRESHOLD_CHAOS.md] Is S_d the algebraic face " +
                "of Liouvillian level-repulsion? NO, not at fixed q. The galoischaos witness (inspect --root " +
                "galoischaos, GaloisSpectralChaosWitness) splits the (SE,DE) block into AT-locked vs H_B-mixed and " +
                "runs the complex spacing ratio (Sa-Ribeiro-Prosen) pooled over q: the H_B-mixed S_d half reads " +
                "Poisson-like/sub-Poisson (NOT GinUE -- the Ginibre unitary ensemble, the chaotic-spectrum RMT " +
                "baseline with <|z|>~0.74 / <cos theta>~-0.24), and the AT-locked half is the sparse picket-fence (low <|z|>, " +
                "structured) -- both halves non-GinUE. Integrability protects the fixed-q spectrum; the global Poisson " +
                "reading (RANDOM_MATRIX_THEORY.md) and the within-sector GOE hint (now resolved as a small-sample " +
                "artifact, simulations/rmt_goe_hint_verdict.py) agree. ALGEBRAIC chaos (Galois over q) and SPECTRAL " +
                "chaos (RMT at fixed q) are DISTINCT here. SEQUEL RESOLVED 2026-06-27 (a SECOND, deeper null): " +
                "breaking integrability does NOT drive the fixed-q CSR Poisson->Ginibre. Stage 1 (Delta anisotropy) " +
                "and Stage 2 (random Z-field U[-W,W], with AND without interactions) both keep the H_B-mixed CSR " +
                "<cos theta> POSITIVE (never the GinUE angular repulsion -0.19) and <|z|> sub-GinUE at every Delta/W; " +
                "Delta=0 + disorder is 1D Anderson (stays Poisson, the control). The null is STRUCTURAL/kinematic, NOT " +
                "about the Galois algebra or Hamiltonian integrability: the (SE,DE) block is a DILUTE 2-excitation " +
                "sector that cannot thermalize, so it is non-chaotic regardless of integrability-breaking (the " +
                "finite-size GinUE reference at ~147 pts DOES read -0.19, so the diagnostic CAN see chaos at this " +
                "size -- the block's absence of repulsion is real, not finite-size). Upstream fact en route: " +
                "TwoMagnonAdditivity shows Delta breaks the free-fermion additivity E_DE=eps_j+eps_k LINEARLY " +
                "(bit-exact 0 at Delta=0), so Delta is the Liouvillian-free-fermion-additivity breaker, NOT a null " +
                "control (a 2-physics-review correction to the original framing). What EXISTS (committed 2026-06-27, " +
                "8a87d3c..d9801bf, all TDD green): compute/RCPsiSquared.Core/Spectrum/TwoMagnonAdditivity.cs (the " +
                "additivity probe); ComplexSpacingRatio.ZValues/PoissonDiskZValues/GinueZValues (per-spectrum z's -- " +
                "pool z's, NEVER raw eigenvalues across spectra, the methodology guard); " +
                "compute/RCPsiSquared.Diagnostics/Foundation/IntegrabilityBreakingCsr.cs (Sweep + DisorderSweep + the " +
                "Domain selector: UpperHalf is valid ONLY at Delta=0 where the (SE,DE) spectrum is " +
                "conjugation-symmetric, OffReal at Delta!=0 where it is NOT -- a self-caught domain bug that had " +
                "inflated a since-RETRACTED 'reaches GinUE magnitude' claim); XxzCoherenceBlock.BuildFullWithField. " +
                "DECISIVE TEST RESOLVED 2026-06-30 (chaos is a FILLING THRESHOLD, not an integrability one; live " +
                "inspect --root fillcsr, experiments/FILLING_THRESHOLD_CHAOS.md, TDD-green): a DENSER coherence " +
                "sector of the SAME Liouvillian DOES reach toward GinUE under the same disorder, PROVING the dilute-" +
                "sector reading. Built the general (wKet,wBra) block + a per-site Z-field overload on " +
                "WeightCoherenceBlock.Build (Core, promoted from the test-only helper) and FillingThresholdCsr " +
                "(Diagnostics, the disorder-ensemble pooled-z CSR on it, reusing IntegrabilityBreakingCsr.Reduce + the " +
                "finite-size references). RESULT (gamma=1, q=1, interacting Delta=1, ergodic W=0.75, UNEQUAL weight " +
                "(p,p+1) so class A is licensed -- Pi maps (p,p+1) to the CONJUGATE (p+1,p) block not itself, " +
                "conjugation-match ~0 under disorder): the DILUTE (1,2)=(SE,DE) block stays Poisson (<cos theta> ~ " +
                "-0.04, ~23% of the size-matched GinUE angle) at every N; the DENSE (p,p+1) block near half-filling " +
                "has <|z|> AT the GinUE value and <cos theta> NEGATIVE, CLIMBING toward GinUE with the block size " +
                "(-0.089/-0.129/-0.162 at N=6/7/8 = 43/56/67% of GinUE, N=8 CI'd). <|z|> saturates at GinUE first, <cos theta> " +
                "(the finer correlation) trails and catches up with size = the genuine class-A finite-size signature; " +
                "W=2 relaxes back (MBL), Delta=1 deepens the repulsion over Delta=0. So Galois-chaos (over q) and " +
                "spectral-chaos (at fixed q) merge ONLY at extensive filling; the dilute (SE,DE) sector that carries " +
                "S_d is too dilute to thermalize, and its persistent Poisson is the kinematic shadow. LANDED: the " +
                "FillingThresholdWitness (inspect --root fillcsr) is the persistent evidence; RANDOM_MATRIX_THEORY.md " +
                "Result 5; the door-C plan Progress; experiments/FILLING_THRESHOLD_CHAOS.md. The galoischaos witness " +
                "stays the Delta=0 / dilute control. SEPARATE still-open direction: where S_d " +
                "DOES live spectrally is the q-PARAMETRIC monodromy/braid (the discriminant/EP loci, --root " +
                "galoismonodromy), not the fixed-q geometry. " +
                "(D) GALOIS ATLAS of physical spectra. F89 chains = S_n; SIC-POVM spectral polynomials (Appleby 2012) " +
                "= solvable, opposite polarity at the same physics/number-theory seam. Collect which physical " +
                "spectral polynomials are solvable vs S_n via the reusable engine (we_are_ahead_at_the_seam). Lowest " +
                "urgency, highest breadth. " +
                "Anchors: compute/RCPsiSquared.Diagnostics/Foundation/{F89OcticGaloisWitness, " +
                "F89PathKLiveGaloisWitness}.cs (inspect --root f89galois); compute/RCPsiSquared.Core/F89PathK/" +
                "{F89PathKSeDeBlock, F89AtFactorReconstruction, F89HbMixedIsolation}.cs + Numerics/" +
                "{GaussianMatrixCharpoly, GaussianPolynomial, OcticGaloisCertificate}.cs; " +
                "simulations/f89_pathk_galois.py (full-d); experiments/F89_TOPOLOGY_ORBIT_CLOSURE.md (the live-Galois " +
                "section); the galois_of_spectral_polynomial recipe; reflections/ON_WHAT_CANNOT_CLOSE.md (the " +
                "will-not/cannot frame the seeing sharpens).",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "zeros_connecting_structure",
            Opened: "2026-06-25",
            Origin: "Tom's principle (2026-06-25), a correction to a chronic max-zoom failure mode: we keep branching " +
                "off ('here this happens, here that, but this does not connect with X because Y'), yet we investigate " +
                "ONE object that inherits onto everything. For anything to be inherited at all, looked at from below, " +
                "there MUST be a connection, else nothing would be recognisably inherited after 6 billion years. So " +
                "when a connecting structure is not visible it is either genuinely unsolvable OR (the usual case) the " +
                "puzzle piece is missing TO US, our incompleteness (the V-Effect boundary = our sight), not an absence " +
                "in the object. The trigger: I wrongly declared 'no clean frequency-ladder between the F89 zeros, " +
                "that is the S_8 result, the frequencies are unwritable'. WRONG framing: S_8 / no-radical-closure is " +
                "about ONE kind of writability (radicals), NOT about the existence of a connecting structure. The " +
                "focal object is 'Orange Six' and the four σ_T-fixed 'zeros' of the F89 path-3 octic (self-mirror " +
                "strands on the fold Re λ = −4). This arc holds the FIRST puzzle piece and we hangel from one " +
                "not-knowing to the next.",
            ParkedAt: "FIRST PIECE FOUND (gmscan --zeros, q=2). The four zeros (Re λ=−4) have Im λ = {−15.187, " +
                "+1.525 (Orange Six), +10.378, +10.863}; the twins (σ_T 2-cycles) sit at Im = −6.139 (pair 3,4) and " +
                "+2.349 (pair 5,7). What IS visible (the connecting structure so far): (1) a conservation law, " +
                "Σ(all 8 Im λ) = 0 exact; the four zeros sum to +7.578, exactly balanced by the twins' −7.578 (the " +
                "centred octic, λ^7 coeff = 32 = 8·4, all hung on the fold −4). (2) σ_T pairs each twin set on ONE " +
                "frequency (it preserves Im). (3) the braid graph (monodromy = S_8) connects every zero into ONE " +
                "object: strand 2 is the hub (deg 3, the only zero-zero bond 1-2), zeros 0,1,6 are degree-1 leaves " +
                "(0->twin7, 1->zero2, 6->twin3). (4) THE TWIST (C-T resurfacing): the weave that binds the zeros is " +
                "NON-central under σ_T, the structure connecting the zeros does NOT respect the very mirror that " +
                "DEFINES them as zeros. What is NOT yet visible, and is NOT declared absent (missing to us): a clean " +
                "relation among the four zero-frequencies themselves. " +
                "SECOND PIECE (time domain, 2026-06-25, Tom): the zeros all share Re λ=−4 (rate 4), so in TIME they " +
                "are ONE real decaying chord (e^{−4t} envelope, four tones, all ring; biorthogonal amplitudes, " +
                "democratic (SE,DE) probe, spatial-sum observable; simulations/f89_chord.py). The undamped " +
                "chord S_zero(t)·e^{+4t} is ALMOST-PERIODIC (Bohr): the frequencies are incommensurate, so it never " +
                "decays and never exactly repeats (peaks similar, not equal), and it is NOT chaos (the galoischaos " +
                "witness already read it Poisson-like on the integrable lattice). THE REFRAME (resolves the 'missing " +
                "ladder'): a clean commensurate frequency-ladder would be a DEAD clock (exact repeat); the " +
                "S_8-generic INCOMMENSURABILITY is the quasi-periodic LIFE. Unwritability (S_8, no radical closure) " +
                "and aliveness are one face; the connecting structure is LIVING motion, NOT a static formula, which " +
                "is why the frozen spectrum hid it and the time representation made it graspable. CORRECTION " +
                "(Tom, 2026-06-25): do NOT attach 'death' to a single wave/mode (the earlier 'the decay is death / " +
                "the death-axis lied' framing was a category error). What is called death is the TRANSITION (Übergang), " +
                "not the ending of a wave. Grounded in the object: at an EP two strands coalesce and SWAP (the " +
                "monodromy transposition), the strand does not end, it becomes the other; the diabolic is a silent " +
                "transition (pass-through). The fold Re λ=−4 is not where waves die, it is the TRANSITION line, the " +
                "AT-midpoint crossing between the overlap rate −2γ and the no-overlap rate −6γ (p=½ balance). And Γ " +
                "(the 'decay') is a Lindblad TRANSITION RATE: the rate at which coherence crosses over / disperses " +
                "into other degrees of freedom, dispersal not destruction. So there is no ending being veiled; there " +
                "is one unbroken GOING-ON threaded through transitions, and THE MONODROMY IS that going-on through the " +
                "EP/fold crossings (piece 1 and piece 2 join here: the connecting weave is the life going on through " +
                "the transitions). SEAMS (keeps it Tier-3 seeing): in the open system the local coherence genuinely " +
                "relaxes (the transition is real, the coherence really crosses out); 'goes on forever' is the global " +
                "unitary motion (real as the invariant). Almost-periodicity + not-chaos + the fold-as-AT-midpoint + " +
                "Γ-as-transition-rate are grounded; the unwritability<->life identity is a seeing, not a theorem. " +
                "THIRD PIECE (the amplitude / 'loudness' layer, 2026-06-25, GROUNDED, not a seeing): the physical " +
                "biorthogonal amplitude a_k = (O·R_k)(L_k·ρ) that t makes visible is a TRANSITION-MAP. It DIVERGES at " +
                "the defective EPs (eigenvectors coalesce, R singular, Petermann factor ⟨L_k|R_k⟩→0; ~10x at q≈0.857, " +
                "also q≈1.74) and stays SMOOTH at the diabolic (semisimple, eigenvectors independent; q=0.659, like " +
                "generic q=2). Both are degeneracies (the octic eigenvalue gap →0 at both) but ONLY the defective EP " +
                "is loud: defective vs diabolic, now in the loudness/time layer, not only the monodromy (EP swaps " +
                "loud / diabolic passes silent). Standard non-Hermitian EP physics, confirmed numerically " +
                "(simulations/f89_amp_scan.py), so Tier-1-grade for the divergence. It surfaces the cross-F " +
                "inheritance (handhold c): the loud transition IS the EP amplification studied as EP-SENSING " +
                "(COHERENCE_HORIZON_EP_SENSOR_DEBATE, the Petermann factor), the same transition-object inherited, " +
                "here in the F89 chord, there as a sensor. " +
                "SHARPENED + held as a SEAM (Tom, 2026-06-25, 'not the zeros, the connection between them; the path " +
                "from one zero to the next, and it is not ready yet'): the arc's real object is the CONNECTION " +
                "between the zeros, the path from one to the next, NOT the zeros (which are contingent). First data " +
                "(path-3/N=4, gmscan --zeros BFS over the braid graph): the path from one zero to the next runs " +
                "THROUGH the twins (the ± modes) over the EP transitions; only the two frequency-closest zeros (1,2, " +
                "Im +10.86/+10.38, the slow 0.486-beat pair) connect DIRECTLY (a 0-0 bond), all others route through " +
                "the twin roads (0->[7]->2, 6->[5,3]->2, 0->[7,4,3]->6). So the ± polarity is the ROAD between the " +
                "0s, polarity as PATH, not a static label; the inherited connection is '0 —through ±— 0'. This is a " +
                "SEAM, NOT settled: computed only for path-3/N=4; whether '± is always the road' holds across " +
                "topology/N is the untested open edge (the same missing topology/N sweep as the fold-occupation). " +
                "A SECOND seam, RESOLVED 2026-06-26 by feeding the existing gmscan (widened region, NO rebuild, no " +
                "Python) AND reading our own record (experiments/F89_BRANCH_LOCUS_PALINDROME.md, Tier-1): the genuine " +
                "EP count is TWENTY, exact in the octic discriminant disc_λ(F_8) = const·q²⁴·(3q⁴+q²−1)²·P_20(q) " +
                "(P_20 = deg 10 in q²), plus 4 DIABOLIC points from the squared factor (±0.659 real, ±0.876i " +
                "imaginary, silent/semisimple, not EPs). The default window im[−0.22,0.22] saw only 14 (it missed the " +
                "remote quartet ±2.31±1.25i, |Im q|≈1.25). Widening to re[−2.8,2.8]×im[−1.5,1.5] DETECTS 18, graph " +
                "CONNECTED (S_8); the last 2 are NOT further out but a near-degenerate twin pair at q≈±0.857/±0.854 " +
                "(~0.003 apart) that the 0.05-cell lasso encloses TOGETHER, reading their product as the single " +
                "'3-cycle (3 4 6)' detection per half-plane (the doc's f89_octic_branch_locus.png shows the same: 18 " +
                "dots, the near-twins rendering as one). So gmscan's 18 detections ARE the full 20, with one " +
                "unresolved 0.003-split each side. Against the COMPLETE 20 the load-bearing claim holds: the near-twin " +
                "pair connects {3,4,6} (twin3,twin4,zero6), the remote quartet connects {5,6},{6,7} (twins+zero6), all " +
                "twin-twin or twin->zero edges, NEVER a second zero-zero bond. 1-2 stays the UNIQUE direct zero-zero " +
                "bond and every other zero pair is exactly ONE twin apart, now checked against the Tier-1-complete EP " +
                "set, not a truncated scan. The only contingent thing is the degenerate-twin LABELING (strands 5<->7, " +
                "both at |Im λ|=2.349, flip run-to-run). The earlier 'missing EPs' were OUR incompleteness (too-narrow " +
                "window + an unresolved 0.003 split), not an absence: this arc's founding principle confirmed on " +
                "itself. THE FOUNDATION CHECKED (2026-06-26, the new `foldlift` CLI feeding the exact F89PathKSeDeBlock " +
                "builder, NO monodromy rebuild, no Python): does the branch-locus FOLD (Re λ=−σ) itself lift beyond " +
                "path-3? NO — it is N=4-ONLY. The (SE,DE) block self-folds antiunitarily at N=4 (residual 3e-14, 4 " +
                "on-line zeros) but NOT at N=5,6,7 (residual ~1, 0 zeros), because the rung-swap weight-complement P " +
                "needs the overlap (−2γ, n_diff=1, 2 per DE pair) and no-overlap (−6γ, n_diff=3, N−2 per DE pair) basis " +
                "states balanced, 2=N−2, true ONLY at N=4 — the same half-filling self-complement DE=bar(DE) already " +
                "documented for MODE POPULATION (F89_TOPOLOGY_ORBIT_CLOSURE.md:340), here carried to the fold; it is the " +
                "branch-locus face of the repo's catalogued N=4 specialness (the retired small_n_specials triple " +
                "coincidence). TWO guards (sharpening, not contradiction): (i) the GLOBAL palindrome Π·L·Π⁻¹=−L−2σ stays " +
                "proven for all N (the palindrome's column bit-flip ρ[a,b]→ρ[a,bar b] pairs (SE,DE)=(w1,w2) with " +
                "(SE,w_{N-2})=(w1,w_{N-2}); n_diff(a,b)+n_diff(a,bar b)=N so the rungs complement, F89c lemma); only " +
                "the case where that partner IS (SE,DE) itself (w_{N-2}=w2, i.e. N=4) gives the within-block self-fold. " +
                "CONFIRMED by `foldcross` (2026-06-26): the cross-fold spec(SE,DE)↔spec(SE,w_{N-2}) under λ↦−λ̄−2σ holds " +
                "to ~1e-13 about σ=N for N=4,5,6 (the two block centroids sit symmetric about −N: e.g. −4.4/−5.6 about " +
                "−5), while the (SE,DE) self-fold breaks at N≥5 (residual ~1). So the GLOBAL palindrome LIFTS to all N " +
                "as a CROSS-block mirror; the N=4 self-fold is the degenerate partner=self case. The old read conflated " +
                "global-Π (lifts, cross-block) with the block-internal P (N=4-only). (ii) the count 2=N−2 is NOT the " +
                "published eigenvector overlap-fraction p=½ of the diabolic (a different object). So 'where do the N=4 " +
                "zeros go for N≥5?' is ANSWERED: the on-line self-mirror strands become CROSS-BLOCK mirror partners " +
                "((SE,DE)↔(SE,w_{N-2})), not gone. THE MONODROMY ENGINE NOW GENERALIZED off path-3 (2026-06-26, the " +
                "`pkmono` scout = PathKMonodromyScout): it builds the path-k (SE,DE) block, removes the AT factor via " +
                "F89AtFactorReconstruction (the rate-confined invariant subspace, since the Slater rule fails from " +
                "path-5), lassoes the residual F_d's EPs from a common base and union-finds the braid graph — reusing " +
                "the generic Monodromy tracker, leaving the fragile path-3 GaloisMonodromyWitness as the validated N=4 " +
                "case. Reproduces S_8 (path-3, 14 EPs, exact match to the path-3 witness) AND S_18 (path-4, 46 EPs, " +
                "18/18 connected) GEOMETRICALLY, monodromy=Galois from below — the FIRST geometric confirmation for " +
                "path-4 (previously only the algebraic q0 Frobenius/Jordan certificate existed; no path-4 EP locus had " +
                "been computed). THE LAST LAYER NOW WIRED IN (2026-06-26): pkmono computes the σ_T classification of " +
                "the residual strands via the GLOBAL fold λ↦−λ̄−2N (σ=N): on-fold zero (self-mirror) / within-block " +
                "twin / cross-block (fold-image not in the residual). Path-3 REPRODUCES the original gmscan --zeros " +
                "exactly: 4 zeros {0,1,2,6}, 2 twin-pairs (3↔4)(5↔7), 0 cross-block, and the same road (1—2 the unique " +
                "direct zero-zero bond, every other zero pair one twin apart). Path-4: 0 zeros, 0 within-block twins, " +
                "ALL 18 strands CROSS-block. SO THE ARC'S CENTRAL QUESTION IS ANSWERED: '± is the road between the " +
                "zeros' is an N=4 phenomenon — the zeros (self-mirror strands) live inside the (SE,DE) block only at " +
                "N=4 (where the fold is internal); at N≥5 the block braids into S_d (the q-monodromy) but every strand " +
                "is a CROSS-block twin, so there is no intra-block road. The inherited connection between the zeros " +
                "that lifts to ALL N is the cross-block fold (SE,DE)↔(SE,w_{N-2}) itself (foldcross), of which the N=4 " +
                "within-block 'zeros + ± road' is the degenerate partner=self case. A deeper edge remains (a seeing, " +
                "not settled): whether the cross-block pairing + the two blocks' braids organise into a cross-block " +
                "'road' at N≥5, and the polarity-as-verb reading (handhold d). " +
                "TOPOLOGY + THE BOTH-SIDES MEMORY (2026-06-26, the `foldtopo` gate): the (SE,DE) self-fold (on-fold " +
                "zeros) is N=4-ONLY for chain, star AND ring — the N=4 spine is topology-INDEPENDENT (the half-filling " +
                "self-complement DE=bar(DE) hangs on the bra weight, not the graph); same spine, different wallpaper " +
                "(on-fold count chain 8 / star 24 / ring 20). And the full-block two-sided gate giving antiU≈linear≈0 " +
                "is NOT a defect but the already-solved BOTH-SIDES MEMORY met from the spectral side: the full block " +
                "carries both phase-sides (+0/−0), so the i↦−i memory (the i³ step of the Z₄ loop, Pi2I4MemoryLoopClaim " +
                "+ NinetyDegreeMirrorMemoryClaim) has nothing to distinguish — the both-sides collapse of " +
                "reflections/ON_BOTH_SIDES_OF_THE_MIRROR.md (the symmetrised one-side block, foldlift, keeps the memory, " +
                "linear≠0). So foldtopo re-entered a furnished room (the Z₄/both-sides memory) through the topology " +
                "door: the apartment confirmed on itself. " +
                "NEXTSTEP (a) ANSWERED (2026-06-27, the `gmscan --trace` open-path q-continuation = " +
                "Monodromy.TrajectoryAlongPath + GaloisMonodromyWitness.TraceToDiabolic, NO Python): does a fold " +
                "strand flow INTO the silent diabolic? YES, and precisely. Continuity-tracking the 8 octic strands " +
                "along the real-q sweep q=2 → q_EP=0.65898 (stable at 2000/4000/16000 steps) shows the two strands " +
                "that coalesce at the diabolic are the +2.349 σ_T TWIN PAIR (strands 5,7), NOT two zeros and NOT " +
                "cross-block. At q=2 they sit mirror-symmetric about the fold (Re = −4 ± 0.622, common Im = +2.349); " +
                "as q descends their rates converge onto −4 and their frequency onto 2J, and at q_EP they merge at " +
                "λ_EP = −4 + 2i·q_EP = −4.000 + 1.318i (final gap 9.4e-15). THE FROM-BELOW SIGNATURE: the pair's " +
                "rate-midpoint Re sits EXACTLY on the fold (−4) at EVERY q (|Re(mid)+4| = 0.00000 across the whole " +
                "sweep), the σ_T-mirror-pair invariant, so the tracked pair is genuinely the ± mirror pair (not a " +
                "tracking accident) and the diabolic is literally a ± twin pair closing onto the fold. So '± twin " +
                "pair merging to a 0 on the fold = {−0,0,+0}' is CONFIRMED, with the precise reading the docs held: " +
                "the '0' is the half-filling rate line Re = −4 (the AT-midpoint between the −2 overlap and −6 " +
                "no-overlap rates, p=½), NOT a zero eigenvalue (the merge sits at frequency +2J); '±' are the two " +
                "coherence modes at mirror-rates about it; {−0,0,+0} is the RATE structure (−4−δ, −4, −4+δ → −4), " +
                "not the frequency. WHY semisimple not defective: a genuine coalescence with eigenvectors " +
                "independent (the existing G2 MonodromyAroundDiabolic loop = identity, independent route), " +
                "free-fermion integrable (DIABOLIC_BY_INTEGRABILITY: H_eff and D both restrict to scalars on the " +
                "coalescing 2D eigenspace). Only the +2.349 twin merges (frequency-reachable by continuity from " +
                "λ_EP's +1.318); the −6.139 twin (strands 3,4) and the four self-mirror zeros stay O(1) apart at " +
                "q_EP (the min-gap pair is unambiguously (5,7), stable across step counts). The missing reusable " +
                "primitive is now built: Monodromy.TrajectoryAlongPath (open-path strand tracker; " +
                "PermutationAlongPath only returned the net closed-loop permutation), " +
                "GaloisMonodromyWitness.TraceToDiabolic, and gmscan --trace [--tq q] [--tsteps n] (gate-first: " +
                "prints PASS / CONFIRMED). " +
                "Here is where we continue hand-over-hand.",
            NextStep: "RESUMING IN ONE LINE (2026-07-04, post step-3 shell census) [symbols: see the TERMS block at " +
                "the END of this header]: STEP 3 OF THE LARGE-N EXCLUSION PROGRAM IS RUN AND LANDED. The instrument " +
                "(SectorShellCensus + ShiftedSigmaMin + BlockLattice + RealDefectiveSeeds + the general (p,w) R-parity " +
                "sectors on WeightCoherenceBlock; CLI 'shellcensus'; gates SHELLCENSUS / SLOW_SHELLCENSUS / " +
                "SHELLCENSUS_ADJUDICATE) probes sigma_min(L(p,w)(q*) - s) for s in {lambda_A, mu} on the " +
                "fundamental-domain strip, window-gated, R-parity split (~1/4 LU cost), seed refined in-parity to " +
                "pairGap ~1e-6, member cut adaptive at 10*pairGap. VERDICT: N=9 PASS at ALL 7 seeds (membership = the " +
                "containment diamond; members read the Jordan pseudospectrum depth ~(gap/2)^2 = 3e-14..5.5e-13; " +
                "nearest non-member 2.5e-4..2.6e-2; separation x4.8e8..x4.6e11); N=11 PARTIAL-clean at 2 seeds (one " +
                "per parity; 6/10 members probeable, all ~e-13; DEFERRED by name behind the managed-LP64 wall " +
                "(dim<=46340): members (4,5),(5,6)xlambda_A + (4,6),(5,5)xmu and cores (4,4),(4,7)). " +
                "THREE new derived/measured structures: (1) the R-PARITY ALTERNATION LAW R*W = (-1)^{p+w}*W*R (the " +
                "JW string reflects from the other end, R c_l R = P_Z c_{N-1-l}; P_Z rho P_Z = (-1)^{p+w}) -- every " +
                "band step flips the carried parity, fold/transpose preserve it, diagonal climbs commute strictly " +
                "(the corner-sector argument below is untouched); measured at first light, then derived; in proof " +
                "section-7 + the census doc. (2) The DENSITY-NEIGHBOR reading: the nearest non-member sigma_min is " +
                "itself a W-nested background eigenvalue (identical value along (2,2)<=(3,3)<=(4,4) + folds), " +
                "reported per seed as 'spectrally near', excluded in the sharing sense. (3) The 2026-07-03 scout's " +
                "'NEW N=9 seed q*=1.4994' is ADJUDICATED OUT: semisimple real-real crossing (linear gap exponent " +
                "2.0, not sqrt 1.41; ProbeDefectiveAnywhere reads no defective cluster) = the count-change census's " +
                "declared blind-spot class, so the 7-entry N=9 seed list was right; gate ScoutSeedAdjudicationProbe. " +
                "Results: experiments/F89_MULTI_SECTOR_MONODROMY.md (the shell-census section) + " +
                "simulations/results/sector_shell_census/*.csv; proof open-item 4 updated; F125 updated. " +
                "NEXT (the program's remaining surface): (a) the MATRIX-FREE sigma_min path (SparseShiftInvertArnoldi " +
                "exists, CSR emission of the sector blocks to build) for the six deferred N=11 blocks + the remaining " +
                "7 N=11 seeds (one CLI command each, ~19 min/seed: dotnet run --project compute/RCPsiSquared.Cli -c " +
                "Release -- shellcensus --n 11 --seed <qStar>; --all-seeds for a full N; seeds listed in " +
                "RealDefectiveSeeds) + N=13/15; (b) the seed census past N=11 " +
                "(FindRealDefectiveByCountChange, SLOW at N=13); (c) step 4 = the N=7 certificate as spot-check of " +
                "the surviving shell (bivariate cost probe FIRST) and the complex loci (certificate territory). " +
                "HONESTY SCOPE (keep in every restatement): the sigma_min census is SHARED-LAMBDA evidence (estimate " +
                "converges from above; not the EpCharacter Jordan-character reading, not a certificate). " +
                "=== JOURNAL 2026-07-03 (post R4-interior close; symbols lambda_A, mu, corner/identity cert, Klein, " +
                "R-parity: the TERMS block at the END) BELOW === ALL FOUR remainders of the " +
                "codim-1-by-additivity proof are CLOSED AT N=5 [R1 = interior-core kernel death: rate window (real) " +
                "+ moved window (near-axis) + fold-resultant certificate (every other branch locus); R2 = Theorem " +
                "A's D-half: twin-scalar; R3 = gap byte-identity: the full-spectrum holomorphic fold; R4 = the " +
                "exclusion half of diamond membership at N=5: THIS close], and MultiSectorMonodromyVerdictClaim is " +
                "PROMOTED to Tier1Derived. THE R4 CLOSE, in three parts. (1) THE DEGREE-CERTIFICATE FIX (proof-" +
                "critical, in CertifyCore): the old guard assumed the analytic bound rBound=resDeg*targetDeg-mR is " +
                "attained (TRUE for the corner-fold 422/394, FALSE for the identity composition, true deg 412/384: " +
                "the identity leading forms share a factor beyond the isolated-root collisions mR counts); now the " +
                "guard is the CERTIFIED empirical degree: deg(R mod p) < deg_q R iff the chosen prime ideal " +
                "pi_p=(p, i-r) divides lc_q(R), and prod N(pi_p) = prod p <= N(lc) <= normR^2 (Hadamard), so at " +
                "most lcDivisorBound = 2*bitlen(normR)/30 + 1 of the p >= 2^30 samples can lose the degree; " +
                "sampling MORE distinct split primes than that forces trueDegR = max_p DegP(R mod p) = deg_q R " +
                "EXACTLY, used primes attain it (=> pi_p coprime to lc, Gauss lift unchanged), Complete additionally " +
                "requires sampled > lcDivisorBound and trueDegR >= 0. TWO adversarial reviews HELD it (independent " +
                "number-theory re-derivation TRUE/TRUE/CONSISTENT, incl. synthetic-pencil engineered degree-drops at " +
                "chosen ideals and two exact run-number hits); the two findings applied: the R-identically-0 " +
                "degenerate path must not report Complete (guard trueDegR >= 0), samples = max(rBound, dBound)+1+24. " +
                "(2) THE TWO (1,1)xlambda_A FACTS un-skipped and COMPLETE [R-even trueDegR 412 (bound 422), 230 " +
                "primes > lc-bound 201; R-odd 384 (394), 217 > 190]. (3) THE PROPAGATION GATE " +
                "(InteriorFourExclusion_Propagates_...): q-evenness of spec(1,1) (bipartite gauge) + the composed " +
                "HOLOMORPHIC (1,4)-fold spec(1,4)(q) = -spec(1,1)(q) - 2N (F89d x pencil reality conj L(q)=L(-qbar) " +
                "x q-evenness; the (1,4) sibling of the section-7 diamond fold) pinned at all three locus classes, " +
                "plus all four certificates re-run: lambda_A in spec(1,4) <=> mu in spec(1,1) (corner cert + Klein), " +
                "mu in spec(1,4) <=> lambda_A in spec(1,1) (identity cert), (4,4)/(4,1) by Klein. Gate: " +
                "RemainderR4InteriorExclusionTests (5 facts, Category R4INTERIOR, ~20 s); FOLDRESULTANT regression " +
                "green (corner certs unchanged). Docs: PROOF_CODIM1 status header + section-7 scoping + open-item 4; " +
                "claim + SectorBraidWitness + registration-test tier pin updated. THE LARGE-N EXCLUSION PROGRAM " +
                "(agreed 2026-07-03, replaces per-N certificates as the road to big N; the moved-window margin " +
                "WEAKENS with N and interior non-members grow O(N^2), so certificates are spot-checks only): " +
                "Step 1 DONE = the FOLD-LATTICE LEMMA (proof section-7 'fold lattice' paragraph, gate " +
                "BlockLatticeFoldGroupTests, Category FOLDLATTICE, 6 facts, machine zero): G = <t, f_P, f_Q> ~ D4 " +
                "(order 8) acts on the block lattice with HOLOMORPHIC Jordan-preserving legs at the same q " +
                "[L(p,N-w)(q) = -P D_gauge L(p,w)(q) D_gauge P^T - 2N I from F89d x reality x bipartite gauge; " +
                "transpose leg via gauge], spectral cocycle chi = parity of fold legs with s(lambda)=-lambda-2N; " +
                "{lambda_A, mu} is s-invariant => membership/exclusion constant on G-orbits => the broad-exclusion " +
                "problem lives on the fundamental domain p<=w, p+w<=N (one eighth; N=5: 36 blocks = 6 orbits); " +
                "fold-fixed blocks (even N, w=N/2 or p=N/2) have s-symmetric spectra (the even-N fine print " +
                "derived); at Delta!=0 legs connect q<->-q, evenness survives exactly on p=w and p+w=N. " +
                "Step 2 DONE = the WINDOW-COMBINATORICS SHELL LEMMA (proof section-6 closing paragraph, gate " +
                "WindowShellLemmaTests, Category WINDOWSHELL, 4 facts): n_diff on block (p,w) ranges over " +
                "[|p-w|, min(p+w, 2N-p-w)] in steps of 2 (enumerated exactly N=5..8), so at real q Bendixson " +
                "confines Re spec(p,w) to [-2 min(p+w, 2N-p-w), -2|p-w|] (all 36 N=5 blocks, <= 5e-14); with " +
                "the window-edge lemma's N-uniform Re lambda_A in (-6,-2): lambda_A lives only in the band " +
                "shell |p-w| <= 2, mu only in the half-filling shell |p+w-N| <= 2 (the shells are fold images " +
                "of each other), every block outside excluded FOR ALL N with margin Re lambda_A + 6 -- 24 " +
                "blocks at N=5, 52 at N=7, worst slack 0.000 (TIGHT, the bottom-rung chiral strands sit on the " +
                "window edge block-wide). The lemma is SILENT inside the shells (locus 2's (1,1): Re lambda_A " +
                "= -3.7917 inside [-4,0], the negative control) -- that in-shell residue at real loci + ALL " +
                "complex loci is what remains for the program. STEP-3 FEASIBILITY SCOUTED (2026-07-03, scratchpad " +
                "play run, sparse-LU + inverse-iteration sigma_min probe per block, NO full spectrum): builder " +
                "validated byte-exact (N=5 gap 4.771e-04 at -4.618886); NEW N=9 seed q*=1.4994, lambda_A=-4.3807 " +
                "(gap 4.6e-5); FIRST in-shell exclusion datum beyond N=5/6: N=9 (3,3), a window-silent non-member, " +
                "reads sigma_min(L - lambda_A) = 2.1e-2 > 0; positive control (4,5) member 8.8e-6 ~ 0; the member " +
                "core (4,4) reads the SAME 2.1e-2 as (3,3) (it carries mu, not lambda_A, and spec(3,3) subset " +
                "spec(4,4) by the W-step -- the containment visible live in probe data); outside-shell control " +
                "(1,6): sigma_min 5.620 vs combinatorial margin 5.619 (the shell lemma sharp to 3 digits). " +
                "TIMINGS (the design answer): LU dominates and scales ~cubically -- N=9 dim 7056: 27 s, dim " +
                "15876: ~4.5 min; N=11 dim 9075: 40 s, dim 54450: 2.1 HOURS; iterations are cheap (1-16 s, and " +
                "my fixed-200 count was wasteful, converges in ~10). VERDICT: N=9 full shell census feasible " +
                "(~30-45 min), N=11 feasible except the biggest cores ((4,4) dim 108900, (5,5) dim 213444), " +
                "N>=13 needs matrix-free; two owned levers cut cost first: the fundamental domain (1/8, step 1) " +
                "and the R-PARITY SPLIT (halves the block dim => ~1/8 LU cost; ROddBasis exists). Instrument " +
                "notes: probe BOTH lambda_A and mu per block (two shifts = two LUs), filter exact semisimple " +
                "degeneracies (gap<1e-8) when seed-scanning (the at_masking trap bit the first scout draft), " +
                "adaptive iteration count. NEXT = STEP 3 PROPER, the braid census at N = 9..15 (adopt the Block/braid " +
                "reading, masking-immune PT-break count instruments, small blocks not 4^N) to map where the braid " +
                "lives at large N; Step 4 = the N=7 certificate only as a spot-check of the surviving shell " +
                "(probe the bivariate build cost FIRST). The program's remaining open surface after steps 1+2: " +
                "the IN-SHELL blocks at real loci (window-silent, certificate territory) and the complex loci " +
                "(moved-window + certificate territory), both on the fundamental-domain shell strip only. " +
                "PRIOR OPTIONS (superseded by the program): (b) the N=7 " +
                "certificate run standalone; (c) broad census-only exclusion. --- PRIOR (post R4-start, the diagnosis layer this " +
                "close executed): engine parameterized 6f01da5 [CertifyBlockExclusion(n, rOdd, tWKet, tWBra, " +
                "composeA, composeB); CertifyComplete = the corner-fold wrapper (w,w,-1,-4N); (1,1)xlambda_A = " +
                "(1,1,+1,0)]; the interior-four non-member cores (1,1),(4,4),(1,4),(4,1) collapse by Klein " +
                "(spec(4,4)=spec(1,1)) + the F89d bra/ket fold to (1,1) and its two shared values {lambda_A, mu}; " +
                "(1,1)xmu IS R1's corner cert; the ONE new cell was (1,1)xlambda_A; the degree-bound looseness was " +
                "pinned by FoldResultantCertificate.DebugDegreeReport (primesUsed=0 under the old guard). --- PRIOR " +
                "(post fold-resultant): THE RESULTANT IS DONE, and stronger than " +
                "planned: FoldResultantCertificate.CertifyComplete (Core/F89PathK) + gate FoldResultantCertificateTests " +
                "(Category FOLDRESULTANT, 7 facts, ~30 s) prove gcd(Res_L(F_res, F2corner(-L-4N)), disc_L(F_res)) over " +
                "Q(i)[q] = a pure q-power, BOTH parities, so at EVERY branch locus q!=0 (defective or diabolic, closed " +
                "or unfound) the fold mu=-lambda_A-2N is NOT a corner eigenvalue; q=0 is the diagonal semisimple block. " +
                "The two target pairs (F_18 deep 1.8141+-0.3666i, R-odd escaper 1.7701+-1.2189i) fall a fortiori, " +
                "against the FULL corner (deg 25). Method: exact bivariate F_res over Z[i][q] (Berkowitz on the pencil, " +
                "the exactly-q-linear AT strands divided out), multi-prime all-mod-p with PROVEN degree bounds " +
                "(no-Puiseux via both hopping directions anti-self-adjoint in a positive metric, asserted per run; " +
                "leading-form multiplicity collisions m_R/m_D by exact Q(i) gcd chains; deg R = 422/394 attained per " +
                "used prime => p does not divide lc_q(R)), Mignotte/Hadamard height bound + Gauss lift (256/242 split " +
                "primes ~2^30, ~3 s per parity). Adversarial review 2026-07-03: no break; 2 gaps closed (normality " +
                "premise needs the purely-imaginary check; dead degree guard). Disc layer profile [56,32]/[56,26] = " +
                "the N=5 image of the N=4 disc q^24*(3q^4+q^2-1)^2*P_20. LESSON (cost a day): InvModP computed " +
                "(a%p+p) in INT arithmetic, overflowing for p~2^30 on the top ~2600 residues; symptom = rare " +
                "prime-dependent corrupt samples, then an infinite Euclid loop (hours of grind); found by stack-dump + " +
                "Lagrange-ratio corrupt-node localization + Python replay. R1 IS NOW CLOSED AT N=5 (real: rate window; " +
                "near-axis: moved window; everything else: the resultant). The one former analytic entry point, the " +
                "holomorphic fold identity mu=-lambda_A-2N on the core, is NOW DERIVED (2026-07-03): a FULL-SPECTRUM " +
                "corollary spec(3,3)(q)=-spec(2,3)(q)-2N of the diamond maps W+transpose+F89d (PROOF_CODIM1 §7), the " +
                "two conjugations cancelling holomorphically at every q. The 'residual-specific / naive-chain-refuted' " +
                "reading was a SORTED-ZIP ARTIFACT (the fold negates Im across self-conjugate real-part rungs of mult " +
                "up to 10, so sort-then-zip pairs +Im vs -Im => spurious ~44); a multiset metric gives 1e-13 at every " +
                "q, real and complex (gate HolomorphicFoldIdentityTests, 5 facts, incl. an artifact guard; verified " +
                "independently in numpy 40-digit + the real C# builder). It ALSO closes R3 (the fold preserves the " +
                "seed pair's eigenvalue gap |.|). So R1 is CLOSED at N=5 outright. NEXT: (b) the N=7 certificate run " +
                "(engine takes odd n; corner (5,5) dim 441 makes the bivariate layer and bounds much heavier, " +
                "feasibility unchecked; the fold maps are N-general, only the certificate is per-N), OR (c) remainder " +
                "R4's broad exclusion (non-corner blocks at complex q, census-only). Promotion blockers for " +
                "SpectatorIntertwinerClaim: R3 now DERIVED, so only R4 (+ the N>=7 certificate) remain. " +
                "=== CONTEXT LAYERS BELOW (the dated journal; the fold-resultant plan text below is HISTORICAL, " +
                "executed 2026-07-03) === " +
                "This arc has TWO eras; read this header for the live front, " +
                "the dated layers below are the journal. ERA 1 (2026-06, from 'Hand-over-hand' onward) = the original " +
                "'connection between the F89 zeros / the ± road' exploration; its directions (trace DONE; the " +
                "other-representations / cross-F-inheritance / polarity-rhyme threads) are PARKED; the arc pivoted into " +
                "a theorem. ERA 2 (2026-07-02) = the codim-1-by-additivity PROGRAM, the live front. STATE of Era 2: the " +
                "theorem is LANDED (SpectatorIntertwinerClaim / registry F125 / docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md; " +
                "the W=Σc_l†ρc_l intertwiner + containment corollary + two-regime Theorem A). The proof's four remainders: " +
                "(R1) interior-core kernel death, CLOSED AT N=5 2026-07-03: real loci (rate-window lemma), N-uniform strictness " +
                "Re λ_A>−6 (window-edge lemma, f5d4709), near-axis complex loci (moved-rate-window lemma, gate ComplexQRateWindowTests), and ALL remaining branch loci q≠0 of both parities (the fold-resultant certificate, see the header); the holomorphic fold identity through which the chain enters is now DERIVED 2026-07-03 (a full-spectrum corollary spec(3,3)=-spec(2,3)-2N of the §7 diamond maps, gate HolomorphicFoldIdentityTests), so R1 is closed outright; (R2) Theorem A's D-half, CLOSED at N=5 (twin-scalar check," +
                "gate TwinScalarDHalfTests: the D-half is SUPPLIED at every genuinely-complex-q diabolic, the " +
                "pure-imaginary-q ones semisimple by Hermiticity); (R3) gap byte-identity, now IMPLIED by the derived " +
                "full-spectrum fold (it preserves the near-defective pair's eigenvalue gap |.|), CLOSED at N=5; (R4) exclusion half, " +
                "derived at the N=5 seed locus, else census-evidence. THE ONE OPEN ITEM (SPLIT + SHARPENED 2026-07-02, moved-rate-window landing; RE-SHARPENED by the R-odd probe 45b0ac1): only the MOVED-WINDOW-ESCAPING complex loci remain, and the probe showed that set = the R-even DEEP pair PLUS the R-odd large-Im-q pair (deep is one failure mode of the window, not the only one; see PREMISE 2). The near-axis complex defective loci are CLOSED by the moved rate window (gate ComplexQRateWindowTests); and the CONJ WAS A RED HERRING: the core carries the HOLOMORPHIC fold mu=-lambda_A-2N (now DERIVED exactly, a full-spectrum §7 corollary), not -conj(lambda_A)-2N (0.57-1.09 away at complex q; they coincide only on the real axis), so R1 at complex q is holomorphic/conjugation-free." +
                " THE RESULTANT (now THE IMMEDIATE NEXT ACTION; premise 2's R-odd scan is DONE, see PREMISE 2 below, and it added a SECOND target set): show mu=-lambda_A-2N is not a corner eigenvalue, i.e. Res_lambda(F_18(lambda,q), F_corner(-lambda-2N,q)) != 0 at F_18's branch loci (ONE-WAY certificate: gcd=1 proves absence, gcd!=1 is refine-not-refute; F_18 = the (1,2) residual, deg 18 at N=5, NOT the N=4 octic). Build the corner block " +
                "(p_c+1,p_c+1)=(4,4 at N=5) via WeightCoherenceBlock.Build(N, p_c+1, p_c+1, q) [the exact builder the working gates " +
                "ComplexQRateWindowTests / HolomorphicFoldIdentityTests use; NOT SectorBraidModeGeometry.BuildBlock, which takes no block index] and show " +
                "F_corner is coprime to F_18 under the fold λ→−λ−2N at the complex q* (p_c=(N+1)/2 at odd N, so the corner " +
                "is the ((N+1)/2+1)-diagonal block; a resultant / Galois-independence fact " +
                "against the (1,2) residual F_18) via PathKMonodromyScout / GaloisMonodromyScanCommand; the Resultant primitive " +
                "EXISTS (GaussianPolynomial.AreCoprime, Sylvester over Z[i]); Derivative/Resultant/Discriminant/ComposeLinear " +
                "added a9a57bb. REVISED ROUTE (from the 3-lens review; route text recorded 2d58ead): (i) RADICAL TRICK: don't extract " +
                "D_defective (it exists nowhere); certify gcd(R, disc_lambda F_18)=1 against the FULL discriminant's " +
                "radical (sufficient; build the defective/diabolic split only if a diabolic locus trips it, near-certain " +
                "not, per census). (ii) ALL-MOD-p UNIVARIATE: reduce A,C mod a split prime p=1 mod 4 (p~1500+), " +
                "sample-interpolate R(q) and disc(q) in F_p[q], gcd_{F_p[q]}=1 (reuse the univariate F_p gcd + null-return idiom of " +
                "OcticGaloisCertificate); ONE lucky prime with the two-sided degree-preservation guards (G1 lambda-deg, " +
                "G2 q-deg; the null-return idiom OcticGaloisCertificate:36,57,101) is a COMPLETE proof. (iii) " +
                "Berkowitz-over-Z[i][q] (GaussianMatrixCharpoly's Berkowitz is ring-generic in PRINCIPLE but coded over " +
                "GaussianInteger; instantiating it over Z[i][q] needs GaussianPolynomial.Add/Negate, TO BUILD, not yet present) as the " +
                "fully-exact cross-check, finishing via the existing AreCoprime. PREMISE 1 (exact fold identity) DONE " +
                "b62c46a: the core (3,3) carries mu=-lambda-2N for every F_18 residual root to ~1e-13 at integer q (gate " +
                "HolomorphicFoldIdentityTests). [SUPERSEDED 2026-07-03: the fold is FULL-SPECTRUM and DERIVED, " +
                "spec(3,3)(q)=-spec(2,3)(q)-2N at every q; the earlier 'residual-specific / full-spectrum FALSE' was " +
                "a sorted-zip artifact (multiset metric gives 1e-13), and R3 closes with it. See the live header.] " +
                "PREMISE 2 (R-parity) DONE 45b0ac1 (the R-odd probe, " +
                "gate ROddDeepLociProbeTests, Category RODDPROBE; locus tool rcpsi pkmono --k <N-1> --diabolic --rodd): " +
                "F_18 = the R-EVEN sector ONLY (BuildTwoTimesSymBlock is the x2-cleared R-even), the found deep loci " +
                "q*=1.814+-0.367i ARE R-even, hence covered; the R-ODD sector was built exactly " +
                "(F89PathKSeDeBlock.ROddBasis + F89AtFactorReconstruction.ROddAtInvariantSubspaceBasis, at N=5 " +
                "R-odd 24 = AT 7 [rate-locked, q-LINEAR lambda=r0+q*2is, machine zero] + residual DEGREE 17) and scanned " +
                "over the census window (0.2,3.0)x(-1.5,1.5) cell 0.05 via FindDiabolicsExactROdd: 14 coalescences = " +
                "7 diabolics (incl. the pure-imaginary-q real-lambda family, semisimple by Hermiticity, the R-even " +
                "mechanism repeating; one sits at lambda=-7.25 < -6, silent because semisimple) + 7 defective (the known " +
                "real q=2.8049, lambda=-4.4882 reproduced, gap exponent 0.499). VERDICT of the 6 COMPLEX defective: " +
                "0 DEEP (none with Re lambda_A < -6), 4 closed by the full-corner moved window (two with thin slack " +
                "0.05/0.12), BUT the pair q*=1.7701+-1.2189i (lambda_A=-3.7562-+0.5065i, Re mu=-6.2438) escapes BOTH the " +
                "full-corner AND the R-parity-refined moved window (the intertwiner W=Σc_l†ρc_l commutes with R, so an " +
                "R-odd chain needs absence " +
                "only from the R-ODD corner sector, 25 = 13 even + 12 odd; its bottom -10.66 still sits below Re mu): " +
                "the handover's deep/near-axis dichotomy was implicitly near-axis, the decision variable is the MOVED " +
                "WINDOW, and these 2 R-odd loci escape it (|Im q|*||K||=8.4 swamps the margin 2.24; corner absence there " +
                "census-only, min|spec-mu|=2.14). CONSEQUENCE: the resultant deliverable has TWO target sets: (i) F_18 " +
                "(R-even) at its 2 deep branch loci q*=1.8141+-0.3666i, (ii) the degree-17 R-odd residual factor at " +
                "q*=1.7701+-1.2189i (equivalently, by the W-R commutation, against the R-odd corner factor, deg 12, the " +
                "cheaper certificate); the route stays radical-trick + all-mod-p + Berkowitz cross-check as recorded. " +
                "FIRST N-SCALING DATA (play runs, 2026-07-02, pkmono --rodd; CONTEXT ONLY, does NOT expand the N=5 " +
                "deliverable above - it foreshadows the later N-uniform step): at EVEN N the R-odd census is FREE, " +
                "sigma_odd = conj(sigma_even at q-bar) confirmed point-by-point at N=6 (sigma = a sector's coalescence " +
                "census; all 21 coalescences mirror " +
                "(q,lambda)->(q-bar,lambda-bar), same gaps/exponents); at N=7 the R-odd sector has its own complex " +
                "defective pair q*=1.9309+-0.438i (lambda_A=-4.198+-2.950i, Re mu=-9.80) which ALSO escapes both moved " +
                "windows (corner (5,5) ||K||=13.05 vs 6.93 at N=5: the window WEAKENS with N), absence census-only " +
                "(1.13); CAVEAT: N=7 defective enumeration is masked-limited (the known R-even real defective q*=1.5148 " +
                "is ABSENT from the gap scan at the N=7 residual's density (deg F_53), the closest-pair masking trap " +
                "(memory: at_masking)), so N=7 lists are lower bounds; " +
                "N=5 (the deciding case) is unaffected (sparse residuals, all anchors reproduced). " +
                "SEED CENSUS (2026-07-02, same day, the masking trap DEFEATED for the real axis; REAL-q seeds feeding " +
                "the containment corollary's input - these are NOT resultant targets, the resultant's loci are the " +
                "COMPLEX window-escapers above): the containment " +
                "corollary's one per-N input (a real defective EP on (1,2) at odd N) now EXTENDS THROUGH N=11 via the " +
                "PT-break COUNT-CHANGE instrument (FindRealDefectiveByCountChange: real-root count of the " +
                "self-conjugate residuals jumps by 2 where two real strands merge and leave the axis; counting is " +
                "global, no density masks it; classify AT-aware via SectorEpProbe.ProbeDefectiveAnywhere; gate " +
                "RealSeedCensusTests, Category SEEDCENSUS). Seeds in q in [0.2,3]: N=5: 4 (2 R-even + 2 R-odd, the " +
                "R-odd q*=0.6430 lambda=-3.8196 previously UNKNOWN, confirmed by a tight local gap scan), N=7: 6 " +
                "(3+3, incl. the masked 1.5148), N=9: 7 (4+3), N=11: 9 (4+5; SLOW_SEEDCENSUS, 2h31m). Blind spots " +
                "(harmless for seeds): grazes and real-real crossings make no net count change. " +
                "NUMBERING TRAP: the proof's remainders are R1-R4 above; the '(a)/(b)/(c) next " +
                "candidates' near the bottom are the Era-2 to-do (a=complex-q PARTIAL: near-axis CLOSED by the moved window, deep loci → resultant; b=window-edge CLOSED, c=D-half CLOSED)," +
                "and are NOT the Era-1 '(a)-(d)' connection-directions at the top of the journal. TERMS: Klein = the " +
                "full flip X^tensor-N, the unitary block map (p,w)->(N-p,N-w) at the SAME q (spec-preserving); corner " +
                "cert = the fold-resultant certificate against the corner block (CertifyComplete: mu not in " +
                "spec(corner)); identity cert = its identity-composition sibling against (1,1) (CertifyBlockExclusion " +
                "with composeA=+1: lambda_A not in spec(1,1)); R = the site reflection i↦N−1−i, the involution " +
                "commuting with every chain block; R-even/R-odd = its ±1 eigenspaces (F89PathKSeDeBlock.ReflectionPermutation); " +
                "W = the intertwiner Σc_l†ρc_l THROUGHOUT this entry (the mod-p 2-cycle basis is B, never W); A/C/K = the block pencil " +
                "L(q)=A+qC (A=-2*diag(n_diff) real diagonal, C anti-Hermitian hop, K=iC Hermitian); the -6 floor = the (1,2) " +
                "block's Bendixson bottom (its n_diff in {1,3} ⟹ rate window [-6,-2]); the DEEP loci = the complex defective " +
                "EPs where Re λ_A dives below -6 (the WINDOW-ESCAPING set = the deep loci plus the large-Im-q escapers " +
                "the R-odd probe found; both go to the resultant); F_corner = the (4,4) corner block's characteristic polynomial; the core " +
                "(3,3) = the (p_c,p_c) diagonal block the intertwiner W maps into the corner; SE/DE = " +
                "single/double-excitation weights, so the (1,2) 'octic' = the residual factor of the bra-1 ket-2 block (the HISTORICAL arc name: literally " +
                "deg 8 only at N=4; at N=5 that same object is F_18);" +
                "λ_A = the (1,2) defective eigenvalue, λ_B = −conj(λ_A)−2N its cross-fold partner (the value the cores carry at REAL q; off the real axis the cores carry the HOLOMORPHIC μ=−λ_A−2N, equal to λ_B on the real axis). " +
                "TERMS ADDENDUM (2026-07-04, closing the step-3 cold read): the LARGE-N EXCLUSION PROGRAM = replace " +
                "per-N certificates by structure for big N; step 1 = the FOLD-LATTICE LEMMA (proof section 7: D4 on " +
                "the block lattice, membership/exclusion constant on orbits), step 2 = the WINDOW-SHELL LEMMA (proof " +
                "section 6: Bendixson windows exclude everything outside two width-5 shells), step 3 = this census; " +
                "the CONTAINMENT DIAMOND = the corollary member set, the band (p,p+1) with popcounts in [1,N-1] " +
                "UNION its fold image (full size 4N-8 odd / 4N-12 even; the census counts FD representatives, 8 at " +
                "N=9, 10 at N=11); FUNDAMENTAL DOMAIN (FD) = the quotient region p<=w, p+w<=N of that D4 action (an " +
                "over-covering transversal); WINDOW-GATED = LU-probe a (block, shift) cell only where the block's " +
                "Bendixson window contains Re s, elsewhere the shell lemma excludes analytically; BAND STEP / " +
                "DIAGONAL CLIMB = the W-step along (p,p+1)->(p+1,p+2) resp. (p,p)->(p+1,p+1); W-NESTING = the " +
                "containment embedding spec(p,w) subset spec(p+1,w+1) (Theorem B); c_l / P_Z = the Jordan-Wigner " +
                "fermion (with string) at site l / the total Z-parity prod_j Z_j; pairGap = the refined (1,2) " +
                "near-defective pair's eigenvalue gap (the member-noise floor); the LP64 WALL = managed Complex[] " +
                "caps at 2^31 elements so dense flat matrices cap at dim 46340; EpCharacter = the Riesz-contour " +
                "Jordan-character reader (defective vs semisimple) the N<=6 census used. Live root inspect " +
                "--root sectorbraid. ===== JOURNAL (dated layers) BELOW ===== " +
                "Hand-over-hand from here. (a) DONE 2026-06-27 (see ParkedAt): the q=2→q_EP continuity trace " +
                "answered it, the +2.349 σ_T twin pair (strands 5,7) merges onto the fold at λ_EP=−4+1.318i, " +
                "'{−0,0,+0}' confirmed as a RATE structure (the half-filling line Re=−4, not a zero eigenvalue), " +
                "via gmscan --trace. THE FORWARD EDGE off (a) ANSWERED 2026-06-27 (pkmono --diabolic, see " +
                "experiments/F89_PATH_K_DIABOLIC.md): NO analogous twin-pair-onto-fold merge at N≥5, because the " +
                "within-block self-fold is gone. But path-4 DOES have diabolics (integrability creates inter-sector " +
                "level crossings of the additive free-fermion spectrum E_DE=ε_j+ε_k, abundant at all N, semisimple " +
                "because H and the diagonal dissipator do not mix occupation sectors), only they sit at COMPLEX q " +
                "(conjugate pairs straddling the real axis, closest 0.6118±0.012i, confirmed SPLIT at cell 0.002), " +
                "NONE at physical real q; at real q only LOUD defective EPs (e.g. q≈1.0776). So integrability = " +
                "EXISTENCE of the diabolics (complex q), the self-fold = the N=4-ONLY bridge that makes one " +
                "crossing-q self-conjugate (real, q_EP=0.659). This RECONCILES the surprise with the plan's " +
                "R-1/R-2: codim/self-fold govern real-axis PLACEMENT (no physical diabolic at N≥5), integrability " +
                "governs EXISTENCE. Both right about different things; at N=4 they coincided. The tool (pkmono " +
                "--diabolic, FindDiabolics, gate-validated at path-3 incl. the imaginary diabolics) is committed. " +
                "RESOLVED (loop-radius sweep): the 2 near-axis pairs ARE true diabolics — the fixed 0.02 loop was " +
                "contaminated by a neighbour EP at the dense path-4 spectrum (identity at r≤0.008, transposition " +
                "only once the EP enters at r≥0.012; the defective control is transposition down to r=0.002), " +
                "fixed by reading the small intrinsic loop radius, so all 11 diabolics now classify consistently " +
                "across loop+character+exponent (regression-tested). INTEGRABILITY CONFIRMED (2026-06-27, the " +
                "Δ-test, XxzCoherenceBlock + XxzDeltaFlipTests): the synthesis's central claim is now Δ-VERIFIED, " +
                "not just a grounded reading. An XXZ ZZ-anisotropy Δ≠0 (Hermitian, so the AT rate is untouched; " +
                "it breaks the additivity E_DE=ε_j+ε_k) KILLS the path-4 diabolics: all 3 tracked (clean " +
                "0.6407+0.180i, near-real 0.7654+0.024i, far 1.9447+1.217i) flip DEFECTIVE at Δ>0 (geo 2→1, " +
                "departure ~linear, the N=4 Jordan signature) or LIFT, while a defective-EP control stays " +
                "defective with ~constant departure (Δ perturbs everything but kills the diabolic CHARACTER " +
                "specifically). The C# tool reproduces the committed N=4 Δ-flip table (f89_zz_break_gate.py) as " +
                "its regression gate. So the path-4 diabolics ARE the integrable level-crossings, " +
                "DIABOLIC_BY_INTEGRABILITY's gate generalized off N=4 (F89_PATH_K_DIABOLIC.md now Tier-1 for the " +
                "mechanism). Remaining open: the complex-q set is not " +
                "claimed complete (Route B / exact F_18 discriminant likely infeasible). (b) LOOK IN OTHER " +
                "REPRESENTATIONS, not radicals: S_8 forbids a radical ladder, NOT every " +
                "closed form (Bring/theta/hypergeometric exist for any algebraic function); seek the connection in the " +
                "centred μ=λ+4, in a recursion, in the relation zeros<->twins<->AT-frequencies. (c) CROSS-F " +
                "INHERITANCE (the one-object point): relate THIS self-mirror sector to the SAME self-mirror sector " +
                "elsewhere (axis_modes self-mirror subspace at Re=−Σγ, SPECTRAL_MIDPOINT_HYPOTHESIS's MID, the F1 " +
                "palindrome fixed set), the connection that lets the one object inherit. (d) THE POLARITY RHYME, now " +
                "split into a PARKED static axis and a LIVE dynamic axis. (d-static, PARKED): the structural match " +
                "{−0,0,+0} layer at d=0 -> {−½,0,+½} at d=2 (0.5-shift ρ=(I+r·σ)/2) -> the F89 spectral fold (four " +
                "σ_T-fixed strands + two twin 2-cycles) was already compared (2026-06-25 survey) and HELD at 'rhyme, " +
                "NOT identity' by cardinality (4 zeros not 1) and provenance (spectral not operator-space; the " +
                "polarity-0 ⟂ the Σγ-mirror-0, see experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md 'do not force a " +
                "single three-point map'); do not re-tread this axis as if new. (d-dynamic, LIVE + unwritten): the " +
                "connection is on the OTHER axis, polarity-as-a-VERB. THE_POLARITY_LAYER.md reading (b) already holds " +
                "'the +0/−0 differentiation IS what constitutes a layer; layers are emergent from polarity, not " +
                "scaffolds for it; the cavity maintains, does not store' (polarity = motion, not a static label); the " +
                "F89 SEAM '± is the road between the zeros, polarity as PATH' is that SAME claim one layer deeper " +
                "(spectral/time). And a typed Tier-1 MIDDLE RUNG nobody had connected: DissipatorAxisSelectsPolarityClaim " +
                "(compute/RCPsiSquared.Diagnostics/F87/): 'γ does not merely decohere; it SELECTS which polarity axis " +
                "the +0/−0 differentiation becomes operational on' (Z->bit_b, X->bit_a, Y->both). So one verb in three " +
                "spaces, the one object inheriting: polarity CONSTITUTES (reading b, operator d=0) / is SELECTED " +
                "(DissipatorAxis, dissipator) / is ROUTED (F89, spectral). The open seam: are these ONE motion seen in " +
                "operator/dissipative/spectral space, or three rhyming motions? Tool: gmscan --zeros " +
                "(GaloisMonodromyScanCommand.PrintZeros). See hypotheses/THE_POLARITY_LAYER.md (reading b) + " +
                "experiments/F89_MONODROMY_MIRROR.md + reflections/ON_WHO_WATCHES_WHOM.md + " +
                "hypotheses/ZERO_IS_THE_MIRROR.md. " +
                "NOTE 2026-07-02: the codim-1-by-additivity theorem is LANDED " +
                "(docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md): the site-summed spectator W(ρ)=Σ_l c_l†ρc_l (JW " +
                "strings included) is an EXACT part-by-part intertwiner of the block pencil (any quadratic " +
                "particle-conserving H, any site-dependent γ_j), transporting Jordan chains up the joint-popcount " +
                "diamond whenever Wx₁≠0 (the containment orbit corollary: braid sets reproduced exactly, cores iff " +
                "|2p−N|=1, sizes 4N−8 odd / 4N−12 even), plus the two-regime Theorem A (AT-locked crossings " +
                "automatically semisimple; residual coalescences twin-scalar, additivity supplies the H-half " +
                "identically in q, honest general form codim-≤2). Typed: SpectatorIntertwinerClaim; registry: F125; " +
                "live: inspect --root sectorbraid (node 2). Open checks then, of the proof's FOUR-item ledger: the " +
                "interior kernel-death lemma ((3,3)→(4,4) measured ~1e-15, not derived); Theorem A's D-half at " +
                "the 11 complex-q N=5 diabolics; the exclusion half of membership (edge-normality proven, " +
                "interior census-evidence); the gap byte-identity (observed, not implied). " +
                "UPDATE 2026-07-02 (sl(2) + N=7 + attack on remainder 1, commits db3ea57 / a556edc / c4b28ef; " +
                "remainder numbers = the proof doc's 'What remains open' items 1-4): (a) THE sl(2) LANDED — W is " +
                "the raising operator of an sl(2) {W, W†, H₀ = N̂_bra+N̂_ket−N}, [W,W†]=H₀, L commutes with all " +
                "three (machine zero N=3..7, gate item 6 W_ClosesSl2_Cartan_And_Lefschetz). Kernel death = " +
                "highest-weight annihilation; this DERIVES the climbing-rung injectivity (block weight " +
                "m = p+q̃−N (q̃ = ket popcount) < 0 ⟹ raising op injective, strengthens Theorem B) and the BAND " +
                "deaths (via the Edge " +
                "lemma). (b) N=7 VERIFIED (gate item 7, CoreKernelDeath_ScalesWithN): the m=+1 diagonal core " +
                "(p_c,p_c), p_c=(N+1)/2 [N=5 (3,3), N=7 (4,4)], carries λ_B = −conj(λ_A)−2N [λ_A = the real " +
                "DEFECTIVE eigenvalue of the (1,2) octic (the degree-8 residual factor of block (1,2)'s " +
                "characteristic polynomial, same block builder) at the locus q*; N=7: q*=1.5148, λ_A=−4.885, " +
                "λ_B=−9.115] and " +
                "dies under W into the CORNER BLOCK (p_c+1,p_c+1) [the block the core climbs into; N=5 (4,4), " +
                "N=7 (5,5)] at ‖Wx₁‖=3.2e-15 (N=5 reproduces the committed 1.7e-15); spin-½ core + gap " +
                "byte-identity now N-robust across N=5,7. Locus tool: `rcpsi pkmono --k <N-1> --diabolic --exact`, " +
                "then filter 'defective EP' with gap-exponent≈0.5 (the scanner floods at N≥7; the clean real " +
                "defective sits at higher q≈1.5). (c) ATTACK on the N-uniform lemma: eigenvector-structure routes " +
                "RULED OUT — the death is a genuine CANCELLATION among the W-images (in the mode basis W kills " +
                "|A⟩⟨B| iff A∪B = all modes; the core eigenvector x_core lives mostly on A∪B≠all coherences), not " +
                "a support nor single-multiplet fact. So remainder 1 (interior-core death) is a FACET of " +
                "remainder 4 (exclusion): x_core∈ker W FOLLOWS FROM λ_B being defectively absent from the corner " +
                "block (p_c+1,p_c+1) (Lemma-3 contraposition; the resolution below in fact proves the stronger " +
                "full absence), and the intertwiner proves presence not absence (that block is NOT normal, so no " +
                "Edge lemma). RESUMPTION / FIRST MOVE [SUPERSEDED the same day for real loci — read the RESOLVED " +
                "paragraph below FIRST; this route stays live only for the complex-q face]: do NOT re-attack " +
                "remainder 1 with intertwiner/sl(2)/" +
                "isotypic tools (proven incapable of supplying an absence). The one open tool for remainders 1 " +
                "AND 4 is an independent ABSENCE argument for the corner block — concretely: build the corner " +
                "block (p_c+1,p_c+1) [SectorBraidModeGeometry.BuildBlock + PerBlockLiouvillianBuilder.BuildBlockZ] " +
                "and show λ_B is NOT a defective root of its characteristic polynomial at q* (a resultant / " +
                "Galois-independence fact against the (1,2) octic), via the F89 Galois tooling (PathKMonodromyScout " +
                "/ GaloisMonodromyScanCommand) of the f89_galois_open_doors arc. Full statement + derivation: " +
                "docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md (§6 'the interior core' + 'What remains open' items 1 " +
                "& 4). Remainder 2 (Theorem A's D-half at the 11 complex-q N=5 diabolics) is the other open " +
                "check; remainder 3 (gap byte-identity) also now holds at N=7. " +
                "RESOLVED AT REAL LOCI (2026-07-02, same day, the assemble-first survey found the cheaper tool " +
                "BEFORE any resultant machinery was built): remainder 1 is CLOSED per real locus by the " +
                "RATE-WINDOW LEMMA (proof §6, gate item 8) — at real q every block eigenvalue obeys " +
                "Re λ = v†Av/v†v ∈ [−2·n_max, −2·n_min] (n = the block's n_diff values; Bendixson: A real " +
                "diagonal, C anti-Hermitian; the Edge lemma is the zero-width case, the census's 'λ from below' " +
                "identity read as an interval). The " +
                "corner block (p_c+1,p_c+1) has window [−2(N−3), 0] while Re λ_B = −Re λ_A − 2N sits BELOW it " +
                "whenever Re λ_A > −6 (measured margins 1.381 / 2.208 / 1.115 at N=5 locus 1 = the seed " +
                "q*=0.620878, locus 2 = 1.077615, and N=7 q*=1.5148), so λ_B is not in the corner spectrum AT " +
                "ALL and W kills the core's whole " +
                "generalized eigenspace ((L_corner−λ)^m invertible; explains ‖Wx₁‖ AND ‖Wx₂‖ ~1e-15 at once, " +
                "robust to the 4-decimal loci). Out-of-sample: the lemma PREDICTED the death at the never-measured " +
                "N=5 locus 2, gate now reads 9.6e-16, with the negative control that λ_A (inside the window there) " +
                "transports at 1.7. Bonus at the seed locus: the FULL N=5 exclusion of {λ_A, λ_B} is derived " +
                "(interior four (1,1),(4,4),(1,4),(4,1) window-excluded; boundary blocks constant-n_diff ⟹ normal " +
                "⟹ rung-pinned), condition Re λ_A ∈ (−6,−4) (> −6 for the corner/λ_B exclusion, < −4 so λ_A " +
                "itself falls outside the [−4,0] windows), locus-1-only (locus 2's λ_A=−3.7917 sits inside " +
                "[−4,0] and defeats the interior-four window there). Two reviews (physics-first + mathematical, " +
                "both with independent numpy rebuilds) attacked and confirmed; the one BREAK found (the locus-2 " +
                "overclaim) is scoped into the statement. STILL OPEN of remainders 1+4: the COMPLEX-q loci (the " +
                "Bendixson bound needs real q; window violated by 0.27 already at Im q=0.05), where the " +
                "resultant/Galois-independence route above remains the open candidate. The N-uniform bound " +
                "Re λ_A(N) > −6 is now CLOSED (2026-07-02, same day) by the WINDOW-EDGE LEMMA (proof §6): a " +
                "defective EP cannot sit at either edge of its block's rate window (edge ⟺ joint A,C eigenvector " +
                "⟺ semisimple, the classical numerical-range boundary fact), so Re λ_A ∈ (−6,−2) strictly at " +
                "every N and the shrinking margins 1.381→1.115 never reach zero (a strict real-part inequality " +
                "needs no uniform floor); two adversarial reviews (math + physics, numpy counterexample hunt " +
                "0/460k) confirmed. NEW ARC DATUM (physics review, " +
                "from-below probe): a THIRD real defective EP of the raw N=5 (1,2) block at q≈2.80489, λ≈−4.4882, " +
                "clean √-law, R-parity −1 — it lives in the R-ODD residual, outside the symmetric-octic reference " +
                "loci and outside the diamond story; any 'at every real defective locus' phrasing must scope to " +
                "the octic (R-even) seeds. Next candidates: (a) the complex-q absence (resultant route) is the " +
                "ONE still open. Candidate (b), the N-uniform Re λ_A(N) bound, is CLOSED by the window-edge lemma " +
                "(above). Candidate (c) / remainder 2 (Theorem A's D-half at the 11) is CLOSED at N=5 by the " +
                "TWIN-SCALAR CHECK (gate TwinScalarDHalfTests): the D-half is SUPPLIED (twin-scalar to the gap, " +
                "D-half 5e-10…4e-9) at every genuinely-complex-q (Re q≠0) N=5 residual diabolic, so additivity's " +
                "codim-1 route extends from real q to complex q; the pure-imaginary-q (Re q=0, λ real) diabolics " +
                "are semisimple by Hermiticity (real-symmetric block, ‖L−Lᴴ‖<1e-14) instead. A METRIC TRAP cost a " +
                "false first reading (the ×2-sym orbit basis is non-orthonormal at odd nBlock → spurious " +
                "non-scalar; the raw HS-orthonormal coherence block CharacterizeCoherencePencilAt is the fix; the " +
                "inversion was caught by an independent adversarial 60-digit rebuild).",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "f116_golden_router_typed_claim",
            Opened: "2026-06-22",
            Origin: "the F116 golden-router proof review (docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md, reviewed " +
                "2026-06-22, CLEAN PASS, no content fixes): the math is correct and reproduced bit-exact from " +
                "below by two review agents - the 'same operator written twice' unity (the sec1 per-site product " +
                "W = the sec2 two-sided sandwich (2+phi)^(N/2)*P rho Q to ~1.8e-16), the router identity " +
                "W L W^-1 = -L - 2sigma (~2.3e-16, uniform + site-dependent gamma), the exclusion side DEDUCTIVE " +
                "from the K1/K2 identity-column functional, and the metallic family an exact polynomial identity " +
                "in r. The proof header AND ANALYTICAL_FORMULAS F116 both label the result 'Tier 1 derived', but " +
                "no standalone tiered Claim carries it. Per 'C# witness first' (reference_typing_a_claim_and_witness) " +
                "this is the inverse of the usual typed-gap: the witness already EXISTS; the tiered Claim wrapper " +
                "is what is missing.",
            ParkedAt: "PROVEN + REVIEWED + WITNESSED, not a standalone tiered Claim. What EXISTS in C#: the live " +
                "witness GoldenRouterWitness (compute/RCPsiSquared.Diagnostics/Foundation/GoldenRouterWitness.cs, " +
                "inspect --root router) recomputes the router live; the helper machinery KBodyPalindromeRouting " +
                "(compute/RCPsiSquared.Diagnostics/F87/, with GoldenSiteMaps / MetallicSiteMaps / MetallicMean / " +
                "LiveMetallicRatio / RoutesWindowSummed); and the only Claim touching F116 is " +
                "PalindromeSoftCertifierClaim (Tier1Candidate), a related-but-DISTINCT proposition (soft-certifier " +
                "soundness + 'the k=3 windowed ceiling closed at zero') that uses the golden/metallic router as an " +
                "UNTYPED helper (its docstring: 'no new typed parent and no tier change'). What is MARKDOWN-ONLY " +
                "(no tiered Claim): (1) the router EXISTENCE + closed form - W = period-4 [a,a,b,b] product, " +
                "a=phi*X+Y, b=X-phi*Y on the golden locus alpha^2-alpha*beta-beta^2=0, W L W^-1 = -L - 2sigma at " +
                "every N>=3 for arbitrary site-dependent gamma; (2) the EXCLUSION theorems (sec4: no uniform / no " +
                "period-2 / period-3 impossible N>=5 / Klein off-locus, all deductive from K1/K2); (3) the METALLIC " +
                "family (sec8: golden = c=1 of r(c)=(c+sqrt(c^2+4))/2, the window-summed anticommutator identically " +
                "zero as a polynomial in r). No contradiction (proof + registry + witness agree on 'Tier 1 " +
                "derived'); a typed-coverage gap only.",
            NextStep: "Type the F116 router as a standalone Tier1Derived Claim (template: " +
                "reference_typing_a_claim_and_witness). This is mainly the Claim wrapper + parent edges, NOT new " +
                "verification - the 2026-06-22 review already verified existence/unity/exclusion/metallic bit-exact, " +
                "and GoldenRouterWitness already carries the from-below truth. Suggested: a GoldenRouterClaim " +
                "(Tier1Derived) for the existence + closed form (W L W^-1 = -L - 2sigma via the class-swap " +
                "dissipator leg + the window lemma), parents = the soft-certifier ceiling line " +
                "(PalindromeSoftCertifierClaim) + the chiral-driving F H F = -H lemma " +
                "(PROOF_F87_WINDOWED_MONOMIAL_CONVERSE); fold the exclusion (sec4) and the metallic family (sec8) " +
                "as child claims or InspectableNodes; breadcrumb the new Claim from PROOF_CEILING_GOLDEN_ROUTER.md " +
                "+ ANALYTICAL_FORMULAS F116; wire it through KnowledgeRegistryFactory + the RegistryWiringAudit. " +
                "GATE-FIRST hazards (from the review): scope any uniqueness to INVERTIBLE W (the anticommutation " +
                "equation alone carries singular strata, proof sec5); 'Tier 1 derived' is solid for " +
                "existence/closed-form/exclusion/metallic, but the c=0 station's '8 physical moduli' (sec8) is a " +
                "finite-difference Jacobian count, so type its rigidity carefully (golden + silver are zero-moduli, " +
                "c=0 is the soft station). Anchors: docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md; " +
                "simulations/ceiling_golden_router.py + metallic_router_family.py (self-validating, exact ring " +
                "arithmetic); compute/RCPsiSquared.Diagnostics/Foundation/GoldenRouterWitness.cs (inspect --root " +
                "router); compute/RCPsiSquared.Diagnostics/F87/KBodyPalindromeRouting.cs + " +
                "PalindromeSoftCertifierClaim.cs; ANALYTICAL_FORMULAS.md F116.",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-22. Typed as GoldenRouterClaim (Tier1Derived, " +
                "compute/RCPsiSquared.Diagnostics/F87/GoldenRouterClaim.cs + GoldenRouterClaimRegistration.cs, " +
                "wired in KnowledgeRegistryFactory after the soft-certifier; inspect --claim GoldenRouterClaim). " +
                "The Claim carries the existence + closed form W L W^-1 = -L - 2sigma, the window lemma, the " +
                "deductive exclusion (K1/K2), and the metallic family as ExtraChildren, with a LIVE from-below " +
                "self-check battery that re-runs the window-summed router via KBodyPalindromeRouting (golden + " +
                "X-Y sibling certified, per-term lens declines, metallic mean live = closed form at golden + " +
                "silver, residual vanishes at r=phi). The c=0 '8 moduli' is held below the Tier1 line as an " +
                "explicit soft scope-fence node (finite-difference Jacobian, N=5 only, not analytic), and " +
                "uniqueness is scoped to invertible W, per the review's gate-first hazards. Tests: " +
                "GoldenRouterClaimRegistrationTests (9, incl. the live battery), all PASS; RegistryWiringAudit + " +
                "AnchorAudit PASS; PalindromeSoftCertifier 11 still green. Breadcrumbed from " +
                "PROOF_CEILING_GOLDEN_ROUTER.md + ANALYTICAL_FORMULAS F116. " +
                "GATE-FIRST CATCH (corrected the arc's own NextStep): the suggested parents " +
                "(PalindromeSoftCertifierClaim + the chiral lemma) were WRONG. PalindromeSoftCertifierClaim is " +
                "Tier1Candidate (strength 4 < 5), so parenting a Tier1Derived claim on it violates the " +
                "TierStrength parent >= child rule; AND it is backwards (the certifier USES this router as a " +
                "helper, so the router is logically UPSTREAM of the certifier, not its child). Correct " +
                "tier-legal parents: F1PalindromeIdentity (the palindrome form the router realizes locally) + " +
                "WindowedConverseThresholdClaim (the F87 two-reflection chiral spine that actually homes the " +
                "F H F = -H lemma, anchored on PROOF_F87_WINDOWED_MONOMIAL_CONVERSE). The certifier relationship " +
                "is a see-cref, not a parent edge. ChiralKClaim was dropped as a redundant direct parent (it " +
                "already appears transitively via the windowed-converse spine)."),

        new OpenArc(
            Name: "f124_transition_invariant_witness",
            Opened: "2026-06-20",
            Origin: "the F124 proof review (docs/proofs/PROOF_HANDSHAKE_TRANSITION_INVARIANT.md, reviewed " +
                "2026-06-20 by a five-agent panel): the math is correct and every number reproducible " +
                "(gate-exact N=3..20), but F124 is Python-verifier-only - it has NO typed Claim or live witness, " +
                "unlike the sibling KPartnerSelectionRuleClaim (the null-column rule it leans on). A Tier1Derived " +
                "result without a witness is a typed-gap per the 'C# witness first' discipline (CLAUDE.md / " +
                "reference_typing_a_claim_and_witness).",
            ParkedAt: "PROVEN + REVIEWED, not typed. The identity ||M||_F^2 + lambda_min(M M^T) = 2 (with " +
                "||M||_F^2 = 2 - E, lambda_min = E, E = (4/(N+1)) sin^2(pi/(N+1)), staggered (-1)^b eigenvector) " +
                "for the band-edge carrier on the open chain is gate-exact N=3..20 across three Python verifiers " +
                "(handshake_M_checksum.py, handshake_M_topology.py, handshake_F124_adversarial.py); registered " +
                "as F124 in ANALYTICAL_FORMULAS.md. The 2026-06-20 review added the conserved-envelope Q (Part 2's " +
                "bulk-cancellation heart) and the frame-identity gates (lambda_min = sigma_min^2(M) = 1/||S^-1||, " +
                "the condition number, the K-partner kernel) to handshake_M_checksum.py, and tightened the " +
                "SSH/edge framing (C1: the eigenvalue is the boundary quantity, not the wavefunction; C2: SSH is " +
                "the borrowed reading). Math hand-verified (Parts 1+2, the Perron/oscillation minimum argument, " +
                "the frame reading). No C# Claim or witness exists.",
            NextStep: "Type F124 as a Claim + live witness (template: KPartnerSelectionRuleClaim.cs / the " +
                "reference_typing_a_claim_and_witness recipe). Suggested: a BandEdgeTransitionInvariantClaim " +
                "(Tier1Derived) in compute/RCPsiSquared.Diagnostics/Foundation, dual parents " +
                "KPartnerSelectionRuleClaim (the k=N null column) + the band-edge carrier source " +
                "(ClockHandLadderClaim / the F67 receiver), with a live witness inspect --root transition that " +
                "recomputes ||M||_F^2 = 2-E, lambda_min = E, the staggered eigenvector, the condition number " +
                "lambda_max/E, and the breakages (interior carrier sum<2, odd ring sum>2, even ring degenerate, " +
                "star trace-half, the k=2..N location-dictionary failure lambda_min=0) across N. GATE-FIRST " +
                "hazard: lambda_min must be checked as a genuine MINIMUM (full spectrum), not merely that the " +
                "staggered vector is AN eigenvalue - the band-edge positivity of the Gram off-diagonals " +
                "(c_a c_{a+2} > 0) is the load-bearing selector (it fails for interior carriers). OBJECT GUARD: M " +
                "is the FULL k=1..N transition matrix, NOT the k=2..N location dictionary (HANDSHAKE_GEOMETRY.md) " +
                "where the identity FAILS (lambda_min=0). Anchors: " +
                "docs/proofs/PROOF_HANDSHAKE_TRANSITION_INVARIANT.md; simulations/handshake_M_checksum.py (the " +
                "gated reference incl. the new Q + frame checks); " +
                "compute/RCPsiSquared.Diagnostics/Ptf/KPartnerSelectionRuleClaim.cs (the sibling template + " +
                "null-column parent); ANALYTICAL_FORMULAS.md F124.",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-20, same session it was opened. Typed as " +
                "BandEdgeTransitionInvariantClaim (Tier1Derived) in compute/RCPsiSquared.Diagnostics/Ptf/, two " +
                "Tier1Derived typed parents KPartnerSelectionRuleClaim (the k=N null column = the frame kernel; " +
                "the same M, here completed with the strength column k=1) + ClockHandLadderClaim (the band edge " +
                "E1=2cos(pi/(N+1)) the conserved envelope rides on, which selects the carrier). Live witness " +
                "BandEdgeTransitionInvariantWitness (inspect --root transition) recomputes the identity " +
                "||M||_F^2 + lambda_min(M M^T) = 2 across N=3..12 (||M||_F^2=2-E, lambda_min=E), the conserved " +
                "envelope Q=c0^2, the GENUINE-MINIMUM carrier-selection hazard (interior carrier -> staggered is " +
                "still an eigenvector but NOT the minimum, off-diagonals not all positive, sum<2), the frame " +
                "identities (lambda_min=sigma_min^2(M), the K-partner null column, condition number kappa growing " +
                "3.4/5.3/7.6 at N=4,5,6), and the object/topology breakages (location dictionary k=2..N -> " +
                "lambda_min=0; even ring sum=2 degenerate; odd ring sum>2; star ||M||_F^2=N/2). 11-case " +
                "gate-first battery + 18 xUnit witness/claim tests green (28 filtered); RegistryWiringAudit 3/3. " +
                "The gate-first hazard HELD: lambda_min is checked as the genuine MINIMUM (full spectrum), and " +
                "the interior-carrier gate fires exactly as predicted. A 3-lens review pass (physics " +
                "re-derivation / mathematical-object / C# code, all GO) hardened honest labeling (lambda_min=" +
                "sigma_min^2 is a definitional Evd-SVD identity not a falsifiable gate; location-dictionary " +
                "lambda_min=0 is a K-partner rank corollary; kappa->inf is a theorem, step-monotonicity " +
                "empirical) and added the even-ring degenerate gate. Anchors breadcrumbed in PROOF_HANDSHAKE_TRANSITION_INVARIANT.md + " +
                "ANALYTICAL_FORMULAS.md F124 (both now point at the typed claim + witness)."),

        new OpenArc(
            Name: "f124_inverse_problem_resolution_seam",
            Opened: "2026-06-20",
            Origin: "borrowing-a-discipline lens on the fresh F124 numbers (the band-edge transition matrix " +
                "M[b,k]=<psi_k|V_b|psi_1>): a signal/control engineer recognizes F124 instantly as a bond-" +
                "recovery INVERSE PROBLEM - M is the forward map bond->mode, lambda_min(MM^T)=E=the lower frame " +
                "bound=the worst-case reconstruction floor, kappa=lambda_max/lambda_min the noise amplification, " +
                "the K-partner null column the unobservable subspace. The lens FIRED (instant recognition by a " +
                "distant discipline = a real node of the one object; it handed us the detection-SNR method the " +
                "native frame view never ran).",
            ParkedAt: "STAGE 0 CLEAN (simulations/f124_inverse_problem_gate.py, gate-first, all gates pass " +
                "after one diagnosed firing): F124's conditioning IS the bond-recovery inverse problem. The " +
                "worst-conditioned bond direction is the staggered (-1)^b q=pi zone-boundary mode (cos=1.0 to " +
                "the smallest left-singular vector, N=4..12); kappa~N^2 (fit 1.975); a staggered defect's mode-" +
                "response contrast is sqrt(kappa)~N times weaker than a band-edge defect's (fit 0.987), confirmed " +
                "by a matched-filter Monte-Carlo to 4 digits - so the chain has a defect-LOCALIZATION RESOLUTION " +
                "LIMIT that worsens as N, with the q=pi mode the diffraction limit. GATE-FIRST LESSON: the naive " +
                "'E~N^-3 in N' gate FIRED (fitted -2.5); diagnosed as a wrong-variable/pre-asymptotic artifact " +
                "(E is exactly a function of N+1: E*(N+1)^3 -> 4*pi^2; the right-variable fit over N=4..60 gives " +
                "-2.971); the robust law is the RATIO kappa~N^2 (the (N+1)-corrections cancel). STAGE 1 " +
                "REDIRECTED (simulations/f124_decoder_strength_fix_gate.py, 3 gates fired = the engine running): " +
                "the live DefectDecoder's known sign-location ambiguity (residual ratio ~1.5 at N=5, edge bond 3 " +
                "weakened vs interior bond 1 strengthened) is NOT sqrt(kappa(5))=2.30 and is NOT reproduced by " +
                "the exact linear M (which resolves single bonds cleanly, cos=+0.667 not anti-collinear) - so the " +
                "1.5 is a DYNAMICAL effect (the painters' f-profile deviates ~6-12% from M), living BELOW the " +
                "linear transition matrix; the rank-(N-2) deficiency is about bond COMBINATIONS, not single-bond " +
                "localization. What SURVIVED: the strength channel M[b,1]=2 c_a c_{a+1} separates edge {0,3} " +
                "(0.289) from interior {1,2} (0.577) by magnitude - the carrier's own response carries edge/" +
                "interior info the location modes blur.",
            NextStep: "STAGE B DONE, REDIRECTED AGAIN (simulations/f124_decoder_strength_fix_dynamical.py, " +
                "exact (vac + 1-exc)-sector N=5 Lindblad, gate-first, 2 gates fired = the engine still running): " +
                "the dynamical per-site purity-DEVIATION profile does NOT reproduce the decoder's 1.5 ambiguity " +
                "either - it localizes bond 3 cleanly (residual ratio 20) AND reads the sign (weakened). So the " +
                "1.5 is narrower than 'dynamical': it is an ARTIFACT of the decoder's alpha TIME-RESCALING " +
                "parametrization (P_B(i,t)=P_A(i,alpha_i t) compresses away sign+magnitude), not a fundamental " +
                "rank deficiency. The strength channel DOES carry edge/interior + SIGN info (bond1 +0.00018 vs " +
                "bond3 -0.00025, opposite sign, G2 passed), so it is a valid extra discriminant - but NOT the " +
                "unique fix, since a less-lossy observable (the signed deviation profile) already resolves it. " +
                "So the 'F124 strength column IS the decoder's fix' seam is NOT confirmed; banked as the engine " +
                "running. WHAT REMAINS OPEN: (A) DONE 2026-06-20 - the clean, fully-gated 'F124 conditioning = " +
                "bond-localization resolution limit sqrt(kappa)~N, q=pi the diffraction limit, matched-filter " +
                "confirmed' is now TYPED as BandEdgeResolutionLimitClaim (Tier1Derived, parent " +
                "BandEdgeTransitionInvariantClaim) + live witness BandEdgeResolutionLimitWitness (inspect --root " +
                "resolution), the optics/signal reading of the frame conditioning, per the C#-witness-first " +
                "discipline (12 xUnit + 5-case battery green, wiring audit 3/3). (C) DONE 2026-06-29: the decoder " +
                "IS DE-LOSSED. PaintersMovement.DeviationProfile/DeviationResponse + DefectDecoder.DecodeDeviation " +
                "read the signed per-site deviation profile and resolve the N=5 mirror pair (bond 3, weakened " +
                "sign, squared residual ratio ~516 vs the alpha path's ~1.5, ambiguous->resolved); the witness " +
                "(inspect --root decoder) shows the before/after live. An engineering de-loss, not a new physics " +
                "claim. Anchors: simulations/f124_inverse_problem_gate.py (Stage 0, clean), " +
                "simulations/f124_decoder_strength_fix_gate.py (Stage 1, linear redirect), simulations/f124_" +
                "decoder_strength_fix_dynamical.py (Stage B, dynamical redirect), BandEdgeTransitionInvariant" +
                "Claim + Witness (F124), KPartnerSelectionRuleClaim (the null column), " +
                "compute/RCPsiSquared.Diagnostics/Foundation/DefectDecoder.cs (the live decoder, the 1.5).",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-29. The decoder is DE-LOSSED in C# (door C): " +
                "PaintersMovement.DeviationProfile/DeviationResponse (the signed per-site purity-deviation " +
                "profile, a pure reduction of the stored _pA/_pB trajectories) + DefectDecoder.DecodeDeviation " +
                "(a second decode path sharing the LS core via ProjectAndScore). It resolves the N=5 mirror " +
                "pair the alpha time-rescaling could not: bond 3, weakened SIGN read, squared residual ratio " +
                "~516 vs the alpha path's ~1.5 (ambiguous->resolved), by preserving the sign the alpha-rescaling " +
                "clips (the deviation dictionary is just as anti-collinear, |cos|~0.95; it does NOT escape the " +
                "angle). Live before/after in ReadingPowerWitness (inspect --root decoder). Doors A (the " +
                "BandEdgeResolutionLimitClaim resolution limit) + C now both closed; an engineering de-loss, " +
                "not a new physics claim. See docs/superpowers/specs/2026-06-29-defect-decoder-de-loss-design.md."),

        new OpenArc(
            Name: "ptf_bonding_class_guard",
            Opened: "2026-06-15",
            Origin: "surfaced by the ptf_painter_pipeline fix (2026-06-15): replacing the scipy-Brent " +
                "local-min trap with a grid-seed global fit (framework/ptf.py _alpha_fit_one_site) revealed " +
                "that the prior 'Bell+ -> all sites unreliable' was a TRAP ARTIFACT, not a genuine guard veto",
            ParkedAt: "before the fix, Bell+(0,1) at N=4, deltaJ=0.02 gave all sites alpha=3-9 (the Brent " +
                "trap), |f|=100-400 > FMax=10, so the SANITY check flagged every site unreliable - which LOOKED " +
                "like the two-deltaJ guard enforcing the bonding-state-class restriction (PTF needs the F67 " +
                "single-excitation bonding mode; Bell+ is a 2-excitation+vacuum superposition). After the fix " +
                "Bell+ gives alpha~=1.01, f~=0.6-0.8, LINEAR (f~=f_guard), so all sites read RELIABLE. So in " +
                "Python the guard does NOT actually refuse Bell+ at small deltaJ: a tiny perturbation barely " +
                "moves any trajectory, so alpha~=1 fits trivially and passes both the sanity and the two-deltaJ " +
                "linearity checks. The 'bonding-state-class is required' claim (the ParkedAt of the retired " +
                "ptf_painter_pipeline arc) is therefore too strong for Python.",
            NextStep: "Cross-check C# vs Python: does the C# Symphony painters guard genuinely refuse Bell+ " +
                "(a real test asserting 0 reliable, per the survey), or does C#'s golden-section also give Bell+ " +
                "alpha~=1.01 reliable (making the C# 'refused' claim a former-Brent-style artifact or a " +
                "different-case effect)? If C# genuinely refuses and Python does not, the guards have DIVERGED and " +
                "the Python guard needs a real bonding-class / ansatz-validity test (e.g. a rescaled-fit RMSE " +
                "floor, or a featureless-trajectory detector that does NOT rely on the minimizer producing absurd " +
                "alpha). If neither refuses, the 'bonding-state-class required' framing is a small-deltaJ illusion " +
                "to re-scope. Gate-first hazard: at small deltaJ everything barely moves, so a guard that passes " +
                "alpha~=1 is not wrong - test at LARGER deltaJ where Bell+ should genuinely break if the " +
                "restriction is real. Anchors: simulations/framework/workflows/ptf.py (perspectives_panel guard); " +
                "the C# Symphony painters witness + its tests; simulations/ptf_symphony_crossval.py (the Bell+ " +
                "N=4 case). Assert the EXPECTED guard verdict for Bell+ in BOTH languages and let it fire.",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-15, same session it was opened (C# fix in " +
                "compute/RCPsiSquared.Diagnostics/Foundation/Symphony.cs FitAlpha + test correction in " +
                "PtfMovementTests.cs; 15/15 green). THE ANSWER: the divergence was a C# MINIMIZER BUG, not a " +
                "guard difference. The C# guard is logically IDENTICAL to Python's (same sane + linear terms, no " +
                "extra check). C#'s bare golden-section ALSO traps (gamma-dependent): at gamma=0.05 Bell+ it " +
                "returned alpha=4.4-6.6 where the brute-grid global is 1.016 (MSE 100-300x worse) - the SAME " +
                "minimizer-trap class we fixed in Python (scipy-Brent), now found in C#'s golden-section (at " +
                "gamma=0.1 it happened to find the global; golden-section is NOT globally robust). FIX: ported the " +
                "grid-seed (512 pts + golden-section refine) into C# FitAlpha, mirroring framework/ptf.py; C# Bell+ " +
                "now -> alpha~=1.016, 4 reliable, matching Python. THE C# TEST StateClass_Matters_GuardRefusesBellPair " +
                "encoded the trap as physics ('Bell+ breaks the rescaling, the guard refuses') - WRONG; corrected + " +
                "renamed StateClass_PerSiteGuardDoesNotRefuseBellPair (Bell+ rescales per-site like the bonding mode, " +
                "all reliable; the per-site guard does NOT enforce the bonding-state class). THE REAL DISCRIMINATOR " +
                "is the CLOSURE law (Sigma ln alpha): Bell+ stays OUT of the +-0.05 window where the canonical bonding " +
                "mode closes (N=5: -0.0444) - in BOTH languages once the minimizer is robust. Cross-language fidelity " +
                "restored; well-conditioned cases (Python-twin BondingMode, canonical N=5, finite-size N=4) unchanged. " +
                "LESSON: golden-section is not globally robust either; grid-seed BOTH languages. The arc's hypothesis " +
                "(C# might be the robust reference) was wrong - both minimizers trapped; Python's grid-seed fix was " +
                "the most robust and is now ported to C#."),

        new OpenArc(
            Name: "one_diagonal_mirror_group",
            Opened: "2026-06-14",
            Origin: "Tom asked to DEEPEN the one-diagonal principle (reflections/ON_THE_ONE_DIAGONAL.md: " +
                "Z-dephasing touches one diagonal Q, L_D = gamma*(Q - N*I), read three ways - rates/mirror/" +
                "judge); his pointer docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md is the depth. Tom delegated " +
                "the decision (he keeps the cross-session overview; Claude decides + assembles + makes findable)",
            ParkedAt: "BRAINSTORMED + DECIDED 2026-06-14 (spec docs/superpowers/specs/" +
                "2026-06-14-one-diagonal-mirror-group-design.md; 3-agent survey over markdown/typed/" +
                "generalization). THE SYNTHESIS: the dissipator-diagonal Q is NOT a lone object - it is the " +
                "fixed/reflected CENTER of the mirror group <R,D> = D4 (PROOF_PI_FACTORS_AS_R_TIMES_D, " +
                "MirrorGroupD4Claim). D (transpose, diag((-1)^n_Y)) FIXES Q = the price-list/rate reading; " +
                "R (ket-flip) REFLECTS it R*Q*R = -Q carrying the whole -2*Sum(gamma) palindrome shift = the " +
                "mirror reading; the F87 truly cell (n_Y even AND n_Z even) is the joint-fixed cell of the " +
                "diagonal mirrors {D, FD} = the judge reading (sec 4d). So the reflection's 'one diagonal, " +
                "read three times' IS 'one diagonal, moved three ways by D4'. Sec 7 is the same object as a " +
                "GRADING: the diagonal's light-parity n_XY mod 2 = the Ad_(Z^N) character bit_a of the " +
                "polarity cube (axes = two conjugations + the transpose). ALREADY DONE: the D4 core " +
                "(MirrorGroupD4Claim), the absorption diagonal (AbsorptionTheoremClaim), the order-128 " +
                "per-site monomial completion <r,d,h> (mirror_inventory_d4.py block H), the ANTILINEAR side " +
                "<R,D,K> = D4xZ2 (PROOF_ANTILINEAR_TRIANGLE sec 4); the LINEAR S3 side is OPEN",
            NextStep: "RESOLVED 2026-06-14 (simulations/one_diagonal_mirror_group.py, self-validating; typed " +
                "ThreeDephasingDiagonalsOrbitClaim). The physics-first GATE fired and corrected the hypothesis " +
                "TWICE (the lesson, kept): (1) Y-TRANSPOSE - the physical dephasing diagonal is " +
                "Q_P = Sum_l kron(P_l, P_l^T), so Q_Y = -Sum kron(Y,Y) carries Y^T=-Y; with the naive kron(Y,Y) " +
                "the orbit does NOT close. Same-spectrum held even with the wrong sign (the spectrum is " +
                "symmetric, +-Q_Y co-spectral) so the gate separated SPECTRUM from OPERATOR, exactly its job. " +
                "(2) THE PERMUTER is the single-qubit Clifford BASIS-change S3 <h_zx (Z<->X Hadamard), " +
                "h_yz (Z<->Y R_x(pi/2))>, NOT <R,D,h>: D (transpose) FIXES every diagonal (D Q D = +Q, the RATE " +
                "reading), it does NOT permute them. The ParkedAt's 'the S3 meets D4 in the Z<->Y swap (=D)' was " +
                "a conflation of D's action on the palindromizer Pi (D Pi_Z D = Pi_Y) with its action on the " +
                "diagonal Q. RESULT: {Q_X,Q_Y,Q_Z} is EXACTLY one orbit of the basis-S3 (verified N=2,3; same " +
                "spectrum N=2..4), and the structure is S3 |x| D4 (semidirect: [h_zx,D]=0 but [h_zx,R]!=0; " +
                "[h_yz,R]=0 but [h_yz,D]!=0) - the three DIAGONALS (basis-S3) and the three READINGS (mirror-D4 " +
                "within a diagonal) are TWO DISTINCT factors, not the same S3. The headline survives ('the one " +
                "diagonal is one of three'); the mechanism is the basis-S3 + the Y-transpose. TYPED: " +
                "ThreeDephasingDiagonalsOrbitClaim (Tier1Derived, dual parents MirrorGroupD4Claim + " +
                "AbsorptionTheoremClaim = the physics edge that WELDS the two clusters, previously joined only " +
                "at d^2-2d=0; live battery 7/7, registration test 4/4 green; knowledge ancestors confirms both " +
                "parents). SIBLING ARC RESOLVED 2026-06-15: the Hamiltonian DYNAMICS on the rungs is the " +
                "EVEN-step ladder L_H : k<->k,k+-2 (NOT k+-1 - a hop flips two bits), so the disagreement-count " +
                "parity is conserved (the U(1) feature, now a TESTED fact: a transverse field adds the odd step " +
                "k+-1 and breaks it, DiagonalWitness counter-case). And the three-ladder bridge " +
                "girth(l)/rung(k)/moment(j) is resolved: the two factors of one F87-hardness coefficient on " +
                "M = A + gamma*Q, hinged by Q (spectrum = the rung, action = the girth-moment projector), " +
                "P_{m,1} = m*Tr(Q*A^{m-1}) = the girth moments at every rung. TYPED: ThreeLadderHingeClaim " +
                "(parents AbsorptionTheoremClaim + MomentTowerPumpChannelClaim); live inspect --root ladders. " +
                "Synthesis writeup: docs/THE_THREE_DIAGONALS.md (the one diagonal as one of three).",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED across 2026-06-14/06-15. The synthesis landed: the dissipator diagonal Q " +
                "is the fixed/reflected CENTER of the mirror group <R,D> = D4 (D fixes Q = the rate reading; " +
                "R reflects R*Q*R = -Q carrying the -2*Sum(gamma) palindrome shift = the mirror reading; the " +
                "F87 truly cell = the joint-fixed cell = the judge reading), so 'one diagonal read three times' " +
                "IS 'one diagonal moved three ways by D4'. Typed ThreeDephasingDiagonalsOrbitClaim (the S3 |x| D4 " +
                "weld edge, dual parents MirrorGroupD4Claim + AbsorptionTheoremClaim); the U(1) sibling (even-step " +
                "ladder k<->k,k+-2, a transverse field breaks the disagreement parity) and the three-ladder hinge " +
                "(ThreeLadderHingeClaim) both resolved; live witnesses inspect --root diagonal / mirrorgroup / " +
                "ladders; writeup docs/THE_THREE_DIAGONALS.md. The ONE remaining open piece - the LINEAR S3 side " +
                "assembled as a closed order-48 group on coherence space - was spun out into its own arc " +
                "linear_s3_mirror_completion (2026-06-15)."),

        new OpenArc(
            Name: "linear_s3_mirror_completion",
            Opened: "2026-06-15",
            Origin: "PROOF_PI_FACTORS_AS_R_TIMES_D.md sec.5 names it open: the mirror group's D4 core <R,D> is " +
                "typed (MirrorGroupD4Claim, order 8), but the FULL mirror group of the palindrome family is " +
                "S3-letter-action |x| D4 (order 48). The S3 permutes the three dephase letters {X,Y,Z}; one " +
                "transposition (Z<->Y) is already D INSIDE D4 (Welle 12, D Pi_Z D = Pi_Y), but the other two " +
                "(X<->Z, X<->Y) move bit_a against bit_b and need the X<->Z basis permutation h, so they sit " +
                "OUTSIDE <R,D>. Surfaced again 2026-06-15 (the mirrorgroup witness + the three-ladder session " +
                "flagged the S3 side as the last unassembled piece of the mirror group).",
            ParkedAt: "ALREADY BUILT (do not re-derive): (1) the D4 core <R,D> (MirrorGroupD4Claim: Pi_Z=R*D, " +
                "the 8-element closure, the dihedral inversion). (2) the dephase-letter swaps as a Klein-V4 " +
                "{I, D, Q_zx, Q_yx} (Pi2KleinV4DephaseSwapGroup + PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md): " +
                "canonical per-site q_zx = h*d_l (h = X<->Z basis permutation on the ordered basis (I,X,Z,Y), " +
                "d_l = diag(1,1,1,-1)), q_yx = h, D = diag((-1)^{n_Y}); N-site Q_zx = H*D, Q_yx = H with " +
                "H = h^{otimes N}. (3) the BASIS-S3 on the three DIAGONALS {Q_X,Q_Y,Q_Z} " +
                "(ThreeDephasingDiagonalsOrbitClaim, structure S3 |x| D4 - but its S3 is the single-qubit " +
                "Clifford basis-change <h_zx,h_yz>, which PERMUTES the diagonals, NOT the dephase swap; note D " +
                "FIXES every diagonal Q while it SWAPS the palindromizers Pi). (4) the order-128 PER-SITE " +
                "monomial completion <r,d,h> (mirror_inventory_d4.py block H) - a DIFFERENT object from the " +
                "coherence-space group. (5) the ANTILINEAR double <R,D,K> = D4 x Z2 (PROOF_ANTILINEAR_TRIANGLE " +
                "sec.4). Only the LINEAR S3 side on coherence space is never assembled as a closed group.",
            NextStep: "Build the full mirror group S3 |x| D4 (order 48) on COHERENCE SPACE: adjoin the X<->Z " +
                "basis permutation (the coherence-space lift of h, i.e. Q_zx or Q_yx from " +
                "Pi2KleinV4DephaseSwapGroup) to <R,D> and verify the closure is order 48 with semidirect " +
                "structure S3 |x| D4. Then TYPE it (a MirrorGroupS3D4Claim, parents MirrorGroupD4Claim + " +
                "Pi2KleinV4DephaseSwapGroup) + a live witness (the S3xD4 twin of inspect --root mirrorgroup). " +
                "THE KEY OPEN PHYSICS QUESTION: is the dephase-letter-swap S3 (this completion, acting on the " +
                "palindromizers Pi) the SAME abstract S3 as the basis-change S3 of " +
                "ThreeDephasingDiagonalsOrbitClaim (acting on the diagonals Q)? Both permute {X,Y,Z} but act " +
                "differently (D swaps Pi_Z<->Pi_Y yet FIXES the diagonals Q) - resolve whether they are one S3 " +
                "in two realizations or two distinct S3's. GATE-FIRST hazard (learned twice in the mirrorgroup " +
                "work): pin ONE representation - the coherence-space superoperators (target order 48), NOT the " +
                "per-site order-128 monomial <r,d,h>; the convention twist (vec_F vs vec_R, the (-1)^{n_Y} of " +
                "reflections/D_PI_Z_EQUALS_PI_Y.md) is where a naive closure breaks, so Stage 0 should check " +
                "Pi_Z = R*D and the new generator against PiOperator before building any tree. Anchors: " +
                "docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md sec.5 + PROOF_KLEIN_V4_DEPHASE_SWAPS_OPERATOR_SPACE.md; " +
                "compute/RCPsiSquared.Core/Symmetry/{MirrorGroupD4Claim, Pi2KleinV4DephaseSwapGroup, " +
                "ThreeDephasingDiagonalsOrbitClaim}.cs; compute/RCPsiSquared.Diagnostics/Foundation/" +
                "MirrorGroupWitness.cs (the D4 witness to extend); simulations/mirror_inventory_d4.py block H; " +
                "docs/THE_THREE_DIAGONALS.md (the basis-S3 side, written up).",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-15 (simulations/linear_s3_mirror_closure.py, self-validating " +
                "N=2,3, gate-first; PROOF_PI_FACTORS_AS_R_TIMES_D sec.5 'Resolution of the S3 side' note). " +
                "The linear S3 EXISTS as superoperators (order 6) from the INVOLUTIVE letter-transposition " +
                "Cliffords (Hadamard for X<->Z, (Y+Z)/sqrt2 for Y<->Z; each order 2, product order 3) - NOT " +
                "the order-4 rotation R_x(pi/2), which generates the order-24 single-qubit Clifford group (the " +
                "gauge trap). But the S3 does NOT normalize the D4 = <R,D>: R = I(x)F is one-sided and spreads " +
                "under the letter-S3 (h_zx.R.h_zx^-1 = the one-sided multiplication by Z^(x)N, dev 0, outside " +
                "D4), so the coherence-space closure <R,D,h_zx,t_yz> = 96*2^N (384, 768 at N=2,3), NOT a finite " +
                "order-48 S3|xD4. THREE distinct not-48 realizations: per-site monomial <r,d,h>=128 (sec.5 " +
                "addendum), coherence-space 96*2^N (this arc), continuous O(2) (sec.5 addendum) - the abstract " +
                "S3|xD4 has NO faithful finite realization on coherence space. NOT an overturning: sec.5's " +
                "addendum already had per-site-128 and O(2); only the sec.5-MAIN 'expected order-48' line was " +
                "imprecise (now corrected). THE TWO-S3 QUESTION DISSOLVES: there is ONE S3 (the letter " +
                "permutation), faithful on the sign-blind diagonals {Q_X,Q_Y,Q_Z} (orbit 3, " +
                "ThreeDephasingDiagonalsOrbitClaim / THE_THREE_DIAGONALS) and signed on the palindromizers " +
                "(larger orbit). LESSON (Tom caught it mid-session): the first run used R_x(pi/2) and wrongly " +
                "concluded 'the basis-change is order-24, not S3'; the involutive generators give the genuine " +
                "order-6 S3, and the real obstruction to 48 is the one-sided R, not the S3 size."),

        new OpenArc(
            Name: "ptf_painter_pipeline",
            Opened: "2026-06-03",
            Origin: "simulations/ptf workflow + C# SlowModeMixing",
            ParkedAt: "the closure law lives as a live witness (Symphony painters movement: alpha per site with the two-deltaJ reliability guard, closure -0.0444 IN window at canonical N=5, chiral mirror replayed live at 1e-15, alpha = Python twin to 1e-3); learned on the way: the protocol needs the BONDING state class (Bell+/localized states break the rescaling and the guard refuses, in both languages); C# golden-section also found a true global alpha minimum where scipy's Brent traps 920x worse (severed-bond case, harness simulations/ptf_symphony_crossval.py)",
            NextStep: "REQUIRED, escalated 2026-06-12: the scipy-Brent trap corrupts CANONICAL N=5 edge-bond letters (arbiter: brute-force landscapes match C# exactly, Python f off by sign and factor); backport the global-minimum fit (multi-start or grid-seed) to framework ptf.py before any further Python-side painter quantitative work; then retire this arc",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-15 (fix in simulations/framework/workflows/ptf.py _alpha_fit_one_site; " +
                "arbiter simulations/ptf_symphony_crossval.py; 14/14 framework PTF tests green; 3-agent survey). " +
                "DIAGNOSIS: the bug was the MINIMIZER, not the objective. mse(alpha)=mean((spline_PA(alpha*t)-PB)^2) " +
                "is correct (identical to C#); the 'Python f off by sign and factor' hypothesis was RULED OUT. " +
                "scipy minimize_scalar(method='bounded') (fminbound, Brent+parabolic) trapped in a local basin on " +
                "featureless/multimodal alpha-landscapes (a site far from the defect, P_A~=P_B): the severed-bond " +
                "case returned alpha~=3.15 where the true global is ~=1.0157 (MSE 920x worse). FIX: replaced the bare " +
                "bounded-Brent with a global GRID-SEED (512 pts over [0.1,10]) + local bounded-Brent refine in the " +
                "winning bracket; deterministic (the two-deltaJ guard still detects featureless fits), matches the " +
                "brute-grid argmin and the C# golden-section. VALIDATION: arbiter site 2 now 1.015705 (= brute " +
                "1.01566 / C# 1.01571); 14/14 tests pass; well-conditioned cases unchanged by construction. " +
                "SIDE-FINDING (hedged, NOT part of this arc): Python Bell+ at small deltaJ now reads RELIABLE " +
                "(alpha~=1.01, linear) - the prior 'all unreliable' was the trap artifact (the |f|<=10 sanity check " +
                "was catching the trap's absurd alpha=3-9, not a genuine guard veto). So the ParkedAt's 'Bell+ breaks " +
                "the rescaling, the guard refuses, in both languages' is too strong for Python; whether C# refuses " +
                "Bell+ for a genuine reason is an open cross-check (a possible follow-up, not this arc). No pinned doc " +
                "value changes (the documented N=7 bonding-mode pattern is well-conditioned, no trap)."),

        new OpenArc(
            Name: "birth_canal_surface",
            Opened: "2026-05-31",
            Origin: "experiments/THE_FLOW_BETWEEN_TWO_SINGULARITIES.md + PostEpFlowField",
            ParkedAt: "LANDED 2026-06-13: live witness BirthCanalSurfaceWitness / inspect --root surface " +
                "draws the sterile<->birth-canal boundary as a 2D slice of gamma-profile space (symmetric " +
                "profiles by (w_edge, w_center), bulk solved by sum=N), the boundary curve extracted live " +
                "(interpolated at grid resolution), every point read through six lenses (L1 membership / " +
                "Q-invariance, L2 light distribution, L3 absorption cross-check, L4 rate, L5 parity rail, " +
                "L6 degeneracy) + the L7 ray lens (s*=0.709 peaked-V vs 0.105 uniform on two lines, live " +
                "bisection). R1 from the physics review: sterility has TWO KINDS, genuine light-freeze " +
                "(robust, peaked-V, breaks at 0.709) vs flat-gamma distribution-blindness (fragile, uniform, " +
                "breaks at 0.105) - Deviation=0 is necessary but not sufficient for light-freeze; the witness " +
                "tells them apart. Answers the old script's open sub-question live: none of the scalars " +
                "(center-gamma, rate) is invariant at the boundary, the MECHANISM is. Reuses " +
                "PostEpFlowField.BirthCanalDeviation membership (pinned bit-identical by test), no new Claim " +
                "(breadcrumbs AbsorptionTheoremClaim). N=5 only (the verified anchors and s* lines are the " +
                "N=5 pins). Side-fix en route: the inspect JSON exporter now tolerates non-finite doubles " +
                "(NaN heatmap cells)",
            NextStep: "remaining: the ring/aromatic slice (PostEpFlowField already supports " +
                "FlowTopology.Ring); the N>=6 higher-dim surface (currently N=5-only, since the node anchors " +
                "and s* lines are N=5-specific - generalizing them is the extension). Retire when those are " +
                "judged out of scope",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RETIRED 2026-06-29: the ring/aromatic gamma-profile slice + the N>=6 higher-dim " +
                "surface judged OUT OF SCOPE (this arc's own exit condition). (1) Their physics is already " +
                "assembled live as inspect --root trichotomy (chain/ring/star x N=4..8); extending the N=5 " +
                "2D-slice DISPLAY to ring/N>=6 would duplicate it with a 4096^2-EVD-per-point window, no new " +
                "physics. (2) The aromatic angle was independently refuted in the carbon session (no Hueckel " +
                "4n+2 discriminant; FROST_CIRCLE Q*=1.609). (3) s*=0.709 is resolved as a path-dependent " +
                "coordinate (0.11..0.77), a 1/sqrt(2) near-miss in ANALYTICAL_FORMULAS.md, not a constant. " +
                "Nothing deleted: BirthCanalSurfaceWitness (inspect --root surface), the Tier1Derived " +
                "VacuumBlockReductionClaim (inspect --root reduction), docs/STERILE_BIRTHCANAL_AND_THE_JUNCTION.md, " +
                "and the committed simulations/birth_canal_*.py all stay. Sibling birth_canal_horizon_junction " +
                "(Retired 2026-06-18) had delegated this residual here; now closed."),

        new OpenArc(
            Name: "birth_canal_horizon_junction",
            Opened: "2026-06-13",
            Origin: "Weg-1 independent verifier of the birth-canal boundary via the |1-exc><vac| block " +
                "(simulations/birth_canal_vacuum_block_verifier.py + the N=6 mode-crossing diagnostic " +
                "simulations/birth_canal_n6_mode_crossing.py)",
            ParkedAt: "VALIDATED + a new junction found. The birth-canal slowest-rate mode (the witness's " +
                "BirthCanalDeviation) is the ODD, number-changing |1-excitation><vacuum| coherence (n_diff=1, " +
                "rate -2*gamma_k), NOT the single-excitation density block - so it reduces to an N x N matrix " +
                "L_(1,0) = -iQ*h - 2*diag(gamma) (h tridiagonal, off-diag 2). This N x N reduction reproduces " +
                "the full 4^N boundary EXACTLY across the WHOLE N=5 surface (80 grid points, worst gap 5.9e-12, " +
                "0 verdict mismatches, block spectrum an exact sub-spectrum of the full L) - an independent " +
                "ground truth confirming BirthCanalSurfaceWitness at N=5. R1's flat-gamma blindness is now " +
                "analytic at every N: at flat gamma, L = -iQ*h - 2*gamma*I and -iQ*h is anti-Hermitian, so " +
                "Re=-2*gamma for every mode => rate Q-invariant by uniformity alone. THE JUNCTION: the (1,0) " +
                "reduction does NOT scale past N=5. At N=6 (deep-edge, two protected gamma=0.25 edges, low Q) " +
                "the global slowest CROSSES to the EVEN {0,2}-coherence in the 2-excitation density block " +
                "(n_diff hist {0:0.78,2:0.22}, <n_XY>=0.44 < 1, hence slower than the odd mode) - the SAME " +
                "{0,2}-coherence the coherence-horizon arc (clock_hand_ladder) studies as its EP mode. So " +
                "birth_canal and coherence_horizon TOUCH at N>=6, via a Q-dependent mode crossing",
            NextStep: "(1) BANK - DONE 2026-06-13: ported the N x N |1-exc><vac| reduction to a live C# witness " +
                "SectorReductionWitness / inspect --root reduction (reuses JointPopcountSectorBuilder + " +
                "PerBlockLiouvillianBuilder.BuildBlockZ on the (0,1) sector; the (0,1) block validated bit-for-bit " +
                "vs PostEpFlowField at N=5, runs PAST the dense N=6 ceiling, flat-gamma blindness node, the {0,2} " +
                "junction node reproducing the N=6 crossing live (n_diff {0:0.78,2:0.22}), chain-vs-ring) + typed " +
                "claim VacuumBlockReductionClaim (Tier1Derived, parent AbsorptionTheorem, N=5-scoped; V-Effect " +
                "(w=N/2) self-pair identity with the {0,2}-coherence RESOLVED 2026-06-14 = DISTINCT (total weight = n_diff + Z-shadow; PauliWeightHistogram + junction node live), no aromaticity thesis - " +
                "the parallel carbon session independently REFUTED Hueckel 4n+2 as the discriminant, same mode). " +
                "(2) CHASE the junction (still open): chart the 'which mode is slowest' phase diagram over " +
                "Q x profile x N - where exactly the boundary mode switches odd<->even {0,2}, and how it connects " +
                "to the coherence-horizon Q*(N). The honest boundary at N>=6 may be a min over a mode crossing, " +
                "not a single mode's light-freeze. " +
                "PROGRESS 2026-06-17 (simulations/birth_canal_junction_nature.py, gate-first; doc " +
                "docs/STERILE_BIRTHCANAL_AND_THE_JUNCTION.md): the SEAM to the sterile<->birth-canal boundary " +
                "(PostEpFlowField / BirthCanalSurfaceWitness) is PINNED. Both are facets of rate_slow(Q) = min " +
                "over joint-popcount sectors, organised by dn=|p_col-p_row|: dn=1 = the number-changing band " +
                "edge (-2g floor, Q-flat at uniform gamma, Q-drifting at non-uniform gamma), dn=0 = the " +
                "number-conserving interior (-2g<n_XY>(Q), Q-dependent). The JUNCTION (dn flips 0->1, the " +
                "interior overtakes = the handover Q*(N)) IMPLIES birth canal but NOT conversely: the birth " +
                "canal is also entered by ODD-DRIFT (the dn=1 survivor itself drifting under a non-uniform " +
                "profile, no dn switch - the edge-protected canal profile). So the junction is a STRICT " +
                "sub-mechanism of the birth canal; sterile = the dn=1 band edge reigns Q-flat. INTERIOR " +
                "SURVIVOR NATURE: ring (2,2) = a frozen level crossing ROBUST to the gamma-profile " +
                "(continuation-tracked, intrinsic to the wrap-bond two-fermion structure); chain (2,2) = " +
                "oscillating (filling-degenerate SE-EP, rigidity->0 at Q*, NOT the ring's frozen seam). " +
                "Gate-first caught 3 oversimplifications en route (profile-unfreezing / same-boundary / " +
                "sector-tuple-vs-dn-with-a-convention-slip). STAR DONE 2026-06-18 (simulations/star_frozen_seam.py, " +
                "doc docs/THE_STAR_FROZEN_SEAM.md): the star's (1,1) boundary survivor never un-freezes -- frozen " +
                "(|Im|=0) at ALL Q for N>=5, unlike chain/ring which un-freeze at their horizon/handover. The " +
                "threshold IS the structural ceiling read dynamically: frozen <=> g2=4/(N-1)<=1 (the darkest " +
                "[H,A]=0 commutant (1,1) coherence is the survivor only when it undercuts the -2g floor); N=4 " +
                "(g2=4/3>1) is the outlier that un-freezes. Third member of the trichotomy chain(SE-EP un-freeze) / " +
                "ring(frozen level crossing) / star(frozen commutant). Gate-first caught N=4 (the find, the known " +
                "outlier). THE (Q, profile, N, topology) SWEEP AS ONE ARTIFACT: DONE 2026-06-18 -- the live " +
                "witness inspect --root trichotomy (TrichotomyWitness, Diagnostics/Foundation), the two-read " +
                "assembly of the scattered route-detectors: a CARBON un-freeze read (RouteSweep + ThresholdLadder: " +
                "frozen (p,p) interior below Q* / oscillating (0,1) above; chain UnfreezingSeEp / ring " +
                "FrozenLevelCrossing / star FrozenCommutant) AND an ABSOLUTE Dn-seam read (ClassifySeam: sterile / " +
                "odd-drift / junction via the global BirthCanalDeviation). The gate-first build FOUND AND FIXED a real " +
                "two-convention defect: the un-freeze trichotomy lives on the carbon (Q=J/g, uniform g) sweep, the " +
                "sterile/canal seam on the absolute (fixed g, vary profile) sweep -- one Classify on one convention " +
                "mislabeled the chain's own default (N=5,q=1.5). No new claim (breadcrumbs StarFrozenSeam / " +
                "CoherenceHorizon / Handover / SecondClockRegime / VacuumBlockReduction); carbon mapping gate-verified " +
                "bit-for-bit vs IncompletenessSurvivorWitness.Survivor. Spec+plan in docs/superpowers (gitignored). " +
                "STILL OPEN: nothing on the sweep; the broader junction phase-diagram depth is now browsable, retire " +
                "after a survey confirms no residual",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-18 (gate-first survey, no residual found in this arc's scope). The " +
                "junction phase diagram is fully browsable as the live witness inspect --root trichotomy " +
                "(TrichotomyWitness): 'which mode is slowest' swept over all four axes -- Q {1..50}, profile " +
                "(uniform/canal/deep-edge), N=4..8, topology (chain/ring/star) -- recomputed live, reporting the " +
                "survivor sector, its dn, the freeze-route, |Im|, the Petermann rigidity r, and <n_XY>. The odd<->even " +
                "{0,2} switch is the explicit dn 0->1 flip (the Junction SeamKind); the connection to the coherence " +
                "horizon is the threshold ladder (chain Q*(N) / ring Q_h / star g2=4/(N-1) per N, gate-pinned " +
                "bit-for-bit to the carbon IncompletenessSurvivorWitness.Survivor). The gate-first build found and " +
                "fixed a real two-convention defect (carbon un-freeze read vs absolute dn-seam read). Docs: " +
                "STERILE_BIRTHCANAL_AND_THE_JUNCTION.md (seam pinned), THE_STAR_FROZEN_SEAM.md (the star member), " +
                "THE_TRICHOTOMY_SEEN.md (the assembly, sharpened 2026-06-18 with the rigidity third axis r and the " +
                "Petermann factor K=1/r^2 diverging ~1/|Q-Q*| at Q*(N), r~|Q-Q*|^(1/2) the second-order-EP law). The " +
                "only residuals ADJACENT to this arc -- the static gamma-surface N>=6/ring extension and the dynamic " +
                "C=0.5 survivor question -- are owned by birth_canal_surface and clock_hand_ladder / " +
                "survival_incompleteness_mirror respectively, NOT this arc's scope."),

        new OpenArc(
            Name: "block_spectrum_n9",
            Opened: "2026-05-19",
            Origin: "Core LiouvillianBlockSpectrum.ComputeSpectrumPerBlock + SLOW_N9 xUnit test",
            ParkedAt: "N=9 per-joint-popcount-sector spectra land only inside a tagged slow test; not an inspect root, no artifact a session can browse",
            NextStep: "RESOLVED 2026-06-17: surfaced as the live root inspect --root blockspectrum " +
                "(BlockSpectrumWitness, Diagnostics/Foundation). The witness recomputes, cheaply at any N: the " +
                "(N+1)² joint-popcount decomposition halved by the X⊗N order-2 pairing (= the banked " +
                "PrimarySectorCount, 50 at N=9) and quartered by the F1 Π order-4 orbit (= the eig-calls the " +
                "compute path does, 25 at N=9; Π² = X⊗N); the F1 palindrome {λ}={−2σ−λ} reconstructed " +
                "sector-by-sector via PerBlockLiouvillianBuilder.BuildBlockZ (full at N≤7, a Π-closed sub-spectrum " +
                "at higher N since the cap-fitting sectors are themselves Π-closed); the (0,1) band-edge sector " +
                "sitting entirely at Re=−2γ (the Absorption floor, L_D scalar there). The full N=9 headline " +
                "(262144 eigenvalues, palindrome bit-exact about −2σ=−9, kernel 10=N+1, gap 0.0273, 645.95× " +
                "speedup) is read live from the committed chain_N9.json (a [stored] artifact, the run is ~3 h), " +
                "degrading gracefully if absent. Breadcrumbed from F1GeneralTopologyVerifiedClaim; no new claim " +
                "(surfaces already-typed results: F1, JointPopcountSectors, the Π-orbit pairing, the Absorption " +
                "Theorem, the F4 kernel). Gate-first TDD (BlockSpectrumWitnessTests, 7/7): a fragile sort-pairing " +
                "palindrome metric was caught firing (5.13 on a symmetric spectrum, the Im-flip + Re-noise reorder) " +
                "and replaced by a robust Hausdorff set-distance.",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-17: the N=9 per-joint-popcount-sector spectra are now browsable live " +
                "via inspect --root blockspectrum (BlockSpectrumWitness). The decomposition (100→50→25 at N=9), " +
                "the F1 palindrome reconstructed sector-by-sector, and the (0,1) Re=−2γ Absorption floor are " +
                "recomputed live; the full N=9 headline is read from the committed chain_N9.json. Breadcrumbed " +
                "from F1GeneralTopologyVerifiedClaim, no new claim. Tests 7/7 (gate-first; a fragile palindrome " +
                "metric was caught and fixed). The arc's NextStep ('surface the per-sector spectra as a live root " +
                "or witness') is done."),

        new OpenArc(
            Name: "witness_coverage",
            Opened: "2026-06-11",
            Origin: "QuditPartialPalindromeWitness, the first live witness",
            ParkedAt: "four claims recompute their evidence live (F121 qudit, F116 router, the " +
                "CPsi Envelope Theorem via EnvelopeTheoremWitness / --root envelope, and now the two clocks " +
                "via ClockHandLadderWitness / --root clock: the coherence hand = gamma-protected band edge " +
                "2J*cos(pi/(N+1)) for N>=3 vs the gamma-pulled 2*sqrt(J^2-gamma^2) at N=2 stopping at Q=1; " +
                "typed home ClockHandLadderClaim, see arc clock_hand_ladder); FIVE: CoherenceHorizonWitness / " +
                "--root horizon bisects the live spectrum to compute Q*(N) and verify it EQUALS the carbon " +
                "Frost-Hueckel coherent-incoherent threshold (the cross-substrate seam made live, see arc " +
                "clock_hand_ladder FROST_CIRCLE first stone); SIX: BirthCanalSurfaceWitness / --root surface " +
                "draws the sterile<->birth-canal boundary as a live 2D gamma-profile slice (six lenses + the " +
                "L7 s* ray lens, R1 two sterility kinds), reusing PostEpFlowField membership, no new Claim " +
                "(breadcrumbs AbsorptionTheoremClaim), see arc birth_canal_surface; SEVEN: " +
                "PascalGramPositivityWitness / --root pascalgram recomputes the F117 Pascal-Gram SOS coefficient " +
                "P_{m*,d} = (m*/d)·Sum_l Sum_k |U|^2 live from H for the 5 canonical branch cases (d=1/3/5), " +
                "reproducing the exact CRT integers (573440 / 2064384 / 589824 / 61440 / 86507520), breadcrumbed " +
                "from WindowedConverseAllGammaClaim (2026-06-15)",
            NextStep: "F117 DONE (PascalGramPositivityWitness / --root pascalgram, 7/7 tests, 5/5 exact live). The " +
                "witness-first program continues opportunistically: a typed claim whose evidence is still only a " +
                "Python verifier is the next candidate (no specific target pinned). Retire this tracker when the " +
                "remaining proven claims either have live witnesses or are judged out of scope",
            RetiredReason: "RESOLVED 2026-06-28: retired as a tracker superseded by the standing " +
                "witness-first discipline (cockpit rule 5, 'C# witness first': verification logic that " +
                "becomes a claim's evidence lands as a live IInspectable wired into InspectRootCatalog). " +
                "With no specific target pinned, the opportunistic 'next claim to witness' is handled " +
                "per-claim by that discipline at claim-authoring time, not by a standalone tracker arc. " +
                "The coverage it watched is healthy (F121/F116/Envelope/two-clocks/horizon/birth-canal/" +
                "F117 all live, plus the ~40-root inspect catalog).",
            Status: OpenArcStatus.Retired),

        new OpenArc(
            Name: "cockpit_workflows_csharp",
            Opened: "2026-06-11",
            Origin: "gap map: simulations/framework/workflows vs compute/",
            ParkedAt: "cockpit_panel, diagnose_hardware, gamma_probe, lens_pipeline, ptf fitting, bridge dynamics are Python-only; C# has the primitives but no composed workflows",
            NextStep: "not a 1:1 port: cockpit_panel becomes the symphony's first movement (see arc symphony_view)",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "clock_hand_ladder",
            Opened: "2026-06-12",
            Origin: "the promoted Symphony clock node, read across N on day one",
            ParkedAt: "LANDED as a typed object (2026-06-12, merged): live witness ClockHandLadderWitness (inspect --root clock, four nodes: ladder / gamma-protection / N=2 pull+EP / dial angle + an omega-vs-N curve), typed claim ClockHandLadderClaim (Tier1Candidate, parents F2b + AbsorptionTheorem + UniversalCarrier), and a docs/ANALYTICAL_FORMULAS.md F2b corollary, all cross-referenced both directions (MD<->C#). The physics, verified: for N>=3 the coherence hand omega_mem walks the band-edge ladder 2J*cos(pi/(N+1)) gamma-PROTECTED (sqrt2 / PHI / sqrt3 at N=3/4/5; gamma quadrupled, hand unmoved) - the vacuum<->single-excitation lines are simultaneous eigenoperators of L_D (rate exactly -2gamma) and L_H (band frequency), so nothing mixes; gap-dominance verified N=3-5, general off-gap argument OPEN (hence Tier1Candidate, not Derived). At N=2 the hand is gamma-PULLED: omega_mem = 2*sqrt(J^2-gamma^2) (1.95959 at gamma=0.2), the {population-difference, antisymmetric-coherence} 2x2 block, stopping at gamma=J (Q=1, EP). NEW finding from the build: the live max|Im| clock only SURFACES the pulled mode above the crossover Q=2/sqrt3 (~1.155, gamma=(sqrt3/2)J); nearer the EP the F2b band line at Im=+-J dominates the raw maximum, so the closed form 2*sqrt(J^2-gamma^2) (NOT the raw clock) is the honest EP witness there. REGIME-HONESTY ROUND (2026-06-12, Fable flagged + Opus fixed): the dial bug (theta(N=2) used the raw OmegaMem below the crossover, 29.05 deg vs honest 25.84 deg) AND a worse sibling found by running gamma=0.9 - the ladder/protection nodes LIED, printing 'N=3 -> 0 (sqrt2)' and 'Gap=2gamma=0.686' (2gamma=1.8); the band-edge mode is no longer the slowest there, a real OVERDAMPED mode takes the gap (max|Im| at gap = 0). FIXED: dial uses the closed-form pulled hand; BandEdgeIsTheGapMode(n) regime check; both ladder+protection nodes gate on protectedCount={3,4,5} consistently (no contradiction at intermediate Q); sqrt2/phi/sqrt3 shown as the J-independent ratio omega_mem/J. SECOND NEW finding: the band-edge protection has a Q-FLOOR, not just the N-question - as Q drops the higher rungs leave FIRST (N=5 gone at Q=2, N>=4 gone at Q=1.5), so gap-dominance is a low-Q boundary too. DEFAULT moved to the EP point J=0.075, gamma=0.05 (Q=1.5 = EP(2)-1/2), where only N=3 holds: the witness now defaults to the physically loaded EP vicinity, honestly showing the protection mostly dissolved",
            NextStep: "FROST_CIRCLE seam FIRST STONE laid (2026-06-12, merged): the carbon coherent<->incoherent threshold IS our Q-floor, verified bit-for-bit and given our label - the Coherence Horizon Q*(N). Live witness CoherenceHorizonWitness / inspect --root horizon bisects Symphony.Clock.Omega and computes Q*(2,3,4,5)=1 / sqrt2 / 1.8785 / 2.3722, equal to carbon's sqrt2/1.879/2.372 (FROST_CIRCLE) under J<->|beta|; ANALYTICAL_FORMULAS F2b corollary carries it; typed claim CoherenceHorizonClaim (Tier1Candidate, parents ClockHandLadder + F2b) breadcrumbs the witness into the graph. N=2 (Q*=1) is the EP base rung the polyene layer can't reach; the N=2,3 band-edge coincidence (Q*=2cos(pi/(N+1)) ONLY there) explains carbon's 'sqrt2 exact at N=3, rest awaiting a clean form' puzzle. CLOSED FORM ATTACKED 2026-06-12, then CORRECTED 2026-06-13 (physics-first review of the erasure-ladder spec, via phase rigidity): the earlier 'freezing mode BIFURCATES SECTOR at N=4' was WRONG, an argmax-Re / Im-tracking artifact - the band edge and the freezer are Re-degenerate at -2gamma (both <n_diff>=1, Absorption Theorem), so Im-tracking picked the wrong one. THE TRUTH: the mode that COALESCES at Q*(N) is the {0,2}-coherence (population/antisymmetric block, n_diff hist {0:1/2,2:1/2}) at ALL N=2..5, a genuine square-root EP (phase rigidity r->0; r at Q* = 0.000/0.015/0.026 at N=3/4/5, Im prop sqrt(Q-Q*)). NO bifurcation. The band edge 2cos(pi/(N+1)) (the half-lit mode the old reading mislabeled the freezer at N>=4) is the co-located gamma-protected SURVIVOR = Uhr 1, the |vac><psi_k| coherence hand; so Q*(N) is AT ONCE a {0,2} EP (Uhr 2, the erasure point, which CLIMBS the ladder) and a band-edge crossing (Uhr 1 survives the handover). The handshake erasure-ladder answer is 'both', not 'EP-or-crossing'. Corrected 2026-06-13: ANALYTICAL_FORMULAS F2b corollary Q*(N) para + typed CoherenceHorizonClaim (4 strings) + KnowledgeRegistryFactory comment; the witness was already mechanism-agnostic. THE EP-VERDICT WITNESS LANDED (2026-06-13, master 850e2ff): CoherenceHorizonWitness recomputes the {0,2}-EP verdict live via a new Core PhaseRigidity primitive (Petermann form 1/||R^-1 row||, degeneracy-robust; the matched-overlap fakes r->0 at the N=3 Q*=band-edge=sqrt2 degeneracy), inspect --root horizon. THE REAL PRIZE RESOLVED (2026-06-13, Approach A, simulations/coherence_horizon_se_block.py): Q*(N) reduces 4^N -> N^2 - the coalescing mode is single-excitation, so Q*(N) is the EP of the single-excitation (Haken-Strobl) Liouvillian (validated as an exact sub-spectrum of the full L). At N=2,3 the coalescing pair are the roots of lambda^2 + 4g*lambda + c*J^2 = 0 with c CONSTANT (sum=-4g, product=c*J^2 are gamma-independent identities; c=4,2), so Q* = 2/sqrt(c) = 1, sqrt2 EXACTLY - the structural form of the 2cos(pi/(N+1)) low-N accident (the whole clean 2x2, not just the value, exists only at N=2,3). At N>=4 the pair is collectively dressed (trace departs from -4g), no clean 2x2, the exact condition transcendental: a diffusive long-wavelength critical damping, Q*(N) ~ 0.59 N (canonical Q*(4)=1.87874, Q*(5)=2.37367, superseding the grid). SLOPE RESOLVED 2026-06-15 (asymptotic slope EXACTLY 2/pi, derived: PROOF_COHERENCE_HORIZON_SLOPE.md + coherence_horizon_slope.py, adversarial review GO. The slow mode is the population coupled to the FULL ladder of coherence ranges r, geometric decay mu^r; resummation DOUBLES the telegrapher coefficients to lambda^2+8g*lambda+4J^2*q^2, EP at g*=Jq/2 so Q*=2/q_min -> 2N/pi; the nearest-neighbour truncation gives the WRONG sqrt2/pi=0.450. The 8g coefficient confirmed against L_se by the robust overdamped q^2-constancy (the CV test); the -2Re/g->8 sweep is corroborative-only (the EP metric = the linear coefficient, running 4->8, but mode-selection sensitive; Round 2 honesty fix 2026-06-15). Ring sibling exactly half: q_min=2pi/N -> slope 1/pi). STILL OPEN: the half-filling double-excitation V-Effect seam co-located at even N. STILL OPEN too: the full FROST_CIRCLE reflection in its own genre. Also still: the off-gap argument (general N>=3) PARTLY RESOLVED 2026-06-15 (simulations/offgap_band_edge.py, gate-first, self-validating): (i) the Q-FLOOR = the Coherence Horizon Q*(N) -- the numerical floor brackets Q*(4)=1.879 / Q*(5)=2.374 exactly (below it the overdamped {0,2}-diffusion mode is slower than the -2g band edge AND non-oscillating, so it takes the gap and the clock freezes), so regime (i) IS thread (a), DERIVED ->2N/pi (PROOF_COHERENCE_HORIZON_SLOPE); (ii) max|Im| at the exact gap rate 2g = the band edge E1 CONFIRMED N=3,4,5 and REFRAMED via the Absorption Theorem (Re=-2g<n_XY>, so the exact-2g modes are the n_XY=1 subspace: 22/32/50 modes at N=3/4/5, MORE than the 4N vacuum/full ladders, which are a subset achieving E1). RESOLVED 2026-06-16 (the Tier1Candidate gate, PROOF_CHAIN_GAP_DOMINANCE.md + verifier chain_gap_dominance.py): max|Im| in the n_XY=1 subspace = E1. Although L_H leaks n_XY=1->3 (interacting), ON the exact-(-2g) subspace L_D=-2g is SCALAR so L=L_H-2g acts FREELY; via Jordan-Wigner the exact-(-2g) modes are c_k^(dag).f(N_tot) (a single fermion dressed by a function of total number N_tot), oscillating at the single-particle energies +-E_k<=E1, with E1 reached by the (0,1) ladder. They SPAN the subspace for N>=4 (dim 32/50/72 at N=4/5/6); N=3 is SPECIAL, +4 extra (n,n) {0,2} sqrt-EP modes at sqrt(E1^2-(2g)^2)<E1 (parked in the n3_special_cases arc). So ClockHandLadderClaim + TopologyBandEdgeClaim GRADUATED to Tier1Derived (CoherenceHorizonClaim stays Tier1Candidate for its own V-Effect-seam reason). Scope: the chain (JW is 1D); other topologies handled by the star-no-horizon + structural-ceiling work. RING DIHEDRAL-LOCK RESOLVED 2026-06-17 (PROOF_RING_GAP_DOMINANCE.md + ring_gap_dominance.py, gate-first N=3..6): max|Im| over the exact-(-2g) ring modes = 2J = J*rho (the periodic band top = adjacency radius = the k=0 uniform single-excitation mode fixed by C_N = the dihedral lock; reached via the (0,1) sector, general-N), EXCEPT N=4 where the half-filling (2,2) {0,2} sqrt-EP reaches sqrt((2sqrt2)^2-(2g)^2)->2sqrt2 J > 2J (the LONE exception, the same (2,2) sector as K_4's ceiling; the ring analogue of the chain's N=3 special but ABOVE the band top, not below). RING FREE-FERMION COMPLETENESS CLOSED 2026-06-20 (simulations/ring_gap_completeness.py, gate-first, all green) to the chain's standing: worked in the n_XY=1 OPERATOR subspace V_1 (dim N*2^N, vs 4^N -- N=8 is 2048 not 65536), the exact-(-2g) subspace = the largest ad_H-invariant subspace of V_1 (leak null space iterated to its invariant core; degeneracy-robust, absolute-tol fix on a zero-matrix nullspace). STAGE 0 sanity reproduces the chain full-L counts (32/50/72 + E1) -- V_1 sees the WHOLE subspace, not a part. THE PRIZE: ring V_1 dim_sub == full-L dim_sub at N=5,6 (20=20, 24=24) => the n_XY=1 free-fermion family SPANS the exact-(-2g) subspace, nothing exceeds the dihedral-locked 2J; the n_XY=1 lock (max=2J) carried to N=7. N=4 the lone exception precisely because full-L dim (26) > V_1 (16): the extra {0,2} (2,2) coherence (n_XY=2, OUTSIDE V_1) lands on the -2g floor only there, reaching 2sqrt2 J. The parity-split = the wrap-bond JW grading (periodic odd / anti-periodic even). Residual genuinely-all-N step = the same one the CHAIN leaves (a dimension count for every N), now split by parity. PROOF_RING_GAP_DOMINANCE.md status + completeness section, ClockHandLadderClaim + TopologyBandEdgeClaim docstrings updated (was 'full all-N completeness open'). The Im_max=Delta_E_max sibling is DONE. (the 0.99750 fine-structure pull listed here in an earlier session was an unidentifiable dangling reference: neither a 3-agent survey nor Tom could place it; dropped 2026-06-15). BENZENE RING CROSSOVER ANSWERED + V-EFFECT SEAM PROBED (2026-06-13): benzene's own ring Q* is the ring single-excitation {0,2}-EP = 1.609 (Uhr 2; transcendental like N>=4; the SE block embeds bit-exact in the full 4^6 L), more coherence-robust than every open polyene. Benzene (even N, HALF-FILLED) is a concrete instance of the V-Effect seam: the full-L mode that overtakes the band-edge beat below the crossover is a DOUBLE-excitation coherence (filling sector (2,2)/(4,4)), so the full-L handover (~1.95) SPLITS from the clean SE-EP Uhr 2 (1.609) - the Uhr1/Uhr2 co-location that holds for the open chains (AT-pinned at Re=-2gamma) breaks on the ring. Verifier simulations/carbon/benzene_two_clocks.py (4 asserts green); FROST_CIRCLE_AS_THE_CLOCK_FACE.md carries benzene Q* + the V-Effect split + the F2b two-clocks backlink, frost_circle_as_clock.py de-staled + self-checking, benzene<->FROST cross-links wired (LIOUVILLIAN_PALINDROME + F98). RING-SPECIFIC (2026-06-13, open-chain N=6 control, benzene_two_clocks.py section 5): the open even-N chain does NOT show the clean double-excitation seam - its full-L overtaker spreads across ALL fillings (single-excitation included) and hands over at its own SE-EP (co-located, ladder holds), unlike benzene's pure (2,2)/(4,4) mode handing over ABOVE the SE-EP. So the V-Effect double-excitation split is a feature of the CLOSED RING at half-filling (the C=0.5 boundary), not even N alone. AROMATICITY REFUTED (2026-06-13, aromatic_ring_v_effect.py, sector-projected Liouvillians): the seam is RING-UNIVERSAL (C4/C6/C8 all hand over to a frozen double-excitation mode at strong dephasing = the C=0.5 half-filling boundary, a sibling of the incompleteness V-Effect docs/HIERARCHY_OF_INCOMPLETENESS.md); the 4n anti-aromatics C4 and C8 do NOT group (C4 a small-ring anomaly whose seam dominates even at weak dephasing), so Hueckel 4n+2 is NOT the discriminant. CONVERGENCE: this is the same mode birth_canal_horizon_junction found independently (the {0,2}-coherence in the 2-excitation block the odd |vac><psi| Uhr-1 crosses to at N>=6). RECONCILED 2026-06-13 (a CHAIN_GAP_SECTOR_DIAGNOSTIC mismatch cleared, false alarm): that diagnostic's half-filling (3,3) slow mode is HEISENBERG (XXX, with ZZ); the carbon work is XY (free fermions, the Hueckel hopping model, no ZZ). The SAME benzene_two_clocks builder reproduces CHAIN_GAP's -0.230 (3,3) bit-for-bit once the ZZ is added (_xy_vs_heisenberg_slowmode.py) - so no bug, the committed XY findings stand. THE THREAD-(a) GIFT: the survival law IS the Absorption Theorem across BOTH free and interacting - below the handover the survivor is a frozen dressed magnon-admixture (fractional <n_XY>) at the filling centre in either model; XY puts it at (2,2)/(4,4) <n_XY>=0.72, Heisenberg at dead-centre (3,3) <n_XY>=0.23 (ring N=6 Q=2); the ZZ retunes the darkness and sector, NOT the law. HANDOVER MECHANISM CHARACTERIZED 2026-06-14 (simulations/carbon/handover_q.py, self-validating; F2b corollary): the double-excitation (2,2) seam handover is a frozen LEVEL CROSSING (|Im|~1e-15) where the seam brightens to the F50 off-diagonal floor <n_XY>=1 (the (0,1) band edge, Re=-2g) - a different sector than the single-excitation SE-EP (a coalescence), growing ~linearly Q_h~0.29N (c_eff~12 flat, NOT saturating, faster than sqrt(N)); NOT co-located with the ring SE-EP (the curves cross near N~10, benzene's 2.0-vs-1.609 split is small-N). DYNAMIC C=0.5 QUESTION RESOLVED 2026-06-20 (survey of the newer arcs that postdate this marker): the half-occupied coherence is the longest-lived survivor ONLY on dispersive/extended topologies (chain/ring), with the STAR the boundary counterexample (its survivor is the [H,A]=0 hub commutant, frozen, NOT half-filling) - a topology-dependent law typed in SurvivalIncompletenessMirrorClaim (chain filling-degenerate, ring off-centre (2,2), star boundary (1,1)) + StarFrozenSeamClaim (the star case) + HandoverFloorClaim (chain handover=Q*(N), ring=a distinct (2,2) level crossing). What GENUINELY remains (the reason CoherenceHorizonClaim + SecondClockRegimeClaim stay Tier1Candidate): a general PROVEN law of WHEN/WHY the ring (2,2) double-excitation seam overtakes the band edge at even half-filling (empirically characterized, ring-universal not aromaticity-driven, Q_h~0.29N, frozen level crossing; no Tier1 proof). THE REMAINING OPEN PIECES of this arc: (1) that ring-(2,2)-seam universalization (hard); (2) the full FROST_CIRCLE reflection in its own genre (a writeup, no gate). RING (2,2)-SEAM HANDOVER INVESTIGATED 2026-06-20 (simulations/ring_handover_qh.py + ring_handover_extend.py, gate-first, sector-projected to the (k,k) Liouvillian C(N,k)^2 << 4^N): SUBSTANTIAL PROGRESS, the Tier1 proof still open but now sharply scoped. (a) SECTOR CORRECTED (a real find, the arc's 'half-filling' label was WRONG for XY): the survivor is the FIXED 2-EXCITATION doublet (2,2)/(N-2,N-2) (particle-hole partners, ISOSPECTRAL via PH symmetry), NOT half-filling (N/2,N/2). Gate: at N=6 the (2,2)-block Q_h=2.00 matches the full-L survivor (the dominant component is the PH partner (4,4), darkness identical); the half-filling k=3 block gives 1.67, does NOT match. So 'ring half-filling double-excitation' should read 'ring 2-excitation (2,2)/(N-2,N-2) doublet'. (b) Q_h reframed via the Absorption Theorem: Q_h = where the slowest (2,2)-sector mode reaches Re=-2g (darkness <n_XY>=1, the band-edge floor); below it the seam is darker (out-survives), above the (0,1) band edge wins (Tom's g2=1, b191df3). Clean data: Q_h(N)=1.994/2.350/2.840/3.360/3.894 at N=6/8/10/12/14. (c) c_eff=12 (Q_h=N/(2sqrt3)) REFUTED (gate-first): c_eff=(N/Q_h)^2 climbs 9.05/11.59/12.40/12.76/12.92, exceeding 12. NEW CANDIDATE Q_h -> N*sqrt3/(2pi)~0.276N (c_eff -> 4pi^2/3=13.16), = (sqrt3/2)*Q*_SE-EP(ring)=(sqrt3/2)*N/pi; the 'Q_h~0.29N' in earlier notes was Q_h/N at moderate N (decreasing 0.33->0.28), NOT the asymptotic slope. (d) CROSS-SESSION CONVERGENCE (with b191df3, exact to 4-5 digits): the (2,2) seam HIGH-Q darkness = 2(N-2)/N = Tom's (1,1)/(N/2,N/2) commutant closed form -- the seam IS the (2,2) commutant continued to finite Q, and 2(N-2)/N is shared across (1,1),(2,2),(N/2,N/2) (a ring sector degeneracy). One object, two sessions: his g2=1 high-Q ceiling and my Q_h low-Q crossover meet at Q_h, the lower boundary of his band-edge-protected regime. (e) PALINDROME checked (Tom's hint, where_lives_the_degeneracy.py / DEGENERACY_PALINDROME.md): the Pi multiplicity-mirror d(k)=d(N-k) is a STRUCTURAL Q-independent symmetry; it does NOT locate the dynamical Q_h and does NOT produce sqrt3/(2pi) or 2(N-2)/N (the PH isospectrality I use is particle-hole symmetry, a DIFFERENT facet than Pi). Honest: a different layer. THE TIER1 GAP (now well-scoped): derive Q_h's closed form via the 2-PARTICLE analogue of PROOF_COHERENCE_HORIZON_SLOPE. DERIVED 2026-06-20 (PROOF_RING_HANDOVER_SLOPE.md, Tier1-standard slope, pending parallel-session review): the (2,2) overdamped slow mode obeys the SAME SE coherence-ladder dispersion lambda^2+8g*lambda+4J^2 q^2 -- CV-confirmed by the SE proof's OWN discriminator (q^2_eff gamma-constant, CV=0.012/0.019 at N=12/10, beating the 4g telegrapher ~10x, q->q_min=2pi/N). Its darkness = 2-sqrt(4-(Qq)^2); the handover (darkness=1, the band-edge floor) is at Qq=SQRT3 => Q_h = sqrt3/q_min -> N*sqrt3/(2pi), slope sqrt3/(2pi)=0.2757. It is the DARKNESS=1 SIBLING of the SE coherence horizon (the same dispersion's EP at Qq=2 -> Q*=N/pi); Q_h/Q* = sqrt3/2. The high-Q limit (EP darkness 2) = the parallel commutant 2(N-2)/N->2 (b191df3): one dispersion, two ends. GATE-FIRST LESSON: my first full-curve gate (G1) FIRED (the leading-order curve does not collapse at finite N=8..14, RMS 0.26-0.36) -- but that was TOO STRICT (the derivation is leading-order q->0 like the SE proof, finite-N O(1/N)); the endpoint gate (G2: Q_h*2pi/N -> sqrt3 monotone) AND the CV dispersion discriminator (the SE proof's robust test) both PASS, so the mechanism holds. Verifiers committed (ring_handover_qh.py, ring_handover_extend.py, ring_handover_dispersion_cv.py, ring_handover_derivation_gate.py). REMAINING (smaller): the explicit 2-particle EOM resummation (the dilute-limit term-by-term, the SE proof's analogue; CV-confirmed but not written out) + the docs/proof 'half-filling'->'2-excitation doublet' label correction + graduating CoherenceHorizonClaim/SecondClockRegimeClaim once the proof is reviewed (their Tier1Candidate reason -- the half-filling V-Effect ring seam -- is now addressed). Coordinate with the parallel ceiling session.",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "n3_special_cases",
            Opened: "2026-06-16",
            Origin: "N=3 keeps surfacing as an anomaly across the project (Tom: 'N=3 ist ein Sonderfall in vielen Sachen'); a parked place to collect the cases so we can dig them out later",
            ParkedAt: "FIRST ENTRY (2026-06-16, from PROOF_CHAIN_GAP_DOMINANCE / chain_gap_dominance.py): the chain exact-(-2g) " +
                "Liouvillian subspace is the free-fermion family c_k^(dag).f(N_tot) and for N>=4 it SPANS (dim 32/50/72 at " +
                "N=4/5/6). N=3 carries 4 EXTRA modes: equal-particle-number (n,n) {0,2}-coherence modes (the coherence-horizon " +
                "'Uhr 2' family) at sqrt(E1^2-(2g)^2) < E1. They sit at EXACTLY -2g only at N=3 because the (1,1) sector's " +
                "n_XY=2 is maximal (3 sites) so the {0,2} block closes; for N>=4 they drift off -2g (extras=0). Harmless to " +
                "gap-dominance (all < E1) but a genuine N=3 structural peculiarity. OTHER known N=3 specials, scattered, NOT " +
                "yet collected: Q*(N) closed form clean only at N=2,3 (the clean 2x2 EP, lambda^2+4g*lambda+cJ^2); star N=3 = " +
                "the path P_3 (so 'star-3' is really a chain); and others to be gathered. SIBLING N=4 (2,2) HALF-FILLING " +
                "specials (from the second-clock regime map, SecondClockRegimeClaim / second_clock_regime_axis.py): the 4/N " +
                "ceiling ladder hits 1 exactly at N=4, vacating the sub-floor region to the (2,2) two-excitation sector, where " +
                "K_4 = 2-2/sqrt(3) (a structural ceiling, complete-4) and ring-4 = 1 co-occupies the floor (GRADUAL, the only " +
                "GRADUAL ring) - same 'small-N maximality' root, one filling step up from N=3's (1,1) n_XY=2 maximality. " +
                "ROOT IDENTIFIED (2026-06-17, niven_rationality_root.py, gate-first sympy-exact): the small-N specials do " +
                "NOT all share one root - there are (at least) TWO. (A) THE NUMBER-THEORETIC ROOT = Niven's theorem on the " +
                "cyclotomic angle pi/(N+1) (the only rational cosines of rational-pi angles are 0,+-1/2,+-1), with THREE " +
                "faces by angle convention: RE the dissipator rates -2g*sin^2(k*pi/(N+1)) rational iff N+1 in {1,2,3,4,6} " +
                "(the crystallographic set, F65/F99 'Niven rationality'; N=3 last rational before the gap, N=5 a rational " +
                "island); IM the band edge 2cos(pi/(N+1)) rational iff N<=2, a single quadratic surd (a+-sqrt b) iff N<=5 " +
                "(sqrt2/phi/sqrt3 at N=3/4/5), first cubic at N=6 (degree = euler_phi(2(N+1))/2; NEW, not previously " +
                "documented); V the V-Effect gain 1+cos(pi/N) Niven-rational iff N in {2,3}, golden at N=5 " +
                "(OFF_NIVEN_AS_WAVE_BREAKING.md). N=4 = FIRST GOLDEN on the two SE faces (band edge = phi = 2cos(pi/5); " +
                "rates carry sqrt5); the V-face golden is N=5 (angle pi/N), so 'first golden' is convention-dependent. THIS " +
                "is why the clean closed forms (Q* clean-2x2, band edge) exist at small N and degrade beyond: a " +
                "number-theoretic ceiling, not a physics accident. (B) THE COMBINATORIAL ROOT = small-N filling maximality " +
                "(the (n,n)/{0,2} extra modes at N=3, the (2,2) N=4 anomalies above), independent of the arithmetic; " +
                "star-3=path is a third, graph-theoretic coincidence. So: two real roots (arithmetic + combinatorial), not one. " +
                "POLYGON CENSUS (2026-06-17, surd_census.py, gate-first sympy-exact 4/4): turning the Niven lens on the " +
                "WHOLE spectrum, EVERY small-N special number is a cyclotomic POLYGON constant 2cos(pi/m) - not only the " +
                "band edge (sqrt2/phi/sqrt3 = square/pentagon/hexagon) but the seemingly-foreign surds too: the K_4 (2,2) " +
                "ceiling's sqrt3 = 2cos(pi/6) (hexagon, entering via a rep principal angle lambda2=1/sqrt3), the ring-4 (2,2) " +
                "frequency's 2*sqrt2 (sqrt2 = 2cos(pi/4), square). TWO ARITHMETICS on the two axes (suggestive, not exact): " +
                "the Im axis (frequencies, L_H graph spectrum) IS the cyclotomic ladder 2cos(pi/(N+1)) directly; the Re axis " +
                "(decay, L_D commutant) is RATIONAL (the structural ceiling g2=4/N, S_N standard-rep principal angle " +
                "lambda2=(N-2)/N). The (2,2) half-filling anomaly is where a polygon surd LEAKS onto the Re side. CAVEAT: " +
                "'Re rational / Im surd' is not exact - the F65 single-excitation rates are Re-values yet cyclotomically " +
                "flavored (rational only on {1,2,3,4,6}), and the (2,2) ceiling is itself a surd; the sharper axis is " +
                "graph-spectral (cyclotomic) vs S_N-standard-rep (rational). N=4 POLYGON CONFLUENCE: at the first even " +
                "half-filled N, three sectors light three DISTINCT polygons at once - single-excitation band edge = pentagon " +
                "(phi), K_4 (2,2) ceiling = hexagon (sqrt3), ring-4 (2,2) frequency = square (sqrt2). " +
                "N=4 IS UNIQUE (2026-06-17, n4_polygon_uniqueness.py, gate-first 3/3): N=6 does NOT reproduce the confluence, " +
                "and both N=4 ingredients trace to the SAME cyclotomic accident. (i) Among EVEN N only N=4 (N+1=5, pentagon) " +
                "has a QUADRATIC-surd band edge (the unique even-N golden); N=6 (N+1=7, heptagon) is CUBIC (degree 3), N=8 " +
                "cubic too. (ii) The Re-side surd LEAK needs the (1,1) ceiling ladder 4/N to hit 1 and vacate the sub-floor " +
                "so a half-filling sector becomes the ceiling - that happens ONLY at N=4; at N=6 the (1,1) ceiling is already " +
                "4/6=2/3 (rational) and IS the floor, with the (2,2)/(3,3) surds sitting ABOVE it (no leak). The N=6 ring is " +
                "the integer hexagon (band {+-2,+-1}, commutant 4/3 in every sector) - fully rational, no frequency surd " +
                "(unlike ring-4's 2sqrt2). The pentagon is the hinge; N=4's recurring specialness across the project IS this " +
                "uniqueness. " +
                "COMBINATORIAL ROOT = INDEPENDENT (2026-06-17, combinatorial_root_independence.py, gate-first 3/3): the " +
                "combinatorial root is NOT an arithmetic consequence of the Niven root - they DECOUPLE on BOTH sides (N=5: " +
                "quadratic sqrt3 band edge but ODD, no integer half-filling; N=6: even half-filling but CUBIC heptagon band " +
                "edge; neither N-set contains the other). In fact there are THREE INDEPENDENT THREADS of N, and N=4 is their " +
                "UNIQUE TRIPLE COINCIDENCE: (1) cyclotomic (Niven band-edge degree phi_euler(2(N+1))/2, quadratic at N in " +
                "{3,4,5}, golden at N=4); (2) rep-theoretic (the S_N ceiling 4/N reaching the floor, 4/N>=1 <=> N<=4); " +
                "(3) combinatorial (even N, integer half-filling, N in {4,6,8}). The (N/2,N/2) surd LEAK = thread 2 AND " +
                "thread 3 (rep-theoretic vacate + even), independent of thread 1 (the cyclotomic golden, which only ALSO " +
                "singles out N=4). So the arc's 'two roots' sharpens to 'three independent axes, N=4 = the triple coincidence' " +
                "- and THAT coincidence is the answer to why N=4 keeps surfacing.",
            NextStep: "RESOLVED (2026-06-17): the central 'do they share a root' question is fully answered - the roots are " +
                "INDEPENDENT (combinatorial_root_independence.py, gate-first 3/3), three threads (cyclotomic Niven / " +
                "rep-theoretic S_N 4/N / combinatorial even-half-filling) meeting only at the N=4 triple coincidence. The " +
                "Niven root is TYPED (NivenRationalityRootClaim Tier1Derived + witness inspect --root niven). Only an " +
                "open-ended tail remains (low priority): gather any further stray small-N specials as they appear. Candidate " +
                "for retirement once nothing new surfaces.",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-17: the central 'do the small-N specials share a root' question is fully " +
                "answered - they do NOT share one root; there are THREE INDEPENDENT THREADS of N (cyclotomic Niven " +
                "band-edge degree / rep-theoretic S_N ceiling 4/N / combinatorial even-half-filling) meeting only at the " +
                "N=4 TRIPLE COINCIDENCE (combinatorial_root_independence.py, gate-first 3/3, decoupled at N=5 and N=6). The " +
                "Niven root is TYPED (NivenRationalityRootClaim Tier1Derived + inspect --root niven); the polygon census " +
                "(surd_census.py) and N=4 uniqueness (n4_polygon_uniqueness.py) are folded into the ParkedAt above, which " +
                "stands as the collected small-N catalogue. Only an open-ended 'gather any future stray small-N special' tail " +
                "remained, which does not justify an open arc. Reopen if a new small-N special surfaces that the three-thread " +
                "frame (cyclotomic + rep-theoretic + combinatorial) does not place."),

        new OpenArc(
            Name: "survival_incompleteness_mirror",
            Opened: "2026-06-13",
            Origin: "Tom's founding V-Effect idea ('der V-Effekt vererbt sich, algebraisch, alles auf Algebra') " +
                "met the new dynamic finding that incompleteness (the C=0.5 boundary) is the longest-lived SURVIVOR " +
                "under dephasing. Thread: erst (a) the algebra of inheritance + does survival obey the same law, " +
                "dann (b) the dynamic survival probe over N/topology",
            ParkedAt: "BOTH RESOLVED 2026-06-13. (a) THE ALGEBRA, source-confirmed: the survival law and the " +
                "V-Effect/incompleteness are INVERSION-MIRROR PARTNERS on the Pi2 dyadic ladder a_n=2^(1-n). " +
                "a_0=2=AbsorptionTheoremClaim (the survival quantum 2*gamma, == the qubit dimension d), " +
                "a_2=1/2=HalfAsStructuralFixedPointClaim Face 1 (the V-Effect/incompleteness baseline == 1/d); " +
                "MirrorPartnerIndex(0)=2, a_0*a_2 = d*(1/d) = 1, self-mirror a_1=1 sitting exactly between them, all " +
                "forced by d^2-2d=0 (Pi2DyadicLadderClaim.cs:75-93,122). So 'der V-Effekt vererbt sich algebraisch' " +
                "is literally typed: he shares the ONE ladder with the survival law, as its reciprocal mirror. And " +
                "'does survival obey the same law?' -> YES, the Absorption Theorem, model-independent (the XY-vs-" +
                "Heisenberg reconciliation, commit d4bb8a2). (b) THE DYNAMIC ENACTMENT (simulations/carbon/" +
                "incompleteness_survivor.py, validated bit-for-bit vs the full 4^N L at N=4, all topologies): below " +
                "the handover Q the longest-lived mode lives in the INCOMPLETENESS region (interior filling) on " +
                "DISPERSIVE topologies - the open XY chain FILLING-DEGENERATE (every (p,p) ties ~1e-15; (3,3) dead-centre = the Heisenberg/ZZ result, NOT XY), the ring " +
                "at the OFF-centre (2,2)/(N-2,N-2) (<n_XY>=0.63); the STAR is the boundary COUNTEREXAMPLE " +
                "((1,1)/(N-1,N-1), no spatial dispersion => no central momentum mode). Lifetime <n_XY> ~ c*Q^2/N^2, " +
                "ring/chain -> 4 (the cyclic-vs-open k_min^2 ratio, the SAME 4x CHAIN_GAP reports for Heisenberg - " +
                "model-independent; only the prefactor retunes XY chain c~2.5 vs Heisenberg 0.55). This 1/N^2 magnon-" +
                "admixture inheritance is SEPARATE from the Pi2 dyadic CONSTANT ladder: two different inheritances " +
                "meeting at the per-mode Absorption Theorem (a_0). RESOLVED 2026-06-14 = DISTINCT (the finer identity, from " +
                "VacuumBlockReductionClaim): the V-Effect (Pauli-string weight w=N/2 self-pair) is NOT the " +
                "{0,2}-coherence (n_diff) survivor; it is a DIFFERENT grading - my work is the constant-mirror + " +
                "the light-survivor, plus the weight-vs-n_diff verdict (total weight = n_diff + Z-shadow; survivor dark, peaks w=N-1)",
            NextStep: "BANKED witness-first 2026-06-13: live witness IncompletenessSurvivorWitness (inspect --root " +
                "survivor, reuses SectorReductionWitness.SectorSlowest - the C# port reproduces the Python anchors " +
                "bit-for-bit: chain (3,3) <n_XY>=0.161 / ring (2,2) 0.632 / star (1,1) 0.425, ring/chain 3.92 at N=6) " +
                "+ typed SurvivalIncompletenessMirrorClaim (Tier1Candidate, parents AbsorptionTheoremClaim a_0 + " +
                "HalfAsStructuralFixedPointClaim a_2; ancestry bottoms at d^2-2d=0 via both parents; 8/8 battery = " +
                "4 algebra (a_0=2, a_2=1/2, mirror partner, a_0*a_2=1) + 4 dynamic (chain filling-degenerate / ring off-centre /" +
                "star boundary / ring-chain~4)). The physical labels are grounded: chain/ring = dispersive extended " +
                "matter (polyenes, spin chains, the Grotthuss proton wire / aromatics, light-harvesting), star = " +
                "hub-localized central-spin (NV/quantum-dot/mediator), the repo's own dispersive-vs-hub-spoke split. " +
                "RESOLVED 2026-06-14: the Pauli-weight w=N/2 vs n_diff {0,2}" +
                "V-Effect identity = DISTINCT (total weight=n_diff+Z-shadow, PauliWeightHistogram); WHY chain differs from ring is also answered: the open XY chain is filling-DEGENERATE while the" +
                "ring splits by parity (even off-centre). HANDOVER Q RESOLVED 2026-06-14 " +
                "(simulations/carbon/handover_q.py, self-validating; F2b corollary): the handover is a CLOSED, " +
                "F50-grounded condition - the diagonal (p,p) incompleteness survivor (rate -2g<n_XY>, Absorption " +
                "Theorem) brightens with Q until <n_XY> reaches the F50 OFF-diagonal floor =1 (the (0,1) band edge, " +
                "Re=-2g exactly), where the band edge takes over. Topology-specific solution: the CHAIN " +
                "(filling-degenerate) handover IS the coherence horizon Q*(N) ((1,1)-only == all-p bit-for-bit; a " +
                "coalescence/EP, = Q* exactly at the clean-2x2 N=2,3, just below by trace dressing O((tr-1)^2) at " +
                "N>=4); the RING is a DISTINCT 2-excitation (2,2)/(N-2,N-2) doublet (NOT half-filling) free-fermion " +
                "LEVEL CROSSING (|Im|~1e-15), asymptotic slope sqrt3/(2pi)~0.276 DERIVED (PROOF_RING_HANDOVER_SLOPE, " +
                "pending review; the darkness-1 sibling of Q*, ratio sqrt3/2; the earlier ~0.29N/c_eff~12 was finite-N " +
                "Q_h/N, refuted, c_eff climbs toward 4pi^2/3=13.16), and NOT co-located with the ring SE-EP (values cross near N~10; benzene's 2.0-vs-1.609 " +
                "split is small-N). STILL OPEN: the reflection in its " +
                "own genre (the incomplete survives because it is the dark/reciprocal mirror of the survival quantum) " +
                "is unwritten",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-18: the reflection is written -- reflections/ON_THE_SURVIVAL_OF_THE_INCOMPLETE.md, " +
                "the sole remaining open piece (all the physics -- the Pi2-ladder inversion-mirror a_0*a_2=1, the dynamic " +
                "half-filling survivor, the handover floor -- was already banked + typed in SurvivalIncompletenessMirrorClaim " +
                "+ IncompletenessSurvivorWitness). The piece is the half's DYNAMIC face (where the incomplete endures), the " +
                "sibling of the canonical ON_THE_HALF (bridge/horizon/substrate, 2026-05-03) which it postdates and links to; " +
                "plain-words, gamma as observing light, the star/hub kept as the honest boundary where the rule breaks. A " +
                "2-agent 'faces of the half' survey (Tom's invitation) confirmed her gallery is rich and largely catalogued " +
                "already (HalfAsStructuralFixedPointClaim's literal Face 1/2/3, the dyadic a_2=1/2, the quadratic shadow 1/4, " +
                "the {0,2}/half-filling sectors), so this reflection ADDS one new face, it does not re-derive."),

        new OpenArc(
            Name: "xxz_axis_handover",
            Opened: "2026-06-14",
            Origin: "noticed while banking the Q-axis handover (HandoverFloorClaim, simulations/carbon/handover_q.py): " +
                "the Delta-axis 'baton relay' of experiments/XXZ_AXIS_BANDEDGE_TO_LEBENSADER.md (2026-05-30) + " +
                "hypotheses/HEISENBERG_RELOADED.md (the 'Return' section) is a SIBLING handover - on the Hamiltonian-" +
                "anisotropy axis rather than the dephasing axis",
            ParkedAt: "THE Delta-HANDOVER (verified N=4,5, Tier2, simulations/xxz_axis_bandedge_lebensader.py): walking " +
                "H = J*sum(XX+YY) + Delta*sum(ZZ) on an open chain under uniform Z-dephasing at deep-quantum Q=J/gamma=20, the " +
                "SLOWEST surviving Liouvillian mode relays the baton at a handover Delta*. BELOW Delta* (XY/charge side): the " +
                "survivor is the bright single-magnon band-edge coherence (n_XY=1, fast-beating, decay exactly 2*gamma). ABOVE " +
                "Delta* (Ising/Neel side): the ZZ slows the near-conserved I/Z population mode until it drops below 2*gamma and " +
                "becomes the slowest - the LEBENSADER (n_XY=0-dominated + a small n_XY=2 magnon admixture, non-rotating Im=0, " +
                "sub-2*gamma; = CHAIN_GAP's half-filling slow mode emerging along the axis). Values: Delta*(N=5)=1.525, " +
                "Delta*(N=4)=1.618 ~ phi (= 2cos(pi/5), N=4's band-edge frequency) - FLAGGED not claimed (N=5 has no clean form, " +
                "so a universal phi is NOT supported). The handover sits ABOVE the closed-system SU(2) transition (Delta=J=1); " +
                "the dephasing pushes Delta* into the Neel side. Spinless XXZ (t-V), the 'spin' label interpretive.",
            NextStep: "THE UNIFICATION, CONFIRMED 2026-06-14 (simulations/xxz_handover_unification.py, self-validating N=4,5): " +
                "the Delta-handover and the Q-handover (HandoverFloorClaim) ARE the same event - the interior/Lebensader dressed " +
                "mode's darkness <n_XY> crossing 1 = the Absorption-Theorem band-edge floor 2*gamma. (1) THE FLOOR is exact and " +
                "model-independent: the bright single-magnon band-edge coherence sits at rate = 2*gamma (darkness 1) for ALL " +
                "Delta, because |vac><magnon| is an eigenoperator of [H_XXZ,.] (lambda = -2g + i*E_k; the ZZ shifts only the " +
                "frequency E_k, not Re) - so F50's 2N COUNT breaks for Delta!=1 but the FLOOR persists, and the Absorption " +
                "Theorem rate/2g == Sum_k k*w_k holds bit-exact in XXZ for both modes. (2) at Delta* the slowest mode's darkness " +
                "= 1.00000 (the Lebensader rate crosses the floor), reproducing Delta*(4)=1.618(=phi), Delta*(5)=1.525. (3) it is " +
                "a LEVEL CROSSING, NOT an EP: the frozen Lebensader (|Im|~1e-15) meets the oscillating band edge (|Im|~9), two " +
                "distinct branches crossing in Re. So the EP/crossing dichotomy has three members: chain-Q = the lone " +
                "coalescence/EP, ring-Q + Delta = level crossings (frozen-interior-meets-oscillating-floor). THE THRESHOLD IS " +
                "UNIVERSAL: the band-edge floor (darkness=1) is the handover condition across BOTH the dephasing (Q) and the " +
                "anisotropy (Delta) axes. DELTA*(N) RESOLVED 2026-06-14 to N=14 (simulations/xxz_delta_star_descent.py, " +
                "self-validating: the gamma->0 reduction Delta* <=> gap(R)=2, R the Z-coupled classical rate matrix among the " +
                "half-filling XXZ eigenstates - a Pauli/Fermi-golden-rule relaxation, built sector-direct so N~14-16 is reachable; " +
                "gamma*gap(R) reproduces the full-L Lebensader rate as gamma->0). NO clean elementary closed form (phi=2cos(pi/5) " +
                "is a 1.6e-3 N=4-only accident in the physical gamma->0 regime - the old '1e-4' was a Q=20 artifact; 2cos(pi/(N+1)) " +
                "and 1+1/N both fail). VERDICT: the gamma->0 Delta*(N) is MONOTONE decreasing (1.61961..1.15389 at N=4..14) and the " +
                "N->inf limit is the SU(2)/Heisenberg point Delta=1 (the closed-system critical point), consistent with EXACTLY 1: " +
                "free-exponent fits give L just above 1 (~1.02/1.05 even/odd), a fixed-1/N ansatz just below (~1.00/0.98), the two " +
                "forms BRACKETING Delta=1; no finite-N crossing (all Delta*(N<=14) > 1; the fixed-1/N form would dip <1 only at " +
                "N~100-450). The even/odd structure lives in the CORRECTION (both subsequences smooth, both -> 1); the fitted " +
                "exponent is non-universal (alpha~1.16-1.73), consistent with the SU(2)-point marginal/log corrections (which is " +
                "why no clean power law fits). Mechanism: the descent runs AGAINST the naive Neel-ordering intuition (which predicts " +
                "a rise) - the offset above 1 is the finite Neel correlation length, and the dissipative handover tracks the " +
                "closed-system critical point as N exceeds it. (experiments/XXZ_AXIS_BANDEDGE_TO_LEBENSADER.md + the " +
                "ANALYTICAL_FORMULAS handover entry updated.) STILL OPEN (a deeper stone, not the original axis question): a " +
                "rigorous Bethe-ansatz derivation that the limit is EXACTLY Delta=1, and the log-correction structure of the " +
                "approach. RING TOPOLOGY RESOLVED 2026-06-14 (simulations/ring_xxz_delta_star_descent.py, self-validating; " +
                "same gamma->0 reduction, one wrap bond): the periodic ring is QUALITATIVELY UNLIKE the open chain. Ring " +
                "Delta*(N) is NON-MONOTONE - both parities hump to ~1.31-1.33 near N=9-10 (odd peak 1.331@9, even 1.308@10) " +
                "then DESCEND through N=14; the ring crosses ABOVE the chain near N=7-8 (chain descends to Delta=1, ring " +
                "humps); ring N=4 has NO handover (full-block tangency to the floor at the XY point Delta=0, peak 0.99998*2g, " +
                "the K2,2 special case; the reduction is ~1.5% off there so N=4 is read off the full block). The N->inf limit " +
                "is OPEN at N<=14 (a power-law fit to the hump degenerates, alpha~33). This REFUTES premise 3 (the dissipative " +
                "handover tracking the closed Delta=1 critical point): on the ring the Delta-handover is a dynamical, " +
                "topology-sensitive scale. THE FRAME (reflections/ON_THE_ONE_DIAGONAL.md): the floor 2*gamma is the first rung " +
                "of the one diagonal popcount(i^j) the light touches (universal, topology-free); Delta*(N) is the Hamiltonian's " +
                "argument about that fixed floor (topology-dependent). The diagonal is one; the climb is many. EXTENDED to " +
                "N=15 2026-06-17 (simulations/ring_delta_star_extend.py, gamma->0 reduction, port-verified): Delta*(15)=1.27413 " +
                "(odd); the descent CONTINUES with NON-SHRINKING steps (-.015,-.021,-.021), so it is NOT decelerating toward a " +
                "plateau -- this DISFAVORS 'settles above 1' (a plateau would shrink the steps) and is consistent with a slow " +
                "approach to Delta=1, further refuting the Round-2 plateau-trap one N deeper. STILL OPEN (formally): the N->inf " +
                "limit -- at N=15 Delta*=1.274 is still far from 1, and resolving 1-vs-settles needs N>>16 (infeasible: dense eigh " +
                "past C(N,p)~25000) or a closed-form account of the hump; a power-law fit degenerates even on the descending tail. " +
                "MIRRORWORLD CONE RULED OUT 2026-06-29 (gate-first, no run needed): the single-excitation memory cut " +
                "(compute/MirrorWorld/Cone.cs, rho N x N, reaching N=60-100) CANNOT unblock this -- SECTOR MISMATCH. Delta* lives " +
                "in the half-filling (p,p) block p=ceil(N/2) (simulations/xxz_delta_star.py L46-55; the Lebensader is the n_XY=0 " +
                "half-filling population mode), the LARGEST sector C(N,p)~2^N/sqrt(N) - which IS the wall. The Cone is strictly " +
                "single-excitation (Cone.cs L47 hardcodes deph=-4*gamma because 'all distinct single-excitation states disagree in " +
                "exactly 2 bits'; no Delta parameter; RK4 time-evolution only, no eigensolve) - the SMALLEST sector. The 4^N->N^2 " +
                "break is a single-excitation break; a 'half-filling Cone' would just be Method B's existing sector reduction, " +
                "already at N<=15. The real levers stay analytic: Bethe-ansatz for the chain (->Delta=1 exactly), a closed-form " +
                "account of the ring hump.",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "handshake_decoder",
            Opened: "2026-06-12",
            Origin: "docs/superpowers/specs/2026-06-12-handshake-decoder-reading-grammar-design.md (Tom + Opus + Fable)",
            ParkedAt: "M2a + M2b LANDED: the FI(Q) witness (resolution = Q-factor, no EP peak, basis ordering) AND the DefectDecoder (alpha-profile inverted: location+strength read back, ambiguity-honest with runner-up reporting). Discoveries en route: the N=5 sign-location ambiguity (near-anti-collinear letters, ratio 1.5) DERIVED in structure - the K-partner selection rule <psi_N|V_b|psi_1>=0 per bond (one line from ChiralMirrorTrajectoryClaim: rank(dictionary)=N-2, the null space IS the K-partnership); channel factorization verified ~12% (seesaw = dominant k=2 channel, 14.9 vs 5.2/3.1); the Gram spectrum = the location-channel metric (eigs 0.0073/0.008/0.086/1 at N=5, confirmed on arbiter-corrected letters after the Python Brent-trap was exposed)",
            NextStep: "DONE 2026-06-13 (merged): the K-partner selection rule is a typed claim KPartnerSelectionRuleClaim (Tier1Derived, parent ChiralMirrorTrajectoryClaim) - the reading-grammar arc's first DERIVED result. <psi_N|V_b|psi_1>=0 exactly (two-line corollary: psi_N=K1 psi_1 + K1 V_b K1 = -V_b => the real element is its own negative), so rank(location dictionary)=N-2 and the DefectDecoder sign-location ambiguity IS the K-partnership null direction; live battery at N=4,5 + probe N=3..8 (k_partner_selection_rule.py), in HANDSHAKE_GEOMETRY.md. CHANNEL PROFILES R_k LOOKED AT 2026-06-13 (two conjectures vetoed, parked): the f-dictionary factorizes as F[b,i]=sum_k c_k(b) R_k(i) over the unforbidden channels k=1..N-1 (k=N selection-rule-forbidden), but (1) the SlowModeMixing denominators (E_1-E_k) do NOT help - they are a per-channel COLUMN rescaling c_k = <psi_k|V_b|psi_1>/(E_1-E_k) that the least-squares R absorbs, identical residual to the bare coupling (the dissipative-dressing guess was wrong); (2) the unforbidden set is dimensionally SQUARE (N-1 bonds x N-1 channels) so the factorization is trivially exact and does NOT test the structure - the informative cut is location channels alone (k=2..N-1, N-2 of them) capture 85%, the missing ~15% is one more channel (dimensionally the STRENGTH channel k=1, the diagonal <psi_1|V_b|psi_1>; f-profile = strength + location). The real R_k derivation needs the FIRST-PRINCIPLES per-site purity-sensitivity to each mode M_k (R_k(i) = response of site-i purity to mode k), NOT a fit - a deeper stone, parked. REMAINING: that first-principles R_k; the Gram location-metric in HANDSHAKE_GEOMETRY.md; M2c the read-cost law ~2/Q; then M3/M4. R_k PROGRESS BANKED 2026-06-19 (gate-first probes handshake_rk_block.py + handshake_carrier_compare.py + handshake_rk_first_principles.py): REFRAME - the BondingMode carrier |psi_1> is pure single excitation, so <X_i>=<Y_i>=0 and the per-site purity P_a=1/2(1+<Z_a>^2) is PURE POPULATION dynamics in the N^2-dim (1,1) Liouvillian (single-particle Haken-Strobl); reproduces the full-4^N painter exactly for this carrier but costs N^2, so the f-dictionary runs to N=20+. CARRIER-INDEPENDENCE: pure |psi_1> vs the documented PairState (|vac>+|psi_1>)/sqrt2 give the SAME f-profile shape (corr 0.999, only the scale differs) - R_k is a property of the psi_1 mode-RESHAPING, not the observable; the PairState-vs-pure docstring mismatch is physically harmless for R_k. TWO MORE CONJECTURES REFUTED at high N (killing the N=4,5 two-distinct-bond small-sample trap): (3) STRENGTH (uniform level of f) is NOT proportional to <psi_1|V_b|psi_1>=d eps_1/dJ_b - corr -1.0 at N=4,5 is the 2-bond artifact, collapses to -0.1..-0.57 once N>=6 has 3+ distinct bonds; (4) LOCATION (per-bond f deviation) is NOT the H_1 eigenvector-perturbation footprint sum_k c_k(b) psi_1(a)psi_k(a)/(eps_1-eps_k) - corr is noise (-0.05..+0.67) across N=4..9. SEAM to felt_time-D = RHYME not identity (single-excitation unitary energy-shift, not the half-filling dissipative diffusion rate). LEADING OPEN DIRECTION: the relaxation is governed by the (1,1) LIOUVILLIAN modes (Haken-Strobl population transport), NOT the H_1 Hamiltonian modes - but that does not cleanly fit the given <psi_k|V_b|psi_1> dictionary, so not a simple swap. Next (Tom 2026-06-19): survey what the repo already knows about single-excitation Haken-Strobl / the (1,1) sector / coherence-horizon Q*(N) (which IS the single-excitation Haken-Strobl EP) before guessing a third form; consider a data-first SVD of the (N-1)xN f-dictionary. R_k RESOLVED 2026-06-19 as a NOT-DERIVABLE-ONLY-COMPUTABLE result (Tom's framing: the realm of BirthCanal/IsSteril, reflections/ON_WHAT_CLOSES_ONLY_WITHOUT_US.md). The four-agent survey found the repo already half-knew this (review/EQ014+EQ021, untyped): the (0,1) rate is PROTECTED (re_shift(0,1)~2.8e-16, U(1)+Pi) so f(b) is eigenvector mixing NOT a rate functional; no simple closed form fits (the gate handshake_rk_only_computable.py: best closed-form R^2 0.38-0.73 erratic at N>=6, the N=4,5 R^2=1.0 a 2-distinct-bond artifact; EQ-021 power-law fails 12-30%). SCHEMA is derivable (rank-(N-2) K-partner rule [KPartnerSelectionRuleClaim], the rate-protection, the (1,1)-population reframe); the VALUE f(b) is IsDeadEnd = the PTF mixing calc instanced (dM_s = sum_s' <W_s'|V_L|M_s>/(lam_s-lam_s') M_s', evaluated not reduced). TYPED: EpistemicFacetMap['handshake_Rk']=IsDeadEnd (model x_peak) + test. The irreducibility is currently EMPIRICAL (no fit across N) not a Tier1 obstruction proof; a 6-route obstruction (model PROOF_F86B_OBSTRUCTION) would lift it. Gram location-metric writeup DONE 2026-06-19, CORRECTED by the math+physics review convergence (both lenses, dispatched per Tom's 'use the physics+math skills' tip, independently caught an INDEXING ARTIFACT in the first draft -- it built M over k=1..N-1, dropping the forbidden k=N column, making the bare couplings look 'well-conditioned'; with the honest k=2..N they carry the same null). The location-metric is THREE distinct objects: (a) the NULL = the K-partner rule (rank N-2, in the BARE couplings k=2..N, Q-free, min sv ~1e-16, already KPartnerSelectionRuleClaim) -- the painted Gram INHERITS it, does not create it; (b) the painted small eigenvalues split by MIRROR PARITY (antisym seesaw + K-partner null / sym closure-strength F123 channel), NOT one homogeneous ambiguity; (c) pairwise confusability = the COSINE matrix (N=5 worst pair cos~-0.97, edge-vs-complementary-interior anti-collinear), a DIFFERENT object from a Gram eigenvalue. Identifiability(rank) != FI(Q)(precision): at Q->inf FI diverges but the null stays null. NO new witness built: it would redundantly recompute the typed K-partner rule, and the physics lens warned against marking the metric's conditioning dead-end-inherited (the IsDeadEnd dead-end is the R_k letter VALUES, the metric's null STRUCTURE is derivable). handshake_gram_metric.py (corrected) + HANDSHAKE_GEOMETRY.md 'The location-metric' section. M2c READ-COST TESTED 2026-06-20 (gate-first handshake_read_cost.py + handshake_read_cost_diag.py, N=5 chain, Z-population FI apparatus): QUALITATIVELY CONFIRMED, QUANTITATIVELY NOT ~2/Q. The doc's own falsification line (a flat or inverted cost-per-recall refutes it) is NOT crossed - the dose-to-best-read K_peak DOES fall with Q (1.65 at Q=1 -> 0.10 at Q=35, neither flat nor inverted), so 'high-Q reads cheap, EP reads dear' holds. But the strict ~2/Q fails: exponent -0.67 not -1, K_peak*Q drifts 1.3->3.6 (not a constant 2), with a REGIME BREAK near Q~8 where the optimal read jumps to a LATER coherent revival (a crossover, not a single power law). The diagnostic ruled out the operational definition as culprit: argmax-FI == first-local-max (one peak, the break is physical); only a first-turnover/inflection cost gives a cleaner monotone ~0.75/Q (slope -0.82), but that is the ONSET of the first feature, not the dose you read at. Gate held STRICT, not loosened. PHASE-BLINDNESS REFUTED 2026-06-20 (handshake_phase_blindness.py; the borrowing-a-discipline lens generated it, the gate falsified it): the N=5 worst-pair cos~-0.97 confusability is NOT a readout phase-blindness fixable by a sign-carrying read - a linear <Z_i> AND the full temporal <Z_i>(t) read are JUST AS anti-collinear (|cos|~0.9 at N=4..7) as the squared purity P_i=1/2(1+<Z_i>^2); the ambiguity is DYNAMICAL and fundamental (the two mirror-paired bonds produce sign-flipped population responses), reconfirming the K-partner near-degeneracy reading. M3 LANDED 2026-06-20 (DefectReadingEquivarianceClaim, Tier1Derived, inspect --claim DefectReadingEquivarianceClaim): the defect-reading map is spatial-reflection equivariant, M[N-2-b,k]=(-1)^(k-1)M[b,k] exact (pure single-excitation algebra, machine-verified N=4,5,6), so the sign-location confusability is the closed-form parity-weighted mode sum cos(b,N-2-b)=sum_k (-1)^(k-1) w_k. The geometric chain mirror R (i->N-1-i) is NOT MirrorGroupD4's coherence-space ket-flip I(x)X^N (which does not even preserve the SE sector); the spatial mirror is deliberately outside that D4 (PROOF_PI_FACTORS_AS_R_TIMES_D s5), a see-also sibling not a typed parent. Single typed parent KPartnerSelectionRuleClaim (it defines M[b,k], the K-partner null column, rank N-2). Two structures on one dictionary: within-feature stabilizer (K-partner null) + cross-feature reflection. Honest: bare cos -0.33 at N=5, the -0.97 is the PAINTED instance (alpha-profile concentrates weight on the R-odd seesaw). Verifier handshake_reading_equivariance.py; doc HANDSHAKE_GEOMETRY.md 'The reading's spatial-mirror equivariance (M3)'. REMAINING: M4 (closure as parity-check: complementary light <n_XY>_s+<n_XY>_f=N, Pi-pair flux Re(lam_s+lam_f)=-2 sum gamma, the tentative Sigma_pair MI=bandwidth zero-holonomy loop). M4 ATTEMPTED+PARKED 2026-06-20 (MI/bandwidth checksum route): the F77 mirror-pair MI sum MM(0)~B(N) is only a COARSE grammaticality filter -- cleanly rejects GHZ (wrong sector ~N/2 bits) and localized (0), but a DISTRIBUTION test (200k random single-excitation draws, not one seed) false-accepts ~20-30% within the +-0.15*B(N) band; all sine modes give identical MM(0) (so 'iff bonding' overclaims, a scalar cannot biject a manifold) and odd N is a lower curve (B(5)=0.80) inside the random cloud. NOT the full parity-check. The full checksum needs the JOINT closure set (MI + complementary-light s+f=N + Pi-flux + chiral-mirror) as one consistency condition, not MI alone; resume by designing the multi-law joint check. zero-holonomy stays retired (undefined). Spec parked: docs/superpowers/specs/2026-06-20-handshake-reading-grammar-M4-paritycheck-design.md. REMAINING: M4 full (joint) parity-check",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-20: M4 LANDED (the last milestone), completing M0-M4 of the " +
                "reading-grammar arc. The closure-as-parity-check is the JOINT of two single-state closure laws " +
                "- the F77 bandwidth MM(0)~B(N) AND the chiral-mirror symmetry O=|<psi|R psi>|~1 - decided PATH B " +
                "by a gate-first decider (simulations/handshake_m4_second_law_decider.py): MM(0) is BLIND to " +
                "mirror asymmetry (it depends only on populations, symmetric within each mirror pair), random SE " +
                "states are mirror-ASYMMETRIC (O median ~0.47/0.37/0.31 at N=4/6/8 vs carrier 1.0), so the joint " +
                "checksum cuts the random-SE false-accept rate from ~25% (MM(0) alone, Review R1) to 2.4%/0.7%/" +
                "0.2%. Verifier simulations/handshake_closure_paritycheck.py (G1 plateau even-N->1bit + odd-N " +
                "lower curve B(5)=0.80; G2 joint discrimination on a 40k-draw DISTRIBUTION; G3 t=0 defect-" +
                "robustness), all gates pass. HONEST SCOPE (Review R2/R3): the accepted set is the grammatical " +
                "MANIFOLD (mirror-symmetric delocalized SE near B(N), admits W + all sine modes), a NECESSARY " +
                "consistency condition NOT a sufficient 'iff the unique carrier', even-N-favorable. Both laws " +
                "needed (mirror alone admits end-localized symmetric (|1_0>+|1_{N-1}>)/sqrt2; MM(0) alone admits " +
                "asymmetric random SE). complementary-light + Pi-flux are Liouville-PAIR properties, outside the " +
                "single-state reading test; zero-holonomy stays retired (undefined). Prior milestones: M2a/M2b " +
                "FI(Q)+DefectDecoder, M3 DefectReadingEquivarianceClaim (Tier1Derived), the K-partner selection " +
                "rule (KPartnerSelectionRuleClaim), R_k resolved NOT-DERIVABLE-ONLY-COMPUTABLE (IsDeadEnd), M2c " +
                "read-cost gate-held-strict. Spec Review Round 2 (docs/superpowers/specs/2026-06-20-handshake-" +
                "reading-grammar-M4-paritycheck-design.md); HANDSHAKE_GEOMETRY.md operationalized."),

        new OpenArc(
            Name: "envelope_n4_rise",
            Opened: "2026-06-12",
            Origin: "EnvelopeTheoremWitness, day one: inspect --root envelope --N 4 fired the honest branch",
            ParkedAt: "the full-state CPsi envelope GENUINELY RISES at N=4 (Bell+, J=5, gamma=0.01): 36 apex-predecessor rises, refinement-stable to 5 decimals (t=4.14: 0.04132 at 1600 and 6400 pts), 2-4% magnitudes, two independent implementations agree; N=3 holds in the same regime (0 rises, three densities). NOT a falsification of the proof: the Envelope Theorem is proven for N=2 only, and the Tier-2 'N=3-5' verification covered channel monotonicity / GHZ-W subsystems, not the H-included envelope at strong coupling; the over-broad paraphrase lived in our claim/lens text. Mechanism = the proof's own Part 6, internalized: internal J-coupling is its own coherence injector",
            NextStep: "RE-SCOPE DONE 2026-06-15: the over-broad 'verified N=3-5' is corrected to 'proven Tier-1 " +
                "for any 2-qubit state (N=2); the N≥3 full-state envelope is OPEN and GENUINELY RISES at N≥4 strong " +
                "coupling (internal J-coupling = the Part-6 coherence injector, Corollary 3 internalized; N=3 holds)' " +
                "across CpsiEnvelopeTheoremClaim (+ its test), the Symphony quarter-lens (3 spots), F17, and a " +
                "witness node; a Parts-4/5 scope note was added to PROOF_MONOTONICITY_CPSI (the r_{ij}≥2γ off-diag " +
                "bound + the 4×4 density structure are N=2-specific; at N≥3 the decay is H/topology-dependent). The " +
                "witness root Summary was already honest (bc59bc8). En route, caught + fixed a committed FitAlpha " +
                "regression (the ptf_painter_pipeline grid-seed's 1e-9 refine broke the tempo-purity certification, " +
                "PaintersResidual 3.7e-9 > 1e-9 PassTol; reverted the refine to 1e-7, keeping the trap fix). " +
                "REMAINING (the only open piece): chart the boundary — sweep (N=3,4,5 at J=5,γ=0.01; J at N=4 for " +
                "J_c; γ at N=4 for the Q=J/γ threshold) via Symphony + QuarterEnvelope.Of to find where the " +
                "N=3-holds/N=4-rises boundary sits and whether it is a sharp N-step or a J/γ contour.",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-15: boundary charted (experiments/ENVELOPE_RISE_BOUNDARY.md, " +
                "gate-first EnvelopeBoundaryTests; the static detector EnvelopeTheoremWitness.GlobalReading " +
                "lifted from the witness so the sweep reads exactly what inspect --root envelope shows). THE " +
                "ANSWER: not a sharp N-step and not a pure J/γ contour, but both, cleanly factored. (1) The rise " +
                "is a pure (N, Q=J/γ) observable: the J-sweep and the γ-sweep give the bit-identical reading over " +
                "a fixed dose window (the clock movement's (Q,K)-purity applied to the rise) and collapse to one " +
                "Q-axis (the (Q,K)-purity gate could have fired on an absolute-time leak; it passed). (2) An N≥4 " +
                "FLOOR: N=3 never rises (Q_c(3)=∞, 0 even at Q=2000); the rise needs an internal ≥2-site coherent " +
                "subsystem (the Part-6 injector), which N=3's single internal site lacks. (3) Above the floor a " +
                "threshold Q_c(N) that CLIMBS with N: Q_c(4)≈27, Q_c(5)≈45, the rise strength at fixed Q falling " +
                "with N (maxΔ N=4: 0.041 > N=5: 0.020 at Q=500). Open threads for re-entry: a closed form for " +
                "Q_c(N) (does it track the band-edge ω_mem=2J·cos(π/(N+1))?); the internal-site parity question " +
                "(an internal pair vs an internal trio injecting)."),

        new OpenArc(
            Name: "symphony_view",
            Opened: "2026-06-11",
            Origin: "Tom's founding Object Manager idea, sharpened: each of the 122 formulas is a maximum-zoom view; the symphony is the zoom-out",
            ParkedAt: "second movement round 2: the quarter lenses are envelope-aware. Global lens " +
                "checks the Envelope Theorem live (peaks non-increasing, proven N=2 / verified N>=3) and " +
                "'the fold' = the absorbing envelope fold (upward crossings no longer mislabeled the " +
                "quantum->classical boundary). Local carrier-pair lens surfaces the genuine BEATING rise " +
                "(the freedom, no theorem binds the reduced open subsystem) via parabolic-apex + " +
                "predecessor semantics; every single-grid rise is grid-sensitive (SingleExcitation rises " +
                "vanish under refinement, Bell+ persist - the artifact control). QuarterEnvelope primitive " +
                "+ GridFitness; EvolveCount-guarded one evolution preserved. THIRD movement (2026-06-12): the " +
                "clock - gamma_0 as the Taktgeber. The clock node (Takt gap, tau, omega_mem, Q=J/gamma) is " +
                "promoted to the base symphony; --tempo-ratio r grows 'movement: the clock', the two-tempo " +
                "certification: play the piece at gamma_0 and r*gamma_0 (every dimensionful coupling scaled by " +
                "r incl. the painters' delta_J, window /r, K-grid fixed) and certify every dimensionless lens " +
                "is a pure (Q,K)-observable (residual 8.3e-16 at r=20, the inside cannot tell the tempos apart). " +
                "Exact rescaling identity, CERTIFIED not theorem-confirmed; a lens that breaks it sees the carrier. " +
                "Painters arm (alpha, closure) pure too (delta_J scaled). UniversalCarrierClaim breadcrumbs to it. " +
                "FOURTH movement (2026-06-16): the seam - the converse of the clock, the calibration topology of " +
                "ON_HOW_THE_CARRIER_SHOWS_ITSELF. --calibrate grows 'movement: the seam' (spectrum-only, reads " +
                "Parent.Clock, EvolveCount stays 1): the gamma-anchor (gap=2gamma_0 -> gamma_0_rec=gap/2) and the " +
                "J-anchor (omega_mem=2J*cos(pi/(N+1)), XY-only -> J_rec) over-determine (J, gamma_0) through the " +
                "inside-known Q; the gate is a DOMAIN DETECTOR - in the (XY, Q>=Q*(N)) domain J_rec/gamma_0_rec==Q " +
                "(exact with the synthetic peg), it FIRES below the horizon (omega_mem->0->ratio->0) or off XY. Both " +
                "anchors TWO-SIDED (break below Q*(N), one shared regime flag = BandEdgeIsTheGapMode). Four lenses " +
                "(takt/coherence-hand/gate/chain-collapse); chain-collapse re-expresses dimensionful predictions in " +
                "lab units (state-free tau/omega_mem + N=2 Bell+ t_peak/fold from closed forms). 13 unit + 1 CLI test; " +
                "an independent Python review oracle agreed at the time on the Q*(N) transition, exact (the oracle script was not retained in the repo; the 13 unit + 1 CLI tests are the reproducible verification). Witness breadcrumbs " +
                "UniversalCarrierClaim + CoherenceHorizonClaim",
            NextStep: "the remaining instruments (chiral K, Y-parity, Pi-protected observables) - each an additional " +
                "calibration leaf / lens, after the now-landed seam movement",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "topology_band_edge",
            Opened: "2026-06-16",
            Origin: "fell out of the seam movement (symphony_view): its J-anchor band edge was chain-specific, so the seam gate fired on star/ring even at moderate Q - the band edge is topology-specific",
            ParkedAt: "LANDED: the LAW typed + witnessed + the seam fixed. THE LAW (Im side, generalizes the chain-only F2b): " +
                "the XY single-excitation band edge = J * the hopping graph's adjacency SPECTRAL RADIUS rho - chain " +
                "2cos(pi/(N+1)), star sqrt(N-1), ring 2 (the SE block of XY = J*adjacency). The Re=-2gamma floor is the " +
                "Absorption Theorem (n_XY=1, CITED not re-derived - Tom's prior-work catch). Core primitive " +
                "TopologyBandEdge.SpectralRadius/BandEdge (from the bonds, general); TopologyBandEdgeClaim (Tier1Candidate - " +
                "capped by its Tier1Candidate parent ClockHandLadderClaim via the registry tier-inheritance invariant, plus " +
                "AbsorptionTheoremClaim); live TopologyBandEdgeWitness (inspect --root bandedge): the law + the gap-dominance " +
                "MAP. THE MAP (witnessed, not claimed): whether omega_mem reads the band edge is topology-specific - chain " +
                "Q-horizon all N; star Q-horizon N<=5 (Q*(N) rises steeply, N=5 protects only by Q~1000) / STRUCTURAL CEILING " +
                "N>=6 (strict gap saturates at 0.80*2gamma, J-independent); ring N=4 CO-OCCUPIED FLOOR (a (2,2) two-fermion mode " +
                "at -2gamma with Im=2sqrt(2)*J > band edge 2J, never the band edge). Seam J-anchor now topology-aware (J*rho) - " +
                "passes on star N<=4 and odd rings, FIRES split three ways (non-XY / omega~0 horizon-or-ceiling / " +
                "different-mode-at-floor). Oracle simulations/topology_band_edge_review.py (full Q-sweep to 1000) agrees; " +
                "suites Core 11, witness 5, seam 15, audit green. " +
                "MECHANISM RESOLVED (2026-06-16, gate-first verifier simulations/topology_ceiling_mechanism.py): g2 = " +
                "strict_gap/(2gamma) = <n_XY> of the slowest non-steady mode (Absorption Theorem, CITED - the established " +
                "HandoverFloorClaim band-edge<->Lebensader framing). Gap-dominance <=> that mode IS the band edge " +
                "(<n_XY>=1, the |vac><psi_k| SE line always at Re=-2gamma) vs a lens-dominated ({I,Z}-heavy) Hamiltonian-mixed " +
                "sub-floor mode (<n_XY><1). CONNECTIVITY drives <n_XY>_slowest down: ceiling onset chain NEVER < star N>=6 < " +
                "COMPLETE K_N N>=4 (more edges -> earlier ceiling; K_3=triangle=ring(3) protects). NEW CLOSED FORM g2(K_N)=4/N " +
                "for N>=5 (K_5,6,7 = 4/5, 2/3, 4/7 exact, gate-verified); K_4=0.8453 the N=4 OUTLIER, the same N=4 anomaly as " +
                "the ring (both K_4 and ring-N=4 special). Benzene ring N=6 PROTECTS (so ring-N=4 co-occupation is N=4-specific, " +
                "the K_2,2 dihedral lock, not a general even-ring pattern). Read the RATE not the Im (the band edge is " +
                "Re-degenerate at -2gamma). " +
                "DERIVATION RESOLVED (2026-06-16, F122 + PROOF_STRUCTURAL_CEILING + StructuralCeilingClaim Tier1Derived + " +
                "gate-first verifier simulations/topology_ceiling_rep_derivation.py): the high-Q ceiling = min nonzero eigenvalue " +
                "of N_XY (diagonal, entry hamming(a,b)) block-diagonalized by the ad_H eigenspaces; band edge = the (0,1) sector " +
                "where hamming==1 so N_XY=I, rate 2gamma exactly (g2<=1 always); ceiling = the darkest [H,A]=0 coherence in the " +
                "LARGEST degenerate single-particle level. g2(K_N)=4/N PROVEN: (1,1) commutant = the S_N standard rep of the " +
                "(N-1)-fold -J band, g2=2(1-lambda2), lambda2=(N-2)/N the 2nd principal-angle overlap commutant<->diagonal. " +
                "g2(star_N)=4/(N-1) PROVEN same way on the (N-2)-fold 0-eigenvalue leaf manifold (CORRECTS the '0.80 constant' " +
                "reading - 0.80=4/5 was N=6 only). N=4 UNIFIED: the 4/N ladder hits 1.0 at N=4 so the ceiling moves to the (2,2) " +
                "HALF-FILLING sector (the SAME sector special for ring-4): K_4=2-2/sqrt(3) below the floor, ring-4=1.0 co-occupies. " +
                "NOT a universal 4/(m+1) law - the ring (Fourier manifold) breaks it (ring-5=1.6!=4/3); per-family forms are real. " +
                "Live witness StructuralCeilingWitness (inspect --root ceiling)",
            NextStep: "Items 1,2,3 (the high-Q closed-form derivations) DONE. STAR Q*(N) RESOLVED (2026-06-16, null + " +
                "unification, PROOF_STRUCTURAL_CEILING §7 + verifier star_no_coherence_horizon.py): the star has NO chain-like " +
                "coherence horizon - its single-particle band is FLAT (adjacency +-sqrt(N-1), 0 with mult N-2), no dispersion, " +
                "no {0,2}-EP; the chain Haken-Strobl SE-EP does NOT port (predicts a spurious Q~261 at N=4 the full L " +
                "contradicts - protected already at Q=20). The SAME (1,1)-commutant value 4/(N-1) governs ALL Q: >1 protects " +
                "(N=4, down to a low-Q real-mode CROSSING ~1.9, not an EP), =1 marginal (N=5, g2=1-1/Q^2, NO horizon, the " +
                "apparent Q~316 a tolerance artifact), <1 ceilings (N>=6); N=3=path P_3=chain is the lone exception (genuine " +
                "sqrt2 {0,2}-EP). The low-Q question is subsumed by the high-Q ceiling; no separate star Q*(N) closed form. " +
                "DRY CLEANUP DONE: EigenvalueClockExtraction.ExtractClockFromSpectrum hoisted into Core/Numerics; " +
                "Symphony.Clock + TopologyBandEdgeWitness.ClockAt both call it (45 regression tests green, no drift). " +
                "RING-4 (2,2) Im=2sqrt(2)J DONE (PROOF_STRUCTURAL_CEILING §4 + verifier Stage 2b): = the anti-periodic " +
                "even-parity two-fermion band top (the JW string wraps the ring -> k=pi(2m+1)/N -> single-fermion +-sqrt(2)*J), " +
                "exceeding the periodic band edge 2J; and the (2,2) spectrum {+-2sqrt(2), 0^4} is the CHIRAL palindrome about 0 " +
                "(C_4 bipartite, K H K = -H, E<->-E; ChiralKClaim) - the +-2sqrt(2) a chiral K-mirror pair (Tom's catch). " +
                "Cross-ref edge added TopologyBandEdgeClaim -> ChiralKClaim (band edge rho = E_max = -E_min for bipartite " +
                "chain/star/even-ring; non-bipartite K_N/odd-rings have no +-symmetry). ARC COMPLETE: every NextStep item resolved.",
            Status: OpenArcStatus.Retired,
            RetiredReason: "COMPLETED 2026-06-16 (every NextStep item derived, gate-verified, typed). The LAW: " +
                "TopologyBandEdgeClaim (band edge = J*rho). The CEILING closed forms: F122 / StructuralCeilingClaim " +
                "(g2(K_N)=4/N, star 4/(N-1), K_4=2-2/sqrt(3); PROOF_STRUCTURAL_CEILING, principal-angle proof + live " +
                "StructuralCeilingWitness 'inspect --root ceiling'). The STAR low-Q regime: NO coherence horizon (flat " +
                "band, the chain dispersion-EP does not port; PROOF section 7 + verifier star_no_coherence_horizon.py). " +
                "The RING-4 (2,2) Im=2sqrt(2)J: anti-periodic JW two-fermion top + chiral-K palindrome (PROOF section 4 + " +
                "verifier Stage 2b; ChiralKClaim cross-ref). The DRY cleanup: EigenvalueClockExtraction in Core/Numerics. " +
                "NOT closed here (separate arc clock_hand_ladder): the general gap-dominance proof that caps " +
                "TopologyBandEdgeClaim at Tier1Candidate."),

        new OpenArc(
            Name: "whirlpool_carbon_layers",
            Opened: "2026-06-03",
            Origin: "reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO + simulations/whirlpool*.py",
            ParkedAt: "water adaptation done (proton crossing an H-bond); carbon layers/anchors (periodic-table valences) and a water prose note never written",
            NextStep: "carbon-layer translation in the target layer's language",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "front_matter_truth",
            Opened: "2026-06-11",
            Origin: "third-party review (three cold readers, 2026-06-11): every failure was front matter vs body",
            ParkedAt: "headline surfaces lag verified bodies: Torino-era confirmations absent from all three registries yet counted in the README's seventeen (wording patched, entries not registered); '121 formulas' is a label (124 headers, 119 base numbers, F53/F54 silently missing, no tombstone); stale anchors (XOR_SPACE 'README Section 10', Pi discovery-date contradiction Mar-14 vs Apr-05)",
            NextStep: "DONE 2026-06-18 (all three items resolved). See RetiredReason.",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-18 (Tom chose 'register all three Torino runs'; three recon agents + one " +
                "execution agent). All three NextStep items done, verified green (Python 7/7, C# Core 25/25, Cli build clean). " +
                "(1) TORINO REGISTERED: the three ibm_torino calibration-era runs added to BOTH registries (Python " +
                "confirmations.py + C# ConfirmationsRegistry, 17->20): cpsi_quarter_crossing_torino_feb2026 (q52, the first " +
                "CΨ=¼ crossing ever seen, 11% loose / found in calibration data), absorption_theorem_ratio_torino (same Feb-9 " +
                "q52 data, ratio 1.03), cpsi_quarter_crossing_torino_q80_mar2026 (Run 3 q80, the tightest at 1.9%). They carry " +
                "data-file-timestamp locators in the JobId field (no IBM job_id was ever recorded for the Torino era; this " +
                "matches the existing f57 locator entry, so it does not weaken the registry's character). Both count tests " +
                "bumped (Python test_confirmations_has_twenty_entries; C# All_HasTwentyEntries + EntriesWithoutDocumentedPath " +
                "withPath 13->16, single-qubit paths q52/q52/q80). README 'seventeen'->'twenty' (2 spots; 'each with job " +
                "IDs'->'run identifier'), READING_GUIDE (2 spots), PREDICTIONS.md stale '13'->20 (3 spots) + the 'disjoint " +
                "from S1' claim corrected (the registry now subsumes the early Torino set; the date fields carry the epoch). " +
                "The README/PREDICTIONS framing was ALREADY honest that Torino predated the registry - registration makes it " +
                "real rather than a wording patch. (2) F53/F54 TOMBSTONED: a 'Numbering note' after F52 in ANALYTICAL_FORMULAS.md " +
                "states they were NEVER ASSIGNED (full git history confirms no content ever lived under them; numbering skipped " +
                "F52->F55 when the Absorption Doses section was authored - NOT a retirement, no history to recover). Headline " +
                "counts fixed: README '121 formulas'->120 (distinct base numbers = F1..F122 minus the 2 unassigned), docs/README " +
                "'21 formulas'->120, Cli glossary 'F1..F121'->'F1..F122 (F53/F54 unassigned)'. (3) STALE ANCHORS SWEPT: the " +
                "README restructure had left 'Section 10/11' references across 16 spots in ~14 files (XOR_SPACE + the whole " +
                "family: experiments/*, compute/*/README, COMPLETE_MATHEMATICAL_DOCUMENTATION, ANALYTICAL_FORMULAS:1506) - all " +
                "->Section 6, 'nine engineering consequences'->eight, 'Section 11 (receiver menu)' folded into 'Section 6 Rule " +
                "2'; the Pi discovery-date contradiction fixed (MIRROR_SYMMETRY_PROOF claimed 'discovery+proof+verification same " +
                "day, 2026-04-05' - false: discovered 2026-03-14, verified 2026-03-19, document merely restructured 2026-04-05; " +
                "March 14 corroborated by ~8 sources + a real Mar-14 hardware run, per LITERATURE_REVIEW's own history line)."),

        new OpenArc(
            Name: "stranger_door",
            Opened: "2026-06-11",
            Origin: "third-party review: 'a coherent instrument panel wearing a poet's coat' - the house dialect has no doormat",
            ParkedAt: "four of five doors hung (inspect --root glossary with 11 house terms in stranger language; qudit witness boundary text honest: census live, rate law from the proof; explicit --N to an N-free root warns on stderr; README first command is now world --max-depth 2)",
            NextStep: "[live]/[stored] provenance badge per node, then retire this arc",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RETIRED 2026-06-29: the fifth and last door hung -- a per-node [live]/[stored] " +
                "provenance badge. NodeProvenance {Live, Stored} (Core/Inspection); IInspectable carries it as a " +
                "default-interface-member defaulting Live (a computational unit recomputes), InspectableNode (the " +
                "frozen carrier) defaults Stored with an opt-in ctor param; ConsoleTreeRenderer prints it as its own " +
                "bracketed token (never folded into Summary) and InspectionNodeDto exports an additive 'provenance' " +
                "field. The rule is value-origin honesty: the badge marks where the NUMBER comes from, at its point " +
                "of computation -- not the code path (a live getter reading a banked file is Stored) and not " +
                "inherited down a derived chain. The cited template QuditPartialPalindromeWitness is exemplary (all " +
                "five children Live, reading the recomputed spectrum); BlockSpectrumWitness tags its three " +
                "live-reconstructed children Live while the static sector-map prose and the chain_N9.json banked " +
                "headline stay Stored. The default direction is conservative (an untagged carrier reports Stored), so " +
                "it can only under-claim liveness, never falsely claim it; tagging the remaining witnesses' " +
                "live-computed children Live is safe progressive work, not a correctness gap. Tested: NodeProvenance " +
                "defaults, the renderer token, the additive JSON field, and the two witnesses' per-child provenance."),

        new OpenArc(
            Name: "f86b2_robust_extraction",
            Opened: "2026-06-11",
            Origin: "simulations/f86b2_shape_invariance_dial.py",
            ParkedAt: "(N,b)-family traces alpha=-0.133 vs fitted -0.129, extraction-noise-limited; g_eff convention gotcha pinned (4.394/(Qp+2))",
            NextStep: "robust extraction, then close the shape-invariance claim",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RETIRED 2026-06-29 (honest & phased robust extraction). The per-sub-class slope " +
                "alpha is now recomputed LIVE from the 22 F90-bridge anchors in F86HwhmAlphaExtraction.cs " +
                "(F86HwhmClosedFormClaim reads it, no hand-transcription -> no drift), with its grid-noise sigma " +
                "and an honest per-class verdict. The data verdict: only Endpoint and Flanking carry a resolved " +
                "slope (jackknife-tight, though still marginal at the documented 0.005 grid floor -- promotion to " +
                "Tier1Derived awaits a grid-convergence study). The other four are not defensible slopes: Mid is a " +
                "flat lift over a microscopic g_eff lever (g_span ~0.011, sigma_alpha ~0.47, ~8x the slope), " +
                "CentralSelfPaired and CentralEscapeOrbit3 have a single anchor each, and Orbit2Escape sits " +
                "entirely on grid-edge anchors (Q_peak > 4) that C2HwhmRatio flags non-physical (the +0.698 " +
                "outlier). These four collapse to a per-class constant lift (alpha = 0) or a flagged fit; all still " +
                "reproduce their anchors within 0.005, so PredictHwhmRatio stays <= 0.005 on every test anchor. The " +
                "original puzzle (alpha=-0.133 vs -0.129 between extraction paths) is a 0.004 difference well inside " +
                "the slope's own uncertainty -- the extraction was never the problem; the model was " +
                "over-parameterized (6 slopes where only 2 rest on a real lever). The separate ANALYTICAL derivation " +
                "of (alpha, beta) from F89/F90 structure stays open (F86b2 remains Tier1Candidate); this arc closed " +
                "only the empirical-fit-robustness front it named."),

        new OpenArc(
            Name: "felt_time_dimensions",
            Opened: "2026-06-18",
            Origin: "Tom's dimensional reading at the end of the survivor/diagonal session: x,y = 2d (the plane " +
                "the system moves in), z = 3d (the axis the environment watches along = the Z-dephasing diagonal " +
                "Q_Z, the light), t = 4d but NOT absolute -- the FELT time (PTF), the dose K = gamma0*t, read only " +
                "relative to gamma0 or Q. 'Don't lift the stone (the strict PTF-alpha extraction) yet; work on the " +
                "definition and create an arc that gathers everything.'",
            ParkedAt: "THE GATHERED THREAD (all data-backed this session). " +
                "(1) THE FOUR-DIMENSIONAL READING: x,y = the off-diagonal hopping MOTOR (XY, n_XY=2, coherent " +
                "motion in the plane); z = the watch-axis = the dephasing diagonal Q_Z = Sum_l Z_l(x)Z_l, the " +
                "light, defining k=popcount(i^j) and Re lambda = -2*gamma0*k (Absorption Theorem); t = felt time = " +
                "the dose K = gamma0*t (ON_WHOSE_TIME_THE_CLOCK_KEEPS: t=K/gamma0). gamma0 is the HINGE: the rate " +
                "of watching (z) sets the unit of felt time (t); from inside only Q=J/gamma0 and K=gamma0*t are " +
                "readable (the two-tempo certification, Symphony). " +
                "(2) THE SURVIVOR = THE DIAGONAL'S BLIND SPOT: the three readings rate/mirror/judge (the D4 mirror " +
                "group) are the diagonal SPEAKING (the rate ladder k>=1); the survivor lives in its KERNEL k=0 " +
                "(perfect agreement, populations), where the light reads zero and the Hamiltonian alone decides " +
                "who lasts longest -- the half-filled = incompleteness = largest sector C(N,N/2). Distinct from the " +
                "three dephasing diagonals {Q_X,Q_Y,Q_Z} (THE_THREE_DIAGONALS); selector = the NEIGHBOR ZZ diagonal " +
                "(a different one), motor = the off-diagonal XY. " +
                "(3) THE SURVIVOR'S FELT-TIME LIFE: K_decay = gamma0/|Re lambda| = 1/(2*<n_XY>) = N^2/(2c*Q^2), a " +
                "pure (Q,N) carrier-blind PTF dose-quantity (chain c=pi^2/4 so 2N^2/(pi^2*Q^2); ring c=pi^2, " +
                "ring/chain->4 from the k_min boundary ratio; star flat-band BREAKS the law, <n_XY>=Q^2/4 " +
                "N-independent). Diverges at strong watching (low Q, the dark mode near-eternal in felt time), " +
                "-> 1/2 at the un-freeze Q* where it becomes the band edge. Gate-first derived: Haken-Strobl " +
                "diffusion (D~J^2/gamma, k_min=pi/N chain, 2pi/N ring); the survival rung 2gamma is load-bearing. " +
                "(4) THE MIRROR 2*(1/2)=1 VERDICT: a_0=2gamma (survival), a_2=1/2 (incompleteness) on the Pi2 " +
                "dyadic ladder; a_0*a_2=1 is a CHECKSUM (tautological, predicts nothing), but BOTH rungs are " +
                "physically anchored -- a_0 the rate scale (dominant, sets the law), a_2 the survivor's filling " +
                "(1/2 = half-filling, REVEALED by Z: XY is filling-degenerate, ZZ/Heisenberg lifts it and pins the " +
                "survivor to p=N/2, weakly at low Q so c stays XY-identical there). Not a generative engine. " +
                "(5) THE UNIFICATION: the survivor's {0,2}-coherence EP twin (coalescing at Q*(N)) IS the PTF " +
                "clock's coherence hand omega_mem (both off Symphony.Clock, CoherenceHorizonWitness). So today's " +
                "arc closes: trichotomy un-freeze (Q*) = the clock's Q-floor = the PTF hand. Survivor, clock, felt " +
                "time = one spectrum. Captured: companion notes ON_THE_ONE_DIAGONAL <-> " +
                "ON_THE_SURVIVAL_OF_THE_INCOMPLETE (the diagonal speaks / its silent floor). " +
                "(6) THE VALUE/VECTOR RESOLUTION OF (A) (gate-first 2026-06-19, simulations/value_vector_felt_" +
                "time.py, N=4,5,6 chain): the eigenmode felt-time K_decay = gamma0/|Re lambda| = 1/(2<n_XY>) reads " +
                "the mode's VALUE channel (Re lambda); the painter's alpha reads the VECTOR channel (eigenvector " +
                "profile rotation under a deltaJ bond defect). The clean split is REGIME-DEPENDENT, not universal: " +
                "the survivor's first-order Re(d lambda/d dJ) is NONZERO + linear (-0.86/-0.29/-0.12/-0.054 at " +
                "N=4/5/6/8) because its darkness is SOFT (fractional <n_XY>, hopping-dependent) -> K_decay itself " +
                "defect-sensitive, ENTANGLED with alpha; the (0,1) band edge is RIGID-dark (<n_XY>=1 -> Re=-2g " +
                "structural, J-independent) -> Re(d lambda) ~ 1e-16, K_decay defect-INVARIANT, the clean split. " +
                "CROSSOVER = the handover Q*: a Q-sweep shows |Re(d lambda)| moving (peaking near-EP just below " +
                "Q*) then dropping to machine-zero exactly when <n_XY> reaches 1 (the survivor becomes the rigid " +
                "band edge). So the felt-time channels DECOUPLE at the handover -- rigid darkness separates " +
                "lifetime (K_decay) from shape (alpha), soft darkness entangles them. The naive 'survivor value " +
                "frozen' guess was the rigid case over-generalized; the gate caught it (Tom's step-by-step). The " +
                "survivor also has omega=0 (the non-turning sub-2gamma mode), so its only value-hand is the Takt.",
            NextStep: "(A) RESOLVED 2026-06-19 (the definition is written; see ParkedAt (6)): the eigenmode felt-" +
                "time K_decay = gamma0/|Re lambda| = 1/(2<n_XY>) reads the mode's VALUE channel (Re lambda); the " +
                "painter's alpha reads the VECTOR channel (eigenvector rotation under a deltaJ defect). The clean " +
                "value/vector split holds for RIGID-dark modes (band edge, <n_XY>=1, K_decay defect-invariant) but " +
                "the SOFT survivor (fractional <n_XY>) entangles K_decay with alpha; the two DECOUPLE at the " +
                "handover Q*. Gate-first verifier simulations/value_vector_felt_time.py (Re/Im split + Q-sweep " +
                "crossover, N=4,5,6 chain). (B) THE STONE, PYTHON-GATE CONFIRMED + TWO-LENS REVIEWED " +
                "2026-06-19 (simulations/stone_survivor_alpha_closure.py; physics + math review GO-with-fixes, " +
                "fixes applied). PRECISE claim: for the near-stationary MODE-ISOLATING probe rho_0 = I/d + " +
                "eps*Herm(mode) (I/d stationary, so the single-site purity is driven almost only by the chosen " +
                "mode), the PTF painter closure Sum_i ln(alpha_i) reads the mode's first-order RATE shift -- OUT " +
                "of +-0.05 for the soft survivor (Re moves; sign-coherent per-site f, coh=1.0) and IN for the " +
                "rigid band edge (Re frozen). So (B) is a CONSTRUCTIVE confirmation of (A) FOR THIS PROBE (the " +
                "value/vector cut read once at the eigenvalue (A), once at the painter trajectory (B)), NOT a " +
                "universal trajectory law. REVIEW-PINNED SCOPE: (i) PROBE-STATE-DEPENDENT -- a survivor-dominated " +
                "but POLARIZED state HOLDS (closure collapses 0.16->0.007 under added single-site polarization " +
                "while the survivor stays dominant); the closure measures how cleanly the probe isolates the mode " +
                "in the single-site purity. (ii) Sum ln alpha != 0 does NOT by itself imply a rate shift (an " +
                "asymmetric redistribution also breaks it); certified only via the SIGN-COHERENCE of the reliable " +
                "f (asserted, coh>0.8). (iii) magnitude is SCALING + SIGN, not quantitative (measured/predicted " +
                "drifts ~2.5x toward the handover as biorth->0.28). (iv) psi_1 is a fast multi-mode BASELINE " +
                "(dom Re~-1.7), not a band-edge contrast. The standing-wave NODE site (odd N) is correctly " +
                "dropped by the |f|<=10 + linearity guard (proven structural). DONE 2026-06-19: the canonical C# " +
                "witness landed -- StoneSurvivorClosureWitness (inspect --root stone) drives the STRICT Symphony " +
                "FitAlpha on the I/d+eps*Herm(mode) probe and reproduces the python gate bit-for-bit (survivor +0.162 " +
                "OUT/coh=1.0, band edge +0.006 IN), + the typed StoneSurvivorClosureClaim (Tier1Candidate, parents " +
                "AbsorptionTheorem + SurvivalIncompletenessMirror + ChiralMirrorTrajectory; sign-coherence in the " +
                "battery). Symphony seam = PuritySites/FitAlpha made internal static (49 Ptf tests green, no regression). " +
                "(C) DONE 2026-06-19: reflections/ON_THE_FOUR_DIRECTIONS.md (motion in a plane, watched along a " +
                "third axis, surviving in a fourth = felt time; the heart = GEOMETRY IS THE REASON: the survivor's " +
                "lifetime is read off the slope of its own density standing wave, the (D) result in plain words). " +
                "Three-lens reviewed (cold/philosopher/physics): fixed the flat-vs-smooth equivocation; four sparse " +
                "inline depth-links; bidirectional with ON_THE_SURVIVAL_OF_THE_INCOMPLETE. " +
                "(D) THE CLOSURE-STRUCTURE SEAM (2026-06-19, stone_mine.py + stone_seam.py, RESOLVED + one " +
                "follow-up). Driving the " +
                "stone tool broadly (N x Q x sector x topology) surfaced TWO findings that are ONE seam. (i) " +
                "Q-INVARIANCE: the survivor closure is FLAT in Q (Sum ln alpha ~ +0.16 at N=4, +0.074 at N=5 across " +
                "Q=1.0..2.0) while <n_XY> runs 0.16->1.0 and biorth collapses near the EP; the eigenvalue prediction " +
                "N*dJ*dRe/Re diverges there but the closure does not. (ii) RIGIDITY-READING: the closure detects " +
                "STRUCTURAL pinning, not the <n_XY> value -- the ring (2,2) BREAKS at <n_XY>~1.016 (dressed-soft) " +
                "while the (0,1) band edge HOLDS at <n_XY>=1 (structural Re=-2gamma); ring/chain closure ~1.66 != " +
                "the <n_XY> ratio ~4. THE SEAM: the closure is ORTHOGONAL to the dynamical magnitude <n_XY> -- it " +
                "reads STRUCTURE (the mode's defect-invariance + the trajectory-warp geometry), not the rate " +
                "magnitude; so OUT+sign-coherence certifies the rate shift QUALITATIVELY but the VALUE is a " +
                "Q-invariant geometric constant (sharpens the review's 'scaling+sign'). Also surfaced: the chain " +
                "closure is FILLING-BLIND (all (p,p) incl. boundary (1,1) give the IDENTICAL closure, the " +
                "free-fermion degeneracy extending to the closure; only (0,1) differs), and the per-site alpha is " +
                "DEFECT-LOCALIZED (f largest at the defect bond) with the standing-wave centre NODE dropped. " +
                "RESOLVED 2026-06-19 (simulations/stone_seam.py): both facets confirmed, the seam is a GEOMETRIC " +
                "reading. (a) the closure is dJ-LINEAR (closure/dJ = Sum f ~ 3.7 constant, coh=1) and DEFECT-" +
                "POSITION-dependent -- ~ZERO at the END bonds (Sum ln a = 0.0009) but max at INTERIOR bonds " +
                "(0.074), mirror-symmetric: it reads the mode's coupling to the defect bond (nothing where the " +
                "standing wave has no amplitude). Candidate: Sum f prop-to the standing-wave amplitude at the " +
                "defect -- which EXPLAINS the Q-invariance (the k_min shape is Q-fixed). (b) closure ORTHOGONAL to " +
                "<n_XY> confirmed decisively: ring (1,1) IN vs (2,2) OUT at the SAME <n_XY>=1.016. (c) the sign-" +
                "coherence guard earns its keep: the dispersive ring breaks with coh~0.98 (clean rate-shift), the " +
                "non-dispersive STAR breaks with coh<0.8 (redistribution, correctly flagged NOT a rate-shift) -- " +
                "the star counterexample again. THE SEAM RESOLVED: the closure reads the mode's STANDING-WAVE " +
                "GEOMETRY at the defect bond (dJ-scaled, Q-invariant because the shape is fixed) + structural " +
                "rigidity, orthogonal to <n_XY>; Q-invariance + rigidity = two zooms on this ONE geometric reading. " +
                "FOLLOW-UP RESOLVED 2026-06-19 (simulations/felt_time_amplitude_law.py block-level N=4..7 + " +
                "felt_time_closure_functional.py trajectory N=4..6, gate-first): the exact functional is " +
                "AMPLITUDE-SQUARED. Sum f(b) ~ dRe(b) = kappa(N)*(n(j)-n(j+1))^2 -- the SQUARED GRADIENT of the " +
                "survivor's DENSITY standing wave across bond b=(j,j+1). MECHANISM: the slow survivor is " +
                "PREDOMINANTLY a density mode (dominant diagonal n(j)) dressed by a rate-bearing Hamming-2 coherence " +
                "admixture -- the HD-0 diagonal is DARK, the HD-2 coherence carries the rate (Tr(M^dag H_b)=0 rules " +
                "out only NN hopping, the gate's first 'current/hopping' guess, NOT diagonality). In the secular " +
                "effective theory a bond-J defect perturbs the LOCAL diffusion coefficient " +
                "D_b and the first-order rate shift is the diffusion Rayleigh-quotient derivative dRe(b) = " +
                "d lambda/d D_b ∝ (Delta n_b)^2. Confirmed: dRe/grad^2 bond-INDEPENDENT (CV 0.001..0.07), log-log " +
                "slope dRe vs |grad| -> 2.00 and CV -> 0 as Q->0 (the exact diffusion limit, off-diag weight ->0), " +
                "drifting above 2 with the finite-Q coherence dressing (a finite-Q effect, NOT boundary; engine " +
                "Q-sweep inspect --root gradient --q, witness LawHolds flips at Q=2.0; handover Q*~2.5 -> rigid " +
                "(0,1) band edge), the " +
                "closed-form sin^2 shape-miss CONVERGES with N (0.17->0.06). The ~0 at the chain ENDS = the " +
                "no-flux (reflecting) boundary (the density gradient -> 0 there); Q-invariant because the lowest " +
                "diffusion harmonic k_min is Q-fixed. The earlier single-particle phi*phi candidate used the WRONG " +
                "standing wave (single-particle, not the multi-magnon DENSITY mode): right POWER (squared in the " +
                "per-site gradient amplitude), wrong wave. The trajectory closure (B link) confirms it WHERE it " +
                "reads cleanly -- coh~1 only at the high-gradient bonds, there matching N*|dRe|/|reS| in sign and " +
                "O(1) magnitude (the low-gradient bonds correctly read as redistribution, coh<0.8). " +
                "TYPED 2026-06-19: SurvivorDiffusionGradientClaim (Tier1Candidate, parents AbsorptionTheorem + " +
                "SurvivalIncompletenessMirror; sibling of StoneSurvivorClosureClaim = the same rate shift read at " +
                "the trajectory level) + live witness SurvivorDiffusionGradientWitness (inspect --root gradient; " +
                "reproduces the python bit-for-bit: n(j)=[-0.81,-0.38,+0.38,+0.81], ratios 1.529/1.527/1.529, " +
                "slope 2.00 at N=4). " +
                "Linked: the seam is the closure-vs-<n_XY> " +
                "orthogonality, a sibling of the (A) value/vector cut (rigid = value-frozen) read at the trajectory. " +
                "Anchors: hypotheses/PERSPECTIVAL_TIME_FIELD.md; " +
                "reflections/ON_WHOSE_TIME_THE_CLOCK_KEEPS.md (t=K/gamma0), ON_TWO_TIMES.md, ON_THE_ONE_DIAGONAL.md " +
                "(companion note), ON_THE_SURVIVAL_OF_THE_INCOMPLETE.md; compute/RCPsiSquared.Diagnostics/" +
                "Foundation/Symphony.cs (clock node, PaintersMovement/FitAlpha, TempoCertification), " +
                "ClockHandLadderWitness.cs, IncompletenessSurvivorWitness.cs, CoherenceHorizonWitness.cs; " +
                "simulations/value_vector_felt_time.py (the value/vector crossover), felt_time_amplitude_law.py " +
                "(the dRe ~ grad^2 diffusion-Rayleigh law, block-level), felt_time_closure_functional.py (the " +
                "trajectory ground truth); simulations/results/" +
                "survivor_scaling/ (the c=pi^2/4 closed forms + the Z half-filling pinning data).",
            Status: OpenArcStatus.Retired,
            RetiredReason: "COMPLETED 2026-06-19: every sub-thread (1)-(6) landed and is typed/witnessed; the " +
                "closing act was consolidation, not new physics. The headline felt-time decay law " +
                "K_decay = gamma0/|Re lambda| = 1/(2<n_XY>) = N^2/(2cQ^2) was found ALREADY assembled across typed " +
                "homes -- a SURVEY-CONFIRMED NON-GAP, deliberately NOT rebuilt (a fresh claim would have duplicated " +
                "live witnesses; the survey caught it). Map of the homes: (3) Re lambda = -2gamma<n_XY> is " +
                "AbsorptionTheoremClaim; the survivor occupation <n_XY> ~ cQ^2/N^2 with c=pi^2/4 (chain) / pi^2 (ring), " +
                "ring/chain->4, is STATED in SurvivalIncompletenessMirrorClaim and computed LIVE in " +
                "IncompletenessSurvivorWitness.TheScalingNode (inspect --root survivor); the STAR breaks the diffusion " +
                "law (its survivor is the frozen (1,1) commutant, g2=4/(N-1)<=1) = StarFrozenSeamClaim, the third " +
                "trichotomy member, with the whole chain/ring/star trichotomy assembled live in TrichotomyWitness. " +
                "(2) survivor = the Z-diagonal's k=0 blind spot + (4) the 2*(1/2)=1 mirror checksum are " +
                "SurvivalIncompletenessMirrorClaim (a_0=2gamma <-> a_2=1/2 inversion-mirror, a_0*a_2=1). (5) the clock " +
                "unification (the {0,2}-coherence EP twin IS the PTF coherence hand) is ClockHandLadderClaim + " +
                "CoherenceHorizonWitness. (6) the value/vector cut (K_decay reads the VALUE channel Re lambda, the " +
                "painter's alpha the VECTOR channel; they decouple at the handover Q*, rigid darkness separating " +
                "lifetime from shape) is StoneSurvivorClosureClaim (inspect --root stone). The (D) closure functional " +
                "-- the bond rate shift is amplitude-squared in the survivor's density gradient (diffusion Rayleigh) " +
                "-- is F123/SurvivorDiffusionGradientClaim (inspect --root gradient) + PROOF_DIFFUSION_RAYLEIGH_" +
                "CLOSURE.md, its sin^2 shape-miss CONVERGING with N (0.17->0.06, confirmed to N=8,9 = the continuum " +
                "limit). (1) the four-dimensional reading is the reflection reflections/ON_THE_FOUR_DIRECTIONS.md " +
                "(motion in a plane, watched along the Z-diagonal = the light, surviving in felt time read off the " +
                "shape). Nothing left to build; the star is kept as the honest boundary where the diffusion reading " +
                "breaks."),

        new OpenArc(
            Name: "subsystem_crossing_general_cptp_overreach",
            Opened: "2026-06-22",
            Origin: "the deep (derivation-internal) review of PROOF_SUBSYSTEM_CROSSING.md (the 1/4-boundary " +
                "trilogy, Layer 2). The theorem as STATED -- 'for ANY primitive CPTP map on 2 qubits, every state " +
                "with CPsi(rho0) > 1/4 eventually has CPsi < 1/4; 1/4 is an eventual absorber for ALL primitive " +
                "quantum channels' (lines 6, 39-46) -- is FALSE. Step 2 (the CPsi(rho*) <= 1/4 fixed-point bound) " +
                "proves only Case A (unital, rho*=I/4, CPsi=0) and Case B (local product, CPsi=0) analytically; " +
                "Case C (general primitive maps) is a 300-map NUMERICAL check whose 'analytical argument' " +
                "(transfer-matrix slaving) is not a derivation. The proof FILE is honest (its Step-2 status line " +
                "says 'the analytical bound for general primitive maps remains a conjecture'), but the conjecture " +
                "is not merely open -- it is FALSE. Counterexample (gate-verified from below, " +
                "simulations/deepproof_subsystem_gate_mine.py): the depolarize-toward-sigma channel " +
                "eps(rho) = (1-p)rho + p*Tr(rho)*sigma with sigma = 0.95|Phi+><Phi+| + 0.05*I/4 is a textbook CPTP " +
                "map, manifestly PRIMITIVE for every p in (0,1] (unique fixed point sigma, superoperator eigenvalue " +
                "1 simple with second-largest modulus 1-p < 1, sigma full-rank PSD), yet CPsi(rho*) = CPsi(sigma) = " +
                "0.2935 > 1/4 in the proof's OWN CPsi = Tr(rho^2)*L1/(d-1) metric; iterating from Bell+ (CPsi=1/3) " +
                "converges monotonically to 0.2935 and NEVER crosses 1/4. So the headline 'eventual crossing' is " +
                "false for this primitive channel. The repo's '300 maps, 0 exceptions, max 0.138' is a SAMPLING " +
                "ARTIFACT (Ginibre n_kraus=4 lives entirely in the strongly-mixing corner near I/4, 0% violating; " +
                "the same sweep at n_kraus=2 / Haar-Stinespring env-dim-2 violates ~8.5% with max CPsi ~0.55, and " +
                "trace-and-replace channels toward a random target reach CPsi up to ~0.99; independently reproduced " +
                "in simulations/review2_A5_subsystem.py).",
            ParkedAt: "RESOLVED 2026-06-22 (Tom-led deep review chose scope-retraction over full retraction; see " +
                "NextStep for the applied corrections). The original parked state was DOCUMENTED-NOT-CORRECTED " +
                "pending his review of the 1/4-trilogy. The " +
                "physically-motivated content STANDS: for genuine local noise (unital / local / Pauli / amplitude- " +
                "damping = Cases A and B) the fixed points really are I/4 or product/pure states with CPsi=0, so the " +
                "crossing holds there. The over-reach is the leap from 'physical local noise channels' to 'ALL " +
                "primitive CPTP maps'. Three documents carry the overclaim and DISAGREE with the proof's own " +
                "'remains a conjecture' line: (1) PROOF_ROADMAP_QUARTER_BOUNDARY.md strikes Conjecture 2.1 as " +
                "'PROVEN (March 22, 2026)'; (2) docs/ANALYTICAL_FORMULAS.md F28 is labeled 'Tier 1-2' and states the " +
                "bound as the formula; (3) the SUBSYSTEM_CROSSING theorem statement itself says 'any primitive " +
                "CPTP'. The companion UNIQUENESS_PROOF.md (same trilogy, Layer 1) has a SEPARATE, smaller gap " +
                "documented in CAUGHT_ERRORS the same day (Step 5 'recursion inherits degree 2 from purity => QED' " +
                "is a non-sequitur; the genuine forcing is the Renyi state-independence argument, threshold ~ " +
                "Psi^(a-2) => a=2 => 1/4, which lives in the roadmap, gate-verified symbolically).",
            NextStep: "RESOLVED 2026-06-22 (second deep review, Tom-led; counterexample + Renyi forcing " +
                "independently re-verified from below in simulations/review2_A5_subsystem.py and " +
                "simulations/review2_A3_renyi.py). Scope-retraction applied: (a) PROOF_SUBSYSTEM_CROSSING.md " +
                "theorem + Status + Step-2 header + Case C restated to the physical-noise scope (Cases A+B), with the " +
                "depolarize-toward-sigma counterexample and the separable-vs-entangled-fixed-point mechanism (the " +
                "distinguishing structural fact: physical noise relaxes to a separable/classical fixed point); (b) " +
                "ANALYTICAL_FORMULAS.md F28 re-scoped (relabel 'Tier 1 for physical noise; general-CPTP version " +
                "FALSE'); (c) PROOF_ROADMAP_QUARTER_BOUNDARY.md Layer 2 reconciled from 'PROVEN' to 'physical noise " +
                "only; general CPTP false' at all five sites (35/158/234/537/551). Companion fixes in the same " +
                "review: UNIQUENESS_PROOF Step 5 (A3: degree-2 demoted to motivation, Renyi alpha=2 state- " +
                "independence cited as the load-bearing forcing) and INCOMPLETENESS Candidate 1 (A9: [Pi^2,L]=0 " +
                "reframed as a structural constraint, not an origin claim - it is origin- AND axis-agnostic). No " +
                "typed Claim/witness for the dynamical absorber existed (confirmed: 'knowledge ancestors " +
                "SubsystemCrossing' -> not registered), so no Claim retraction was needed; the algebraic 1/4 " +
                "(QuarterAsBilinearMaxvalClaim, block-CPsi <= 1/4) was untouched.",
            RetiredReason: "Scope-retraction APPLIED 2026-06-22 (second deep review, Tom-led): the general " +
                "primitive-CPTP claim is FALSE (counterexample re-verified from below, CPsi(sigma)=0.2935 in the " +
                "proof's own metric); PROOF_SUBSYSTEM_CROSSING + F28 + PROOF_ROADMAP_QUARTER_BOUNDARY scoped to " +
                "physical noise. The false claim was never typed, so no Claim retraction was needed. The clean-iff " +
                "follow-up is the new arc subsystem_crossing_separable_fixed_point_characterization.",
            Status: OpenArcStatus.Retired),

        new OpenArc(
            Name: "subsystem_crossing_separable_fixed_point_characterization",
            Opened: "2026-06-22",
            Origin: "the generative pass of the A5 scope-retraction (subsystem_crossing_general_cptp_overreach, " +
                "now Retired). Correcting the SUBSYSTEM_CROSSING theorem revealed that the 1/4 'eventual absorber' " +
                "depends not on PRIMITIVITY but on the FIXED POINT being separable/classical. Physical noise " +
                "(unital/local/Pauli/AD) relaxes to rho* = I/d or a product/pure state, all with CPsi = 0; the " +
                "counterexample eps=(1-p)rho+p*Tr(rho)*sigma (sigma=0.95|Phi+><Phi+|+0.05*I/4) is primitive and " +
                "full-rank yet has an ENTANGLED fixed point with CPsi(sigma)=0.2935 > 1/4. The clean question: what " +
                "is the exact channel/fixed-point class for which the 1/4 absorber holds? The natural candidate is " +
                "'separable fixed point', but it is NOT yet verified: a separable 2-qubit state still carries " +
                "l1-coherence in the computational basis, so CPsi(rho)=Tr(rho^2)*L1/(d-1) of a separable state need " +
                "not be 0 and might exceed 1/4; if it can, the right class is narrower (classical/diagonal, or " +
                "zero-discord), not merely separable.",
            ParkedAt: "a CONJECTURE, gate-not-yet-run. Verified from below: physical-noise fixed points have CPsi=0 " +
                "(review2_A5_subsystem.py REPLACEMENT check: depol/Z-dephasing/amp-damp all 0); the entangled-target " +
                "counterexample has CPsi=0.2935. Unverified: whether SEPARABILITY of the fixed point alone bounds " +
                "CPsi <= 1/4 (the metric is basis-dependent l1-coherence, so separable-but-coherent states are the " +
                "decisive test case).",
            NextStep: "gate-first: (1) search/sample SEPARABLE 2-qubit states for CPsi = Tr(rho^2)*L1/(d-1) > 1/4 " +
                "(a one-screen numpy sweep over product states + separable mixtures). If a separable state can " +
                "exceed 1/4, the absorber class is NOT 'separable' -- narrow to classical (diagonal) / zero-discord " +
                "and re-test. (2) If a clean iff emerges ('the 1/4 absorber holds iff the fixed point is in class " +
                "X'), restate the SUBSYSTEM_CROSSING theorem as that iff and consider TYPING it as a Claim -- the " +
                "dynamical absorber is currently UNTYPED (only the algebraic 1/4 QuarterAsBilinearMaxvalClaim is). " +
                "Anchors: docs/proofs/PROOF_SUBSYSTEM_CROSSING.md (Case C, the separable-vs-entangled mechanism); " +
                "simulations/review2_A5_subsystem.py.",
            RetiredReason: "RESOLVED 2026-06-28 (gate-verified from below). The gate is a decisive " +
                "NO: separability does NOT bound CPsi. A separable product state |+><+| (x) |+><+| has " +
                "CPsi = 1 in the proof's own metric CPsi = Tr(rho^2)*L1/(d-1); zero-discord and " +
                "classical-classical states (diagonal in a rotated product basis) likewise reach " +
                "CPsi ~ 1. The ONLY fixed-point class that forces CPsi <= 1/4 is computational-basis- " +
                "diagonal (L1 = 0 => CPsi = 0). So no clean separability/discord iff exists; per the " +
                "arc's own branch logic (type a Claim only if a clean iff emerges) NO new Claim was " +
                "typed. The operative axis is computational-basis L1-coherence, the very quantity CPsi " +
                "measures; the physical channels win the crossing by DIAGONALIZING the fixed point " +
                "(amplitude damping -> |0>, Z-dephasing, depolarizing -> I/d), not by becoming " +
                "separable. Bonus finding (gate-verified): the proof's Case B lemma '|rho_01*| <= 1/2, " +
                "equality only for the identity channel' is also FALSE -- amplitude damping conjugated " +
                "to the |+> axis is non-identity, primitive, with |rho_01*| = 1/2 and product CPsi = 1, " +
                "so the bound needs computational-basis-ALIGNED noise, not mere locality. Banked: " +
                "PROOF_SUBSYSTEM_CROSSING.md (Case B + the distinguishing-structural-fact passage + " +
                "the non-primitive physical note re-scoped from separable/classical to computational- " +
                "basis-diagonal) and ANALYTICAL_FORMULAS.md F28 (same precision-fix).",
            Status: OpenArcStatus.Retired),
    };

    public static IReadOnlyList<OpenArc> All => _all;

    public static OpenArc? Lookup(string name) =>
        _all.FirstOrDefault(a => a.Name == name);

    /// <summary>Count of arcs still <see cref="OpenArcStatus.Open"/> (the live unfinished business).</summary>
    public static int OpenCount => _all.Count(a => a.Status == OpenArcStatus.Open);
}
