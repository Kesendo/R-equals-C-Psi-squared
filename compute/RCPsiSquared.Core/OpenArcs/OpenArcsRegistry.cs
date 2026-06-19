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
                "the C# Symphony painters witness + its tests; simulations/_ptf_symphony_crossval.py (the Bell+ " +
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
            ParkedAt: "the closure law lives as a live witness (Symphony painters movement: alpha per site with the two-deltaJ reliability guard, closure -0.0444 IN window at canonical N=5, chiral mirror replayed live at 1e-15, alpha = Python twin to 1e-3); learned on the way: the protocol needs the BONDING state class (Bell+/localized states break the rescaling and the guard refuses, in both languages); C# golden-section also found a true global alpha minimum where scipy's Brent traps 920x worse (severed-bond case, harness simulations/_ptf_symphony_crossval.py)",
            NextStep: "REQUIRED, escalated 2026-06-12: the scipy-Brent trap corrupts CANONICAL N=5 edge-bond letters (arbiter: brute-force landscapes match C# exactly, Python f off by sign and factor); backport the global-minimum fit (multi-start or grid-seed) to framework ptf.py before any further Python-side painter quantitative work; then retire this arc",
            Status: OpenArcStatus.Retired,
            RetiredReason: "RESOLVED 2026-06-15 (fix in simulations/framework/workflows/ptf.py _alpha_fit_one_site; " +
                "arbiter simulations/_ptf_symphony_crossval.py; 14/14 framework PTF tests green; 3-agent survey). " +
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
            Status: OpenArcStatus.Open),

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
            Status: OpenArcStatus.Open),

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
            NextStep: "FROST_CIRCLE seam FIRST STONE laid (2026-06-12, merged): the carbon coherent<->incoherent threshold IS our Q-floor, verified bit-for-bit and given our label - the Coherence Horizon Q*(N). Live witness CoherenceHorizonWitness / inspect --root horizon bisects Symphony.Clock.Omega and computes Q*(2,3,4,5)=1 / sqrt2 / 1.8785 / 2.3722, equal to carbon's sqrt2/1.879/2.372 (FROST_CIRCLE) under J<->|beta|; ANALYTICAL_FORMULAS F2b corollary carries it; typed claim CoherenceHorizonClaim (Tier1Candidate, parents ClockHandLadder + F2b) breadcrumbs the witness into the graph. N=2 (Q*=1) is the EP base rung the polyene layer can't reach; the N=2,3 band-edge coincidence (Q*=2cos(pi/(N+1)) ONLY there) explains carbon's 'sqrt2 exact at N=3, rest awaiting a clean form' puzzle. CLOSED FORM ATTACKED 2026-06-12, then CORRECTED 2026-06-13 (physics-first review of the erasure-ladder spec, via phase rigidity): the earlier 'freezing mode BIFURCATES SECTOR at N=4' was WRONG, an argmax-Re / Im-tracking artifact - the band edge and the freezer are Re-degenerate at -2gamma (both <n_diff>=1, Absorption Theorem), so Im-tracking picked the wrong one. THE TRUTH: the mode that COALESCES at Q*(N) is the {0,2}-coherence (population/antisymmetric block, n_diff hist {0:1/2,2:1/2}) at ALL N=2..5, a genuine square-root EP (phase rigidity r->0; r at Q* = 0.000/0.015/0.026 at N=3/4/5, Im prop sqrt(Q-Q*)). NO bifurcation. The band edge 2cos(pi/(N+1)) (the half-lit mode the old reading mislabeled the freezer at N>=4) is the co-located gamma-protected SURVIVOR = Uhr 1, the |vac><psi_k| coherence hand; so Q*(N) is AT ONCE a {0,2} EP (Uhr 2, the erasure point, which CLIMBS the ladder) and a band-edge crossing (Uhr 1 survives the handover). The handshake erasure-ladder answer is 'both', not 'EP-or-crossing'. Corrected 2026-06-13: ANALYTICAL_FORMULAS F2b corollary Q*(N) para + typed CoherenceHorizonClaim (4 strings) + KnowledgeRegistryFactory comment; the witness was already mechanism-agnostic. THE EP-VERDICT WITNESS LANDED (2026-06-13, master 850e2ff): CoherenceHorizonWitness recomputes the {0,2}-EP verdict live via a new Core PhaseRigidity primitive (Petermann form 1/||R^-1 row||, degeneracy-robust; the matched-overlap fakes r->0 at the N=3 Q*=band-edge=sqrt2 degeneracy), inspect --root horizon. THE REAL PRIZE RESOLVED (2026-06-13, Approach A, simulations/coherence_horizon_se_block.py): Q*(N) reduces 4^N -> N^2 - the coalescing mode is single-excitation, so Q*(N) is the EP of the single-excitation (Haken-Strobl) Liouvillian (validated as an exact sub-spectrum of the full L). At N=2,3 the coalescing pair are the roots of lambda^2 + 4g*lambda + c*J^2 = 0 with c CONSTANT (sum=-4g, product=c*J^2 are gamma-independent identities; c=4,2), so Q* = 2/sqrt(c) = 1, sqrt2 EXACTLY - the structural form of the 2cos(pi/(N+1)) low-N accident (the whole clean 2x2, not just the value, exists only at N=2,3). At N>=4 the pair is collectively dressed (trace departs from -4g), no clean 2x2, the exact condition transcendental: a diffusive long-wavelength critical damping, Q*(N) ~ 0.59 N (canonical Q*(4)=1.87874, Q*(5)=2.37367, superseding the grid). SLOPE RESOLVED 2026-06-15 (asymptotic slope EXACTLY 2/pi, derived: PROOF_COHERENCE_HORIZON_SLOPE.md + coherence_horizon_slope.py, adversarial review GO. The slow mode is the population coupled to the FULL ladder of coherence ranges r, geometric decay mu^r; resummation DOUBLES the telegrapher coefficients to lambda^2+8g*lambda+4J^2*q^2, EP at g*=Jq/2 so Q*=2/q_min -> 2N/pi; the nearest-neighbour truncation gives the WRONG sqrt2/pi=0.450. The 8g coefficient confirmed against L_se by the robust overdamped q^2-constancy (the CV test); the -2Re/g->8 sweep is corroborative-only (the EP metric = the linear coefficient, running 4->8, but mode-selection sensitive; Round 2 honesty fix 2026-06-15). Ring sibling exactly half: q_min=2pi/N -> slope 1/pi). STILL OPEN: the half-filling double-excitation V-Effect seam co-located at even N. STILL OPEN too: the full FROST_CIRCLE reflection in its own genre. Also still: the off-gap argument (general N>=3) PARTLY RESOLVED 2026-06-15 (simulations/offgap_band_edge.py, gate-first, self-validating): (i) the Q-FLOOR = the Coherence Horizon Q*(N) -- the numerical floor brackets Q*(4)=1.879 / Q*(5)=2.374 exactly (below it the overdamped {0,2}-diffusion mode is slower than the -2g band edge AND non-oscillating, so it takes the gap and the clock freezes), so regime (i) IS thread (a), DERIVED ->2N/pi (PROOF_COHERENCE_HORIZON_SLOPE); (ii) max|Im| at the exact gap rate 2g = the band edge E1 CONFIRMED N=3,4,5 and REFRAMED via the Absorption Theorem (Re=-2g<n_XY>, so the exact-2g modes are the n_XY=1 subspace: 22/32/50 modes at N=3/4/5, MORE than the 4N vacuum/full ladders, which are a subset achieving E1). RESOLVED 2026-06-16 (the Tier1Candidate gate, PROOF_CHAIN_GAP_DOMINANCE.md + verifier chain_gap_dominance.py): max|Im| in the n_XY=1 subspace = E1. Although L_H leaks n_XY=1->3 (interacting), ON the exact-(-2g) subspace L_D=-2g is SCALAR so L=L_H-2g acts FREELY; via Jordan-Wigner the exact-(-2g) modes are c_k^(dag).f(N_tot) (a single fermion dressed by a function of total number N_tot), oscillating at the single-particle energies +-E_k<=E1, with E1 reached by the (0,1) ladder. They SPAN the subspace for N>=4 (dim 32/50/72 at N=4/5/6); N=3 is SPECIAL, +4 extra (n,n) {0,2} sqrt-EP modes at sqrt(E1^2-(2g)^2)<E1 (parked in the n3_special_cases arc). So ClockHandLadderClaim + TopologyBandEdgeClaim GRADUATED to Tier1Derived (CoherenceHorizonClaim stays Tier1Candidate for its own V-Effect-seam reason). Scope: the chain (JW is 1D); other topologies handled by the star-no-horizon + structural-ceiling work. RING DIHEDRAL-LOCK RESOLVED 2026-06-17 (PROOF_RING_GAP_DOMINANCE.md + ring_gap_dominance.py, gate-first N=3..6): max|Im| over the exact-(-2g) ring modes = 2J = J*rho (the periodic band top = adjacency radius = the k=0 uniform single-excitation mode fixed by C_N = the dihedral lock; reached via the (0,1) sector, general-N), EXCEPT N=4 where the half-filling (2,2) {0,2} sqrt-EP reaches sqrt((2sqrt2)^2-(2g)^2)->2sqrt2 J > 2J (the LONE exception, the same (2,2) sector as K_4's ceiling; the ring analogue of the chain's N=3 special but ABOVE the band top, not below). Open extension: the all-N ring free-fermion completeness (parity-split span). The Im_max=Delta_E_max sibling is DONE. (the 0.99750 fine-structure pull listed here in an earlier session was an unidentifiable dangling reference: neither a 3-agent survey nor Tom could place it; dropped 2026-06-15). BENZENE RING CROSSOVER ANSWERED + V-EFFECT SEAM PROBED (2026-06-13): benzene's own ring Q* is the ring single-excitation {0,2}-EP = 1.609 (Uhr 2; transcendental like N>=4; the SE block embeds bit-exact in the full 4^6 L), more coherence-robust than every open polyene. Benzene (even N, HALF-FILLED) is a concrete instance of the V-Effect seam: the full-L mode that overtakes the band-edge beat below the crossover is a DOUBLE-excitation coherence (filling sector (2,2)/(4,4)), so the full-L handover (~1.95) SPLITS from the clean SE-EP Uhr 2 (1.609) - the Uhr1/Uhr2 co-location that holds for the open chains (AT-pinned at Re=-2gamma) breaks on the ring. Verifier simulations/carbon/benzene_two_clocks.py (4 asserts green); FROST_CIRCLE_AS_THE_CLOCK_FACE.md carries benzene Q* + the V-Effect split + the F2b two-clocks backlink, frost_circle_as_clock.py de-staled + self-checking, benzene<->FROST cross-links wired (LIOUVILLIAN_PALINDROME + F98). RING-SPECIFIC (2026-06-13, open-chain N=6 control, benzene_two_clocks.py section 5): the open even-N chain does NOT show the clean double-excitation seam - its full-L overtaker spreads across ALL fillings (single-excitation included) and hands over at its own SE-EP (co-located, ladder holds), unlike benzene's pure (2,2)/(4,4) mode handing over ABOVE the SE-EP. So the V-Effect double-excitation split is a feature of the CLOSED RING at half-filling (the C=0.5 boundary), not even N alone. AROMATICITY REFUTED (2026-06-13, aromatic_ring_v_effect.py, sector-projected Liouvillians): the seam is RING-UNIVERSAL (C4/C6/C8 all hand over to a frozen double-excitation mode at strong dephasing = the C=0.5 half-filling boundary, a sibling of the incompleteness V-Effect docs/HIERARCHY_OF_INCOMPLETENESS.md); the 4n anti-aromatics C4 and C8 do NOT group (C4 a small-ring anomaly whose seam dominates even at weak dephasing), so Hueckel 4n+2 is NOT the discriminant. CONVERGENCE: this is the same mode birth_canal_horizon_junction found independently (the {0,2}-coherence in the 2-excitation block the odd |vac><psi| Uhr-1 crosses to at N>=6). RECONCILED 2026-06-13 (a CHAIN_GAP_SECTOR_DIAGNOSTIC mismatch cleared, false alarm): that diagnostic's half-filling (3,3) slow mode is HEISENBERG (XXX, with ZZ); the carbon work is XY (free fermions, the Hueckel hopping model, no ZZ). The SAME benzene_two_clocks builder reproduces CHAIN_GAP's -0.230 (3,3) bit-for-bit once the ZZ is added (_xy_vs_heisenberg_slowmode.py) - so no bug, the committed XY findings stand. THE THREAD-(a) GIFT: the survival law IS the Absorption Theorem across BOTH free and interacting - below the handover the survivor is a frozen dressed magnon-admixture (fractional <n_XY>) at the filling centre in either model; XY puts it at (2,2)/(4,4) <n_XY>=0.72, Heisenberg at dead-centre (3,3) <n_XY>=0.23 (ring N=6 Q=2); the ZZ retunes the darkness and sector, NOT the law. HANDOVER MECHANISM CHARACTERIZED 2026-06-14 (simulations/carbon/handover_q.py, self-validating; F2b corollary): the double-excitation (2,2) seam handover is a frozen LEVEL CROSSING (|Im|~1e-15) where the seam brightens to the F50 off-diagonal floor <n_XY>=1 (the (0,1) band edge, Re=-2g) - a different sector than the single-excitation SE-EP (a coalescence), growing ~linearly Q_h~0.29N (c_eff~12 flat, NOT saturating, faster than sqrt(N)); NOT co-located with the ring SE-EP (the curves cross near N~10, benzene's 2.0-vs-1.609 split is small-N). STILL OPEN: the dynamic C=0.5 question (is the half-occupied coherence always the longest-lived survivor over N/topology?)",
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
                "N>=4); the RING is a DISTINCT (2,2) free-fermion LEVEL CROSSING (|Im|~1e-15), growing ~linearly " +
                "Q_h~0.29N (c_eff~12 flat, ~4x chain darkness so ~half the slope), NOT saturating, faster than " +
                "sqrt(N), and NOT co-located with the ring SE-EP (values cross near N~10; benzene's 2.0-vs-1.609 " +
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
                "past C(N,p)~25000) or a closed-form account of the hump; a power-law fit degenerates even on the descending tail.",
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "handshake_decoder",
            Opened: "2026-06-12",
            Origin: "docs/superpowers/specs/2026-06-12-handshake-decoder-reading-grammar-design.md (Tom + Opus + Fable)",
            ParkedAt: "M2a + M2b LANDED: the FI(Q) witness (resolution = Q-factor, no EP peak, basis ordering) AND the DefectDecoder (alpha-profile inverted: location+strength read back, ambiguity-honest with runner-up reporting). Discoveries en route: the N=5 sign-location ambiguity (near-anti-collinear letters, ratio 1.5) DERIVED in structure - the K-partner selection rule <psi_N|V_b|psi_1>=0 per bond (one line from ChiralMirrorTrajectoryClaim: rank(dictionary)=N-2, the null space IS the K-partnership); channel factorization verified ~12% (seesaw = dominant k=2 channel, 14.9 vs 5.2/3.1); the Gram spectrum = the location-channel metric (eigs 0.0073/0.008/0.086/1 at N=5, confirmed on arbiter-corrected letters after the Python Brent-trap was exposed)",
            NextStep: "DONE 2026-06-13 (merged): the K-partner selection rule is a typed claim KPartnerSelectionRuleClaim (Tier1Derived, parent ChiralMirrorTrajectoryClaim) - the reading-grammar arc's first DERIVED result. <psi_N|V_b|psi_1>=0 exactly (two-line corollary: psi_N=K1 psi_1 + K1 V_b K1 = -V_b => the real element is its own negative), so rank(location dictionary)=N-2 and the DefectDecoder sign-location ambiguity IS the K-partnership null direction; live battery at N=4,5 + probe N=3..8 (_k_partner_selection_rule.py), in HANDSHAKE_GEOMETRY.md. CHANNEL PROFILES R_k LOOKED AT 2026-06-13 (two conjectures vetoed, parked): the f-dictionary factorizes as F[b,i]=sum_k c_k(b) R_k(i) over the unforbidden channels k=1..N-1 (k=N selection-rule-forbidden), but (1) the SlowModeMixing denominators (E_1-E_k) do NOT help - they are a per-channel COLUMN rescaling c_k = <psi_k|V_b|psi_1>/(E_1-E_k) that the least-squares R absorbs, identical residual to the bare coupling (the dissipative-dressing guess was wrong); (2) the unforbidden set is dimensionally SQUARE (N-1 bonds x N-1 channels) so the factorization is trivially exact and does NOT test the structure - the informative cut is location channels alone (k=2..N-1, N-2 of them) capture 85%, the missing ~15% is one more channel (dimensionally the STRENGTH channel k=1, the diagonal <psi_1|V_b|psi_1>; f-profile = strength + location). The real R_k derivation needs the FIRST-PRINCIPLES per-site purity-sensitivity to each mode M_k (R_k(i) = response of site-i purity to mode k), NOT a fit - a deeper stone, parked. REMAINING: that first-principles R_k; the Gram location-metric in HANDSHAKE_GEOMETRY.md; M2c the read-cost law ~2/Q; then M3/M4. R_k PROGRESS BANKED 2026-06-19 (gate-first probes _handshake_rk_block.py + _handshake_carrier_compare.py + _handshake_rk_first_principles.py): REFRAME - the BondingMode carrier |psi_1> is pure single excitation, so <X_i>=<Y_i>=0 and the per-site purity P_a=1/2(1+<Z_a>^2) is PURE POPULATION dynamics in the N^2-dim (1,1) Liouvillian (single-particle Haken-Strobl); reproduces the full-4^N painter exactly for this carrier but costs N^2, so the f-dictionary runs to N=20+. CARRIER-INDEPENDENCE: pure |psi_1> vs the documented PairState (|vac>+|psi_1>)/sqrt2 give the SAME f-profile shape (corr 0.999, only the scale differs) - R_k is a property of the psi_1 mode-RESHAPING, not the observable; the PairState-vs-pure docstring mismatch is physically harmless for R_k. TWO MORE CONJECTURES REFUTED at high N (killing the N=4,5 two-distinct-bond small-sample trap): (3) STRENGTH (uniform level of f) is NOT proportional to <psi_1|V_b|psi_1>=d eps_1/dJ_b - corr -1.0 at N=4,5 is the 2-bond artifact, collapses to -0.1..-0.57 once N>=6 has 3+ distinct bonds; (4) LOCATION (per-bond f deviation) is NOT the H_1 eigenvector-perturbation footprint sum_k c_k(b) psi_1(a)psi_k(a)/(eps_1-eps_k) - corr is noise (-0.05..+0.67) across N=4..9. SEAM to felt_time-D = RHYME not identity (single-excitation unitary energy-shift, not the half-filling dissipative diffusion rate). LEADING OPEN DIRECTION: the relaxation is governed by the (1,1) LIOUVILLIAN modes (Haken-Strobl population transport), NOT the H_1 Hamiltonian modes - but that does not cleanly fit the given <psi_k|V_b|psi_1> dictionary, so not a simple swap. Next (Tom 2026-06-19): survey what the repo already knows about single-excitation Haken-Strobl / the (1,1) sector / coherence-horizon Q*(N) (which IS the single-excitation Haken-Strobl EP) before guessing a third form; consider a data-first SVD of the (N-1)xN f-dictionary. R_k RESOLVED 2026-06-19 as a NOT-DERIVABLE-ONLY-COMPUTABLE result (Tom's framing: the realm of BirthCanal/IsSteril, reflections/ON_WHAT_CLOSES_ONLY_WITHOUT_US.md). The four-agent survey found the repo already half-knew this (review/EQ014+EQ021, untyped): the (0,1) rate is PROTECTED (re_shift(0,1)~2.8e-16, U(1)+Pi) so f(b) is eigenvector mixing NOT a rate functional; no simple closed form fits (the gate _handshake_rk_only_computable.py: best closed-form R^2 0.38-0.73 erratic at N>=6, the N=4,5 R^2=1.0 a 2-distinct-bond artifact; EQ-021 power-law fails 12-30%). SCHEMA is derivable (rank-(N-2) K-partner rule [KPartnerSelectionRuleClaim], the rate-protection, the (1,1)-population reframe); the VALUE f(b) is IsDeadEnd = the PTF mixing calc instanced (dM_s = sum_s' <W_s'|V_L|M_s>/(lam_s-lam_s') M_s', evaluated not reduced). TYPED: EpistemicFacetMap['handshake_Rk']=IsDeadEnd (model x_peak) + test. The irreducibility is currently EMPIRICAL (no fit across N) not a Tier1 obstruction proof; a 6-route obstruction (model PROOF_F86B_OBSTRUCTION) would lift it. Gram location-metric writeup DONE 2026-06-19, CORRECTED by the math+physics review convergence (both lenses, dispatched per Tom's 'use the physics+math skills' tip, independently caught an INDEXING ARTIFACT in the first draft -- it built M over k=1..N-1, dropping the forbidden k=N column, making the bare couplings look 'well-conditioned'; with the honest k=2..N they carry the same null). The location-metric is THREE distinct objects: (a) the NULL = the K-partner rule (rank N-2, in the BARE couplings k=2..N, Q-free, min sv ~1e-16, already KPartnerSelectionRuleClaim) -- the painted Gram INHERITS it, does not create it; (b) the painted small eigenvalues split by MIRROR PARITY (antisym seesaw + K-partner null / sym closure-strength F123 channel), NOT one homogeneous ambiguity; (c) pairwise confusability = the COSINE matrix (N=5 worst pair cos~-0.97, edge-vs-complementary-interior anti-collinear), a DIFFERENT object from a Gram eigenvalue. Identifiability(rank) != FI(Q)(precision): at Q->inf FI diverges but the null stays null. NO new witness built: it would redundantly recompute the typed K-partner rule, and the physics lens warned against marking the metric's conditioning dead-end-inherited (the IsDeadEnd dead-end is the R_k letter VALUES, the metric's null STRUCTURE is derivable). _handshake_gram_metric.py (corrected) + HANDSHAKE_GEOMETRY.md 'The location-metric' section. REMAINING: M2c read-cost ~2/Q; M3/M4",
            Status: OpenArcStatus.Open),

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
                "oracle simulations/_seam_movement_review.py agrees (Q*(N) transition exact). Witness breadcrumbs " +
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
            Status: OpenArcStatus.Open),

        new OpenArc(
            Name: "f86b2_robust_extraction",
            Opened: "2026-06-11",
            Origin: "simulations/f86b2_shape_invariance_dial.py",
            ParkedAt: "(N,b)-family traces alpha=-0.133 vs fitted -0.129, extraction-noise-limited; g_eff convention gotcha pinned (4.394/(Qp+2))",
            NextStep: "robust extraction, then close the shape-invariance claim",
            Status: OpenArcStatus.Open),

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
                "2026-06-19 (simulations/_stone_survivor_alpha_closure.py; physics + math review GO-with-fixes, " +
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
                "(D) THE CLOSURE-STRUCTURE SEAM (2026-06-19, _stone_mine.py + _stone_seam.py, RESOLVED + one " +
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
                "RESOLVED 2026-06-19 (simulations/_stone_seam.py): both facets confirmed, the seam is a GEOMETRIC " +
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
                "FOLLOW-UP RESOLVED 2026-06-19 (simulations/_felt_time_amplitude_law.py block-level N=4..7 + " +
                "_felt_time_closure_functional.py trajectory N=4..6, gate-first): the exact functional is " +
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
                "simulations/value_vector_felt_time.py (the value/vector crossover), _felt_time_amplitude_law.py " +
                "(the dRe ~ grad^2 diffusion-Rayleigh law, block-level), _felt_time_closure_functional.py (the " +
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
    };

    public static IReadOnlyList<OpenArc> All => _all;

    public static OpenArc? Lookup(string name) =>
        _all.FirstOrDefault(a => a.Name == name);

    /// <summary>Count of arcs still <see cref="OpenArcStatus.Open"/> (the live unfinished business).</summary>
    public static int OpenCount => _all.Count(a => a.Status == OpenArcStatus.Open);
}
