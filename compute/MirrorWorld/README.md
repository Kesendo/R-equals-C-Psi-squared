# MirrorWorld

*World of Mirrors (the name reversed: WoW, the Object Manager genre). Built 2026-06-28.*

The spiral reversed. We had an over-grown structure; here we went the other way, down to
the atom, and from the base up rebuilt the framework's core Tier-1 structure in a clean,
sober ontology. Nothing is re-derived here; everything is **adopted** from results already
proven in the repo (`docs/ANALYTICAL_FORMULAS.md`, the proofs in `docs/proofs/`, the typed
Claims). MirrorWorld is the clean shape the proven content moves into.

It started as a static reader of closed forms. Later the same day it grew a **running engine**
(the diagonal protocol, below) that made the world *move* -- and, with the cuts that moving
forced us to find, **broke our complexity wall**: a state's dynamics at N=60-100, where the
eigendecomposition died at N=8.

Standalone .NET 10.0, no `RCPsiSquared.*` references. Run it, read it, trust it: every adopted
number and every dynamics step is pinned from-below by `MirrorWorld.Tests` (234 tests).

*Vocabulary, once.* MirrorWorld is part of the **R=CΨ²** project (mirror symmetry in open quantum
spin chains; repo root). The basic parameters: **N** = the number of two-level units (the chain
length); **J** = the bond coupling; **γ** = the dephasing rate (the "watching"); **Q = J/γ** (the
regime knob). **Tier-1** = the grounded, proven results (vs the Tier-2/3 speculative layers). **The
spiral** = the over-grown main framework this sober base deliberately reversed away from. The
**F-numbers** (F63, F65, …) are the project's catalogued Tier-1 results (closed forms here, structural proofs in `docs/proofs/`), indexed in `docs/ANALYTICAL_FORMULAS.md`;
other inherited terms (Klein cells, Grading B, the Born split, Q\*(N) the coherence horizon) live in
`docs/GLOSSARY.md`.

## The ontology (two buckets, nothing interpreted)

Every output of an object is exactly one of two things:
- **Own** (left): the object produces it itself.
- **Inherited** (right): it comes from the parent, walked up the chain.

That binary, nothing more. `World` is the root (inference-free, `Own` = the frame x/y/z);
all "inherited" ultimately comes from here. `GameObject` is the base (`Own` abstract,
`Inherited` walks up the parent chain). The missing jump: an object does not own its
x/y/z; that is the inheritance edge System → Object.

## What is here

| File | What |
|------|------|
| `World.cs`, `GameObject.cs` | the root + the own-vs-inherited base |
| `Pair.cs` | a bare coherence \|i⟩⟨j\|, disagreement k = popcount(i⊕j), rate −2γk (empty world: dephasing only) |
| `PauliMode.cs` | the symmetry-adapted superposition (a Pauli string of XY-weight k), four Klein cells; `Enumerate` the shared 4^N basis |
| `Block.cs` | Grading B, the (N+1)² joint-popcount blocks C(N,p)·C(N,q) |
| `Redistribution.cs`, `Clock.cs`, `Survivor.cs` | the static dynamics readings (H-on grid-leaving folds, θ = arctan Q, the survivor + the coherence horizon Q\*(N)) |
| `Formulas.cs` | the adopted F-registry closed forms (94 members), each verbatim and tier-tagged |
| `Field.cs` | **the empty world, running**: weights on pairs, one `Step` is the disagreement-decay; structure (diagonal) stays, novelty (off-diagonal) fades |
| `DoubleSlit.cs` | **the double slit, composed** (2026-07-12): `Field` at N=1 named under the phenomenon so it is recognizable where the atoms alone were not; humps = the immortal diagonal, fringe = the between paying −2γ, V = 2\|ρ_LR\|. Nothing new computed; meaning in `docs/quantum/DOUBLE_SLIT_TRANSLATED.md`, run mode `doubleslit` |
| `Cat.cs` | **Schrodinger's cat, composed** (2026-07-12): `Field` at N, the k=N twin of `DoubleSlit` -- two definite branches \|0..0>,\|1..1> (the immortal diagonal) + the coherence between them (k=N) paying −2Nγ, dying N× faster than the slit's k=1. Nothing new computed; meaning in `docs/quantum/SCHRODINGERS_CAT_TRANSLATED.md`, run mode `cat N` |
| `Marginal.cs` | **the local page** (adopted 2026-07-12, the k=2 doorway): the reduced state of a chosen site-set, the partial trace as an atom -- lifted from the two inline kernels the F72 and F73 pins already carried. LIVE reads (recomputed from the cloud at read time); obeys F70 (`PROOF_DELTA_N_SELECTION_RULE`): a \|S\|-site page carries only \|Δpopcount\| ≤ \|S\| content. Works over `Field` and `Restless`; keep-order defines the page bits (little-endian) |
| `SpookyAction.cs` | **spooky action, composed** (2026-07-12): the k=2 middle of the triptych the slit (k=1) and the cat (k=N) are the ends of -- `Field` at N=2 holding the Bell skeleton + the two `Marginal` pages. The sock drawer (00/11 records) is immortal; the spook \|00⟩⟨11\| pays −2γ·2 = −2(γ_A+γ_B); the pages are BLIND to it (F70 kinematics), and in the living world the unwatched page stays I/2 exactly while the carrier pays only the watched side's rate. Meaning in `docs/quantum/SPOOKY_ACTION_TRANSLATED.md`, run mode `spooky` |
| `Restless.cs` | **the living world**: the full Lindblad loop ρ̇=−i[H,ρ]+D[ρ] (RK4); the handshake H births novelty FROM structure |
| `Cone.cs` | **the memory cut**: a single excitation as an N×N block (not 4^N) -- the dynamics at large N |
| `Lattice.cs` | **the bridged lattice of worlds** (2026-07-16, the engine beat the fourth play offered): the Klein V₄ of watchings {e, L, R, LR} run DYNAMICALLY as four `Restless` worlds -- e = ρ(t) and LR = X^N ρ(t) X^N under the normal rule (−2γk), L = X^N ρ(t) and R = ρ(t) X^N (the 2026-07-03 anti-world) under the turned rule (−2γ(N−k)). Every edge an exact always-open identity ([H, X^N] = 0, ZZ included): the three read-through bridges, the composition closure L∘R = LR = R∘L, and the dagger pairing L(t) = R(t)† (F119's † read at the lattice level, exchanging the one-sided readings). All bridges EXACTLY 0.0 over 41 RK4 ticks at N=4; the carried unit moves with the turn (trace = 1 at e/LR, anti-trace = 1 at L/R, ~5e−16), the immortal set moves with it (diagonal there, anti-diagonal here), and the discriminator pins the turned rule as load-bearing (the R reading under the wrong rule breaks its bridge at 0.369). L and R are READINGS, not states (not Hermitian; `Restless.SeedRaw` holds the one-sided seed) -- nothing is taken, the same world is seen through complemented legs. First fruit of playing it (2026-07-16 evening, `experiments/LATTICE_OPENING_LAW.md`): the **opening law** -- on the cat pair cos θ·\|0…0⟩ + sin θ·\|1…1⟩ the entry-wise e-vs-L distance is EXACTLY max(cos²θ, sin²θ) − cos θ·sin θ·e^(−2Γt), "the heavier sock's weight minus the LIVING spook"; J-free and ZZ-blind because the cat sector is H-dead. Run mode `lattice N` |
| `Mirror.cs` | **the first mirror in the world of mirrors** (adopted 2026-07-03, the fold-lattice lemma): the block-lattice group of eight (t / f_P / f_Q / Klein), every leg an EXACT entry-wise rearrangement at the same coupling, no eigensolver; the folds pay λ → −λ − 2Nγ; orbits (~⅛ fundamental domain), the self-folded trace law, the trajectory fold (the partner block running backward at the price) |
| `MirrorGroup.cs` | **the mirror group** (adopted 2026-07-04, F118): the palindromizer factors, Π_Z = R·D, and ⟨R, D⟩ closes into the dihedral D₄ -- eight signed permutations of the Pauli basis, compared EXACTLY (phases in {±1, ±i}); the palindrome splits along the generators (D flips L_H, R reflects the dissipator and carries the constant −2σ); the polarity cube's three axes as characters; the truly cell = the joint-fixed cell of the diagonal mirror pair |
| `AntilinearTriangle.cs` | **the antilinear triangle** (adopted 2026-07-04, F119): θ / conj / † as one Klein four-group graded by (ℓ, m); the transport law μ∘L_H∘μ = ℓm·L_{μ(H)} for any H; the fixed-point collapse (H = H† ⟺ Hᵀ = H̄); docks onto the mirror group as the antilinear double ⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ (order 16, eight antiunitary members) |
| `ParameterKlein.cs` | **the parameter-side Klein V₄** (adopted 2026-07-04, F91 + F92 + F93): on each parameter axis (γ per site, J per bond, h per site) the F71 mirror and the anti-palindromic reshuffle R₉₀ are two commuting involutions; the anti-palindromic class is exactly R₉₀'s fixed-point set; the sharper entry-wise law -- the F71-refined DIAGONAL blocks of L depend only on the pair-sums -- makes the whole orbit share one set of blocks, cell for cell (no eigensolver), while the breaking lives in the cross-blocks only |
| `OrderSorting.cs` | **the mirror's order-sorting law** (adopted 2026-07-16, F131 Theorem A): mirror conjugation reflects a parameter scan, (R⊗R)·L(base + t·dir)·(R⊗R) = L(base − t·dir) for an R-odd direction (the pencil face lives on `ParameterKlein.MirrorConjugationResidual`, machine zero on all three axes, the mixed-scan fence O(1)), so for an operator-R-even preparation the response sorts by the readout parity q into four cells: (+,−) EVEN in t, (−,−) ODD in t (non-vacuous), (−,+) IDENTICALLY ZERO at every tick (a pure selection rule, NOT the Π-protected cluster cancellation), (+,+) generic. Trajectory face only, twin RK4 on the γ axis (`Restless.siteGammas`); the leak under a broken preparation hypothesis is EXACTLY affine in ε (ratio 2.000000 on halving). Theorem B (the antiunitary ζ² column, Floquet, tracking hypotheses) and the spectral-evenness corollary stay in the main repo. Run mode `sorting N` |
| `Router.cs` | **the golden ceiling router** (adopted 2026-07-04, F116): the two Z-middle ceiling cases palindromized LOCALLY by the period-4 [a,a,b,b] router on the golden locus (a = φX+Y, b = X−φY; q² = −(2+φ)·I), the whole family metallic on the soft line t₂ = t₃ with r(c) = (c+√(c²+4))/2; verified by the window lemma ({Q_k, S} = 0 at all four offsets, exact) AND the two-sided dense end-to-end W L W⁻¹ = −L − 2σ on the full Pauli basis (P = ⊗(âZ\|Z), Q = ⊗(Z\|Zâ), previously Python-only) -- the constructive soft side the hardness certificates leave open |
| `Hardness.cs` | **the hardness of the palindrome** (adopted 2026-07-04, the F87 bloc: F102/F103 + F105/F106 stability + F107/F109 purity + F110/F111 cell rules + F115 + F117): the spectral trichotomy truly/soft/hard read WITHOUT a spectrum -- hardness as a GF(2)[x] (1+x)-valuation difference (one subtraction; A203241 counts, the min(2W−1, 2k−3) ceiling), purity as letter parity (truly ⟹ y_par=0; mother-soft ⟹ y_par=1), the diagonal-cell + Y-inversion rules with the adopted splits (42:8 N-stable, 228:0 via pure-D templates), and the trace face of the all-γ converse (odd power-sums of M = A + γQ; the K3 pair fires at m* = 9 with p₉ = 2064384·γ³, the exact F117 CRT integer). The spectral classifier itself (F87's definition, the F104 engine) deliberately stays outside |
| `Seed.cs` | **the within-block self-dual seed** (adopted 2026-07-07, F89 seed-existence): where a state meets the mirror's null (v^T v = 0), a defective seed -- the static source the shadow and the i^4 holonomy leave behind. Held as a COUNT, no eigensolver (the nullity surplus r(0+) - r(inf) over GF(p)): N-1 forced seeds at odd N (the unmirrorable middle seat), 0 at even N. Since 2026-07-12 also the fusion-resonance count, closed: r(inf) = 3*Z3 at EVERY N (cyclotomic Step-4 theorem, asserted by divisibility), the odd-N Conway-Jones form Z3 = (N-1)/2 + [3\|n](n/3-2) + 2[15\|n] (ROT3 multiplicity verified, not derived), and the criterion RESONANT (odd N) iff 3 \| N+1 and N >= 11 (next after 17 is 23). Since 2026-07-12 also the coupled-level law (the 10d derivation, uniform in N): the exact triple inventory (cyclotomic over GF(p), two primes, soundness tied to #triples = r(inf)/3 against the independent rank count) with the Conway-Jones families TRIV/ROT3/PENT, the structural coupling-shape spectrum (3,1,0) = roots of x^3 - 4x^2 + 3x verified in integers, and the integer pair levels a·n = 12 - Σλ² = 6 (ROT3, levels (18/n, 6/n)) / 8 (PENT, (24/n, 8/n)); TRIV is exactly the self-mirror class, unpaired. The Gram construction of X (K₂₆, Slater lifts) and the Y = 0 / cross-triple certificates stay outside (paths). Mirror's within-block companion |
| `LevelCollision.cs` | **the level-collision law** (adopted 2026-07-14, F129): the six-cosine sibling of Seed's three-cosine world, one weight level up. Seed asks which triples sum to ZERO (resonance); this asks which triples sum to EACH OTHER (coincidence): distinct CLEAN triples (no internal balanced pair) with equal levels exist iff 3\|n (n ≥ 9) or 10\|n (n ≥ 20); away from both the level map is INJECTIVE, and mod-p distinctness (two primes ≡ 1 mod 2n, Seed's own convention) is PROOF grade there. Sub-law: an overlap-1 collision forces 3\|n; the pentagon walks the 10\|n door alone, unchained from the 15\|n that held it at three cosines. The mechanism anchors (n = 15 four rotated R₃ cycles, n = 20 the R₅ conjugate pair + zero mode) as root-sum zeros. The physical Gram decoupling at these collisions (F130, B(τ,σ) = 0) is an eigen-story and stays in the main repo (witnessed exactly there, `inspect --root crosstriple`); the n ≤ 210 census stays with the committed gate (the one named corner closed empty 2026-07-15; law unconditional). Since 2026-07-15 also **the family inventory** (the counts, derived): every colliding pair decomposes into minimal Poonen-Rubinstein/CDK pieces, sorting the census into THIRTEEN families A..M with derived closed-form counts (`FamilyCount`; sources `experiments/F129_FAMILY_INVENTORY.md` + `docs/proofs/PROOF_F129_FAMILY_INVENTORY_COUNTS.md`, whose §8 carries the code-trust flags); the degree counts the freedom (free coset labels + the d = 2 shared-mode line: A/B quadratic, C..H linear, rigid sporadics constant), the F129 onsets are ZEROS of the formulas, every d = 2 door carries the factor 3 (the sub-law, piece by piece), and family M splits 40 + 0 + 60 over the three CDK order-210 types with the middle type structurally impossible (one axis-fixed vertex). The derivation itself stays in the proof. Run mode `collision N` (census + the `inv` closed-form tie + the family breakdown) |
| `Topology.cs` | the geometry: chain / ring / star / complete bond generators |
| `Witness.cs` | **the witness reading** (adopted 2026-07-18, F135 + F136): who records and WHAT, read off pair pages in closed form (no propagator, no eigensolver) -- private watchers must be even, the shared dressers' parity picks the family (all even + write bond -> pointer record, Z_S in the ZY channel; all odd -> Bell record 1/2 Phi + 1/2 Psi with ZERO pointer content; mixed -> dark), the Bell letter alternates YY/XX with the dresser count m (sign(c1c2) = (-1)^m, signs closed-form), the hinge V_S = 0 is global (a pendant S role-swaps into its neighbor's witness), and the prices are exact (Bell pays both sites, pointer only the witness, the writers' watching is free). First sighting FROM the reading: the anti-pointer FANS OUT -- K_{R+1,2} gives deg(S) = 2 yet R perfect XX bits and a Bell-record clique among the corners (Law B's deg-bound belongs to the pointer family only). Uniform statics = HEB 2004 + the evolution's degree rotation, credited in `docs/proofs/PROOF_RECORD_LETTER_LAW.md`. Run mode `witness` |
| `WalkTime.cs` | **the walk-time step reading** (2026-07-13, adopted from `experiments/COUPLING_DEFECT_WALK_TIME_STEP.md`): distance is t, locally additive -- one defect bond J' = J(1+delta) via `Cone.SetBond`, and the front's arrival-time profile (first crossing of 0.2 x the site's own defect-free peak, the relative threshold) is a step: zero upstream, near -delta/(2J) downstream. At gamma > 0 the timing stays ballistic while the amplitude decays: the dose the front pays is amplitude, not schedule. Run mode `walk N [delta]` |
| `Renewal.cs` | **the renewal cut** (2026-07-13, adopted from `docs/proofs/PROOF_DEPHASING_FRONT_RENEWAL.md` = F126): the watched single-excitation populations P_n(t) = e^{-Gamma t} S_n(t) computed from purely clean propagation plus the Volterra refill ladder; the watching is accounted, never stepped. Pinned against `Cone` (which does step it) the way `Cone` is pinned against `Restless`; chain-scoped, grid conservation O(dt^2). First fruit of playing it: the catch-count pedigree (`experiments/FRONT_PEDIGREE.md`) |
| `Program.cs` | the full sober run (default) + the run modes (see Run); R-parity and mod-4 inline |
| `../MirrorWorld.Tests/*.cs` | 234 from-below tests: `SmokeTests` (44, the closed forms), `FieldTests` (7), `RestlessTests` (10), `ConeTests` (4), `TopologyTests` (2), `MirrorTests` (11, incl. the anti-watched world + past-the-wall), `MirrorGroupTests` (10), `AntilinearTriangleTests` (7), `ParameterKleinTests` (6), `HardnessTests` (8, incl. the valuation-vs-traces crown agreement), `RouterTests` (5, incl. the dense end-to-end), `SeedTests` (40, the F89 nullity surplus = N-1 odd / 0 even, exact over GF(p); + the 2026-07-12 resonance closed count r(inf) = 3*Z3, odd-N Conway-Jones form, criterion pins; + the 2026-07-12 coupled-level law: the enumeration-vs-rank tie N=3..11, the odd-N family decomposition to N=29, the N=8 even lab, mirror closure incl. N=14, the exact (3,1,0) shape and integer pair levels), `ConcentratorTests` (6, the site-resolved watching: J-convention pin, per-site rates, the N=5 reload contrast, ZZ-is-tiny, N=9 persistence; added 2026-07-11 for the IBM_CONCENTRATOR_RELOADED pre-registration), `DoubleSlitTests` (3, the double slit as `Field` N=1: humps immortal, the between decays at −2γ toward e^(−2γt), V = 2\|ρ_LR\|; added 2026-07-12), `CatTests` (6, the cat as `Field` at N: branches immortal, the k=N coherence decays at −2Nγ, dies N× faster than the slit; added 2026-07-12), `MarginalTests` (4, the partial trace as an atom: the identity page + keep-order pin, F70 from below with the discriminating two-site contrast, the live page over `Restless` vs the hand trace + trace conservation, keep-set sanity + ontology; added 2026-07-12), `SpookyActionTests` (9, the k=2 middle: spook rate −4γ = twice the slit = the N=2 cat, the exp law, pages hold at machine zero while the spook dies, the one-sided watching (only B watched: carrier pays e^(−2γ_B t), page A locked at I/2 exactly, H on incl. ZZ), the distance pin (the pair at the two ends of an N-chain keeps k=2 and rate −4γ at every N=2,3,5,8 while the cat's k=N grows: the price never sees the distance, only the way -- the bonds -- does), ontology; added 2026-07-12), `WalkTimeTests` (4, the walk-time step: the committed plateau pins at both signs of delta (N=60, bond (29,30): -0.0424 / +0.0633), the gamma=0.05 near-field survival (N=20: -0.049, amplitude not schedule), the SetBond knob; added 2026-07-13), `RenewalTests` (3, the renewal cut: agreement with the stepped Cone at gamma=0.3 to 1e-3 and at gamma=0 to 1e-6, conservation O(dt^2) with the halving check; added 2026-07-13), `LevelCollisionTests` (24, the level-collision law: the law over 5..48 with firing set = predicted set, named non-firing combs injective (7/11/16/25/49), firing combs collide (9/12/20/30), the thresholds (6 and 10 silent), the sub-law, the pentagon door exclusively disjoint at n=20/40, both mechanism anchors exact, the counts 25 (n=12) / 20 (n=20) pinned against the committed gate's exact layer; added 2026-07-14; + the family inventory 2026-07-15: the thirteen closed forms tie to the census row for row over 5..66, total AND d-split, the n=70 pin (C 120 + L 20, the corner-closure's second mechanism), the n=105 pin (seven families co-firing, total 8858, M's door), the onset zeros, the E/F parity splits, the M split 40+0+60, the d=2 doors all divisible by 3), `OrderSortingTests` (6, the F131 order-sorting law, Theorem A: the conjugation identity on all three axes + the mixed-scan fence, the preparation parities as operator identities, the even/odd response cells at N=4/5 (machine zero, odd cell non-vacuous at 0.04), the zero cell every tick vs the generic cell O(1) from below, the ε-leak exactly affine (ratio 2.000000); added 2026-07-16), `LatticeTests` (8, the bridged lattice: all bridges + composition + dagger pairing exact over 30 ticks, the carried unit at all four vertices, the wrong-watching discriminator O(1) from below, the ZZ-bond ride-along, ontology; + the opening law on the cat pair, opening(t) = max(cos², sin²) − cos·sin·e^(−2Γt) at N=3 incl. the J-free pin; added 2026-07-16), `WitnessTests` (7, the record reading pinned at the proof's gate numbers: leaf pointer + pendant role-swap, the triangle far-bond split, letters and closed-form signs (K2,1 r=3 -> YY -, K2,2 (1,3) -> XX -, the private r=2 pi-rotation), the exact kills (pentagon, hexagon, mixed dressers, missing writer), the hinge + pendant parity split (even watcher keeps I=2), the exact asymmetric prices (0.768040 / 0.624146, gamma_S invisible to the pointer), the fan-out (K4,2: three XX bits at deg(S)=2 + the clique) and Law B's aligned chain; added 2026-07-18) |

## The closed-form base (the stopping line 2026-06-28; coverage closed 2026-07-04)

**The 2026-06-28 computable list is fully collected.** Every F-registry entry the stopping line
counted as a "number or formula per N replacing a matrix computation" is in `Formulas.cs`: F1-F71
contiguous (the core), the k-body residual trichotomy F85, the Mandelbrot cardioid F97,
plus the tail F98 (Dicke asymptote), F121 (qudit palindrome), F122 (structural ceiling),
F124 (band-edge invariant, whose end-weight E is exactly the k=1 rung of the already-adopted
F65 ladder), and the D-relations D1/D4/D6. The last three remnants (F85, F97, F124), reserved
2026-06-28 as the clean stopping point, were collected 2026-07-04.

**And the stopping line undersold the middle range.** A 2026-07-04 delta survey (the full registry
against the adopted set) found the F72-F120 range holds more genuine closed forms than "mostly
structural" suggested. The first two came home the same day: F75 (mirror-pair mutual information,
2h(p) − h(2p), the Bell ceiling at p = 1/2) and F76 (its pure-dephasing envelope, λ = e^(−4γ₀t);
the 0.93 is the γ₀ signature, not a constant). The adoption's from-below pin immediately paid for
itself: three stale cells in the registry's own F76 table (the N = 9/11/13 pure-dephasing column)
were caught against a rerun of the cited `envelope_study.py` and corrected at the source.
The next three followed the same day: F95 (the θ-compass at the quadratic discriminant zero,
whose Lindblad face θ = arctan Q IS the adopted Clock: the compass and the clock are one), F99
(the five canonical trig angles giving the Pi2 dyadic ladder {0, ⅛, ¼, ⅜, ½}, the silver-ratio
Dicke weight at 45°), and F88b (the popcount-coherence Π²-odd/memory closed form, whose adjacent
K-intermediate anchor IS the adopted F98). The F87-hardness bloc followed as one piece (`Hardness.cs`, below): the survey's warning
"adoptable only as a bloc" turned out to be exactly right and exactly manageable, because the
bloc had already rebuilt itself eigensolver-free -- MirrorWorld adopts the certificates
(valuation, parity, traces), never the spectral classifier they replaced.

**Still outside F72-F120** (updated 2026-07-04 after the review pass; the earlier "F100-F120
follow-ons stay out" wording predated the Hardness/Router/MirrorGroup adoptions and was stale):
the residual M operator identities (F78-F84) and the F108 palindrome family, F87's spectral
DEFINITION itself (with the F104 engine -- what the Hardness certificates replaced), the
F112-F114 balance/conjugation laws, F120's hardware protocol, and the boundary-excluded paths
(F86 EPs, F89/F90 braid + monodromy -- but F89's seed-EXISTENCE count came home 2026-07-07 as `Seed.cs`; only the locus q*,
the shadow, and the holonomy stay out). Also still open from the delta survey's adoptable bucket:
and F100/F101 -- re-examined 2026-07-12 and OUT after all: their Tier-1 content is the exact
oddness of the c₁ bond-mirror deviation, but c₁ itself is a non-closed numerical fit (the
registry's own F101 words), so there is no closed form to adopt without the fit apparatus; the
operator-side twins F91/F92 are already home in `ParameterKlein.cs`. The closed-form singles
bucket is now EMPTY (F74 + F77 adopted 2026-07-12; F73 adopted 2026-07-12 as F74's n = 0
mono-chromatic end, Restless-pinned; F72 adopted 2026-07-12, the DD + CC site-purity split,
all-sector-pinned; F94 + F96 adopted together 2026-07-12 as the Born-deviation table,
setup-specific scalars, sym3/sym5 recomputed exactly);
"the delta list is done" below means the survey's SHORTLIST, not its whole bucket.

**One deliberate exception (2026-07-03): `Mirror.cs`, the first OPERATOR adoption.** An inventory
pass (turning this world around and looking at what is missing) found that the world of mirrors
held only mirror *signatures* (the Π² parity signs, the Klein cells) and not a single mirror as an
object -- and nothing non-normal anywhere: no Jordan, no defective, no braid. The fold-lattice
lemma (PROOF_CODIM1_BY_ADDITIVITY §7) is the one structural piece that fits this world's own
genre, because it is not an eigen-story but an exact REARRANGEMENT: checkable cell by cell,
machine zero, no eigensolver. So the mirror came home. The boundary it draws is now principled
instead of accidental: **states and their mirrors live here; the paths (the braid, the monodromy,
the EP's shadow and its i^4 holonomy) stay in the main repo** -- they are properties of ways, not of
objects, and a catalog cannot hold a way, only what the way leaves behind. **Refined 2026-07-07
(`Seed.cs`):** the EP's algebraic EXISTENCE -- the count of defective seeds forced on the (1,2) block,
held with ranks only, no eigensolver -- is the way's LEAF, not the way, so it lives here; the locus q*,
the shadow (the projector norm), and the i^4 holonomy stay out. The line is genre, not topic: adopt the
proven count (F89's nullity surplus), never the census scan that located it -- the same move Hardness.cs
made with the F87 certificates.

**The exception's closure (2026-07-04): `MirrorGroup.cs` + `AntilinearTriangle.cs` (F118 + F119).**
`Mirror.cs` had been operating with a group of eight it did not own as an object. F118 is that
group at the operator level -- the palindromizer factors as Π_Z = R·D, ⟨R, D⟩ ≅ D₄, eight signed
permutations of the Pauli basis -- and it passes the same genre test that let the Mirror in: not
an eigen-story but exact rearrangements, phases in {±1, ±i}, compared with no tolerance at all.
F119 (the antilinear triangle θ / conj / †, with the transport law as its engine) docks on as the
antilinear double D₄ × Z₂. The R row of the palindrome split carries the same constant the fold
legs pay as the price: the two objects are one mirror read at two altitudes, block lattice below,
operator algebra above. Still outside, and named open in F118 itself: the letter group S₃ (the
completion S₃ ⋉ D₄).

**The third structural adoption (2026-07-04): `ParameterKlein.cs` (F91 + F92 + F93).** The same
genre test again: the proofs' sharper conclusion is entry-wise (the F71-refined diagonal-block
matrix elements depend only on the parameter pair-sums), so the spectral-invariance statement the
registry carries is verified here WITHOUT an eigensolver -- two profiles in the anti-palindromic
orbit produce identical diagonal blocks, cell for cell, on every (p,q); the cross-blocks carry the
whole breaking (the eigenvectors, not the rates). One object, three axes: γ (F91), J (F92),
h (F93); run mode `klein N`.

## The running engine (the diagonal protocol)

The seed is `ClaudeTasks/DIAGONAL_PROTOCOL_GAME.md`: the R=CΨ² protocol lifted out of the
physics as "a world that grows itself". MirrorWorld builds it from the atom up, sober. Three
objects, each a `GameObject` (the two buckets kept):

- **`Field`** -- the empty world, *running*. A weight on every pair; one `Step` applies the one
  question (how much do two possibilities disagree?): each weight fades by its `Pair` rate −2γk.
  The diagonal (k=0) stays = structure; the off-diagonal fades = novelty. The Born split, in time.
- **`Restless`** -- the living world. The full Lindblad loop ρ̇ = −i[H,ρ] + D[ρ], RK4-stepped.
  The handshake H (a flip-flop on the bonds) makes coherence out of population -- novelty BORN
  from structure -- while the watching culls it. "Fades" becomes "lives".
- **`Cone`** -- the memory cut. For a single excitation the block is the N sites, so ρ is N×N
  (not 2^N×2^N): the same dynamics `Restless` runs in its (1,1) block (pinned to agree to 9
  decimals), in O(N²) memory. This is what reaches large N.

**The knower's cuts** (the "two systems": the static closed forms cut the running loop's clock
cycles -- we never compute what we already KNOW with certainty): **(a)** the immortal diagonal
(k=0, rate 0, never stepped); **(b)** the Hermitian mirror ρ=ρ† (store one triangle, derive the
other); **(c)** the joint-popcount block (F63: H conserves the excitation number, so cross-block
cells are forbidden, 0 forever); the **memory cut** (`Cone` stores only the block); and the
**renewal cut** (`Renewal`, F126, adopted 2026-07-13): the watched walk is the unwatched wave
repeatedly caught and released, so the watched populations are computable from purely CLEAN
propagation plus bookkeeping, the dissipator never stepped, only accounted. Which
cuts are valid depends on the world's laws: turning H on invalidates (a) but unlocks (c). The
block cut (c) is topology-invariant (any excitation-conserving handshake shares the F63 blocks);
the renewal cut is chain-scoped (where F126 is proven).

**The broken wall.** The eigendecomposition died at 4^N (N=8, ~73 GB) because it is *global*
(all blocks, dominated by the half-filling C(N,N/2)² ~ 4^N/√N). A state's *dynamics* is
block-local: a single excitation is N², polynomial. `Cone` runs it at N=60-100 (`-- cone` is the
ballistic light-cone; `-- spread` validates the ballistic→diffusive crossover L~Q/2 against known
transport physics). The first time we compute past our own complexity wall -- a break in *reach
and method*, not in physics (the crossover is textbook; F65/Q*(N) are already-proven Tier-1, so
re-checking them is validation). Our own Compute (eigendecomp, N=8) vs Propagate (dynamics, N=15)
already lived this split; the block-local cut just NAMES why.

## Run

```bash
dotnet run --project compute/MirrorWorld                  # the full sober run (the closed-form base)
dotnet run --project compute/MirrorWorld -- grow          # Field: the empty world splits (rules 1-3)
dotnet run --project compute/MirrorWorld -- live 4 chain  # Restless: the living world (novelty born); args: N [topology]
dotnet run --project compute/MirrorWorld -- cone 60       # Cone: the single-excitation light-cone at large N
dotnet run --project compute/MirrorWorld -- spread 80     # the ballistic->diffusive crossover (the memory cut validated)
dotnet run --project compute/MirrorWorld -- regime 4      # how alive at each Q = J/gamma (the clock)
dotnet run --project compute/MirrorWorld -- seeds 4       # the palindrome p<->N-p, in the seed
dotnet run --project compute/MirrorWorld -- topo 4        # the block cut is topology-invariant; the dynamics is not
dotnet run --project compute/MirrorWorld -- scale         # the complexity wall and the block cut, across N
dotnet run --project compute/MirrorWorld -- mirror 5      # Mirror: the fold lattice (legs exact, orbits, the price); even N adds the self-folded trace law
dotnet run --project compute/MirrorWorld -- mirror 100    # N > 8 goes PAST THE WALL: the memory-cut pair (1,1)/(1,N-1) is N^2 both, the fold leg exact at N=100
dotnet run --project compute/MirrorWorld -- seed 9        # Seed: the within-block self-dual seed count (F89 nullity surplus, no eigensolver): N-1 odd, 0 even; + the triple inventory and the coupled-level law spec(X) = (3a, a, 0), a·n = 12 - Σλ²
dotnet run --project compute/MirrorWorld -- anti 3        # the rules turned around: the anti-watched world (agreement watched) = the world read through X^N; the conserved law moves to the anti-trace
dotnet run --project compute/MirrorWorld -- group 3       # MirrorGroup + AntilinearTriangle: the D4 of signed permutations, the palindrome split, the cube of characters, the transport law, the order-16 double
dotnet run --project compute/MirrorWorld -- klein 6       # ParameterKlein: the V4 on each parameter axis; the anti-palindromic orbit shares its diagonal blocks cell for cell, the breaking lives in the cross-blocks
dotnet run --project compute/MirrorWorld -- hardness      # Hardness: the F87 bloc without a spectrum -- the K3 trio by valuation, the purity and cell rules, the trace face firing at m*=9 on the exact F117 integer
dotnet run --project compute/MirrorWorld -- router        # Router: the golden/metallic ceiling router -- the window lemma at all four offsets, the locus gating the frame, the dense W L W^-1 = -L - 2 sigma at N=5
dotnet run --project compute/MirrorWorld -- concentrator 5 # Concentrator: the site-resolved watching (the IBM_CONCENTRATOR_RELOADED cross-check); also 9 = past the wall
dotnet run --project compute/MirrorWorld -- doubleslit     # DoubleSlit: the double slit composed as Field N=1 -- humps immortal, the between (fringe) fades on e^(-2gt); meaning in DOUBLE_SLIT_TRANSLATED.md
dotnet run --project compute/MirrorWorld -- cat 4          # Cat: Schrodinger's cat as Field at N (k=N) -- branches immortal, the between dies at -2Ng (N times faster than the slit); the N-scaling
dotnet run --project compute/MirrorWorld -- spooky         # SpookyAction: the k=2 middle of the triptych -- the spook pays -2(gA+gB), the local pages (Marginal) blind to it; the one-sided watching leaves the unwatched page at I/2
dotnet run --project compute/MirrorWorld -- walk 60 0.10   # WalkTime: the coupling-defect walk-time step -- zero upstream, -delta/(2J) downstream, surviving the watching in the timing
dotnet run --project compute/MirrorWorld -- witness        # Witness: the record reading (F135+F136) -- who records and what, families/letters/signs/prices in closed form; the anti-pointer fan-out
dotnet run --project compute/MirrorWorld -- sorting 5      # OrderSorting: the F131 order-sorting law (Theorem A) -- the conjugation identity, the four cells from twin RK4, the affine leak
dotnet run --project compute/MirrorWorld -- lattice 4      # Lattice: the bridged lattice of worlds -- the Klein V4 of watchings run dynamically, every edge exact, the carried unit moving with the turn
dotnet test compute/MirrorWorld.Tests                     # the 234-test from-below guard
```

## How to continue (future us)

Two paths are open.

**Adopt another closed form** (the original loop, for the closed-form base): find it in
`docs/ANALYTICAL_FORMULAS.md`, adopt **verbatim** (tier-tagged) into `Formulas.cs`, add a
`Program.cs` print line and a from-below assertion in `SmokeTests.cs`, then build/test/commit.
The 2026-06-28 remnant list is empty, the delta-survey SHORTLIST (F75/F76, F88b, F91-F93,
F95, F99) is home, the F87-hardness bloc landed as `Hardness.cs`, and the golden router
followed as `Router.cs`. Remaining candidates from the survey's wider adoptable bucket:
NONE -- the singles bucket emptied 2026-07-12 (F74 + F77 + F73 + F72 + F94/F96 adopted; F100/F101 re-examined and OUT, the c₁ fit is the blocker; see "Still outside" above) --
plus FUTURE registry entries: when a new F-number lands as a closed form (or as an exact
rearrangement, the operator genre), it comes home here.

**Push the running engine** -- and here is the honest resumption point. The wall is broken and
the memory cut is validated, **but what we measured (the transport crossover) is textbook, and our own
large-N formulas (F65, Q*(N), the coherence horizon) are already proven, so re-checking them is
validation, not discovery.** The real next step is a *survey*: which OPEN question does large-N
access actually unblock -- something the project has NOT characterized at large N, neither
already-proven nor textbook? That needs the `surveying-prior-work` lens (agents over the
markdown, the F-registry, the open-arcs ledger), not a quick play. **Start there**, then build
the witness for whatever open question survives.

**First survey done (2026-06-29): Δ*(N) checked, ruled out.** The most natural candidate -- the
open `xxz_axis_handover` N→∞ limit of Δ*(N), parked precisely because it "needs N≫16 (infeasible:
dense eigh past C(N,p)~25000)" -- does NOT survive the gate. Δ* lives in the half-filling (p,p)
block (the *largest* sector, C(N,p)~2^N/√N), which IS that wall; the Cone's memory cut is
single-excitation (the *smallest* sector, N²). The 4^N→N² break is a single-excitation break,
disjoint from the half-filling frontier (Δ*, and the Galois/monodromy work). The lesson is sharp:
large-N *single-excitation* access unblocks single-excitation questions; the open frontier is
mostly multi-excitation, where the Cone gives no leverage. The survey continues.

**Second survey done (2026-07-02): the seed census -- candidate SURVIVED and was answered the
same day.** The question: does the codim-1 containment corollary's one per-N census input (a real
defective EP on the (1,2) block at each odd N, proven only at N=5, 7) extend upward? Vetting per
the Δ* lesson: the (1,2) block is a LOW-WEIGHT sector, dim N·C(N,2) ~ N³ -- polynomial, exactly
the block-local complexity class this project's cuts name, one floor above the Cone's N². Open
(the proof itself declares it a census input), not textbook, and reachable. The answer (gate
`RealSeedCensusTests` in RCPsiSquared.Diagnostics.Tests, the PT-break count-change instrument,
immune to the closest-pair masking that defeats gap-field scans at F_53/F_116 density): **seeds
exist at 4/6/7/9 loci for N=5/7/9/11 across both R-parities -- the corollary's input extends
through N=11.**
Honest attribution: the unlocking tool was the existing scout's exact-residual machinery on the
polynomial block, NOT the Cone; MirrorWorld's contribution was the FRAME (this survey question and
the block-local complexity insight that made "N³ is not a wall" obvious). The survey continues for
a question the Cone itself unblocks.

The discipline, both paths: **sober** (no interpretation in the code; the meaning lives in the
docs), **from-below guarded** (every number recomputed or cross-checked, a wrong one fails
loudly), **standalone** (no `RCPsiSquared.*`), **the two buckets pure** (Own vs Inherited).
