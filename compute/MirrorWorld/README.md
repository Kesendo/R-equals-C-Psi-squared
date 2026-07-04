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
number and every dynamics step is pinned from-below by `MirrorWorld.Tests` (93 tests).

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
| `Formulas.cs` | the adopted F-registry closed forms (80 members), each verbatim and tier-tagged |
| `Field.cs` | **the empty world, running**: weights on pairs, one `Step` is the disagreement-decay; structure (diagonal) stays, novelty (off-diagonal) fades |
| `Restless.cs` | **the living world**: the full Lindblad loop ρ̇=−i[H,ρ]+D[ρ] (RK4); the handshake H births novelty FROM structure |
| `Cone.cs` | **the memory cut**: a single excitation as an N×N block (not 4^N) -- the dynamics at large N |
| `Mirror.cs` | **the first mirror in the world of mirrors** (adopted 2026-07-03, the fold-lattice lemma): the block-lattice group of eight (t / f_P / f_Q / Klein), every leg an EXACT entry-wise rearrangement at the same coupling, no eigensolver; the folds pay λ → −λ − 2Nγ; orbits (~⅛ fundamental domain), the self-folded trace law, the trajectory fold (the partner block running backward at the price) |
| `MirrorGroup.cs` | **the mirror group** (adopted 2026-07-04, F118): the palindromizer factors, Π_Z = R·D, and ⟨R, D⟩ closes into the dihedral D₄ -- eight signed permutations of the Pauli basis, compared EXACTLY (phases in {±1, ±i}); the palindrome splits along the generators (D flips L_H, R reflects the dissipator and carries the constant −2σ); the polarity cube's three axes as characters; the truly cell = the joint-fixed cell of the diagonal mirror pair |
| `AntilinearTriangle.cs` | **the antilinear triangle** (adopted 2026-07-04, F119): θ / conj / † as one Klein four-group graded by (ℓ, m); the transport law μ∘L_H∘μ = ℓm·L_{μ(H)} for any H; the fixed-point collapse (H = H† ⟺ Hᵀ = H̄); docks onto the mirror group as the antilinear double ⟨R, D, 𝒦⟩ ≅ D₄ × Z₂ (order 16, eight antiunitary members) |
| `ParameterKlein.cs` | **the parameter-side Klein V₄** (adopted 2026-07-04, F91 + F92 + F93): on each parameter axis (γ per site, J per bond, h per site) the F71 mirror and the anti-palindromic reshuffle R₉₀ are two commuting involutions; the anti-palindromic class is exactly R₉₀'s fixed-point set; the sharper entry-wise law -- the F71-refined DIAGONAL blocks of L depend only on the pair-sums -- makes the whole orbit share one set of blocks, cell for cell (no eigensolver), while the breaking lives in the cross-blocks only |
| `Topology.cs` | the geometry: chain / ring / star / complete bond generators |
| `Program.cs` | the full sober run (default) + the run modes (see Run); R-parity and mod-4 inline |
| `../MirrorWorld.Tests/*.cs` | 93 from-below tests: `SmokeTests` (36, the closed forms), `FieldTests` (7), `RestlessTests` (10), `ConeTests` (4), `TopologyTests` (2), `MirrorTests` (11, incl. the anti-watched world + past-the-wall), `MirrorGroupTests` (10), `AntilinearTriangleTests` (7), `ParameterKleinTests` (6) |

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
Surveyed candidates still open: F88b (popcount-coherence closed form), F95 (θ-compass), F99, and
the rest of the delta list.

**Deliberately left** (this is where adoption stops): the genuinely structural rest of F72-F120:
the residual M identities (F80-F84), the F87 trichotomy, the F100-F120 follow-ons. These are
proofs and operator identities, not closed forms; they stay in `docs/proofs/` where they belong,
not as formula lines here.

**One deliberate exception (2026-07-03): `Mirror.cs`, the first OPERATOR adoption.** An inventory
pass (turning this world around and looking at what is missing) found that the world of mirrors
held only mirror *signatures* (the Π² parity signs, the Klein cells) and not a single mirror as an
object -- and nothing non-normal anywhere: no Jordan, no defective, no braid. The fold-lattice
lemma (PROOF_CODIM1_BY_ADDITIVITY §7) is the one structural piece that fits this world's own
genre, because it is not an eigen-story but an exact REARRANGEMENT: checkable cell by cell,
machine zero, no eigensolver. So the mirror came home. The boundary it draws is now principled
instead of accidental: **states and their mirrors live here; the paths (the braid, the monodromy,
the exceptional points) stay in the main repo** -- they are properties of ways, not of objects,
and a catalog cannot hold a way, only what the way leaves behind.

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
cells are forbidden, 0 forever); and the **memory cut** (`Cone` stores only the block). Which
cuts are valid depends on the world's laws: turning H on invalidates (a) but unlocks (c). The
block cut (c) is topology-invariant (any excitation-conserving handshake shares the F63 blocks).

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
dotnet run --project compute/MirrorWorld -- anti 3        # the rules turned around: the anti-watched world (agreement watched) = the world read through X^N; the conserved law moves to the anti-trace
dotnet run --project compute/MirrorWorld -- group 3       # MirrorGroup + AntilinearTriangle: the D4 of signed permutations, the palindrome split, the cube of characters, the transport law, the order-16 double
dotnet run --project compute/MirrorWorld -- klein 6       # ParameterKlein: the V4 on each parameter axis; the anti-palindromic orbit shares its diagonal blocks cell for cell, the breaking lives in the cross-blocks
dotnet test compute/MirrorWorld.Tests                     # the 93-test from-below guard
```

## How to continue (future us)

Two paths are open.

**Adopt another closed form** (the original loop, for the closed-form base): find it in
`docs/ANALYTICAL_FORMULAS.md`, adopt **verbatim** (tier-tagged) into `Formulas.cs`, add a
`Program.cs` print line and a from-below assertion in `SmokeTests.cs`, then build/test/commit.
The 2026-06-28 remnant list is empty and the first of the delta-survey candidates (F75/F76,
F91-F93) are home; the loop continues with the surveyed rest (F88b, F95, F99, ...) and with
FUTURE registry entries -- when a new F-number lands as a closed form, it comes home here.

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
