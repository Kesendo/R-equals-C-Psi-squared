# MirrorWorld

*World of Mirrors (the name reversed: WoW, the Object Manager genre). Built 2026-06-28.*

The spiral reversed. We had an over-grown structure; here we went the other way, down to
the atom, and from the base up rebuilt the framework's core Tier-1 structure in a clean,
sober ontology. Nothing is re-derived here; everything is **adopted** from results already
proven in the repo (`docs/ANALYTICAL_FORMULAS.md`, the proofs in `docs/proofs/`, the typed
Claims). MirrorWorld is the clean shape the proven content moves into.

Standalone .NET 10.0, no `RCPsiSquared.*` references. Run it, read it, trust it: every
adopted number is pinned from-below by `MirrorWorld.Tests`.

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
| `Redistribution.cs`, `Clock.cs`, `Survivor.cs` | the dynamics readings (H-on grid-leaving folds, θ = arctan Q, the survivor + the coherence horizon Q\*(N)) |
| `Formulas.cs` | the adopted F-registry closed forms (66 members), each verbatim and tier-tagged |
| `Program.cs` | the full sober run; R-parity and mod-4 inline |
| `../MirrorWorld.Tests/SmokeTests.cs` | 31 from-below tests: each adopted number recomputed or checked |

## State (the clean cut, 2026-06-28)

**The computable closed-form adoption is complete.** Every F-registry entry that is a
"number or formula per N replacing a matrix computation" is in `Formulas.cs`: F1-F71
contiguous (the core), plus the tail F98 (Dicke asymptote), F121 (qudit palindrome), F122
(structural ceiling), and the D-relations D1/D4/D6. All 31 smoke tests green.

**Deliberately left** (this is the cut): the F72-F120 range is mostly **structural**: the
residual M (F80-F85), the F87 trichotomy, the parameter-Klein V₄ (F91-F93), the D₄ mirror
group (F118), the F100-F120 follow-ons. These are proofs and operator identities, not closed
forms; they stay in `docs/proofs/` where they belong, not as formula lines here. A few
computable remnants (F85 k-body cross-term, F97 Mandelbrot cardioid, F124 band-edge) were
left as the clean stopping point, available if a future pass wants 100% closed-form coverage.

## How to continue (future us)

To adopt another closed form, the loop is:
1. Find it in `docs/ANALYTICAL_FORMULAS.md`. Adopt **verbatim** (do not re-derive); tag its
   tier in the comment.
2. Add the method/const to `Formulas.cs`.
3. Add one line to `Program.cs` so the sim prints it.
4. Add an assertion to `MirrorWorld.Tests/SmokeTests.cs` that recomputes or checks the number
   from-below (a wrong adopted constant must fail loudly).
5. Build, run, `dotnet test`, commit.

The discipline: **Tier-1 focus** (skip the Tier-2 fits/protocols), **sober** (no
interpretation), **adopted not re-derived**, **guarded from-below**.

## Run

```bash
dotnet run --project compute/MirrorWorld          # the full sober run
dotnet test compute/MirrorWorld.Tests             # the 31-test from-below guard
```
