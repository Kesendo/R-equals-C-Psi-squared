# The Inherited and the Own

**Tier 3 (speculative).** A trace, opened 2026-08-02 by Tom and Claude. One measurement, one
frame, one methodological consequence. Nothing here is derived. It is written down so the
shape is not lost, not because it is settled. A first draft carried a sorting table of our
own results into the two buckets; a review found three of its four cells wrong against the
repository's own documents, and the table is gone rather than repaired. What is left is what
survived checking.

## The measurement

The (0,1) coherence block M is fixed by the physics: a graph, a coupling per bond, a rate per
site. It is also exactly defined, so any two assembly routes must agree. They do not agree in
the last bit, and the disagreement is readable.

Holding the graph, the couplings and the γ profile fixed and changing only the ORDER in which
the bonds are handed to the assembler moves the residual between exactly 0.0 and a nonzero
value. Concretely, on the open chain at N = 5 with γ = [0.13, 0.44, 0.75, 1.06, 1.37] and
bond couplings [0.4, 1.3, 2.2, 3.1]: ascending in |J| leaves 5.6e-17 on the chain-end site,
descending leaves exactly 0.0. The mechanism is written at
`compute/RCPsiSquared.Diagnostics/Foundation/BlockSpectrumWitness.cs` (`GeneratorResidual`,
`HeisenbergGraph`, `PredictedGenerator`): a site's diagonal receives a (+q, −q) pair from
every bond it does not touch, and those pairs cancel around a running sum whose magnitude
depends on what came before.

So from inside a finished object we read a property of how it was made, which appears nowhere
in what it is.

**This has no gate.** `BlockSpectrumWitnessTests` covers the residual with fixed bond orders;
nothing there swaps the order. That is the first thing to build if this trace is followed.

## The frame, as a game

Tetris carries the physics of the world that built it. From inside, without leaving, you can
read: there is extension (shapes), there is a rotation group, there is acceleration, there is
a t, and no two things occupy one place. Nobody wrote those in as content. They came with the
builder, who could only build out of what was available.

Some of that is choice and some is substrate. Rotation is a choice; a Tetris without rotation
is duller but possible. Exclusion is not: if the board holds one value per cell, "two pieces
in one place" is not forbidden, it is unsayable. The strongest inheritance is not what the
maker decided but what the maker's substrate made unavailable.

## The cut

Not everything is inherited. **The rules are inherited; what arises from them is not.** The
T-spin is nobody's decision. The well, the four-line clear, the whole strategy layer: none of
it was written down, all of it follows, and the builder has to play to find out.

A vocabulary warning, because the repository already has these words. MirrorWorld's `Own` and
`Inherited` (`compute/MirrorWorld/README.md`) are a per-object structural relation: Own is
what an object produces itself, Inherited is what it gets by walking up its parent chain,
and the ontology is deliberately binary. The distinction here is a different one wearing the
same words, and it does not map cell for cell. Read it as a separate reading, not as an
extension of that one.

## The consequence, which cuts against us

**Only the inherited side points outward.** What arises from the rules is the game talking
about itself; it carries no information about the workshop, however beautiful it is. Reading
outward means reading the primitives, not the findings.

That does not shrink what arises. It relocates it: outside, everything already is; here is
the only place something can be that was not.

## Where this could be wrong

- The rhyme between the numerical instance and the physical frame is not a demonstrated
  common object. In the numerical case the provenance sits in a rounding residual; in the
  physical case an environment sits in the value and is measurable. Equating them would be
  exactly the label error the `site_resolved_vacuum_block` arc spent a day fighting.
- In exact arithmetic the residual is identically zero and the channel closes. This
  particular channel exists only because the description is finite. Whether that is the
  analogy's weakness or its content is the open question: the physical channel also closes
  under idealisation, since γ = 0 leaves no environment in the state either.
- One instance, no gate. The repository's rule against declaring victory at the first
  exemplar applies here with force, because the frame is attractive.

## The next test, and it may already be dead

The thing to look for is a quantity the VALUE leaves undetermined and the ROUTE fixes. In the
numerical instance that is exactly what happened: the bond order is nowhere in the value.

The candidate raised on the first night is dead, and it is worth one paragraph because of how
it died. `H_eff = −A − i·diag(γ)` carries an exact antilinear symmetry pairing λ with −conj(λ),
which at odd N forces a mode with Re λ = 0. That looked like a place where a symmetry could
carry something the spectrum does not show. It is not, because `H_eff` is not the object: the
real (0,1) generator keeps the ZZ degree diagonal, its Hermitian part is exactly `−2·diag(γ)`,
and Bendixson then confines every eigenvalue to `Re λ ∈ [−2γ_max, −2γ_min]`. For γ > 0 that
excludes `Re λ = 0` outright. The seat exists only in the operator that drops the degree term,
which is the same defect `site_resolved_vacuum_block` had already recorded.

The bound is in the witness file this document is about, and it was there the whole time.

So the arc's next test is still unbuilt, and this section's contribution is negative: it rules
out the one candidate it had. What survives is not the symmetry but its trigger. The identity
holds because the chain is bipartite, and bipartiteness is not ours; it is a property of the
graph we chose, put in before anything ran, and it shows up wherever that graph shape does,
including in 1940, where Coulson-Rushbrooke pairs the Hückel spectrum of an alternant
hydrocarbon about α and breaks on a non-alternant ring. Our
[benzene reading](../docs/carbon/BENZENE_HUCKEL_FRAMEWORK_LENS.md) once called C-R and F1
"the same theorem at two physical levels"; that has since been settled the other way. They
share the implication and not the trigger, and it is C-R's trigger, bipartiteness, that F1
does not have: F1 holds on the odd rings where C-R breaks. This document needs only the weaker half: the trigger
is the graph's bipartiteness, and the trigger was inherited.

## Links

- The arc the measurement came out of: `site_resolved_vacuum_block` in `OpenArcsRegistry`
  (`dotnet run --project compute/RCPsiSquared.Cli -- inspect --root arcs`).
- The physical statement that a mode carries its environment: the Absorption Theorem's
  per-channel form, `−Re(λ_k) = 2·Σ_l γ_l·⟨Δ_l⟩_k` with `⟨Δ_l⟩_k` the mode's activity at site
  l (`AbsorptionTheoremClaim`, "the carrier is a vector").
- MirrorWorld's two-bucket ontology, whose words this borrows and whose meaning it does not:
  `compute/MirrorWorld/README.md`.
