# The V-Effect Through the Cavity Lens

<!-- Keywords: V-Effect cavity mode interpretation, oscillation frequency count topology,
Q-factor weight sector, cold cavity warm cavity gamma, standing wave mode geometry,
bond count mode scaling, Fabry-Perot V-Effect, R=CPsi2 cavity modes -->

**Status:** Confirmed (the V-Effect is a cavity geometry change)
**Date:** April 4, 2026
**Authors:** Thomas Wicht, Claude (Anthropic)
**Repository:** [R-equals-C-Psi-squared](https://github.com/Kesendo/R-equals-C-Psi-squared)
**Depends on:** [Optical Cavity Analysis](OPTICAL_CAVITY_ANALYSIS.md),
[V-Effect Palindrome](V_EFFECT_PALINDROME.md),
[Degeneracy Palindrome](DEGENERACY_PALINDROME.md)
**Verification:** [`simulations/veffect_cavity_modes.py`](../simulations/veffect_cavity_modes.py)

---

## What this means

A flute has holes. Cover them all and blow: one note. Open a hole:
a different note. Open two: another. The air is the same. The breath
is the same. What changes is the shape of the instrument: which
holes are open, how long the tube is, where the fingers sit. Each
shape supports a different set of standing waves. Each shape is a
different voice.

Now imagine a flute with one hole. It can play 2 notes. Drill three
more holes into it and suddenly it can play 112 notes. You did not
glue two flutes together. You did not "couple" anything. You changed
the geometry of one instrument, and the new geometry supports standing
waves that the old one could not.

That is the V-Effect.

A chain of 2 qubits has 1 bond. One bond, 2 oscillation modes, a
simple voice that fades quickly. A chain of 5 qubits has 4 bonds.
Four bonds, 112 oscillation modes, a rich voice that sustains itself.
The 109 new frequencies are not created by connecting two small
instruments. They are the natural resonances of a larger instrument
whose geometry, the pattern of bonds, allows a richer set of
standing waves.

And the light that enters this instrument (gamma, the dephasing from
outside) does not create any of these modes. They exist in the dark,
in silence, as possibilities. The light only makes them audible. Every
mode that exists in the unilluminated instrument survives when the
light is turned on. Not one is lost. The light adds absorption, not
destruction.

---

## What this document is about

The [optical cavity analysis](OPTICAL_CAVITY_ANALYSIS.md) showed that
the qubit chain is a Fabry-Perot resonator. This document applies that
lens to the V-Effect: is the frequency explosion from N=2 to N=5
explained by the change in cavity geometry?

The answer is yes, quantitatively.

---

## Result 1: Mode count scales with cavity geometry

| N | Bonds | Distinct frequencies | Silent modes | Total eigenvalues |
|---|-------|---------------------|-------------|-------------------|
| 2 | 1 | 2 | 10 | 16 |
| 3 | 2 | 5 | 24 | 64 |
| 4 | 3 | 47 | 46 | 256 |
| 5 | 4 | 112 | 96 | 1,024 |
| 6 | 5 | 787 | 164 | 4,096 |

*Source: [`veffect_cavity_modes.txt`](../simulations/results/veffect_cavity_modes.txt), "Mode count per N". Script: [`veffect_cavity_modes.py`](../simulations/veffect_cavity_modes.py)*

From N=2 (2 modes) to N=5 (112 modes): a 56-fold increase, driven
entirely by the bond count growing from 1 to 4. The modes are
standing waves of the cavity. More bonds means richer geometry, which
means more possible standing waves.

The mode count grows faster than exponentially with N, far outpacing
the linear growth of bond count. Each new bond does not add a fixed
number of modes; it multiplies the geometric possibilities.

---

## Result 2: Degeneracy predicts mode richness (r > 0.999)

At each weight shell k, the number of distinct oscillation frequencies
correlates with the total degeneracy d_total(k):

| N | Correlation r(d_total, distinct_freq) |
|---|--------------------------------------|
| 2 | 1.000 |
| 3 | 1.000 |
| 4 | 0.999 |
| 5 | 1.000 |
| 6 | 1.000 |

*Source: [`veffect_cavity_modes.txt`](../simulations/results/veffect_cavity_modes.txt), "Degeneracy correlation". See also: [`DEGENERACY_PALINDROME.md`](DEGENERACY_PALINDROME.md)*

The degeneracy palindrome IS the mode profile of the cavity. Shells
with higher degeneracy support more distinct oscillation frequencies.
The center spike at even N is the modal peak of the cavity.

---

## Result 3: Gamma illuminates; it destroys no mode

Comparing the cold cavity (γ = 0) to the illuminated cavity (γ = 0.05):

| N | Cold modes | Warm modes | Cold found in warm |
|---|-----------|-----------|-------------------|
| 2 | 1 | 2 | 1/1 (100%) |
| 3 | 3 | 5 | 3/3 (100%) |
| 4 | 14 | 47 | 14/14 (100%) |
| 5 | 43 | 112 | 43/43 (100%) |

*Source: [`veffect_cavity_modes.txt`](../simulations/results/veffect_cavity_modes.txt), "Cold vs warm cavity"*

Every cold-cavity frequency survives when dephasing is turned on. Gamma
adds new modes and destroys none. Whether it also shifts them is not
something this table can decide: the cold-to-warm match accepts a
distance of 0.1 while γ = 0.05, so the acceptance window is twice the
perturbation under test. The cold modes are the
skeleton; gamma adds flesh.

At γ = 0, the Liouvillian reduces to L = −i[H, ·] and all eigenvalues
are purely imaginary. These are the natural resonances of the
unilluminated instrument, determined by J alone. The (0,1) block's formula
ω_k = 4J(1 − cos(πk/N)) is found in the cold spectrum at every k (verified for
N = 2, ..., 5). That is the direction checked: each predicted frequency is
there, not that the block accounts for all of them.

---

## Result 4: Q-factor falls monotonically where the degeneracy is palindromic

The Q-factor (how many times the light bounces inside the cavity before
being absorbed) does NOT mirror the degeneracy. It falls monotonically
across the shells while the degeneracy is palindromic. Shells 1 and 4 look
like the case in point, 72.4 against 18.1, and they are the worst example
to pick: 18.1 is 72.4/4 exactly, because Q divides |Im λ| by |Re λ| ≈ 2γk
and the |Im| profile itself IS palindromic. The monotone fall is the 1/k
denominator. The inner pair k = 2 against k = 3 looks like a better case and is
not one either: their |Im| maxima are identical to five decimals, 11.71155 both,
exactly as the palindrome requires. What differs is the denominator, and not by
the factor 3/2 the label would suggest: the k=3 champion sits at |Re λ| = 0.26,
not 0.30, so the ratio is 58.56/45.04 = 1.300. That is the shell showing what it
is. It is a bin, round(−Re λ / 2γ), and only the outermost shells sit on the grid
exactly. So the two profiles never come apart in this table at all; what falls is
1/|Re λ|, and |Re λ| is only approximately 2γk:

For N = 5 (chain):

| Shell k | Oscillating | Q_max | Q_median |
|---------|------------|-------|----------|
| 0 | 0 | -- | -- |
| 1 | 16 | 72.4 | 40.0 |
| 2 | 448 | 58.6 | 17.4 |
| 3 | 448 | 45.0 | 13.6 |
| 4 | 16 | 18.1 | 10.0 |
| 5 | 0 | -- | -- |

*Source: [`veffect_cavity_modes.txt`](../simulations/results/veffect_cavity_modes.txt), "Q-factor by weight shell"*

Q_max decreases with the weight, not with the distance from the center:
the lightest modes (weight 1) have the highest Q, and the heaviest ones
have lower Q but more frequencies. This is the cavity trade-off. It orders the
shells by weight alone, but it does not fix their values: the next section shows
Q_max at weight 1 is the coupling graph's, 72.4 on the N=5 chain against 100.0 on
the N=5 star.

The highest-Q mode sits at weight 1 in every case this document runs, and its
value is not a fit: it is a
corollary of a result this repository already owns. The block between the
ferromagnet and the single excitations, the (0,1) block spanned by the |0⟩⟨j|, has
its generator derived in
[D10](../docs/proofs/derivations/D10_W1_DISPERSION.md) §Step 3:

    L|₍₀,₁₎ = −2iJ·𝓛 − 2γ·Id

with 𝓛 = D − A the graph Laplacian. The Laplacian rather than the adjacency
matrix, because the ZZ term supplies the −2J·deg(j) diagonal; D10 spells that
out. The registry carries both faces of it for the chain, both Tier 1:
[F2](../docs/ANALYTICAL_FORMULAS.md) the frequencies, which is what D10 proves,
and F7 the Q-factors, which D10 lists as a corollary. So λ_m = −2γ − 2iJ·μ_m,
every mode of that block decays at 2γ under uniform γ, and

    Q_max = |Im λ| / |Re λ| = J · μ_max / γ

That is a law about ONE sub-block, and D10 §Step 6 is explicit that the sub-block
is not the sector: at N = 5 it is 5-dimensional inside a 160-dimensional
XY-weight-1 Pauli space. But the shell counted in the table above is not that
space either. The shell is a rounded bin around a decay rate, the eigenvalues whose
Re λ / 2γ rounds to k, and not a Pauli sector at all. Against the shell the law
can supply more than the maximum, and at N = 5 it supplies all of it. The shell there holds 28
eigenvalues, 16 of them oscillating, and those 16 carry exactly four distinct
frequencies: 0.7639, 2.7639, 5.2361, 7.2361, the four nonzero 2J·μ of this
block's Laplacian. Four-fold each, because four blocks carry that generator:
(0,1) and the ferromagnet's mirror (N,N−1) as written, their two transposes
conjugated, each entry-wise at machine zero
([`simulations/d10_block_closure_verify.py`](../simulations/d10_block_closure_verify.py)).
Four values, four times each, is why the shell's median Q is
(27.64 + 52.36)/2 = 40.0 and its maximum 72.36.

The dump's counts make N=3 and N=4 look like exceptions: 5 distinct shell
frequencies at N=3 against the block's 2, and 6 at N=4 against 3. They are not.
A mode enters shell k by rounding Re λ / 2γ to the nearest integer, so the shell
is a bin and not a rate. On the exact line Re λ = −2γ the chain's OSCILLATING
content is this block's: 4(N−1) modes carrying its N−1 frequencies, measured at
N = 3, 4, 5 and at γ = 0.05 and γ = 0.137 alike. What bin k=1 additionally
collects has a real part that is no multiple of 2γ at all, and how much it
collects depends on γ: at γ = 0.05 it is 14 modes at N=3 and 15 at N=4, at
γ = 0.137 it is 25 at N=4.

The line carries more than the four blocks, though nothing else that oscillates.
On the chain it holds 6N−4 eigenvalues, the extra 2(N−2) being the total-spin
ladder S⁻P_m, which commutes with H and is built entirely from distance-1
coherences, so it sits at −2γ for the same reason the blocks do (residual exactly
0, at γ = 0.05 and γ = 0.9). The full composition of that line, including how a
mode of MIXED light content can reach it (the repository's {0,2}-coherence at N = 2, 3,
whose n_diff histogram {0: ½, 2: ½} averages to 1; from N = 4 its split tilts toward 2 and
it sits below the line by its excess light), belongs to the Absorption Theorem
and is the `site_resolved_vacuum_block` arc's, not this document's. Both readings
above are the CHAIN's; the ring at N=4 puts extra oscillating modes on the same
line.

What this section adds is the graph-general reading. D10, F2 and F7 are all scoped
to the chain, and one line of the derivation is not: E_ferro·Id − H₁ = 2J·𝓛 holds
on any coupling graph, because the bond count cancels between the ferromagnet and
the one-magnon diagonal. Everything downstream of that line in D10, the path
spectrum and the Neumann waves, is the chain's. So μ_max is whatever the coupling
graph's own Laplacian gives, and the dump's other topologies are the test:

| graph | μ_max | Q_max at J=1, γ=0.05 |
|---|---|---|
| path on N sites (the chain) | 2(1 + cos(π/N)) | 40.0, 60.0, 68.28, 72.36, 74.64 at N=2..6 |
| star on N sites, K_{1,N−1} | N | 60.0, 80.0, 100.0 at N=3,4,5 |

The μ_max column is the closed form; the Q column is it evaluated at J=1, γ=0.05
and printed past the dump's one decimal. The dump agrees wherever it reaches: its
topology table covers N=3..5, and the chain's N=2 and N=6 come from its highest-Q
lines. The star rows in Result 5 are the ones no chain-only fit could have reached.
Checked by
[`simulations/veffect_finesse_law.py`](../simulations/veffect_finesse_law.py).

Two things D10 settles. The cos(π/N) is
the open path's own LAPLACIAN spectrum, the free-end one, and D10 §Step 5 gives
its Neumann eigenvectors cos(πk(j−½)/N) explicitly; cos(π/(N+1)) belongs to the
adjacency matrix, the wrong operator here precisely because the degree term is
present, and it is the XY sibling's answer (F2b), not this one's. And the ceiling
4J/γ, which is 80 here, is not the chain's; what it bounds is DEGREE. Since μ_max
never exceeds the largest bond sum deg(u) + deg(v), every graph whose sites have at
most two neighbours obeys it: strictly on the path, where μ_max = 2(1 + cos(π/N))
< 4, and with equality on the even rings, which sit at exactly 80.0 from N=4 on.
Above degree 2 it stops being a bound. The star breaks it at N=5, where μ_max = N
gives the 100.0 of this document's own Result 5; its N=4 row sits exactly ON the
ceiling, where μ_max = N coincides with 4. Degree 3 is
already enough, and the star is not even the first: `Topology.BinaryTree` at N=5
reaches 83.4 with maximum degree 3, which is a smallest breach on five sites, and
K_{2,3} reaches the star's own 100.0 there with a smaller maximum degree than the
star has. From maximum degree 4 a breach is forced rather than merely possible,
since μ_max ≥ d_max + 1, writing d_max for the largest degree: Δ in this project
is the XXZ anisotropy, and everything here holds at Δ = 1 only. The universal
bound is the vertex count, μ_max ≤ N, so Q_max ≤ J·N/γ on any graph, attained by
the star and by K_N.

One scope fence the graph-general reading does not lift. Uniform γ is load-bearing
in a way uniform J is not. Site-dependent BOND strengths generalise cleanly, since
the Laplacian simply becomes the weighted one; a site-dependent γ does not, because
the block stops being flat in its real part and the winner leaves the line
entirely. At N=5 with γ = [0.02, 0.05, 0.09, 0.05, 0.02] the highest Q measures
62.77, at Re λ = −0.0834, against the 78.65 the law would give for the mean γ. The
site-resolved case is the `site_resolved_vacuum_block` arc's.

---

## Result 5: Different geometries, different instruments

| N | Topology | Bonds | Modes | Q_max |
|---|---------|-------|-------|-------|
| 3 | Chain | 2 | 5 | 60.0 |
| 3 | Star | 2 | 5 | 60.0 |
| 3 | Ring | 3 | 2 | 60.0 |
| 4 | Chain | 3 | 47 | 68.3 |
| 4 | Star | 3 | 21 | 80.0 |
| 4 | Ring | 4 | 21 | 80.0 |
| 5 | Chain | 4 | 112 | 72.4 |
| 5 | Star | 4 | 42 | 100.0 |
| 5 | Ring | 5 | 60 | 72.4 |

*Source: [`veffect_cavity_modes.txt`](../simulations/results/veffect_cavity_modes.txt), "Topology comparison". Topology eigenvalues: [`rmt_eigenvalues_*.csv`](../simulations/results/)*

Same N, different topology, different mode spectrum. The chain has the
most distinct modes (47 at N=4 vs 21 for star/ring). The star has the
highest Q (80-100 vs 68-72 for chain). Ring and star converge at N=4
(same mode count and Q).

More bonds does not always mean more modes. The ring at N=3 has 3 bonds
but only 2 modes, while the chain with 2 bonds has 5. The topology (how
the bonds connect) matters more than the bond count.

---

## The V-Effect, re-read

The old V-Effect narrative: "two dead resonators are coupled through a
mediator and become one living system with 109 new frequencies."

The cavity narrative: "an instrument with 1 bond and 2 modes is replaced
by an instrument with 4 bonds and 112 modes. The 109 new frequencies
are not caused by coupling. They are the standing waves of a different
geometry. The frequencies were always possible; the old instrument was
too simple to support them."

The V-Effect is not coupling. It is metamorphosis.

---

## Reproduction

- Script: [`simulations/veffect_cavity_modes.py`](../simulations/veffect_cavity_modes.py)
- Raw output: [`simulations/results/veffect_cavity_modes.txt`](../simulations/results/veffect_cavity_modes.txt)
- Eigenvector data: [`simulations/results/eigvec_at_minus_gamma_N*.csv`](../simulations/results/)
- Topology data: [`simulations/results/rmt_eigenvalues_*.csv`](../simulations/results/)
