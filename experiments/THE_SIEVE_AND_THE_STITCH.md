# The Sieve and the Stitch

<!-- Keywords: sieve, stitch, every N, every length, family census, residue, far element, flat section, quadratic
phase, binomial system, local holonomy, free chain, ring law, F158 end count, complement connection Theorem 2,
F173, MirrorWorld instrument, how to run the sieve -->

**Status:** two instruments of [MirrorWorld](../compute/MirrorWorld/README.md), `Sieve.cs` and `Stitch.cs`, built
2026-11-10 (the stitch's free chain 2026-11-11), pinned from below in `compute/MirrorWorld.Tests/SieveTests.cs` and
`StitchTests.cs`. The census runs quoted here were made on 2026-11-20 with the commands shown beside them; the read ring
cells a law decides are printed by the family runs, and which of the rest have no law comes from calling
`Stitch.RingLaw` on each listed cell in code, the same day. The first census's
counts and several letter lines are also pinned in the tests.

## What this is about

Take a chain of qubits, or a ring, and let every qubit stand in light: the repository's word, here, for dephasing along Z ([what the light does](../docs/quantum/DEPHASING_TRANSLATED.md)). The qubits are coupled by a Hamiltonian built from one small pattern, a **cell**: a few Pauli strings, such as XZ + ZY (X on one site and Z on the next, plus Z on one site and Y on the next), laid on every pair of neighbouring sites, or on every run of sites as long as the pattern. Every way the state can fade, an eigenmode of the generator built from H and the light together, has a rate and a frequency, and the whole spectrum can be mirror-symmetric about one centre, the sum of the light's rates over the sites: for each way of fading there is a partner that oscillates at the same frequency and fades as much faster than the centre as the first fades slower. That is the palindrome. [F158](../docs/proofs/PROOF_PALINDROME_TWO_END_COUNT.md) decides it at one length by a count at the two ends of the spectrum (the stitch's paragraph below says what is counted), with 2^N Pauli strings at each end, a number that doubles with every site. The question one wants answered is different: does this cell keep the palindrome at **every** length, on the chain and on the ring?

**The sieve** takes a whole family of cells, for instance every cell of one to three of the nine two-letter bond terms with signs ± (417 cells), and sorts each one, on the chain and on the ring, with the cheapest exact tests first. First the stitch's fast no (below). Then the readers, tests whose yes covers many lengths at once: a map of operators built one site at a time, the mirror [the Theodolite](THE_CUBE_SITE_BY_SITE.md) builds from the Pauli letters, and the stitch. Only the lengths that no yes covers and no no breaks go to F158's count, one length at a time, up to ten sites, and on the ring two laws in N from [F173](../docs/proofs/PROOF_SIGNED_MASK_BALANCE.md) are held against the answers. Each cell ends, on each geometry, with one of five verdicts, read, broken, mixed, decided by a law, or residue, defined exactly further down. A residue is a cell that pairs at every length examined that no reader covers: something carries the palindrome there that the instruments do not yet know. Every residue these families have shown has, once understood, become a new reader.

**The stitch** is the reader the first residue taught. F158 counts two kinds of operators: those that never change at all, and those that fade at exactly twice the centre without oscillating. The palindrome holds exactly when there are as many independent operators of the one kind as of the other, and one operator of the second kind that has an inverse, a carrier in the complement connection's sense, makes it hold. Under light along Z the second kind are the operators X^⊗N·D, where X^⊗N flips every bit and D is diagonal, that commute with H. Write H_xy for the entry of H between bitstrings x and y, x̄ for x with every bit flipped, and d_x for the entries of D. Theorem 2 of [the complement connection](../docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md) then asks H_x̄ȳ·d_y = d_x·H_xy for every pair x, y. Where such a D exists with no zero entry the palindrome holds; where every such D must vanish somewhere it is broken. A [colouring](THE_PALINDROME_AS_A_COLOURING.md) is the case where the entries of D are signs set by the parity of chosen bits, a linear rule. The stitch tries a rule one order up: D a product of one-site phases and of two-site phases on the pairs of sites a cell reaches, the same factors on every site (quadratic where the factors are ±1 and ±i). Each equation then reads "a product of the unknown factors equals a ratio of Gaussian integers" (numbers a + bi with whole a and b; the cell's coefficients are whole numbers), and such equations are solved exactly. Each equation looks only at a few sites around the bits it flips, so from a fixed length on, the stitch's window, the set of distinct equations no longer changes. The window is 3·span sites on the ring and 2·(3·span − 2), at least 4, on the chain, span being the number of sites the pattern covers. A solution at the window reads every longer chain or ring; an edge (a pair x, y with H_xy ≠ 0) whose two entries H_xy and H_x̄ȳ differ in size, or a diagonal entry unequal to its complement's, breaks them all; a system with no solution of this form leaves them silent. Below the window each length is solved on its own. The fast no is the stitch's cheapest part: an edge whose two entries differ in size, or a loop of bitstrings around which the ratios H_xy/H_x̄ȳ do not multiply to 1, found inside a stretch of up to twelve sites. Such a stretch recurs on every longer chain, and in the bulk on every long enough ring, so the palindrome breaks at every length from there on.

**What they found.** On the 417 bond cells the chain splits cleanly from three sites on: 149 read at every length, 268 broken at every length, nothing between. On the ring the other tests and F158's count leave a residue, XZ ± ZY and its twin YZ ± ZX, with or without ZZ: broken on every chain examined, paired on every ring examined. The stitch reads it on every ring, and F173's Corollary E2 decides all twelve cells there too (`Stitch.RingLaw` returns it for each). For XZ + ZY a carrier is the Clifford operator W = X^⊗N·S^⊗N·∏ CZ, CZ on every bond of the ring (S† in place of S for XZ − ZY): it turns each XZ term into the ZY term one bond back and each ZY term into the XZ term one bond forward, a reshuffling of the ring's terms that the chain's two ends cut. CZ on every bond carries the twin onto this family, up to the sign of H. With the stitch the ring's residue is empty. A second family, cells of one or two three-site patterns with a letter at each end, left a residue on the chain, X·(X ± Y)·X with its relatives, and the stitch's second reader, the free chain, was built from it (for cells whose flipping terms all flip translates of one set of bits, no edge vanishing; it checks Theorem 2's conditions on the sizes and the diagonal, and that every square of two flips closes). The cell pairs on every chain: the stitch reads three sites, F158's count four and five, the free chain every length from six on. On the ring F173's Corollary E2 decides it from five sites on, pairing exactly where 3 does not divide the length, and the rings of three and four pair as well.

**Where the reading stops.** A yes of a reader that speaks for many lengths is a proof at every length it covers, on the chain usually all of them, on the ring sometimes only some; a law decides every length from its start on; a no holds for every length only where the output says so. Everything else holds at the lengths examined and says nothing past them, and a silence, where no test answers, is not a no. The stitch, its fast no, its free chain and the ring laws read light along Z on every site, and the sieve turns light along one other letter on every site to Z for them; it reads only jumps that are single Pauli letters. F158's count and Theorem 2 do not depend on the light's rates as long as all are positive, and the command line sets every rate to 1.

*Where the arc goes on: this page is step 5 of [Story 10](../docs/READING_GUIDE.md#story-10-what-carries-the-palindrome-september-to-november-2026) of the reading guide, where its two roads meet: the sieve's readers come from [the Theodolite](THE_CUBE_SITE_BY_SITE.md) (step 4) and from Theorem 2 of [the complement connection](../docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md) (step 2); step 6 is [F173](../docs/proofs/PROOF_SIGNED_MASK_BALANCE.md), whose ring laws the sieve asks.*

<!-- STAGE 0 SWEEP RECORD (for us, not for the reader: hidden on GitHub, read by every session, agent and hook)

What the repo already held (sweep, 2026-11-20). The instruments themselves: compute/MirrorWorld/Sieve.cs and
Stitch.cs with their header comments, the README rows of both, the run modes in Program.cs (sieve, sieve cells),
SieveTests.cs and StitchTests.cs. The pages they read: PROOF_PALINDROME_COMPLEMENT_CONNECTION (Theorem 2, the base of
the stitch; Theorem 8, the absent edges, live as PalindromeComplementConnectionClaim.PlanarHoppingPalindrome with PlanarHoppingWitness at the inspect root
planarhopping and not wired into the sieve), PROOF_SIGNED_MASK_BALANCE (F173: Theorem D and Corollary E2 adopted as
Stitch.RingLaw; Corollary E1 and Theorem E name the free chain; Corollary E2 decides XZ ± ZY on every ring and Theorem
E names its carrier; its 24 law cells and 864 two-template cells of the window family; its record names two gaps,
SOFTNESS_IS_N_DEPENDENT against Stitch.FreeEveryN, cited in the body here, and Corollary E3 against the sieve's field
family, left open here), THE_CUBE_SITE_BY_SITE (the Theodolite's every-N reading, the sieve's first yes),
PROOF_PALINDROME_TWO_END_COUNT (F158, the end count asked on the rest). docs/READING_GUIDE.md Story 10 item 5
describes the pair and links this page. docs/ANALYTICAL_FORMULAS.md: F173's entry names the sieve's law cells, its
field family and Stitch.FreeEveryN, with SOFTNESS_IS_N_DEPENDENT and F116's golden router adjacent; its line on the 8
Corollary E2 cells says "break exactly at 3 | N" without the guard (the ring of three pairs), as do the comment in
SieveTests.cs and CLAUDE.md's Stitch bullet, all outside this change;
no F-number belongs to the sieve or the stitch. OpenArcs registry: no arc names the sieve or the stitch; an arc that asks an
every-N question about another object is f138_converse_failures (whether a palindrome without a colouring needs a
blind dephased seat at every N), and the settled f115_valuation_gate_width (verdicts on windowed pairs, its window
floor settled by F173 Theorem B). fw.Confirmations: nothing (no hardware reading). docs/GLOSSARY.md: no
entry for either word. Typed layer: no Claim or witness for the sieve or the stitch; SignedMaskBalanceClaim holds
second implementations of what the sieve asks (RingBreaks for Theorem D, RingOneShape for Corollary E2) and of what it
does not (Chain for Theorems B and B′, DecideCharacter for Theorem E, ChainTwoShapes and ChainShapes for Corollaries
E3 and E3′), and SignedMaskBalanceWitness walks the same 36 three-site templates as the window family; the two never
meet. The word "the stitch" is taken elsewhere for another object: SecondClockRegimeClaim and SecondClockRegimeWitness
(TheStitch, two regimes of the second clock), KnowledgeRegistryFactory's comment on it,
SecondClockRegimeClaimRegistration in Runtime, BlockSpectrumWitness, the
secondclock root's description in Cli's InspectCommand, and FOLD_SHADOW_IN_EXISTING_HARDWARE's cross-run stitch (reused in
ABSORPTION_RUNG_LADDER_HARDWARE_PREDICTION);
MirrorWorld's README calls its stitch a working name, and neither word is in check_spent_words. docs/CAUGHT_ERRORS.md:
entries on F173's sweep naming the stitch, one on Sieve.Count and negative rates (fixed in the code: a negative rate is
refused), and one on the golden window's certificate said to hold on any topology (the router exists on the ring only
at 4 | N). experiments/: no other page on the instruments; experiments/SYMMETRY_CENSUS.md is another census
(symmetries of the Heisenberg Liouvillian), not this one.

Adjacency asked: SOFTNESS_IS_N_DEPENDENT argues by hand what the free chain checks (components are cubes, their
cycles generated by squares), and its crossing XXXX + XYYY + YYYX is pinned in SieveTests.cs as a Mixed chain cell;
cited in the body. NON_HEISENBERG_PALINDROME and TWO_TERM_PALINDROME_KLEIN_ROUTING read the 36 all-plus two-term cells
of the bond family on the chain (TWO_TERM at N = 3, 4, 5; NON_HEISENBERG's counts), and V_EFFECT_FINE_STRUCTURE ("Longer chains") holds the line at every
N ≥ 3 by F158's two certificates; the sieve's 22/14 split on that slice agrees (checked 2026-11-20, cited in the
body). The F87 / F115 classifier line (Diagnostics PalindromeMaskClassifier, WindowedHardnessClaim,
WindowedObstructionScan, THE_PALINDROME_CLASSIFIER, MirrorWorld's Hardness) decides the same object at fixed N or
above a window floor, not every N on the ring. EndCount.CliffordSymmetry and THE_PALINDROME_AS_A_COLOURING's "A
Clifford that permutes the bonds" search at one N the Cliffords that permute H's terms, of which the stitch's W for
XZ + ZY is one read at every N (Stitch.cs's header says so; linked in the body). SieveTests.cs holds the sieve against
Theorem 1 and F138's clause 2 on hopping bonds with fields (128 agreements; all 64 cells inside Theorem 8). F116's golden window XZX + XZY + YZX: palindromic on the open
chain at every N ≥ 3, on the ring only at 4 | N (PROOF_CEILING_GOLDEN_ROUTER §6), named under the open items. F169
decides Heisenberg with a field on every connected graph, the control in the paragraph on light along X or Y.
-->

## How to run them

Both run from the repository root. The first `dotnet run` builds; later runs can skip the build with `--no-build`
placed before the `--`, as in `dotnet run --project compute/MirrorWorld -c Release --no-build -- sieve 3 6 Z`.
The grammar: `sieve [terms 1..9] [maxN 3..10] [jumps] [bonds | fields | windows]` for a family (bonds by default),
`sieve cells MAXN TEMPLATES [JUMPS] [NAME]` for one cell, and `sieve cells MAXN @FILE` for a list.

| What | Command | Time on 2026-11-20 |
|------|---------|--------------------|
| the first census: every cell of one to three bond terms, signs ±, Z light, N up to 6 | `dotnet run --project compute/MirrorWorld -c Release -- sieve 3 6 Z` | 8 s |
| the window family: cells of one or two three-site templates with letters at both ends | `... -- sieve 2 6 Z windows` | 90 s |
| the field family: bond cells of one or two terms with every uniform field from X, Y, Z | `... -- sieve 2 6 Z fields` | 35 s |
| one cell of one's own | `... -- sieve cells 9 XXX:1,XYX:1 Z` | seconds |
| a list of cells, one per line (a line starting with `#` is skipped) | `... -- sieve cells 9 @cells.txt` | |

A family run prints, per geometry, the count of each verdict, how many broken cells are broken for every N and how
many mixed ones are broken past their last pairing, one line for each Mixed, Residue or Law cell (a Law line with the law's text),
and, where there are any, the count of read ring cells with broken letters whose every N a ring law decides, with one
`Read+law` line each (letters and the law's text); its tail "broken from N = b on" is the cells mode's "broken for every N from b". Family names such as
`+XZ-ZY` or `+XZ-ZY|+X-Z` (bond terms, then fields) are the cells mode's `XZ:1,ZY:-1` and `XZ:1,ZY:-1,X:1,Z:-1`.

A cell is written as templates with whole-number coefficients, `XZ:1,ZY:-1`, each template laid on every run of
consecutive sites as long as itself, on the chain and on the ring; a one-letter template is a field on every site,
and an inner `I` is a site the template skips (`XIY:1` reaches three sites). A template padded with `I` at an end is
read without the padding, so on the chain its end placements differ from the padded template's; the program prints a
note when that happens. Each letter after the cell is a jump on every site at rate 1 (`Z`, light along Z, by default;
`XZ` would put two jumps on every site); text after the jumps names the cell. In code, `Sieve.Run(name, cell, maxN)`
returns one `Entry` with a `Pass` per geometry, and `Stitch.EveryN(cell, ring)`, `Stitch.Decide(cell, n, ring)`,
`Stitch.LocalHolonomy(cell)`, `Stitch.FreeEveryN(cell, ring)` and `Stitch.RingLaw(cell)` ask the stitch's readers one
by one; these read light along Z only, and `Sieve.TurnedToZ(cell)` gives the turned cell that `Sieve.Run` hands them
under X or Y light.

## Reading a line

Each geometry prints a verdict and one letter per length, the chain from max(3, the cell's span) sites, the ring from
the smallest ring the cell fits on:

| Letter | Who answered | Holds |
|--------|--------------|-------|
| E | the Theodolite's every-N reading (a map built one site at a time) | at every length its construction covers |
| S | the stitch reads this length, and its solve at the window reads every length from the window on (below the window, by a solve at this length) | as said |
| s | the stitch's solve at this length, with no reading at the window | at this length |
| F | the free chain | from its first length on (on the ring, where no listed modulus divides N) |
| C | the stitch's no: a condition of Theorem 2 fails | exact; with "broken for every N from" it holds from there on |
| P / B | F158's count, palindromic / broken, exact | at this length |
| p / b | F158's count by its ranks alone: the exact check did not close, the bounds point this way | at this length, as evidence |

E, S and F count as covered. Among the uncovered letters s, P and p count as palindromic, and B, b and C as broken.
The suffix ", theodolite" says that the Theodolite's reading exists for the cell, "free from n" gives the free chain's
first length and "off d | N" the lengths it leaves out on the ring; on the ring the text of a law of F173, where one
exists, follows, with "every N decided" where it decides every length and "not every N decided" where it does not (`Pass.DecidedByLaw`: its guard lies among the
lengths examined or one past them, and every letter below the guard is exact).

The verdicts, each about the lengths examined. **Read**: a reader covers some length and every uncovered length is
broken; or the free chain covers every length from its first on and every length below it is an exact yes.
**Broken**: nothing covers and every length is broken. **Mixed**: some uncovered lengths palindromic, others broken.
**Residue**: every uncovered length palindromic (and not read by the free chain's tail). **Law**: a ring cell that would be Mixed or a Residue, whose tail
the stitch's no has not already decided (otherwise it keeps its verdict), and whose every length from the law's guard
on is decided by a ring law of F173, every length below the guard examined exactly. Wherever a law exists the sieve
checks every exact letter from its guard on against it, and a disagreement throws. Past the lengths examined a verdict
says only what its readers and laws say.

One line decoded, the cell X·(X + Y)·X, that is `XXX:1,XYX:1` under Z light, N up to 9:

```
  XXX:1,XYX:1 (Z)
    chain Read     N from 3: sPPFFFF, free from 6
    ring  Law      N from 3: PssBPFB, free from 8 off 3 | N, F173 Corollary E2: one flip shape, read through gcd(f, x^N + 1), f = 1 + x + x^2, from N = 5, every N decided
```

On the chain the stitch solves N = 3 alone (s), F158 finds N = 4 and 5 palindromic (P), and the free chain covers
every length from 6 on (F); every length below 6 being an exact yes, the cell is read at every length. On the ring:
N = 3 palindromic by the count, 4 and 5 by the stitch at that length, 6 broken, 7 palindromic, 8 covered by the free
chain (from 8 on, except where 3 divides N), 9 broken. The verdict Law comes from F173's Corollary E2, which
`Stitch.RingLaw` returns for this cell: one flip shape, f = 1 + x + x², read through gcd(f, x^N + 1), from N = 5 on,
and it pairs exactly where 3 does not divide N. The free chain's argument, that every component is a cube of flips
whose cycles are generated by squares, is the one [SOFTNESS_IS_N_DEPENDENT](SOFTNESS_IS_N_DEPENDENT.md) makes by hand
for its crossing.

## What the censuses hold

**The first census** (`sieve 3 6 Z`, 417 cells, all pinned in `SieveTests.cs`). Chain: Read 149, Broken 268, every
broken cell broken for every N from 3, nothing between. Below three sites the sieve does not look, and there a cell
can differ: XX + XY pairs on two sites (its two components are single edges of equal moduli), though it is broken
from three on. Ring: Read 161, Broken 224, all broken for every N, Mixed 32, no residue. The 32 mixed ring cells pair
on the triangle alone (letters `PCCC` and `sCCC`, 16 each) and are broken for every N from 4. Without the stitch as
second pass (`Sieve.Run(..., stitch: false)` in `SieveTests.cs`) the ring's residue is XZ ± ZY and YZ ± ZX with or
without ZZ, twelve cells: broken on every chain examined, paired on every ring examined by a far element that is no
product (pinned for all twelve at N = 3 to 5 and for +XZ+ZY at 6 too; XXX + YYY at three sites for +XZ+ZY), the two
families carried onto each other by CZ on every ring bond up to the sign of H. With the stitch they read `SSSS`. The
carrier W is one of the term-permuting Cliffords that [the colouring hub](THE_PALINDROME_AS_A_COLOURING.md) searches
for at one length (`EndCount.CliffordSymmetry`), read here at every length. On the 36 two-term cells with both signs +
the chain splits 22 read and 14 broken, the counts of [NON_HEISENBERG_PALINDROME](NON_HEISENBERG_PALINDROME.md), the split that
[TWO_TERM_PALINDROME_KLEIN_ROUTING](TWO_TERM_PALINDROME_KLEIN_ROUTING.md) reads at N = 3, 4, 5 and
[V_EFFECT_FINE_STRUCTURE](V_EFFECT_FINE_STRUCTURE.md#longer-chains) holds at every N ≥ 3 by F158's two certificates.

Of the 161 read ring cells, 13 are read on even rings only and broken at the odd lengths examined, N = 3 and 5,
without "broken for every N" (all 417 cells through `sieve cells 6`): XY, YX, XY ± ZZ, YX ± ZZ, XY − YX,
XY − YX ± ZZ, XX ± XY − YY and XX ± YX − YY. For the six with a single string per flip set, `Stitch.RingLaw` returns
F173's Theorem D from N = 3 on, which breaks every odd ring; the verdict stays Read and the pass says so (`DecidedByLaw`), and the family mode counts
these cells and lists each as `Read+law`. For the other seven it
returns none. The three pure Dzyaloshinskii-Moriya cells XY − YX (± ZZ) are decided by Theorem 8 of the complement
connection, live as `PlanarHoppingWitness` (the inspect root `planarhopping`) but not asked by the sieve: pure DM
bonds pair exactly on bipartite graphs, so every odd ring breaks. For XX ± XY − YY and XX ± YX − YY this page names no
law.

**The window family** (`sieve 2 6 Z windows`, 1296 cells: one or two of the 36 three-site templates with a letter at
each end and any middle, signs ±). Chain: Read 673, Broken 95, Mixed 528, no residue; every broken cell is broken for
every N from 3, and every mixed cell pairs at N = 3, or at N = 3 and 4, and is broken for every N past its last
pairing. Of the 673 read, 8 are X·(X ± Y)·X and its relatives, read by the free chain (pinned in `SieveTests.cs` as
`sPPF`). Ring: Read 657, Broken 567, Mixed 48, Law 24, no residue; the broken and the 48 mixed past their last pairing
are broken for every N. The 24 law cells are decided by F173: 16 with one string per flip set by Theorem D (pairing
exactly at even N), 8 relatives of X·(X ± Y)·X by Corollary E2 (from N = 5 on breaking exactly where 3 divides N; all
8 pair on the ring of three).

As in the first census, 56 of the 657 read ring cells carry broken letters at lengths no reader covers, without
"broken for every N" (all 1296 cells through `sieve cells 6`). `Stitch.RingLaw` returns a law for 48 of them (Theorem D
for 40, Corollary E2 with f = 1 + x² for 8), each from N = 5 on, which decides those lengths from 5 on, and the family
mode lists them as `Read+law`; N = 3 and 4 are decided by their exact letters. For 8 it returns none. Of those, XIY − YIX is pure DM on the graph of
bonds two sites apart, and Theorem 8 decides it: it pairs exactly where 4 divides N (at N = 4 its two placements on
each pair cancel and H vanishes), the letters reading `BEBB` at N = 3 to 6. For XZX + YZY and the six pairs with sign +
of XXZ, YYZ, ZXX and ZYY this page names no law.

**The field family** (`sieve 2 6 Z fields`, 2187 cells: every bond cell of one or two terms with every uniform field
pattern from X, Y and Z, coefficients 0 and ±1). Chain: Read 173, Broken 2014, nothing between, every broken cell
broken for every N from 3. Ring: Read 177, Broken 2010, nothing between, every broken cell broken for every N. Of the
177 read ring cells, 7 carry broken letters (all 2187 cells through `sieve cells 6`), and they are the field-free
cells XY, YX, XY ± ZZ, YX ± ZZ and XY − YX of the first census's list, with the same laws (six decided by Theorem D,
XY − YX by none); no cell with a field is among them.

**Light along X or Y** (`sieve 2 6 X` and `sieve 2 6 Y`, 81 cells each: every bond cell of one or two terms). Light
along one letter other than Z on every site is turned to Z for the stitch's readers by a proper rotation of the
letters (`Sieve.TurnedToZ`: X ↔ Z with Y → −Y, or Y ↔ Z with X → −X, a rotation by π that leaves the spectrum as it
is), while the Theodolite and F158's count read the cell as it is. Under either light: chain Read 53, Broken 28; ring
Read 57, Broken 24; no residue, every broken cell broken for every N, and six read ring cells decided at every N by
Theorem D. Every cell answers at every length up to N = 5 as its turned cell under Z (pinned in `SieveTests.cs` for both
lights). Among them, ten cells the Theodolite reads only turned (XY ± XZ, YX ± ZX and YZ + ZY under X light, and their
five images under Y light); the sieve reads them through the stitch on the turned cell, and with the stitch's readers left out each
is a residue on the chain or the ring. XZ + ZY under X light, which is ZX − XY under Z, reads `CCCCCCC` on chain and ring up to N = 9
(`sieve cells 9 XZ:1,ZY:1 X`), broken for every N from 3.

## Where the reading stops

- The stitch, its fast no, the free chain and the ring laws read light along Z on every site; light along one other
  letter on every site is turned to Z for them first. Under jumps of several letters on every site the sieve asks only
  the Theodolite and F158's count, so a cell broken there is broken at the lengths examined and nothing more.
- The sieve reads jumps that are single Pauli letters (X, Y or Z) at rates of zero or more; F158's count does not
  depend on positive rates, a rate of zero drops the jump, and a negative rate is refused.
- F158's count is asked up to ten sites (`maxN` from 3 to 10), and the chain from three sites on. The stitch's window
  and one site more must fit in 16 sites (spans up to 3 on the chain, up to 5 on the ring); past that the window's solve and
  the solves below it are not asked, while the fast no (up to span 9), the free chain (up to span 7) and the ring laws still are. The fast no reads two of Theorem 2's three conditions, the sizes and the loops (the diagonal condition is
  read in the window's equations), on stretches of up to 12 sites, carrying the ratios modulo a prime, so a loop it
  misses is a silence, never a yes. A ring law holds from its guard on; the code evaluates Corollary E2 up to N = 63
  and Theorem D at any N.
- A cycle running round the whole ring is not local, and neither the fast no nor the window's solve sees it; F173's
  ring laws read those where they apply.
- A ring cell marked Read is a yes at the lengths its reader covers; at the other lengths it is what the letters say,
  and past them what its law says where the pass says every N decided.
- A silence is not a no, and a residue is not yet a statement about every length.

## Where they can grow

None of these is built.

1. **The Theodolite on the turned cell.** Under X or Y light the Theodolite reads the cell as it is, and ten cells of
   the two-term family it reads only turned to Z; the stitch reads them instead. Asking the Theodolite on the turned
   cell as well would let its every-length reading cover them.
2. **Theorem 8 as a reader.** Hopping bonds with a DM part and planar fields have absent edges (both entries of an
   edge zero), where the free chain and the ring laws are silent, and the carriers Theorem 8 builds for DM bonds turn
   from site to site (or close only round the whole ring), which the stitch's uniform factors cannot give. Theorem 8 of the complement connection
   decides them at every length, and its live witness already exists. The three pure DM cells named under the first
   census and XIY − YIX in the window family would be decided by it.
3. **F173's chain laws.** Theorems B and B′ of F173 decide cells with one string per flip set on the open chain at
   every length, three templates and more included; `SignedMaskBalanceClaim` holds them and second implementations of
   Theorem D and Corollary E2 beside `Stitch.RingLaw`. The sieve asks none of them.
4. **The next residue.** Every residue these families have shown became a reader: XZ ± ZY the stitch, X·(X ± Y)·X
   the free chain. Wider families (three templates, longer windows, fields beside windows) are where the next one
   would show. F116's golden window XZX + XZY + YZX, palindromic on the open chain at every N ≥ 3 and on the ring only
   where 4 divides N ([PROOF_CEILING_GOLDEN_ROUTER](../docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md) §6), is a
   three-template cell whose two verdicts are known in advance.
