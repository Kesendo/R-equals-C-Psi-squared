# The channel both ways on a path: what a flight could test (a design record, not flown)

**Status: designed, not flown (2026-10-13).** The design sends a Bell pair from either end of a five-qubit line and compares what arrives, to take F119 §7 onto hardware. On this geometry every outcome it could register is an identity of the circuit, a named breach of it, or a reading of the line's T1 profile, and Tom decided not to fly it. This page records why.

**What stands.**
- **The tied entries.** F119's corollary ties two entries of the far pair between the two ends, the |11⟩ population p11 and the |00⟩⟨11| coherence c.
- **The transposed arm keeps their moduli exactly equal** under any static number-conserving coherent error that stays on the bond its gate acts on. Such a gate is its own transpose up to a relative phase of its two sites. On a path, which is a tree, those phases glue into one phase per site.
- **Z dephasing** is a mixture of diagonal unitaries, each its own transpose. It keeps the identity where the way back is the exact reverse of the way out: the layers and the bonds within each layer in reverse order, each error before its gate. The gate's idealized placement (noise around whole gates) is that reverse. Its native placement (idle as late as possible, pulse sequences that are not palindromes) is not, and dephasing placed after every gate in both directions, as `--limits` places it, breaks the identity outright.
- **Amplitude damping** is the device model's one non-unital noise, and the law's condition is linear in F84's net-cooling vector Δγ_l = γ↓_l − γ↑_l: L − Lᵀ and L(I) both are.
  - Its jumps cannot feed p11 and c (F137's block-triangularity). Under limit 3's conditions its no-jump part damps them alike from both ends.
  - In the silence it reaches the one reading the corollary ties there, the silence-read Φ⁺ fidelity in its unnormalized form, which the identity below lets break only where Δγ ≠ 0.
    - **Uniform T1.** Every element number-conserving, every site decaying at one rate over the same intervals, the way back the exact reverse: then a uniform cold T1 leaves the fidelity tied, by counting.
    - **Non-uniform T1.** It moves the fidelity. On the line's profile that is one to two standard deviations with 10 µs delays between the layers, and about three to six with 20 to 40 µs. The shot noise here carries no readout error and not the design's dispersion margin of 1.3, which would lower these by that factor.
    - **Phase errors** of the assumed sizes move it seven to seventeen times less there.
  - The other silence readings, the probability of the silence and the vacuum population alone, are not readings the corollary ties. The probability of the silence moves even under a uniform T1. The vacuum population is tied there by counting but moves under a non-uniform T1, and under unital noise that does not conserve the number.

So any departure the flight could register belongs in advance to a named class:
- the non-uniform part of the day's T1 profile, in the silence-read fidelity, readable with long delays; that profile comes from a properties pull and drifts;
- a breach of the circuit identity:
  - an error that does not conserve the number;
  - an error the way back does not run in exactly reversed order;
  - a bond carrying two different gates;
  - noise between native pulses;
  - asymmetric preparation or readout;
  - drift.

The one class that carries the law's condition could be read cleanly with long enough delays. But what it reads, at a cold bath, is the non-uniform part of a drifting calibration, a population-decay asymmetry present whether or not the way back realizes the transpose. So a flight would check the T1 model, not the transposition. Within the device model's family the identity closes the rest: every departure that carries the law's condition comes from the Δγ profile. What stays open is a design outside that family whose outcome bears on the condition. We do not know one.

*Thomas Wicht and Claude.*

**Gate:** [`simulations/channel_both_ways_gate.py`](../simulations/channel_both_ways_gate.py), the single source of the circuits, the native noise model, the estimator and every number below. Its outputs are in [`simulations/results/channel_both_ways/`](../simulations/results/channel_both_ways/):

- `gate_output.txt`, the design's predictions under the device model (Parts 1 to 3);
- `limits_output.txt`, what the design can and cannot see (`--limits`);
- `search_output.txt`, the seeded parameter search;
- `kingston_lines_2026-10-13T0550Z.json`, the device inputs: an ibm_kingston properties pull, properties of 2026-10-13 07:30 +02:00, line [74, 73, 79, 93, 94], T1 = 295, 185, 280, 305 and 137 µs.

The hardware runner stays in the external pipeline and is not part of this record.

**Theory:** F119 §7, readings 3 and 4, in [PROOF_ANTILINEAR_TRIANGLE.md](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md); F84 and F137 in [ANALYTICAL_FORMULAS.md](../docs/ANALYTICAL_FORMULAS.md).

---

## What the repo already holds

Swept on 2026-10-13 by store. One store-sweep agent went over the markdown stores, the registries and the flight tooling. Review rounds rechecked them, and the load-bearing lines were re-read at source.

- **[`docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md`](../docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md) §7** holds the law and both readings.
  - **The transpose law:** the transposed generator is the transposed system's generator exactly when L is unital (line 129). For Pauli-string jumps, Lᵀ = L_{Hᵀ} + D (line 135).
  - **The breakers** (line 189): "an added σ⁻ breaks the law; σ⁻ with σ⁺ at equal rates keeps it, at unequal rates breaks it".
  - **Reading 3 (line 155, its Bell pair at line 159).** The Bell pair's p11 and c are transfers the corollary fixes. Its |00⟩ population holds a return probability, which agrees on the XX chain (free fermions) and differs on uneven Heisenberg and XXZ chains. Where the readout restricted to what the dynamics reaches is the mirrored preparation, the operator swap covers it.
  - **Reading 4.** "The silence is what makes the readout the mirrored preparation". Read with the rest projected on |0…0⟩, the pair and the probability of the silence are the same from either end, under a number-conserving real H and Z dephasing.
  - **The rigidity list closes on a path:** the edge block's bipartite conjugation and the open arc `site_resolved_vacuum_block`, "the reason that a path is a tree and its hop phases carry no flux".
  - **The flux ring:** "the swap is exact against the reversed ring and broken against the same ring".
- **[`docs/ANALYTICAL_FORMULAS.md`](../docs/ANALYTICAL_FORMULAS.md)**
  - F119 carries its Verified line and no hardware row.
  - **F84, thermal amplitude damping.**
    - Its correction has the closed form √(Σ_l (γ↓_l − γ↑_l)²)·2^(N−1), where "Δγ_l = γ_↓_l − γ_↑_l is the *net* cooling rate at site l".
    - Every Pauli-string dissipator contributes zero (the Pauli-Channel Cancellation Lemma).
    - [NOISE_ASYMMETRY_SYMMETRY_SCALAR.md](../docs/outbound/NOISE_ASYMMETRY_SYMMETRY_SCALAR.md) (lines 125 to 143) lists that scalar's properties: zero for every Pauli channel, independent of the Hamiltonian, blind to T2, degree one in Γ↓ − Γ↑, independent of temperature for a linear bath.
  - **F137, the T1 half-shift.**
    - "σ⁻ρσ⁺ lowers bra and ket together, so L is block-triangular in the joint popcount". Its diagonal blocks act through H_eff = H − (i/2)Σγᵢnᵢ, and the triangularity is pinned in `compute/RCPsiSquared.Core.Tests/BlockSpectrum/T1BreakingInformationalTests.cs`.
    - Its thermal extension adds the heating channel σ⁺, which raises from below.
- **The typed layer.**
  - MirrorWorld's `AntilinearTriangle.SilenceReadingResiduals` pins reading 4 exactly through order t¹² at N = 5.
  - `AntilinearTriangleClaim` (Tier 1 derived) holds the triangle, with no node for §7.
  - `compute/RCPsiSquared.Core/Lindblad/NoJumpGenerator.cs` holds G_H for H = A + iB, the object F137's diagonal blocks are.
  - `F84ThermalAmplitudeDampingPi2Inheritance` (Tier 1 derived, `compute/RCPsiSquared.Core/Symmetry/`) is the typed home of F84's net-cooling vector.
  - The Diagnostics witnesses hold nothing on reciprocity.
- **`experiments/`, the prior flights.**
  - No flight has sent a payload from one end of a chain and then from the other.
  - [IBM_CONCENTRATOR_RELOADED.md](IBM_CONCENTRATOR_RELOADED.md) (ibm_kingston, N = 5) names "the mirrored-end repeat" as its open follow-up (lines 398, 408), as does [CONCENTRATOR_AB_MECHANISM_TEST.md](CONCENTRATOR_AB_MECHANISM_TEST.md) (line 76).
  - [IBM_RECEIVER_ENGINEERING_SKETCH.md](IBM_RECEIVER_ENGINEERING_SKETCH.md) read the end sites (0, 4) of a five-qubit path by tomography, the nearest instrument.
  - The post-selecting designs ([the staircase null test](STAIRCASE_NULLTEST_HARDWARE_PREDICTION.md), [the corner beat](CORNER_BEAT_HARDWARE_PREDICTION.md), [record parity](RECORD_PARITY_HARDWARE_PREDICTION.md)) flew invalid or not at all.
  - The staircase's first flight is the direct precedent for a T1 profile that moves within one session. Its two spectators' T1 ran 160.7 → 172.3 µs and 132.9 → 115.1 µs between the batch's start and end (line 147).
  - [PRICE_PAIR_HARDWARE_PREDICTION.md](PRICE_PAIR_HARDWARE_PREDICTION.md) run 4 measured the ZZ crosstalk the gate uses as a planning size (−3.9 kHz, ibm_marrakesh).
- **The registries.**
  - `fw.Confirmations` and `ConfirmationsRegistry.cs` hold 24 entries, none a two-direction reading.
  - Among them is `f84_heating_leg_attribution_kingston_july2026` ([F81_VIOLATION_HARDWARE_BRIDGE.md](F81_VIOLATION_HARDWARE_BRIDGE.md)).
    - **What it flew:** a two-legged idle T1 protocol on ibm_kingston qubits [82, 83, 13], repeated on ibm_marrakesh [93, 94, 95], none of them on this design's line.
    - **What it found:** the up-legs flat, so the bath is cold. On the Marrakesh repeat, the in-situ γ↓ sat off the same-day calibration by +7, −28 and −22 %, and one qubit's in-situ T1 stood at 129 µs in an earlier fit and at 316 µs on the day of the valid run.
    - **How the document reads the net cooling:** v_l = d⟨Z_l⟩/dt at I/2^N equals γ↓,l − γ↑,l, and its root-sum over sites is "the RMS identity-escape velocity" (line 39).
  - The OpenArcs registry holds `the_gate_that_does_not_gate` (the March work's bidirectional bridge, named in §7 line 121) and `site_resolved_vacuum_block` (above). No arc covers reading 4.
- **[`docs/CAUGHT_ERRORS.md`](../docs/CAUGHT_ERRORS.md)**
  - **2026-08-04** (a Liouvillian that ran H backwards): the gate pins the direction by comparing complex density matrices.
  - **2026-08-05 and 2026-08-18** (the T₂ → γ trap): the gate's dephasing rate is γ = max(0, 1/T₂* − 1/(2T₁))/2.
  - This design's own entry is dated 2026-10-13.
- **[`docs/GLOSSARY.md`](../docs/GLOSSARY.md)** has no "silence". Here, as in reading 4, the word names the projection only.

Adjacency, asked as its own question:

- **Reading 3 and this design are one object.**
- **F137's block-triangularity and limit 3 are one object.** §7 names F137 only in its sweep (line 121), for the words "detailed balance"; that F137 decides where σ⁻ can reach a payload is not drawn there.
- **F84's net-cooling vector and the law's condition are one object:** L(I) = Σ_l Δγ_l Z_l. §7 names F84 only in the same sweep line.
- **The counting law of limit 4 has the shape of a common mode the repository already names.** One is the absorption theorem's collapse to the count ([PROOF_ABSORPTION_THEOREM.md](../docs/proofs/PROOF_ABSORPTION_THEOREM.md) lines 133 to 137: "a function of the count alone … The collapse is where the uniform hypothesis is spent"). The other is [CMRR_BREAK_NONUNIFORM_GAMMA.md](CMRR_BREAK_NONUNIFORM_GAMMA.md) (line 149). There it is dephasing; here it is T1.
- **The no-flux sentence and limit 1 are one object.** It closes §7's rigidity list, which reads it for the one-excitation block, not for a transposed circuit.

## The design

**Payload and geometry.**
- **Payload:** Φ⁺ on one end pair of a five-qubit line, the rest in |0⟩. The receiving pair is read at the far end.
- **Pairs:** (0,1) and (3,4), site 2 between them.
- **Bond gate:** exp(−ia·h_b), with h_b = (XX+YY)/2 + Δ·ZZ/2 and a = J_b·dt·f, as one synthesized two-qubit unitary. J = (1.3, 2.2, 1.4, 1.4), dt = 0.5 (J_b·dt from 0.65 to 1.1), steps m ∈ {3, 5}.
- **Layers:** E = bonds (0,1), (2,3); O = bonds (1,2), (3,4).

**The far pair's entries:** p00, p11 and c = ⟨00|R|11⟩. Under number conservation the habit reading and the silence reading coincide on p11 and c, since two excitations on the pair force the rest empty. They differ on p00: read with the silence it is the vacuum's, read by habit it also holds a return probability.

**Arms:**

| Arm | Out | Back | Δ | What it was to show |
|---|---|---|---|---|
| Main | E/2 (O E)^(m−1) O E/2 | the same | 1.5 | transfer agrees; habit p00 differs |
| Main0 | the same | the same | 0 | transfer agrees; habit p00 agrees (free fermions) |
| Transposed | (E O)^m | (O E)^m | 1.5 | transfer agrees |
| Control | (E O)^m | (E O)^m | 1.5 | transfer differs |

A replica of the main arm checked dispersion, and zero-depth end arms checked preparation and readout. The gated statistics were the transfer D = √(Δp11² + 2Δ|c|²) and the habit |Δp00| between the two directions.

**Under the native device model** (noise between pulses, `gate_output.txt` Part 1):

| Step | Transfer D: main / main0 / transposed / control | Habit \|Δp00\|: main / main0 |
|---|---|---|
| 3 | 6.5·10⁻⁴ / 2.7·10⁻⁴ / 7.2·10⁻⁴ / 0.110 | 0.0587 / 3.7·10⁻³ |
| 5 | 9.0·10⁻⁴ / 2.4·10⁻⁴ / 9.7·10⁻⁴ / 0.062 | 0.1089 / 1.7·10⁻³ |

With the counts-level decision rules of Part 2, the device model passed every planned test with probability 0.97 at the planned dispersion margin. The limits below show how much of that pass is fixed in advance.

## What the design can and cannot see

**Two conditions have to be kept apart.**
- **F119's transpose law** says the transposed generator is the transposed system's generator. It needs a unital L. A coherent error is unital, so it keeps the law.
- **The two ends agreeing** needs more: Lᵀ = L, that is a real H and a dissipator that is its own transpose. A coherent error that is not its own transpose, such as an XY − YX part, breaks this at the level of the generator, because the way back is then no longer the transposed circuit. On the transferred moduli it shows only where a bond carries two different gates (limit 1).
- **In the device model's family** (a real H, Pauli dephasing, thermal amplitude damping) the two conditions coincide, by one identity:

  L − Lᵀ = Σ_l Δγ_l (σ⁻_l·σ⁺_l − σ⁺_l·σ⁻_l),  L(I) = Σ_l Δγ_l Z_l,  Δγ_l = γ↓_l − γ↑_l,

  where A·B is the map ρ ↦ AρB. Both vanish only at Δγ = 0.
- **The check.** `--limits` checks the identity at five sites on a random integer symmetric H with random integer rates and Z dephasing: both residuals are exactly 0. For X and Y dephasing it holds algebraically, every Pauli string being its own transpose up to a sign that squares away. `--limits` prints and asserts nothing; its exact zeros are read off the printout.

The numbers below are in `limits_output.txt` unless a line names `gate_output.txt`.

### 1. The tree gauge: coherent errors that conserve the number and stay on their bond vanish from the transposed arm, exactly

- **The gauge.** A two-site unitary G that conserves the excitation number has off-diagonal entries of equal modulus in its one-excitation block. So Gᵀ = D·G·D⁻¹, with D a relative phase of its two sites. The per-gate residual is at most 4.3·10⁻¹⁶.
- **The gluing.** On a path, which is a tree, the relative phases glue into one global diagonal D, as long as each bond carries one gate, however often it is repeated. Then U_outᵀ = D·U_back·D⁻¹, and the moduli of the transferred entries agree.
- **The transposed arm** carries one gate per bond. Under a fixed random number-conserving error after every bond gate:
  - its transfer D is at most 4.8·10⁻¹⁶ at every size tried (ε = 0.05, 0.15, 0.3);
  - its complex c difference reaches 0.055 to 0.22.
- **An XY − YX part** after every gate is not its own transpose, but it stays on its bond: the transposed arm reads it as 4 to 6·10⁻¹⁶ at θ = 0.05 to 0.2.
- **The main arm** sees the random errors (2.7·10⁻³, 2.8·10⁻², 0.10) and the XY − YX part (1.1·10⁻² to 5.0·10⁻²) only because its E bonds carry two different gates, the half step and the full step.
  - Take the graph whose vertices are the sites and whose edges are the gates. Two gates on one bond make a double edge, the one cycle a path can carry, and its two phases do not glue.
  - An error of the same size, drawn afresh and placed on the O bonds only, reads 4.3·10⁻¹⁶ on the main arm. On the E bonds only, it reads 7.1·10⁻³.
- **Outside the gauge, the way back has to be the exact reverse.** The transpose of a product is the reversed product of the transposes. So the way back the transposition asks for runs the layers and the bonds within each layer in reverse order, each error before its gate.
  - **What the design's transposed arm does.** It reverses only the layers and keeps each error after its gate. For the gates alone that is the same circuit, because the gates within a layer act on disjoint bonds. It differs for errors that the gauge does not cover:

    | Error | Design's way back | Exact reverse |
    |---|---|---|
    | ZZ of 0.05 rad on the neighbouring bonds, static, past the bond | 2.3·10⁻³ | 2.9·10⁻¹⁶ |
    | Z dephasing p = 0.01 on the bond's left site, stochastic, on the bond | 2.4·10⁻² | 4.5·10⁻¹⁶ |

    A stochastic error breaks the gauge although it stays on its bond, because each of its branches is a different gate on that bond, a double edge again. A static ZZ profile after every gate reads 2.2·10⁻³ on the design's transposed arm (`gate_output.txt`).
  - **A ZZ term acting through a whole layer**, as part of its generator, keeps the layer its own transpose, since that generator is real symmetric.
- **On a ring** the hop phases carry a flux that no diagonal removes. That is the flux reversal of §7.

### 2. The free-chain contrast is reading 3, not the transposition

- **What the contrast is.** The main0 arm (Δ = 0) predicts that the habit-read p00 agrees from both ends, because the free fermions tie it. The main arm (Δ = 1.5) predicts that it differs.
- **It hardly moves when the way back stops being the transposed circuit.** Under XY − YX, main0's habit |Δp00| is quadratic in θ and stays more than an order of magnitude below the main arm's:

  | θ | 0.05 | 0.1 | 0.2 |
  |---|---|---|---|
  | main0 habit \|Δp00\|, max over both steps | 1.9·10⁻⁴ | 7.3·10⁻⁴ | 2.5·10⁻³ |
  | main habit \|Δp00\|, m = 3 | 0.0708 | 0.0693 | 0.0630 |
  | main habit \|Δp00\|, m = 5 | 0.1548 | 0.1674 | 0.1832 |

- **Random number-conserving errors** of one seeded draw move it by 9.5·10⁻⁴, 0.019 and 0.040 at ε = 0.05, 0.15 and 0.3. Only at the largest size does it pass its null of about 0.02.
- **σ⁻ loosens the free fermions' tie a little.** T1 alone gives main0's habit |Δp00| 1.1·10⁻⁴ and 9.3·10⁻⁵ in the native placement (`gate_output.txt` Part 1), far below its null of about 0.02.
- So the contrast separates an interacting chain from a free one, which is reading 3's content. It does not separate a channel whose way back is the transposed circuit from one whose way back is not.

### 3. σ⁻'s jumps do not feed the transferred entries

- **The mechanism is F137.** Amplitude damping lowers bra and ket together, so L is block-triangular in the joint popcount. In the identity above, the part L − Lᵀ shifts the joint popcount by one and has no part inside a diagonal block.
- **This payload starts at the top.** It holds at most two excitations, so the (2,2) and (0,2) blocks that hold p11 and c have nothing above them to feed them.
- **What carries them** is H_eff = H − (i/2)Σγᵢnᵢ, with H real symmetric and the damping diagonal, together with the dephasing jump term in the same blocks. Both are their own transpose. T1 does damp p11 and c through H_eff, alike from both ends.
- **The conditions.** The unnormalized p11 and c then agree from both ends at any T1 profile, as long as:
  - every other element conserves the number;
  - the way back is the exact reverse of the way out;
  - the bath is cold. A heating channel σ⁺ raises from below and is not covered. F84's flight found the bath cold on other qubits, and this line's own excited population is not measured.
- **The check.** `gate_output.txt`'s stress row puts T1 = 12 to 34 µs and T2* = 5.6 to 11.2 µs, unevenly per qubit, at layer boundaries with no gate error, and reads at most 2.8·10⁻¹⁶.
- **Where it leaks.** In the native placement the sx pulses do not conserve the number, and the idle placed as late as possible is not a palindrome. There T1 alone reads 1.3·10⁻⁴ and 7.1·10⁻⁵ (`gate_output.txt` Part 1).

So the law's own breaker cannot be seen in the two entries the design gated.

### 4. Where σ⁻ shows: the silence

Below the top blocks σ⁻ feeds the one- and zero-excitation parts, so the silence carries what limit 3 keeps out of p11 and c. A delay of τ between the main arm's layers gives σ⁻ time to act, and in the idealized placement it keeps the sequence a palindrome; the native placement is not one to begin with. Of the silence's readings, one is a reading the corollary ties.

**The silence-read Φ⁺ fidelity, F_sil = (p00_sil + p11 + 2 Re c)/2, is the corollary's own observable.** It is unnormalized, the silence in the numerator. The renormalized form that reading 4 uses divides by the probability of the silence, which is not tied once σ⁻ acts, so everything here concerns the unnormalized form.
- **Why it is tied.** Its readout, Φ⁺ on the far pair with the rest silent, is the mirrored preparation. So the operator swap ties it wherever L = Lᵀ and the way back is the exact reverse, number conserved or not.
- **Where it can break.** By the identity above, only where Δγ ≠ 0.
- **What the design gated instead.** Its estimator holds all of F_sil's parts, but the design gated the transferred moduli, the habit-read p00 and the control, not F_sil.

| Model | ΔF_sil |
|---|---|
| Depolarizing and Z dephasing, unital but not number-conserving (idealized placement, m = 3) | ≤ 2.2·10⁻¹⁶ |
| Uniform T1, at both coupling sets (idealized placement) | ≤ 3.9·10⁻¹⁶ |
| The line's T1 on the native device model, τ = 10 µs | −5.5·10⁻³ (m = 3), −8.0·10⁻³ (m = 5) |

**Why a uniform cold T1 leaves it tied, although Δγ ≠ 0.**
- **The counting.** Take every element number-conserving, and one decay rate γ on every site acting over the same intervals, as continuous T1 does. Each excitation then decays at γ wherever it sits, so the populations P_k of the number sectors obey

  Ṗ_k = −kγ·P_k + (k+1)γ·P_{k+1},

  whatever H is and whichever end sent.
- **The payload.** Its number distribution (½ at zero, ½ at two) is the same at both ends, so the vacuum population P₀ = p00_sil agrees.
- **With p11 and c tied (limit 3),** ΔF_sil = ΔP₀/2 = 0.
- **P₀'s tie holds for any number-conserving circuit in which every site decays at the same rate over the same intervals.** A schedule that lets the sites of a bond decay at different moments breaks it. The control arm, whose way back is not the transpose, keeps P₀ equal at uniform T1 as well: ΔP₀ ≤ 1.3·10⁻¹⁵ on the main and control arms at m = 3 and 5, with 3 µs of idle per layer.
- **A non-uniform T1 on a cold bath** moves F_sil through the vacuum population, ΔF_sil = ΔP₀/2. With q4 at 137 µs, ΔF_sil is −6.7·10⁻³ against Δp00_sil = −1.3·10⁻².

**What else moves it, and how far it can be read.**
- **Phases.** F_sil contains Re c, which the tree gauge does not tie. The device's phase errors are not measured. At assumed sizes, with no other noise and no delay, ΔF_sil reads:
  - **Rz(0.02) on the bond's left qubit after every gate:** −4.8·10⁻³ and +1.1·10⁻² on the main arm, −7.3·10⁻³ and +1.9·10⁻² on the transposed arm (m = 3, 5);
  - **a static Z profile (sd 0.05 rad):** up to 1.0·10⁻² on the transposed arm.
- **At the signal's operating point** (τ = 10 µs, the line's T1 and T2*, idealized placement) the delays have damped the coherence the phases act on.
  - On the uniform-T1 background the same errors move ΔF_sil by at most 7.6·10⁻⁴.
  - The T1 signal there is −5.8·10⁻³ and −8.8·10⁻³ (m = 3, 5).
- **Shot noise and delay.** Without readout error, on the native device model:

  | τ | ΔF_sil, m = 3 | ΔF_sil, m = 5 | Shot noise |
  |---|---|---|---|
  | 10 µs | −5.5·10⁻³ | −8.0·10⁻³ | 4.0 to 4.1·10⁻³ |
  | 20 µs | −1.4·10⁻² | −1.6·10⁻² | 4.1 to 4.3·10⁻³ |
  | 40 µs | −2.4·10⁻² | −1.8·10⁻² | 4.2 to 4.3·10⁻³ |

  - With T1 switched off it stays below 3·10⁻⁴ in these rows. By 20 µs the coherence has nearly died away, and F_sil carries the vacuum and |11⟩ parts alone.
  - At the design's τ = 0, ΔF_sil is +5.1·10⁻⁴ with T1 and +6.7·10⁻⁴ without it (m = 3): a T1 part of −1.6·10⁻⁴, far below the shot noise.
- **The model's delays carry no ZZ crosstalk.** At the planning size of −3.9 kHz that is about 0.25 rad of conditional phase per 10 µs. It is real symmetric and palindromic, so it keeps c tied. By changing where the payload dwells it would change the size of the T1 signal in the vacuum population; the model does not compute by how much.

**The other silence readings are not tied by the corollary.** The probability of the silence, P_silence, and the vacuum population p00_sil on its own read the far pair with a readout that is not the mirrored preparation, so the corollary does not pair them between the ends.

**The vacuum population.**
- **What ties it and what moves it.** The counting law ties it under a uniform T1. A non-uniform T1 moves it, and so does unital noise that does not conserve the number: under depolarizing and Z dephasing (idealized placement), Δp00_sil is +1.8·10⁻⁴ and +2.6·10⁻⁴ at τ = 0 and 10 µs (m = 3), while ΔF_sil stays at 10⁻¹⁶. That tie is what F_sil has over it.
- **In the cold, number-conserving regime** ΔF_sil = Δp00_sil/2. Read in the ZZ setting alone, at a shot noise of 5.5·10⁻³, p00_sil carries twice the signal of F_sil, whose shot noise is 4.0·10⁻³.

**The probability of the silence.** Under number conservation it is ½ + p11 and inherits p11's tie. Once σ⁻ acts, it moves even under a uniform T1:
- uniform T1 moves P_silence by +2.6·10⁻³ (m = 3, τ = 10 µs), all of it in the silent |01⟩ and |10⟩ parts, while p00_sil and F_sil stay tied;
- on the native device model, P_silence forward minus backward is:

  | τ | m = 3 | m = 5 |
  |---|---|---|
  | 0 | −1.1·10⁻³ | +1.2·10⁻³ |
  | 5 µs | −1.4·10⁻² | −1.4·10⁻² |
  | 10 µs | −2.6·10⁻² | −2.8·10⁻² |

  With T1 switched off it stays at most 2.0·10⁻³. Its shot noise is 2.5·10⁻³ (binomial at P = ½, pooled over the five settings, without readout error).
- its sign and size follow the T1 profile and the couplings. On a model with only T1 and dephasing (m = 3, τ = 10 µs):

  | T1 profile | ΔP_silence |
  |---|---|
  | Uniform 240 µs | +2.6·10⁻³ |
  | 137 µs at q4 | −2.8·10⁻² |
  | 137 µs at q0 | +3.1·10⁻² |
  | 137 µs at q2 | +1.0·10⁻² |
  | Palindromic (137, 240, 240, 240, 137) | +6·10⁻⁴ |

  Uniform T1 on the palindromic couplings J = (1.4, 2.2, 2.2, 1.4) reads −5.4·10⁻³.

### What follows

- **On a path with number-conserving gates, every agreement F119 predicts that this flight could measure is an identity of the circuit**, as long as the way back is the exact reverse of the way out. The design's transposed arm comes close: it ties the moduli under static number-conserving coherent errors on the bond, by the gauge.
- **A flight can break the identity only through what it does not cover:**
  - errors that do not conserve the number: depolarizing, and coherent errors inside the native synthesis, which need not conserve it;
  - errors the way back does not run in exactly reversed order: stochastic ones such as dephasing, and ones that reach past their bond;
  - a bond carrying two different gates;
  - noise placed between native pulses that are not palindromes;
  - asymmetric preparation or readout;
  - drift.
- **Phase errors** reach F_sil through its coherence even where the moduli gauge out.
- **The law's own condition shows only through σ⁻**, and through σ⁺ were the bath warm. Among the tied readings it shows in F_sil, and in the free-chain contrast far below its null. A uniform cold T1 leaves F_sil tied by counting. A non-uniform one moves it by what the day's T1 profile sets, readable with long delays.
- **A loop would bring back a flux that no gauge removes.** But reciprocity against the reversed flux is again an identity of the circuit, and the shortest loop on heavy-hex has twelve qubits.

Each departure this flight could register therefore belongs in advance to a named class. The one that carries the law's condition could be read cleanly with 20 to 40 µs delays. But the readable signal is ΔP₀/2, a population-decay asymmetry that the non-uniform part of a cold, drifting T1 profile sets, present whether or not the way back realizes the transpose. So a flight would check the T1 model, not the transposition. Within the device model's family every departure that carries the law's condition comes from the Δγ profile, by the identity. A design outside that family whose outcome bears on the condition is what stays open. We do not know one.

## What remains usable

The gate runs without QPU access:
- **its circuits;**
- **its native model:** each bond gate as its synthesized pulse sequence, idle as late as possible, decay and depolarizing per pulse;
- **its estimator:** per-qubit confusion pooled from CAL0 and CAL1, inverted as a tensor product, the silence in the numerator;
- **`fidelity_silence`;**
- **`--limits`.**
