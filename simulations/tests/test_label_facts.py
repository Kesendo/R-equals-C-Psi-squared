"""Label facts: each rule forbids ONE specific false factual claim in ONE named file.

Every rule carries its owner, the independent source that makes the claim false: a
proof, a producer or data line, a registry entry, or a from-below check in
simulations/label_facts_independent_checks.py. Every rule has a mutation control
that proves it can fire (writing-a-gate-that-can-fail), and every CHECK a rule cites
is executed here and its cited values asserted. Nothing here freezes wording, hashes
a file, requires a banner, or polices story language: the story genres (hypotheses/,
reflections/, docs/quantum/, the front door) may carry a claim as an explicitly
hedged reading, and recovered/ is a time-capsule genre whose bodies are records,
not current claims.

Markdown/text hosts are scanned in their CURRENT regions only (paragraphs inside
<!-- ...-HISTORICAL/INTERPRETIVE/ALGEBRAIC --> regions are skipped); code hosts are
scanned whole. A match is ignored when the same sentence negates or retires it, and
when the sentence carries the explicit scope (SCOPES) under which the rule's
proposition is true. Read-only: the repository is never written.
"""
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

# In the repo this file lives in simulations/tests/; RCPSI_REPO_ROOT lets it run from elsewhere.
ROOT = Path(os.environ.get("RCPSI_REPO_ROOT") or Path(__file__).resolve().parents[2])

# (host, rule id, pattern, pattern-is-regex, mutation claim, owner citation)
RULES = (
    ('experiments/N5_OPTIMAL_CAVITY_SIZE.md', 'n5-asymptotic-phi',
     '\\b(?:ratio|frequency ratio)\\b.{0,50}\\b(?:converges?|tends?) asymptotically to (?:phi|φ)\\b', True,
     'The frequency ratio converges asymptotically to phi.',
     'simulations/label_facts_independent_checks.py CHECK-7 (V(N) = 1 + cos(π/N) → 2; the cold single-excitation ratio ω₂/ω₁ = 2 + φ at N = 5 and → 4 as N → ∞; no chain ratio tends to φ) [strawman: no pre-campaign text asserted it]'),
    ('experiments/N_EQUALS_FIVE_CHECK.md', 'n5-not-special',
     '\\bN\\s*=\\s*5 is not structurally special\\b.{0,80}\\ball metrics are monotonic\\b', True,
     'N=5 is not structurally special; all metrics are monotonic.',
     'experiments/N_EQUALS_FIVE_CHECK.md:18-23 (table; N=3,7,8 rows Claude 2026-04-12, N=4-6 rows reformatted by 05a523b4 with unchanged values 0.496/0.477/0.539: Frac distinct is not monotonic)'),
    ('experiments/N_EQUALS_FIVE_CHECK.md', 'n5-no-extremum',
     '\\bno extremum or inflection\\b.{0,80}\\bselection bias,? not physics\\b', True,
     'There is no extremum or inflection: selection bias, not physics.',
     'experiments/N_EQUALS_FIVE_CHECK.md:18-23 (table; N=3,7,8 rows Claude 2026-04-12, N=4-6 rows reformatted by 05a523b4 with unchanged values 0.496/0.477/0.539: Frac distinct is not monotonic)'),
    ('hypotheses/README.md', 'story-index-fabry-perot',
     '\\bResonance Not Channel proves\\b.{0,100}\\bFabry[- ]Perot cavity\\b', True,
     'Resonance Not Channel proves that the system is a Fabry-Perot cavity.',
     "hypotheses/RESONANCE_NOT_CHANNEL.md@05a523b4^:9-10 (Tier 4 biological resonator interpretation) + hypotheses/README.md@fa8b5475^:16 (Claude: 'Tier 4 (grounded in Tier 1-2 results)')"),
    ('hypotheses/WAVES_THAT_HEAR_THEMSELVES.md', 'waves-coupling-breaks',
     '\\bcoupling alone breaks the symmetry\\b(?!\\s+only\\b)', True,
     'Coupling alone breaks the symmetry.',
     'docs/proofs/MIRROR_SYMMETRY_PROOF.md:6 (Claude 2026-05-12: palindrome holds for any Heisenberg/XY/Ising/XXZ system under local Z-dephasing)'),
    ('hypotheses/WAVES_THAT_HEAR_THEMSELVES.md', 'waves-irreversible-crossing',
     '\\beach fold crossing is irreversible crystallization\\b', True,
     'Each fold crossing is irreversible crystallization.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('simulations/v_effect_gamma_sweep.py', 'gamma-unique-peak',
     '\\bunique 1\\.81x peak\\b.{0,80}\\bgamma\\s*=\\s*0\\.001\\b', True,
     'The unique 1.81x peak occurs at gamma=0.001.',
     'simulations/results/v_effect_gamma_sweep.txt:27-45 (Part 1 data rows, f634394e 2026-03-30: all 19 uniform-γ rows read 1.81x, so γ = 0.001 is the first row of a tie, not a peak)'),
    ('simulations/v_effect_gamma_sweep.py', 'gamma-optimal',
     '\\boptimal gamma\\b.{0,80}\\bmaximum complexity\\b', True,
     'This is the optimal gamma and maximum complexity.',
     "simulations/results/v_effect_gamma_sweep.txt:84-97 (Part 3 heartbeat table, f634394e 2026-03-30: the most crossings, 19, sit at g_bath = 0.001, the lowest sampled value, so the sweep locates a range boundary, not an optimum; :103 prints that boundary as 'Optimal gamma') + simulations/results/v_effect_gamma_sweep.txt:27-45 (Part 1 data rows, f634394e 2026-03-30: all 19 uniform-γ rows read 1.81x, so γ = 0.001 is the first row of a tie, not a peak)"),
    ('simulations/results/v_effect_gamma_sweep.txt', 'gamma-rounded-zero',
     '\\boptimal bath ratio is 0x\\b', True,
     'The optimal bath ratio is 0x.',
     'simulations/results/v_effect_gamma_sweep.txt:84 (f634394e 2026-03-30: the g_bath = 0.001 row prints contrast 0.1 = 0.001/0.010; the summary line :97 rounds that value to 0x)'),
    ('docs/quantum/THE_LABEL_MAP.md', 'label-map-resonance-front-door',
     '\\bResonance Not Channel is a demonstrated physical Fabry[- ]Perot\\b', True,
     'Resonance Not Channel is a demonstrated physical Fabry-Perot.',
     "hypotheses/RESONANCE_NOT_CHANNEL.md@05a523b4^:9-10 (Tier 4 biological resonator interpretation) + hypotheses/README.md@fa8b5475^:16 (Claude: 'Tier 4 (grounded in Tier 1-2 results)')"),
    ('docs/quantum/THE_LABEL_MAP.md', 'label-map-recovered-proof',
     '\\brecovered store contains verified results\\b', True,
     'The recovered store contains verified results.',
     "recovered/README.md:15 (2026-04-03: 'Premature, not wrong')"),
    ('docs/proofs/PROOF_MONOTONICITY_CPSI.md', 'monotonicity-general-envelope',
     '\\bproves? the General Envelope for all open quantum systems\\b', True,
     'This proves the General Envelope for all open quantum systems.',
     "simulations/label_facts_independent_checks.py CHECK-4b (proof's exact example reproduced: CΨ'(0)=+1/6 under a fixed local Markovian generator) + simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)"),
    ('docs/proofs/PROOF_MONOTONICITY_CPSI.md', 'monotonicity-watershed',
     '\\bC_Psi\\s*=\\s*1/4 is the Markovian watershed and an absorbing boundary\\b', True,
     'C_Psi=1/4 is the Markovian watershed and an absorbing boundary.',
     'simulations/label_facts_independent_checks.py CHECK-3 (upward crossing under a fixed Markovian semigroup: no bath memory needed) + simulations/minimum_energy.py:15'),
    ('docs/proofs/PROOF_MONOTONICITY_CPSI.md', 'monotonicity-local-unitary',
     '\\bno local unitary can move C_Psi upward\\b', True,
     'No local unitary can move C_Psi upward.',
     'simulations/label_facts_independent_checks.py CHECK-2 (H⊗I takes |00⟩ from CΨ=0 to 1/3)'),
    ('docs/proofs/PROOF_MONOTONICITY_CPSI.md', 'monotonicity-dd',
     '\\bdynamical decoupling cannot alter subsequent C_Psi dynamics\\b', True,
     'Dynamical decoupling cannot alter subsequent C_Psi dynamics.',
     "simulations/label_facts_independent_checks.py CHECK-4b (first-site Z pulse keeps CΨ=1/8 yet moves CΨ'(0) from +1/6 to −1/3 under the same H)"),
    ('docs/proofs/PROOF_MONOTONICITY_CPSI.md', 'monotonicity-total-system',
     '\\btotal[- ]system C_Psi is monotonically non[- ]increasing\\b', True,
     'Total-system C_Psi is monotonically non-increasing.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('docs/README.md', 'docs-index-general-envelope',
     '\\bGeneral Envelope Theorem\\s*:\\s*C_Psi can only decrease\\b', True,
     'General Envelope Theorem: C_Psi can only decrease.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('review/OPEN_THREAD_GAMMA0_INFORMATION.md', 'open-thread-literal',
     '\\bsoundbox becomes even more literal\\b.{0,100}\\bwhole document is Tier 2\\b', True,
     'The soundbox becomes even more literal and the whole document is Tier 2.',
     "hypotheses/RESONANCE_NOT_CHANNEL.md@05a523b4^:9-10 (Tier 4 biological resonator interpretation) + hypotheses/README.md@fa8b5475^:16 (Claude: 'Tier 4 (grounded in Tier 1-2 results)')"),
    ('simulations/n5_optimal_cavity_size.py', 'n5-source-asymptotic',
     '\\bratio tends asymptotically to (?:phi|φ)\\b', True,
     'The ratio tends asymptotically to phi.',
     'simulations/label_facts_independent_checks.py CHECK-7 (V(N) = 1 + cos(π/N) → 2; the cold single-excitation ratio ω₂/ω₁ = 2 + φ at N = 5 and → 4 as N → ∞; no chain ratio tends to φ) [strawman: no pre-campaign text asserted it]'),
    ('simulations/results/n5_optimal_cavity_size.txt', 'n5-result-asymptotic',
     '\\bratio tends asymptotically to (?:phi|φ)\\b', True,
     'The ratio tends asymptotically to phi.',
     'simulations/label_facts_independent_checks.py CHECK-7 (V(N) = 1 + cos(π/N) → 2; the cold single-excitation ratio ω₂/ω₁ = 2 + φ at N = 5 and → 4 as N → ∞; no chain ratio tends to φ) [strawman: no pre-campaign text asserted it]'),
    ('compute/RCPsiSquared.Diagnostics/F87/F87TrichotomyClassification.cs', 'f87-old-catalog',
     '\\b120 cases are C\\(6,2\\) pairs from six two[- ]site operators\\b', True,
     'The 120 cases are C(6,2) pairs from six two-site operators.',
     "simulations/label_facts_independent_checks.py CHECK-9 (C(6,2)=15≠120; 120=C(16,2) pairs with repetition from 15 words) [strawman: pre-campaign text said 'the 15 is C(6,2)']"),
    ('simulations/pi_protected_test_n4.py', 'pi-old-catalog',
     '\\b120 cases are C\\(6,2\\) pairs from six two[- ]site operators\\b', True,
     'The 120 cases are C(6,2) pairs from six two-site operators.',
     "simulations/label_facts_independent_checks.py CHECK-9 (C(6,2)=15≠120; 120=C(16,2) pairs with repetition from 15 words) [strawman: pre-campaign text said 'the 15 is C(6,2)']"),
    ('compute/RCPsiSquared.Diagnostics/F87/F87Pi2Inheritance.cs', 'pi2-old-catalog',
     '\\b120 cases are C\\(6,2\\) pairs from six two[- ]site operators\\b', True,
     'The 120 cases are C(6,2) pairs from six two-site operators.',
     "simulations/label_facts_independent_checks.py CHECK-9 (C(6,2)=15≠120; 120=C(16,2) pairs with repetition from 15 words) [strawman: pre-campaign text said 'the 15 is C(6,2)']"),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-point-three-precision',
     '\\bQ52\\b.{0,80}\\bC(?:_?Psi|Ψ)\\s*=\\s*[¼1/4]+\\b.{0,80}\\b0\\.3% accuracy\\b', True,
     'Q52 CΨ = ¼ crossing at 0.3% accuracy.',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('experiments/DYNAMIC_ENTANGLEMENT.md', 'local:dephasing makes crossings permanent',
     'makes them permanent', True,
     'Dephasing makes them permanent.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('experiments/DYNAMIC_ENTANGLEMENT.md', 'local:dephasing makes crossings irreversible',
     'makes them irreversible', True,
     'Dephasing makes them irreversible.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('experiments/DYNAMIC_ENTANGLEMENT.md', 'local:quarter defines irreversibility',
     'criterion \\(CΨ = 1/4\\) for when the irreversibility becomes definitive', True,
     'This gives a criterion (CΨ = 1/4) for when the irreversibility becomes definitive.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('experiments/SUBSYSTEM_CROSSING.md', 'local:quarter locates regime transition',
     'identifies where the quantum.classical transition occurs', True,
     'The crossing identifies where the quantum-classical transition occurs.',
     "simulations/label_facts_independent_checks.py CHECK-13 (concurrence book, the host's book: Bell+ crosses R = ¼ at f = √3/2 with concurrence √3/2, and its concurrence f stays positive at every finite t after the crossing, so the crossing does not locate an entanglement transition)"),
    ('hypotheses/SPECTRAL_MIDPOINT_HYPOTHESIS.md', 'local:quarter is quantum-classical transition',
     'CΨ = ¼ boundary\\s+(?:separates|divides|marks|is)\\s+(?:the\\s+)?(?:(?:boundary|border|transition)\\s+(?:between|from)\\s+)?quantum[^.;]{0,30}classical', True,
     'The CΨ = ¼ boundary separates quantum from classical behaviour.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO.md', 'local:projection equals EP',
     '90°[^.;]{0,60}\\b(?:pinch|projection|flat|shadow|rotation)\\b[^.;]{0,40}\\b(?:is|becomes)\\b[^.;]{0,20}\\bEP\\b', True,
     'At 90° the projection becomes the EP.',
     'experiments/F86_EP_THROUGH_THE_CLOCK.md:97 (Claude 2026-06-21: the EP is eigenvector coalescence, computed from the toy Liouvillian, not a drawn 90° projection)'),
    ('reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO.md', 'local:cusp equals EP',
     "\\bcusp\\s+(?:is|equals|=)\\s+(?:the\\s+)?EP\\b(?!['’]s\\b)", True,
     'The cusp is the EP.',
     'experiments/F86_EP_THROUGH_THE_CLOCK.md:97 (Claude 2026-06-21: the EP is eigenvector coalescence) + docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md:74 (the cusp is a fixed-point recurrence object)'),
    ('visualizations/README.md', 'local:universal quarter crossings',
     'crossing itself is universal', True,
     'The crossing itself is universal.',
     'simulations/label_facts_independent_checks.py CHECK-14 (|++⟩ under Heisenberg + local Z-dephasing stays a product state: its purity-book CΨ starts at 1 and crosses ¼ while its concurrence book stays at 0, so whether a trajectory crosses depends on the book) + CHECK-13 (one Bell+ trajectory crosses at K = 0.03596, 0.03735 and 0.07192 in the concurrence, purity and constant books; simulations/subsystem_crossing_pairs.py:18-23, d9061803 2026-07-20)'),
    ('visualizations/README.md', 'local:quarter physical boundary',
     'boundary is a property of the physics, not of the measurement', True,
     'The boundary is a property of the physics, not of the measurement.',
     'simulations/label_facts_independent_checks.py CHECK-14 (|++⟩ under Heisenberg + local Z-dephasing stays a product state: its purity-book CΨ starts at 1 and crosses ¼ while its concurrence book stays at 0, so whether a trajectory crosses depends on the book) + CHECK-13 (one Bell+ trajectory crosses at K = 0.03596, 0.03735 and 0.07192 in the concurrence, purity and constant books; simulations/subsystem_crossing_pairs.py:18-23, d9061803 2026-07-20)'),
    ('data/ibm_cusp_precision_april2026/README.md', 'local:exact quarter rows',
     'many data points exactly at 1/4', True,
     'There are many data points exactly at 1/4.',
     'data/ibm_cusp_precision_april2026/*.json via simulations/label_facts_independent_checks.py CHECK-10 (the three payloads hold 15, 17 and 19 cpsi rows and none equals 0.25; the only literal 0.25 in them is a delay factor)'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-0',
     'Above ¼, the system is quantum. Below ¼, it is classical.', False,
     'Above ¼, the system is quantum. Below ¼, it is classical.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-11',
     'the change is sudden and irreversible.', False,
     'the change is sudden and irreversible.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-13',
     'Quantum regime. Coherent oscillation', False,
     'Quantum regime. Coherent oscillation',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-15',
     'a traveling boundary: a wave of quantum-to-classical conversion', False,
     'a traveling boundary: a wave of quantum-to-classical conversion',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-20',
     'A chain can sustain exactly one clean quantum-classical boundary', False,
     'A chain can sustain exactly one clean quantum-classical boundary',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-31',
     'Initial entanglement is required.', False,
     'Initial entanglement is required.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: the product state |0,1> crosses; no initial entanglement needed)'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-32',
     'Low noise keeps CΨ above 1/4 permanently.', False,
     'Low noise keeps CΨ above 1/4 permanently.',
     'docs/ANALYTICAL_FORMULAS.md:887 (F25, Claude 2026-04-06) + simulations/label_facts_independent_checks.py CHECK-6 (finite K_fold=0.03735 for every γ>0)'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-33',
     'The system is so coherent that it stays in the quantum regime indefinitely.', False,
     'The system is so coherent that it stays in the quantum regime indefinitely.',
     'docs/ANALYTICAL_FORMULAS.md:887 (F25, Claude 2026-04-06) + simulations/label_facts_independent_checks.py CHECK-6 (finite K_fold=0.03735 for every γ>0)'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-36',
     'holds for the TOTAL system but not for subsystems.', False,
     'holds for the TOTAL system but not for subsystems.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-38',
     '**Why DD fails (proven):**', False,
     '**Why DD fails (proven):**',
     "simulations/label_facts_independent_checks.py CHECK-4b (first-site Z pulse keeps CΨ=1/8 yet moves CΨ'(0) from +1/6 to −1/3 under the same H)"),
    ('experiments/TEMPORAL_SACRIFICE.md', 'temporal-39',
     'Therefore DD cannot change CΨ -- not approximately, but algebraically.', False,
     'Therefore DD cannot change CΨ -- not approximately, but algebraically.',
     "simulations/label_facts_independent_checks.py CHECK-4b (first-site Z pulse keeps CΨ=1/8 yet moves CΨ'(0) from +1/6 to −1/3 under the same H)"),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-0',
     'CΨ = ¼ crossing at 0.3% accuracy', False,
     'CΨ = ¼ crossing at 0.3% accuracy',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-0',
     'CΨ = ¼ crossing at 0.3% accuracy', False,
     'CΨ = ¼ crossing at 0.3% accuracy',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/GLOSSARY.md', 'q52-0',
     'CΨ = ¼ crossing at 0.3% accuracy', False,
     'CΨ = ¼ crossing at 0.3% accuracy',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-0',
     'CΨ = ¼ crossing at 0.3% accuracy', False,
     'CΨ = ¼ crossing at 0.3% accuracy',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-0',
     'CΨ = ¼ crossing at 0.3% accuracy', False,
     'CΨ = ¼ crossing at 0.3% accuracy',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-1',
     'CΨ = 1/4 crossing at 0.3% accuracy', False,
     'CΨ = 1/4 crossing at 0.3% accuracy',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-1',
     'CΨ = 1/4 crossing at 0.3% accuracy', False,
     'CΨ = 1/4 crossing at 0.3% accuracy',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/GLOSSARY.md', 'q52-1',
     'CΨ = 1/4 crossing at 0.3% accuracy', False,
     'CΨ = 1/4 crossing at 0.3% accuracy',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-1',
     'CΨ = 1/4 crossing at 0.3% accuracy', False,
     'CΨ = 1/4 crossing at 0.3% accuracy',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-1',
     'CΨ = 1/4 crossing at 0.3% accuracy', False,
     'CΨ = 1/4 crossing at 0.3% accuracy',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-2',
     '115.0 us, predicted 114.7 us', False,
     '115.0 us, predicted 114.7 us',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-2',
     '115.0 us, predicted 114.7 us', False,
     '115.0 us, predicted 114.7 us',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/GLOSSARY.md', 'q52-2',
     '115.0 us, predicted 114.7 us', False,
     '115.0 us, predicted 114.7 us',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-2',
     '115.0 us, predicted 114.7 us', False,
     '115.0 us, predicted 114.7 us',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-2',
     '115.0 us, predicted 114.7 us', False,
     '115.0 us, predicted 114.7 us',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-3',
     'predicted 114.7 us (0.3% error)', False,
     'predicted 114.7 us (0.3% error)',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-3',
     'predicted 114.7 us (0.3% error)', False,
     'predicted 114.7 us (0.3% error)',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/GLOSSARY.md', 'q52-3',
     'predicted 114.7 us (0.3% error)', False,
     'predicted 114.7 us (0.3% error)',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-3',
     'predicted 114.7 us (0.3% error)', False,
     'predicted 114.7 us (0.3% error)',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-3',
     'predicted 114.7 us (0.3% error)', False,
     'predicted 114.7 us (0.3% error)',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-4',
     'predicted at 114.7 μs and measured at 115.0 μs (0.3% difference)', False,
     'predicted at 114.7 μs and measured at 115.0 μs (0.3% difference)',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-4',
     'predicted at 114.7 μs and measured at 115.0 μs (0.3% difference)', False,
     'predicted at 114.7 μs and measured at 115.0 μs (0.3% difference)',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/GLOSSARY.md', 'q52-4',
     'predicted at 114.7 μs and measured at 115.0 μs (0.3% difference)', False,
     'predicted at 114.7 μs and measured at 115.0 μs (0.3% difference)',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-4',
     'predicted at 114.7 μs and measured at 115.0 μs (0.3% difference)', False,
     'predicted at 114.7 μs and measured at 115.0 μs (0.3% difference)',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-4',
     'predicted at 114.7 μs and measured at 115.0 μs (0.3% difference)', False,
     'predicted at 114.7 μs and measured at 115.0 μs (0.3% difference)',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-5',
     'three competing explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved', False,
     'three competing explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-5',
     'three competing explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved', False,
     'three competing explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('docs/GLOSSARY.md', 'q52-5',
     'three competing explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved', False,
     'three competing explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-5',
     'three competing explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved', False,
     'three competing explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-5',
     'three competing explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved', False,
     'three competing explanations (SPAM (State Preparation And Measurement errors), TLS (two-level system defects in the chip substrate), boundary structure), unresolved',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-6',
     'three competing explanations (SPAM, TLS, boundary structure), unresolved', False,
     'three competing explanations (SPAM, TLS, boundary structure), unresolved',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-6',
     'three competing explanations (SPAM, TLS, boundary structure), unresolved', False,
     'three competing explanations (SPAM, TLS, boundary structure), unresolved',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('docs/GLOSSARY.md', 'q52-6',
     'three competing explanations (SPAM, TLS, boundary structure), unresolved', False,
     'three competing explanations (SPAM, TLS, boundary structure), unresolved',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-6',
     'three competing explanations (SPAM, TLS, boundary structure), unresolved', False,
     'three competing explanations (SPAM, TLS, boundary structure), unresolved',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-6',
     'three competing explanations (SPAM, TLS, boundary structure), unresolved', False,
     'three competing explanations (SPAM, TLS, boundary structure), unresolved',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-7',
     'non-Markovian revivals (observed as excess late-time coherence in the Q52 data)', False,
     'non-Markovian revivals (observed as excess late-time coherence in the Q52 data)',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-7',
     'non-Markovian revivals (observed as excess late-time coherence in the Q52 data)', False,
     'non-Markovian revivals (observed as excess late-time coherence in the Q52 data)',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('docs/GLOSSARY.md', 'q52-7',
     'non-Markovian revivals (observed as excess late-time coherence in the Q52 data)', False,
     'non-Markovian revivals (observed as excess late-time coherence in the Q52 data)',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-7',
     'non-Markovian revivals (observed as excess late-time coherence in the Q52 data)', False,
     'non-Markovian revivals (observed as excess late-time coherence in the Q52 data)',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-7',
     'non-Markovian revivals (observed as excess late-time coherence in the Q52 data)', False,
     'non-Markovian revivals (observed as excess late-time coherence in the Q52 data)',
     'experiments/FIXED_POINT_SHADOW.md@ed7f8a7f^:10 (ab15ba8e 2026-03-23, the March hardware verdict: the Q52 late-time shadow is qubit-specific detuning, not a universal boundary effect)'),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-8',
     'CPsi crosses within sub-1% of prediction on a good qubit', False,
     'CPsi crosses within sub-1% of prediction on a good qubit',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record) + simulations/framework/confirmations.py:346 (q80 Run 3, the tightest Torino record, still deviates 1.9%)"),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-8',
     'CPsi crosses within sub-1% of prediction on a good qubit', False,
     'CPsi crosses within sub-1% of prediction on a good qubit',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record) + simulations/framework/confirmations.py:346 (q80 Run 3, the tightest Torino record, still deviates 1.9%)"),
    ('docs/GLOSSARY.md', 'q52-8',
     'CPsi crosses within sub-1% of prediction on a good qubit', False,
     'CPsi crosses within sub-1% of prediction on a good qubit',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record) + simulations/framework/confirmations.py:346 (q80 Run 3, the tightest Torino record, still deviates 1.9%)"),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-8',
     'CPsi crosses within sub-1% of prediction on a good qubit', False,
     'CPsi crosses within sub-1% of prediction on a good qubit',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record) + simulations/framework/confirmations.py:346 (q80 Run 3, the tightest Torino record, still deviates 1.9%)"),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-8',
     'CPsi crosses within sub-1% of prediction on a good qubit', False,
     'CPsi crosses within sub-1% of prediction on a good qubit',
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record) + simulations/framework/confirmations.py:346 (q80 Run 3, the tightest Torino record, still deviates 1.9%)"),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-9',
     "Nothing contradicts the cockpit's predictions", False,
     "Nothing contradicts the cockpit's predictions",
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-9',
     "Nothing contradicts the cockpit's predictions", False,
     "Nothing contradicts the cockpit's predictions",
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('docs/GLOSSARY.md', 'q52-9',
     "Nothing contradicts the cockpit's predictions", False,
     "Nothing contradicts the cockpit's predictions",
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-9',
     "Nothing contradicts the cockpit's predictions", False,
     "Nothing contradicts the cockpit's predictions",
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-9',
     "Nothing contradicts the cockpit's predictions", False,
     "Nothing contradicts the cockpit's predictions",
     "simulations/framework/confirmations.py:321-322 (registry entry cpsi_quarter_crossing_torino_feb2026, b474863b 2026-06-18 and 91f50b5f 2026-07-28: measured t* = 114.7 μs, t*/T₂* = 1.036, 10.7% above the predicted 0.936) + simulations/label_facts_independent_checks.py CHECK-12 (114.7 μs is the record's stored measured crossing_us and 115.0 μs its re-interpolation from the same density matrices, printed as 'Predicted' and 'Measured' by simulations/cockpit_ibm_hardware.py:319-324; the 0.3% is the gap between two readings of one record)"),
    ('experiments/COCKPIT_UNIVERSALITY.md', 'q52-10',
     'synthesis reports 1.9% deviation', False,
     'Q52 crossing: synthesis reports 1.9% deviation',
     'simulations/framework/confirmations.py:346 (registry entry cpsi_quarter_crossing_torino_q80_mar2026, b474863b 2026-06-18: the 1.9% deviation is q80 in Run 3 with a same-day T2*) + simulations/framework/confirmations.py:322 (91f50b5f 2026-07-28: Q52 sits 10.7% above its prediction)'),
    ('docs/ANALYTICAL_FORMULAS.md', 'q52-10',
     'synthesis reports 1.9% deviation', False,
     'Q52 crossing: synthesis reports 1.9% deviation',
     'simulations/framework/confirmations.py:346 (registry entry cpsi_quarter_crossing_torino_q80_mar2026, b474863b 2026-06-18: the 1.9% deviation is q80 in Run 3 with a same-day T2*) + simulations/framework/confirmations.py:322 (91f50b5f 2026-07-28: Q52 sits 10.7% above its prediction)'),
    ('docs/GLOSSARY.md', 'q52-10',
     'synthesis reports 1.9% deviation', False,
     'Q52 crossing: synthesis reports 1.9% deviation',
     'simulations/framework/confirmations.py:346 (registry entry cpsi_quarter_crossing_torino_q80_mar2026, b474863b 2026-06-18: the 1.9% deviation is q80 in Run 3 with a same-day T2*) + simulations/framework/confirmations.py:322 (91f50b5f 2026-07-28: Q52 sits 10.7% above its prediction)'),
    ('review/OPEN_QUESTIONS_INDEX.md', 'q52-10',
     'synthesis reports 1.9% deviation', False,
     'Q52 crossing: synthesis reports 1.9% deviation',
     'simulations/framework/confirmations.py:346 (registry entry cpsi_quarter_crossing_torino_q80_mar2026, b474863b 2026-06-18: the 1.9% deviation is q80 in Run 3 with a same-day T2*) + simulations/framework/confirmations.py:322 (91f50b5f 2026-07-28: Q52 sits 10.7% above its prediction)'),
    ('review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md', 'q52-10',
     'synthesis reports 1.9% deviation', False,
     'Q52 crossing: synthesis reports 1.9% deviation',
     'simulations/framework/confirmations.py:346 (registry entry cpsi_quarter_crossing_torino_q80_mar2026, b474863b 2026-06-18: the 1.9% deviation is q80 in Run 3 with a same-day T2*) + simulations/framework/confirmations.py:322 (91f50b5f 2026-07-28: Q52 sits 10.7% above its prediction)'),
    ('compute/RCPsiSquared.Diagnostics/Foundation/MirrorSystem.cs', 'mirrorsystem:EVERY_POSITIVE_GAP',
     '/// Rotation selects the slow shelf for every Gap > 0.', False,
     '/// Rotation selects the slow shelf for every Gap > 0.',
     'compute/RCPsiSquared.Diagnostics/Foundation/MirrorSystem.cs:236-238 (Claude 2026-05-29: const tol = 1e-9; if (gap > tol))'),
    ('compute/RCPsiSquared.Diagnostics/Foundation/MirrorSystem.cs', 'mirrorsystem:EXACT_GAMMA_ZERO_BRANCH',
     '/// The no-decay-style branch is exactly gamma=0.', False,
     '/// The no-decay-style branch is exactly gamma=0.',
     'compute/RCPsiSquared.Diagnostics/Foundation/MirrorSystem.cs:236-238 (Claude 2026-05-29: const tol = 1e-9; if (gap > tol))'),
    ('experiments/TEMPORAL_SACRIFICE.md', 'quantum-classical-regime',
     '\\b(?:above|below)\\s+(?:¼|1/4)\\b[^.;]{0,30}\\b(?:it is|the system is|is)\\s+(?:quantum|classical)\\b|\\b(?:enters?|reaches|stays in|is in) the (?:quantum|classical) (?:regime|zone)\\b|\\bquantum operational window\\b', True,
     'Above ¼ the system is quantum, below ¼ it is classical.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('experiments/COHERENCE_DENSITY.md', 'quantum-classical-regime',
     '\\b(?:above|below)\\s+(?:¼|1/4)\\b[^.;]{0,30}\\b(?:it is|the system is|is)\\s+(?:quantum|classical)\\b|\\b(?:enters?|reaches|stays in|is in) the (?:quantum|classical) (?:regime|zone)\\b|\\bquantum operational window\\b', True,
     'Above ¼ the system is quantum, below ¼ it is classical.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('experiments/SUBSYSTEM_CROSSING.md', 'quantum-classical-regime',
     '\\b(?:above|below)\\s+(?:¼|1/4)\\b[^.;]{0,30}\\b(?:it is|the system is|is)\\s+(?:quantum|classical)\\b|\\b(?:enters?|reaches|stays in|is in) the (?:quantum|classical) (?:regime|zone)\\b|\\bquantum operational window\\b', True,
     'Above ¼ the system is quantum, below ¼ it is classical.',
     "simulations/label_facts_independent_checks.py CHECK-13 (concurrence book, the host's book: Bell+ crosses R = ¼ at f = √3/2 with concurrence √3/2, and its concurrence f stays positive at every finite t after the crossing, so the crossing does not locate an entanglement transition) + simulations/label_facts_independent_checks.py CHECK-1 (full-state book: separable |+++⟩ has CΨ = 1; entangled GHZ₃ has CΨ = 1/7)"),
    ('experiments/DYNAMIC_ENTANGLEMENT.md', 'quantum-classical-regime',
     '\\b(?:above|below)\\s+(?:¼|1/4)\\b[^.;]{0,30}\\b(?:it is|the system is|is)\\s+(?:quantum|classical)\\b|\\b(?:enters?|reaches|stays in|is in) the (?:quantum|classical) (?:regime|zone)\\b|\\bquantum operational window\\b', True,
     'Above ¼ the system is quantum, below ¼ it is classical.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('experiments/N_SCALING_BARRIER.md', 'quantum-classical-regime',
     '\\b(?:above|below)\\s+(?:¼|1/4)\\b[^.;]{0,30}\\b(?:it is|the system is|is)\\s+(?:quantum|classical)\\b|\\b(?:enters?|reaches|stays in|is in) the (?:quantum|classical) (?:regime|zone)\\b|\\bquantum operational window\\b', True,
     'Above ¼ the system is quantum, below ¼ it is classical.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('experiments/README.md', 'quantum-classical-regime',
     '\\b(?:above|below)\\s+(?:¼|1/4)\\b[^.;]{0,30}\\b(?:it is|the system is|is)\\s+(?:quantum|classical)\\b|\\b(?:enters?|reaches|stays in|is in) the (?:quantum|classical) (?:regime|zone)\\b|\\bquantum operational window\\b', True,
     'Above ¼ the system is quantum, below ¼ it is classical.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('docs/GLOSSARY.md', 'quantum-classical-regime',
     '\\b(?:above|below)\\s+(?:¼|1/4)\\b[^.;]{0,30}\\b(?:it is|the system is|is)\\s+(?:quantum|classical)\\b|\\b(?:enters?|reaches|stays in|is in) the (?:quantum|classical) (?:regime|zone)\\b|\\bquantum operational window\\b', True,
     'Above ¼ the system is quantum, below ¼ it is classical.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('README.md', 'quantum-classical-regime',
     '\\b(?:above|below)\\s+(?:¼|1/4)\\b[^.;]{0,30}\\b(?:it is|the system is|is)\\s+(?:quantum|classical)\\b|\\b(?:enters?|reaches|stays in|is in) the (?:quantum|classical) (?:regime|zone)\\b|\\bquantum operational window\\b', True,
     'Above ¼ the system is quantum, below ¼ it is classical.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('docs/proofs/PROOF_MONOTONICITY_CPSI.md', 'absorbing-quarter',
     '\\b(?:¼|1/4|quarter)\\b[^.;]{0,40}\\b(?:is|as) an? (?:absorbing|one-way) (?:boundary|barrier|fold)\\b|\\b(?:absorbing|one-way) boundary at (?:¼|1/4)\\b', True,
     'CΨ = ¼ is an absorbing boundary.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('experiments/DYNAMIC_ENTANGLEMENT.md', 'absorbing-quarter',
     '\\b(?:¼|1/4|quarter)\\b[^.;]{0,40}\\b(?:is|as) an? (?:absorbing|one-way) (?:boundary|barrier|fold)\\b|\\b(?:absorbing|one-way) boundary at (?:¼|1/4)\\b', True,
     'CΨ = ¼ is an absorbing boundary.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('docs/README.md', 'absorbing-quarter',
     '\\b(?:¼|1/4|quarter)\\b[^.;]{0,40}\\b(?:is|as) an? (?:absorbing|one-way) (?:boundary|barrier|fold)\\b|\\b(?:absorbing|one-way) boundary at (?:¼|1/4)\\b', True,
     'CΨ = ¼ is an absorbing boundary.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('README.md', 'absorbing-quarter',
     '\\b(?:¼|1/4|quarter)\\b[^.;]{0,40}\\b(?:is|as) an? (?:absorbing|one-way) (?:boundary|barrier|fold)\\b|\\b(?:absorbing|one-way) boundary at (?:¼|1/4)\\b', True,
     'CΨ = ¼ is an absorbing boundary.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('docs/READING_GUIDE.md', 'absorbing-quarter',
     '\\b(?:¼|1/4|quarter)\\b[^.;]{0,40}\\b(?:is|as) an? (?:absorbing|one-way) (?:boundary|barrier|fold)\\b|\\b(?:absorbing|one-way) boundary at (?:¼|1/4)\\b', True,
     'CΨ = ¼ is an absorbing boundary.',
     'simulations/minimum_energy.py:15 (Claude 2026-04-05: |0,1> crosses upward) + simulations/label_facts_independent_checks.py CHECK-3 (Markovian |01> reaches CΨ=0.309 from 0)'),
    ('experiments/BORN_RULE_SHADOW.md', 'born-derivation',
     '\\b(?:derives?|derived|derivation of)\\s+(?:the\\s+|a\\s+)?Born rule\\b', True,
     'This derives the Born rule.',
     'simulations/born_rule_tier1_derivation.py:97-99 (P_u0 = ⟨00|ρ|00⟩: the Born probabilities are the input of the derivation, not its output)'),
    ('docs/proofs/PROOF_F95_ANGLE_AT_QUADRATIC_ZERO.md', 'born-derivation',
     '\\b(?:derives?|derived|derivation of)\\s+(?:the\\s+|a\\s+)?Born rule\\b', True,
     'This derives the Born rule.',
     'simulations/born_rule_tier1_derivation.py:97-99 (P_u0 = ⟨00|ρ|00⟩: the Born probabilities are the input of the derivation, not its output)'),
    ('README.md', 'born-derivation',
     '\\b(?:derives?|derived|derivation of)\\s+(?:the\\s+|a\\s+)?Born rule\\b', True,
     'This derives the Born rule.',
     'simulations/born_rule_tier1_derivation.py:97-99 (P_u0 = ⟨00|ρ|00⟩: the Born probabilities are the input of the derivation, not its output)'),
    ('docs/READING_GUIDE.md', 'born-derivation',
     '\\b(?:derives?|derived|derivation of)\\s+(?:the\\s+|a\\s+)?Born rule\\b', True,
     'This derives the Born rule.',
     'simulations/born_rule_tier1_derivation.py:97-99 (P_u0 = ⟨00|ρ|00⟩: the Born probabilities are the input of the derivation, not its output)'),
    ('docs/ANALYTICAL_FORMULAS.md', 'born-derivation',
     '\\b(?:derives?|derived|derivation of)\\s+(?:the\\s+|a\\s+)?Born rule\\b', True,
     'This derives the Born rule.',
     'simulations/born_rule_tier1_derivation.py:97-99 (P_u0 = ⟨00|ρ|00⟩: the Born probabilities are the input of the derivation, not its output)'),
    ('experiments/UNIVERSAL_QUANTUM_LIFETIME.md', 'old-lifetime',
     '\\b(?:is|gives|predicts|establishes|yields)\\s+(?:a|the)\\s+universal quantum lifetime\\b|\\btwo quantum/classical exits\\b|\\bsimulator.as.measurement\\b', True,
     'This establishes a universal quantum lifetime.',
     'simulations/label_facts_independent_checks.py CHECK-8 (free-|+⟩ root x=0.4239 vs Bell+ K_fold=0.03735: distinct crossings) + docs/ANALYTICAL_FORMULAS.md:887'),
    ('experiments/README.md', 'old-lifetime',
     '\\b(?:is|gives|predicts|establishes|yields)\\s+(?:a|the)\\s+universal quantum lifetime\\b|\\btwo quantum/classical exits\\b|\\bsimulator.as.measurement\\b', True,
     'This establishes a universal quantum lifetime.',
     'simulations/label_facts_independent_checks.py CHECK-8 (free-|+⟩ root x=0.4239 vs Bell+ K_fold=0.03735: distinct crossings) + docs/ANALYTICAL_FORMULAS.md:887'),
    ('docs/ANALYTICAL_FORMULAS.md', 'old-lifetime',
     '\\b(?:is|gives|predicts|establishes|yields)\\s+(?:a|the)\\s+universal quantum lifetime\\b|\\btwo quantum/classical exits\\b|\\bsimulator.as.measurement\\b', True,
     'This establishes a universal quantum lifetime.',
     'simulations/label_facts_independent_checks.py CHECK-8 (free-|+⟩ root x=0.4239 vs Bell+ K_fold=0.03735: distinct crossings) + docs/ANALYTICAL_FORMULAS.md:887'),
    ('docs/GLOSSARY.md', 'old-lifetime',
     '\\b(?:is|gives|predicts|establishes|yields)\\s+(?:a|the)\\s+universal quantum lifetime\\b|\\btwo quantum/classical exits\\b|\\bsimulator.as.measurement\\b', True,
     'This establishes a universal quantum lifetime.',
     'simulations/label_facts_independent_checks.py CHECK-8 (free-|+⟩ root x=0.4239 vs Bell+ K_fold=0.03735: distinct crossings) + docs/ANALYTICAL_FORMULAS.md:887'),
    ('experiments/COHERENCE_DENSITY.md', 'separable-quarter',
     '\\bseparab(?:le|ility)\\b.{0,65}(?:quarter|r\\s*(?:<=|=|<)\\s*1/4|r\\s*≤\\s*1/4)|(?:quarter|r\\s*(?:<=|=|<)\\s*1/4).{0,65}\\bseparab(?:le|ility)\\b', True,
     'Below the quarter the state is separable.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('experiments/N_SCALING_BARRIER.md', 'separable-quarter',
     '\\bseparab(?:le|ility)\\b.{0,65}(?:quarter|r\\s*(?:<=|=|<)\\s*1/4|r\\s*≤\\s*1/4)|(?:quarter|r\\s*(?:<=|=|<)\\s*1/4).{0,65}\\bseparab(?:le|ility)\\b', True,
     'Below the quarter the state is separable.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('docs/GLOSSARY.md', 'separable-quarter',
     '\\bseparab(?:le|ility)\\b.{0,65}(?:quarter|r\\s*(?:<=|=|<)\\s*1/4|r\\s*≤\\s*1/4)|(?:quarter|r\\s*(?:<=|=|<)\\s*1/4).{0,65}\\bseparab(?:le|ility)\\b', True,
     'Below the quarter the state is separable.',
     'simulations/label_facts_independent_checks.py CHECK-1 (separable |+++⟩ has CΨ=1; entangled GHZ₃ has CΨ=1/7)'),
    ('docs/ANALYTICAL_FORMULAS.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('docs/GLOSSARY.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('hypotheses/GAMMA_IS_LIGHT.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('hypotheses/BRIDGE_PROTOCOL.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('hypotheses/GRAVITY_FROM_WAVE_DEATH.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('hypotheses/PERSPECTIVAL_TIME_FIELD.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('hypotheses/RESONANT_RETURN.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('reflections/ON_WHOSE_TIME_THE_CLOCK_KEEPS.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('experiments/CROSSING_TAXONOMY.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('experiments/K_DOSIMETRY.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('experiments/TRAPPED_LIGHT_LOCALIZATION.md', 'f14-universal',
     '(?:gamma.time|K.invariance|decoherence envelope).{0,40}model.independent|model.independent.{0,35}(?:gamma.time|K.invariance|decoherence envelope)|universal (?:gamma.t|γ\\s*t|K) curve|universal progress coordinate|K\\s*=\\s*γ\\s*t.{0,50}invariant for every Hamiltonian', True,
     'K = γ t is invariant for every Hamiltonian.',
     'simulations/gamma_unit_scaling_gate.py:31,45 (Claude 2026-08-29: gamma-alone scaling is not a symmetry; two states, opposite answers)'),
    ('docs/THE_DOUBLE_ROOT.md', 'cusp-ep-identity',
     "\\bcusp\\s+(?:is|equals|=)\\s+(?:the\\s+)?(?:EP|exceptional point)\\b(?!['’]s\\b)|\\b(?:cusp|exceptional point)\\b[^.;]{0,80}\\b(?:same whirlpool|same vortex|one object)\\b", True,
     'The cusp is the EP.',
     'experiments/F86_EP_THROUGH_THE_CLOCK.md:97 (Claude 2026-06-21: the EP is eigenvector coalescence) + docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md:74 (the cusp is a fixed-point recurrence object)'),
    ('reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO.md', 'cusp-ep-identity',
     "\\bcusp\\s+(?:is|equals|=)\\s+(?:the\\s+)?(?:EP|exceptional point)\\b(?!['’]s\\b)|\\b(?:cusp|exceptional point)\\b[^.;]{0,80}\\b(?:same whirlpool|same vortex|one object)\\b", True,
     'The cusp is the EP.',
     'experiments/F86_EP_THROUGH_THE_CLOCK.md:97 (Claude 2026-06-21: the EP is eigenvector coalescence) + docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md:74 (the cusp is a fixed-point recurrence object)'),
    ('experiments/F86_EP_THROUGH_THE_CLOCK.md', 'cusp-ep-identity',
     "\\bcusp\\s+(?:is|equals|=)\\s+(?:the\\s+)?(?:EP|exceptional point)\\b(?!['’]s\\b)|\\b(?:cusp|exceptional point)\\b[^.;]{0,80}\\b(?:same whirlpool|same vortex|one object)\\b", True,
     'The cusp is the EP.',
     'experiments/F86_EP_THROUGH_THE_CLOCK.md:97 (Claude 2026-06-21: the EP is eigenvector coalescence) + docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md:74 (the cusp is a fixed-point recurrence object)'),
    ('experiments/CRITICAL_SLOWING_AT_THE_CUSP.md', 'cusp-ep-identity',
     "\\bcusp\\s+(?:is|equals|=)\\s+(?:the\\s+)?(?:EP|exceptional point)\\b(?!['’]s\\b)|\\b(?:cusp|exceptional point)\\b[^.;]{0,80}\\b(?:same whirlpool|same vortex|one object)\\b", True,
     'The cusp is the EP.',
     'experiments/F86_EP_THROUGH_THE_CLOCK.md:97 (Claude 2026-06-21: the EP is eigenvector coalescence) + docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md:74 (the cusp is a fixed-point recurrence object)'),
    ('reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO.md', 'projection-ep',
     '(?:plotted\\s+marker|projection|pinch)[^.;]{0,40}\\(EP\\)|(?:90°|90 degrees)[^.;]{0,60}\\b(?:pinch|projection|flat|shadow|rotation)\\b[^.;]{0,40}\\b(?:is|becomes|identifies)\\s+(?:the\\s+)?(?:EP|exceptional point)\\b|\\bsingle line\\s+(?:is|identifies|=|becomes)\\s+(?:the\\s+)?(?:exceptional point|EP|pinch)\\b', True,
     'At 90 degrees the flat pinch is the EP.',
     'experiments/F86_EP_THROUGH_THE_CLOCK.md:97 (Claude 2026-06-21: the EP is eigenvector coalescence, computed from the toy Liouvillian, not a drawn 90° projection)'),
    ('simulations/whirlpool_4d.py', 'projection-ep',
     '(?:plotted\\s+marker|projection|pinch)[^.;]{0,40}\\(EP\\)|(?:90°|90 degrees)[^.;]{0,60}\\b(?:pinch|projection|flat|shadow|rotation)\\b[^.;]{0,40}\\b(?:is|becomes|identifies)\\s+(?:the\\s+)?(?:EP|exceptional point)\\b|\\bsingle line\\s+(?:is|identifies|=|becomes)\\s+(?:the\\s+)?(?:exceptional point|EP|pinch)\\b', True,
     'At 90 degrees the flat pinch is the EP.',
     'experiments/F86_EP_THROUGH_THE_CLOCK.md:97 (Claude 2026-06-21: the EP is eigenvector coalescence, computed from the toy Liouvillian, not a drawn 90° projection)'),
    ('experiments/F86_EP_THROUGH_THE_CLOCK.md', 'projection-ep',
     '(?:plotted\\s+marker|projection|pinch)[^.;]{0,40}\\(EP\\)|(?:90°|90 degrees)[^.;]{0,60}\\b(?:pinch|projection|flat|shadow|rotation)\\b[^.;]{0,40}\\b(?:is|becomes|identifies)\\s+(?:the\\s+)?(?:EP|exceptional point)\\b|\\bsingle line\\s+(?:is|identifies|=|becomes)\\s+(?:the\\s+)?(?:exceptional point|EP|pinch)\\b', True,
     'At 90 degrees the flat pinch is the EP.',
     'experiments/F86_EP_THROUGH_THE_CLOCK.md:97 (Claude 2026-06-21: the EP is eigenvector coalescence, computed from the toy Liouvillian, not a drawn 90° projection)'),
    ('simulations/cusp_spiral_2d.py', 'circle-cardioid',
     '\\b(?:quarter|radial|cusp) circle\\b(?:\\s+[\\w-]+){0,3}?\\s+(?:is|equals|=|becomes)\\s+(?:the\\s+)?(?:cardioid|double[- ]root|cusp)\\b', True,
     'The quarter circle is the cardioid.',
     'docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md:74 (Claude 2026-05-17: cardioid c=μ/2−μ²/4 vs the radial |c|=1/4 circle)'),
    ('docs/THE_DOUBLE_ROOT.md', 'circle-cardioid',
     '\\b(?:quarter|radial|cusp) circle\\b(?:\\s+[\\w-]+){0,3}?\\s+(?:is|equals|=|becomes)\\s+(?:the\\s+)?(?:cardioid|double[- ]root|cusp)\\b', True,
     'The quarter circle is the cardioid.',
     'docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md:74 (Claude 2026-05-17: cardioid c=μ/2−μ²/4 vs the radial |c|=1/4 circle)'),
    ('docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md', 'circle-cardioid',
     '\\b(?:quarter|radial|cusp) circle\\b(?:\\s+[\\w-]+){0,3}?\\s+(?:is|equals|=|becomes)\\s+(?:the\\s+)?(?:cardioid|double[- ]root|cusp)\\b', True,
     'The quarter circle is the cardioid.',
     'docs/proofs/PROOF_F97_CARDIOID_HALF_FIXED_POINT.md:74 (Claude 2026-05-17: cardioid c=μ/2−μ²/4 vs the radial |c|=1/4 circle)'),
    ('docs/proofs/PROOF_MONOTONICITY_CPSI.md', 'additive-collective-equivalence',
     'Z\\s*⊗\\s*I\\s*\\+\\s*I\\s*⊗\\s*Z\\s+is\\s+(?:mathematically\\s+)?equivalent\\s+to\\s+two\\s+independent\\s+local\\s+Z\\s+jumps', True,
     'Z ⊗ I + I ⊗ Z is mathematically equivalent to two independent local Z jumps.',
     'simulations/label_facts_independent_checks.py CHECK-5 (||D[Z1]+D[Z2]−D[Z1+Z2]||=8)'),
    ('experiments/NOISE_ROBUSTNESS.md', 'additive-collective-equivalence',
     'Z\\s*⊗\\s*I\\s*\\+\\s*I\\s*⊗\\s*Z\\s+is\\s+(?:mathematically\\s+)?equivalent\\s+to\\s+two\\s+independent\\s+local\\s+Z\\s+jumps', True,
     'Z ⊗ I + I ⊗ Z is mathematically equivalent to two independent local Z jumps.',
     'simulations/label_facts_independent_checks.py CHECK-5 (||D[Z1]+D[Z2]−D[Z1+Z2]||=8)'),
)

# A rule's proposition can be TRUE under an explicit scope. A unit (sentence for a regex
# rule, paragraph for a literal) that carries the scope is exempt for that rule only.
# When `match_may_be` is given, only matches of that shape are exempt: the scope makes
# one direction of the claim true, not every reading of it.
# rule id -> (scope regex, match_may_be regex or None, why the scoped statement is true)
# The H = 0 scope does not reach a sentence that also quantifies over couplings (H ≠ 0,
# every H, all couplings): there the claim is the general one again.
H_ZERO_SCOPE = (r"^(?!.*\b(?:every|all|nonzero|non-zero)\s+(?:H|Hamiltonians?|couplings?)\b)(?!.*\bH\s*≠\s*0)"
                r".*(?:\b[HJ]\s*=\s*0\b(?![.,]\d)|\bwithout (?:a |any )?Hamiltonian\b|\bno Hamiltonian\b"
                r"|\bHamiltonian-free\b)")
SCOPES = {
    "absorbing-quarter": (
        H_ZERO_SCOPE,
        None,
        "With H = 0, local Z-dephasing multiplies each coherence rho_ij by exp(-2 sum_l gamma_l [i_l != j_l] t) "
        "and leaves the diagonal fixed, so purity and l1 both fall and CΨ can only decrease.",
    ),
    "monotonicity-total-system": (
        H_ZERO_SCOPE,
        None,
        "Same H = 0 argument as absorbing-quarter; CHECK-3 needs the Hamiltonian to rise.",
    ),
    "quantum-classical-regime": (
        r"\bconcurrence book\b|\bWootters book\b|\bR_pair\b",
        r"^above\b[\s\S]*\bquantum$",
        "In the concurrence book R = C·l1/3 with l1 <= 3 on two qubits, so R >= 1/4 forces concurrence >= 1/4 > 0; "
        "only the 'above 1/4 is quantum' direction is true there, 'below 1/4 is classical' stays false (CHECK-13).",
    ),
    "separable-quarter": (
        r"\b(?:a|every|any|the) product state\b[^.;]{0,40}\bis separable\b",
        None,
        "A product state is separable by definition, wherever its CΨ sits.",
    ),
    "old-lifetime": (
        r"^(?!.*\b(?:every|all|any)\s+(?:states?|preparations?)\b)"
        r".*(?:(?:free[ -]?)?\|\+(?:⟩|>)\s*(?:model|state|qubit)\b|\bsingle[- ]qubit\b)",
        r"universal quantum lifetime$",
        "Within the free single-qubit |+> model the lifetime t*/T2 = 0.858367 is platform-independent "
        "(docs/ANALYTICAL_FORMULAS.md:770 F12, ee7e4c7b 2026-04-06; CHECK-8); the false claim is a lifetime across states.",
    ),
    "q52-2": (
        r"\bsame (?:\w+ ){0,2}record\b",
        None,
        "Quoting the old pair and naming both numbers as one record is the correction itself (CHECK-12).",
    ),
    "q52-3": (
        r"\bsame (?:\w+ ){0,2}record\b",
        None,
        "Quoting the old pair and naming both numbers as one record is the correction itself (CHECK-12).",
    ),
    "q52-4": (
        r"\bsame (?:\w+ ){0,2}record\b",
        None,
        "Quoting the old pair and naming both numbers as one record is the correction itself (CHECK-12).",
    ),
    "q52-10": (
        r"\bq80\b",
        None,
        "The synthesis does report 1.9%, for q80 (simulations/framework/confirmations.py:346).",
    ),
}
# rule id -> regex the unit must carry for the rule to apply at all: the false claim is
# only the attribution, so the words alone are not the claim.
REQUIRED_CONTEXT = {
    "q52-10": r"\bQ52\b",
}

MARKER = re.compile(
    r"<!--\s*(/?)(?:QUARTER|VEFFECT|CROSSING|F14)-(CURRENT|HISTORICAL|INTERPRETIVE|ALGEBRAIC)\s*-->"
)
NEGATED_BEFORE = re.compile(
    r"\b(?:not|never|no|neither|nor|without|cannot|can't|does not|do not|did not|was not|is not|"
    r"rather than|instead of|refut\w*|disprov\w*|retired|withdrawn|former|historical|old|dropped|"
    r"denied|false|wrong|mislabel\w*|no longer|abandon\w*)\b[^.;:]{0,70}$",
    re.I,
)
NEGATED_AFTER = re.compile(
    r"^[^.;]{0,60}\b(?:is|are|was|were|remains?)\s+(?:not|false|an? (?:historical|interpretive|retired|old|"
    r"withdrawn|superseded)|no longer)\b"
    r"|^[^.;]{0,80}\b(?:not (?:a )?(?:result|claim|theorem|prediction|conclusion|current|established)|"
    r"interpretive|model-specific|retired|withdrawn|refuted|superseded|historical (?:reading|label|status|"
    r"interpretation|nomenclature)|was (?:wrong|false)|does not hold|is not established|is rejected)\b",
    re.I,
)
NEGATED_INSIDE = re.compile(r"\b(?:not|never|no|neither|nor|cannot|without)\b", re.I)
# Story genres (hypotheses/, reflections/, the front door) may carry a claim as an
# explicitly hedged reading; the hedge words below exempt the sentence there only.
HEDGE = re.compile(
    r"\b(?:we read (?:it |this |that )?as|read as|one may read|reads? like|our picture of|as a picture|"
    r"as an image|the image of|metaphor|as if|interpretive invitation|invitation|whether|might|may|could|"
    r"perhaps|imagine|suppose)\b",
    re.I,
)
QUOTED_PROPOSITION = re.compile(r"(?m)^(\s*[-*]\s*)[\"“][^\"”]+[\"”]\s*:")


def read(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8-sig").replace("\r\n", "\n").replace("\r", "\n")


def current_text(path: str, text: str) -> str:
    if not path.endswith((".md", ".txt")):
        return text
    parts = MARKER.split(text)
    kept, mode = [parts[0]], "CURRENT"
    for i in range(1, len(parts), 3):
        closing, region, body = parts[i], parts[i + 1], parts[i + 2]
        mode = "CURRENT" if closing else region
        if mode == "CURRENT":
            kept.append(body)
    return "\n".join(kept)


def paragraphs(text: str):
    text = QUOTED_PROPOSITION.sub(r"\1[quoted proposition]:", text)
    for paragraph in re.split(r"\n[ \t]*\n", text):
        flat = re.sub(r"\s+", " ", paragraph).strip()
        if flat:
            yield flat


def sentences(text: str):
    for flat in paragraphs(text):
        yield from re.split(r"(?<=[.!?])\s+(?=[A-Z0-9*#`\[(\"“])", flat)


def findings(path: str, text: str) -> list[str]:
    out = []
    scan = current_text(path, text)
    story = path.startswith(("hypotheses/", "reflections/", "docs/quantum/")) or path == "MIRROR_THEORY.md"
    for host, rule, pattern, is_regex, claim, cite in (r for r in RULES if r[0] == path):
        expr = pattern if is_regex else re.escape(re.sub(r"\s+", " ", pattern).strip())
        # A false claim can itself be a negative sentence ("Not one old mode survives");
        # the inside-the-match negation guard applies only to rules that do not carry one.
        negative_rule = bool(NEGATED_INSIDE.search(re.sub(r"\\[bsw]|\\s[*+]", " ", pattern)))
        scope, match_may_be, _ = SCOPES.get(rule, (None, None, None))
        context = REQUIRED_CONTEXT.get(rule)
        # a regex names a proposition (sentence-local); a literal phrase may span a full stop
        for unit in (sentences(scan) if is_regex else paragraphs(scan)):
            if context and not re.search(context, unit, re.I):
                continue
            hit = None
            for match in re.finditer(expr, unit, re.I | re.S):
                if NEGATED_BEFORE.search(unit[:match.start()]) or NEGATED_AFTER.search(unit[match.end():]):
                    continue
                if is_regex and not negative_rule and NEGATED_INSIDE.search(match.group()):
                    continue
                if story and HEDGE.search(unit):
                    continue
                if scope and re.search(scope, unit, re.I) and (
                        match_may_be is None or re.search(match_may_be, match.group(), re.I)):
                    continue
                hit = unit
                break
            if hit:
                out.append(f"{path}:{rule}: {hit[:160]}")
                break
    return out


def mutate(path: str, text: str, claim: str) -> str:
    """Insert the claim as current prose: after the first CURRENT marker if any, else after line 1."""
    if path.endswith(".py"):
        return text.rstrip("\n") + "\n# " + claim.replace("\n", " ") + "\n"
    if path.endswith(".cs"):
        return text.rstrip("\n") + "\n// " + claim.replace("\n", " ") + "\n"
    for m in MARKER.finditer(text):
        if m.group(1) or m.group(2) == "CURRENT":
            at = text.find("\n", m.end()) + 1
            return text[:at] + "\n" + claim + "\n\n" + text[at:]
    at = text.find("\n") + 1
    return text[:at] + "\n" + claim + "\n\n" + text[at:]


HOSTS = tuple(dict.fromkeys(r[0] for r in RULES))


def test_rule_table_is_well_formed():
    keys = [(r[0], r[1]) for r in RULES]
    assert len(keys) == len(set(keys)), "duplicate (host, rule)"
    for host, rule, pattern, is_regex, claim, cite in RULES:
        assert (ROOT / host).is_file(), host
        assert claim and cite, (host, rule)
        if is_regex:
            re.compile(pattern)


@pytest.mark.parametrize("host", HOSTS)
def test_host_carries_no_false_claim(host):
    assert findings(host, read(host)) == []


@pytest.mark.parametrize("host,rule,pattern,is_regex,claim,cite", RULES, ids=[f"{r[0]}::{r[1]}" for r in RULES])
def test_rule_fires_on_its_mutation(host, rule, pattern, is_regex, claim, cite):
    text = read(host)
    mutant = mutate(host, text, claim)
    assert mutant != text
    hits = findings(host, mutant)
    assert any(f"{host}:{rule}:" in item for item in hits), (rule, claim, hits[:3])


def test_historical_region_is_not_scanned():
    host, rule, pattern, is_regex, claim, cite = next(r for r in RULES if r[0].endswith(".md"))
    inside = "<!-- QUARTER-CURRENT -->\n# T\n\n<!-- QUARTER-HISTORICAL -->\nHistorical record:\n\n" + claim + "\n"
    assert findings(host, inside) == []
    current = "<!-- QUARTER-CURRENT -->\n# T\n\n" + claim + "\n"
    assert findings(host, current)
    closed = "<!-- VEFFECT-INTERPRETIVE -->\nx\n<!-- /VEFFECT-INTERPRETIVE -->\n\n" + claim + "\n"
    assert findings(host, closed), "closing a bounded region must return to CURRENT"


def test_negated_or_retired_mention_is_not_a_finding():
    host = "experiments/N_EQUALS_FIVE_CHECK.md"
    claim = "N=5 is not structurally special; all metrics are monotonic"
    assert findings(host, claim + ".\n")
    assert findings(host, "The old row said " + claim + ", and that reading was withdrawn.\n") == []
    assert findings(host, claim + ", and this is not a current claim of this document.\n") == []
    proof = "docs/proofs/PROOF_MONOTONICITY_CPSI.md"
    assert findings(proof, "No local unitary can move C_Psi upward.\n")
    assert findings(proof, "The claim that no local unitary can move C_Psi upward is refuted by the Hadamard example.\n") == []


def test_story_hosts_accept_explicitly_hedged_readings_only():
    story = "hypotheses/README.md"
    assert findings(story, "Resonance Not Channel proves that the system is a Fabry-Perot cavity.\n")
    assert findings(story, "One may read it as an image in which Resonance Not Channel proves that the system is a Fabry-Perot cavity.\n") == []
    experiment = "experiments/N_EQUALS_FIVE_CHECK.md"
    hedged = "One may read this as: N=5 is not structurally special; all metrics are monotonic.\n"
    assert findings(experiment, hedged), "a hedge does not exempt an experiment document"


# TRUE sentences that an over-broad pattern used to flag. Each must stay quiet on its host.
TRUE_SCOPED_SENTENCES = (
    ("README.md", "Under H = 0 and local Z-dephasing the quarter is a one-way barrier: CΨ can only fall."),
    ("docs/proofs/PROOF_MONOTONICITY_CPSI.md",
     "Under H = 0 and local Z-dephasing the quarter is a one-way barrier: CΨ can only fall."),
    ("docs/proofs/PROOF_MONOTONICITY_CPSI.md",
     "Under H = 0 and local Z-dephasing, total-system C_Psi is monotonically non-increasing, since the channel "
     "is unital and every off-diagonal element shrinks."),
    ("experiments/SUBSYSTEM_CROSSING.md",
     "In the concurrence book a pair above 1/4 is quantum in the strict sense: its concurrence is positive."),
    ("experiments/SUBSYSTEM_CROSSING.md",
     "In the concurrence book the crossing happens where entanglement lives: R >= 1/4 forces concurrence >= 1/4."),
    ("experiments/COHERENCE_DENSITY.md", "A product state below the quarter is separable, but so is |+++> far above it."),
    ("README.md", "The decoherence readout builds on the generalized Born rule R_i = C_i·Ψ_i² of Born Rule Mirror."),
    ("experiments/BORN_RULE_SHADOW.md",
     "**Status:** Confirmed (Born rule has zero interference between the past and future parts, by linearity of the trace)"),
    ("docs/THE_DOUBLE_ROOT.md", "The quarter circle meets the cardioid only at c = 1/4, which is the cusp."),
    ("docs/THE_DOUBLE_ROOT.md",
     "The recurrence cusp and the toy exceptional point obey the same condition, a vanishing discriminant, "
     "in two different quadratics."),
    ("docs/THE_DOUBLE_ROOT.md", "In the typed graph the cusp is the EP's sibling, not its identity."),
    ("reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO.md", "In the typed graph the cusp is the EP's sibling, not its identity."),
    ("reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO.md",
     "At Q_EP the angle between the two right eigenvectors collapses from 90° to 0, and that collapse is the EP "
     "of the toy Liouvillian."),
    ("simulations/whirlpool_4d.py",
     "# the eigenvector angle collapses from 90 degrees to zero at Q_EP, and that point is the exceptional point"),
    ("experiments/F86_EP_THROUGH_THE_CLOCK.md",
     "The two right eigenvectors coalesce as their angle goes from 90° to 0, and the point where it reaches 0 is the EP."),
    ("docs/GLOSSARY.md", "Eigenstates of H carry no dynamics under the unitary part alone; dephasing then moves them."),
    ("docs/GLOSSARY.md", "F14 extended to the Wootters book gives K = ln(4/3)/8."),
    ("docs/THE_BRIDGE_WAS_ALWAYS_OPEN.md",
     "Coupling breaks symmetry for 14 of the 36 two-bond combinations at N=3 (the V-Effect census)."),
    ("experiments/UNIVERSAL_QUANTUM_LIFETIME.md",
     "The cubic gives a universal quantum lifetime in units of T2 for the free |+> model: t*/T2 = 0.858367 on every platform."),
    ("hypotheses/SPECTRAL_MIDPOINT_HYPOTHESIS.md",
     "February documents took the CΨ = ¼ boundary to separate quantum from classical behaviour; GHZ₃ at CΨ = 1/7 refutes that."),
    ("docs/GLOSSARY.md", "For q80 (Run 3, same-day Ramsey T2*), the synthesis reports 1.9% deviation."),
    ("review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md",
     "- Q80 entry: `experiments/IBM_HARDWARE_SYNTHESIS.md` (line 142): synthesis reports 1.9% deviation for Run 3 on q80."),
    ("experiments/COCKPIT_UNIVERSALITY.md",
     "The cockpit printout reads 115.0 us, predicted 114.7 us, but both numbers come from the same record."),
    ("hypotheses/WAVES_THAT_HEAR_THEMSELVES.md",
     "A Heisenberg bond keeps the palindrome, so coupling alone breaks the symmetry only when the bond itself is "
     "one of the 14 breaking combinations."),
)

# The false side of each narrowing: a scope or a tighter pattern must not swallow it.
STILL_FALSE_SENTENCES = (
    ("README.md", "CΨ = ¼ is an absorbing boundary even when H = 0.5 couples the pair.", "absorbing-quarter"),
    ("README.md", "Under H = 0 and for every coupling the quarter is an absorbing boundary.", "absorbing-quarter"),
    ("docs/proofs/PROOF_MONOTONICITY_CPSI.md", "Total-system C_Psi is monotonically non-increasing for every H.",
     "monotonicity-total-system"),
    ("experiments/SUBSYSTEM_CROSSING.md", "In the concurrence book a pair below 1/4 is classical.",
     "quantum-classical-regime"),
    ("experiments/SUBSYSTEM_CROSSING.md", "Above 1/4 the pair is quantum, whatever its concurrence.",
     "quantum-classical-regime"),
    ("experiments/COHERENCE_DENSITY.md", "Below the quarter the state is separable.", "separable-quarter"),
    ("experiments/COHERENCE_DENSITY.md", "Below the quarter every state is separable, like a product state.",
     "separable-quarter"),
    ("experiments/UNIVERSAL_QUANTUM_LIFETIME.md", "For the free |+> model there are two quantum/classical exits.",
     "old-lifetime"),
    ("experiments/UNIVERSAL_QUANTUM_LIFETIME.md", "This establishes a universal quantum lifetime for every state.",
     "old-lifetime"),
    ("experiments/UNIVERSAL_QUANTUM_LIFETIME.md",
     "The single-qubit cubic gives a universal quantum lifetime for every state.", "old-lifetime"),
    ("review/OPEN_QUESTIONS_INDEX_PROPOSAL_hardware-test.md",
     "Q52 crossing, resolving document: synthesis reports 1.9% deviation.", "q52-10"),
    ("experiments/COCKPIT_UNIVERSALITY.md", "Q52: 115.0 us, predicted 114.7 us.", "q52-2"),
    ("docs/THE_DOUBLE_ROOT.md", "The cusp is the exceptional point.", "cusp-ep-identity"),
    ("docs/THE_DOUBLE_ROOT.md", "The quarter circle is the cardioid boundary.", "circle-cardioid"),
    ("reflections/ON_THE_WHIRLPOOL_YOU_STEER_TO.md", "At 90 degrees the flat pinch is the EP.", "projection-ep"),
    ("README.md", "This framework derives the Born rule.", "born-derivation"),
    ("docs/GLOSSARY.md", "The gamma-time curve is model-independent.", "f14-universal"),
    ("hypotheses/SPECTRAL_MIDPOINT_HYPOTHESIS.md", "The CΨ = ¼ boundary marks quantum versus classical behaviour.",
     "local:quarter is quantum-classical transition"),
)


@pytest.mark.parametrize("host,sentence", TRUE_SCOPED_SENTENCES)
def test_true_scoped_sentence_stays_quiet(host, sentence):
    assert host in HOSTS, host
    assert findings(host, sentence + "\n") == []


@pytest.mark.parametrize("host,sentence,rule", STILL_FALSE_SENTENCES)
def test_narrowed_rule_still_fires_on_its_false_side(host, sentence, rule):
    assert any(f"{host}:{rule}:" in item for item in findings(host, sentence + "\n")), (rule, sentence)


def test_scopes_and_contexts_name_live_rules():
    rule_ids = {r[1] for r in RULES}
    for rule, (scope, match_may_be, why) in SCOPES.items():
        assert rule in rule_ids, rule
        re.compile(scope)
        if match_may_be:
            re.compile(match_may_be)
        assert why, rule
    for rule, context in REQUIRED_CONTEXT.items():
        assert rule in rule_ids, rule
        re.compile(context)


# ---------------------------------------------------------------- the CHECK owners
CHECKS_SCRIPT = ROOT / "simulations/label_facts_independent_checks.py"


@pytest.fixture(scope="module")
def check_values(tmp_path_factory):
    """Run the owner script once, in a scratch working directory, and read its values."""
    completed = subprocess.run(
        [sys.executable, str(CHECKS_SCRIPT), "--json"],
        cwd=tmp_path_factory.mktemp("label_fact_checks"),
        capture_output=True, text=True, encoding="utf-8", timeout=600, check=False,
        env={**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONDONTWRITEBYTECODE": "1"},
    )
    assert completed.returncode == 0, completed.stderr[-3000:]
    return json.loads(completed.stdout)


def _rational(value: str) -> Fraction:
    return Fraction(value)


# CHECK id -> the values the rules cite, asserted exactly where the check is exact and as
# the cited display otherwise (the cited number IS a rounded display of the float).
def _assert_check_1(v):
    assert _rational(v["cpsi_plus3_separable"]) == 1
    assert _rational(v["cpsi_ghz3_entangled"]) == Fraction(1, 7)


def _assert_check_2(v):
    assert _rational(v["cpsi_00"]) == 0
    assert _rational(v["cpsi_after_hadamard_x_identity"]) == Fraction(1, 3)


def _assert_check_3(v):
    assert v["cpsi_at_0"] == 0.0  # |01> carries no coherence: exact
    assert v["cpsi_max"] > 0.25
    assert f"{v['cpsi_max']:.3f}" == "0.309"
    assert 0.0 < v["t_up"] < v["t_max"] < v["t_down"]


def _assert_check_4b(v):
    assert _rational(v["cpsi_rho0"]) == Fraction(1, 8)
    assert _rational(v["cpsi_after_pulse"]) == Fraction(1, 8)
    assert _rational(v["slope_rho0"]) == Fraction(1, 6)
    assert _rational(v["slope_rho0_h_zero"]) == Fraction(-1, 12)
    assert _rational(v["slope_after_pulse_same_h"]) == Fraction(-1, 3)
    assert _rational(v["slope_after_pulse_passive_frame"]) == Fraction(1, 6)


def _assert_check_5(v):
    assert v["norm_squared"] == "64"
    assert v["norm"] == "8"


def _assert_check_6(v):
    assert f"{v['f_star']:.3f}" == "0.861"
    assert f"{v['k_fold']:.5f}" == "0.03735"
    assert v["concurrence_formula"] == "f"
    assert v["concurrence_at_crossing"] == v["f_star"] > 0


def _assert_check_7(v):
    assert v["V_limit"] == "2"
    assert v["V5_minus_1_plus_phi_half"] == "0"
    assert v["omega21_N5_minus_2_plus_phi"] == "0"
    assert v["omega21_limit"] == "4"
    assert f"{v['phi']:.4f}" == "1.6180"
    assert f"{v['V5']:.4f}" == "1.8090"


def _assert_check_8(v):
    assert f"{v['x_root']:.4f}" == "0.4239"
    assert f"{v['t_star_over_T2']:.6f}" == "0.858367"
    assert f"{v['bell_k_fold']:.5f}" == "0.03735"


def _assert_check_9(v):
    assert (v["C_6_2"], v["C_16_2"]) == (15, 120)


def _assert_check_10(v):
    assert v["files"] == 3
    assert v["exactly_quarter_total"] == 0
    per_file = v["per_file"].values()
    assert sorted(entry["rows"] for entry in per_file) == [15, 17, 19]
    assert sorted(entry["within_0_01_of_quarter"] for entry in per_file) == [0, 0, 6]
    assert all(entry["exactly_quarter"] == 0 and entry["phase_keys"] == 0 for entry in per_file)


def _assert_check_12(v):
    assert f"{v['stored_crossing_us']:.1f}" == "114.7"
    assert f"{v['interpolated_crossing_us']:.1f}" == "115.0"
    assert f"{100 * v['relative_gap']:.1f}" == "0.3"
    assert v["stored_pure_dephasing_prediction_t_over_T2"] == 0.8584


def _assert_check_13(v):
    assert v["f_concurrence_book"] == "sqrt(3)/2"
    assert f"{v['k_concurrence_book_value']:.5f}" == "0.03596"
    assert f"{v['k_purity_book_value']:.5f}" == "0.03735"
    assert f"{v['k_constant_book_value']:.5f}" == "0.07192"
    assert v["t_at_gamma_0_05"] == [0.7192, 0.747, 1.4384]
    assert v["concurrence_formula"] == "f"
    assert v["concurrence_at_concurrence_crossing"] == "sqrt(3)/2"


def _assert_check_14(v):
    assert v["commutator_is_zero"] is True
    assert v["product_family_solves_lindblad"] is True
    assert v["spin_flip_product_is_scalar"] is True
    assert v["concurrence_at_a_0_1_to_1"] == ["0"] * 10
    assert v["purity_book_at_a_1"] == "1"
    assert v["purity_book_crossing_count_in_unit_interval"] == 1
    assert 0.0 < v["purity_book_crossing_a"] < 1.0


CHECK_ASSERTIONS = {
    "CHECK-1": _assert_check_1, "CHECK-2": _assert_check_2, "CHECK-3": _assert_check_3,
    "CHECK-4b": _assert_check_4b, "CHECK-5": _assert_check_5, "CHECK-6": _assert_check_6,
    "CHECK-7": _assert_check_7, "CHECK-8": _assert_check_8, "CHECK-9": _assert_check_9,
    "CHECK-10": _assert_check_10, "CHECK-12": _assert_check_12, "CHECK-13": _assert_check_13,
    "CHECK-14": _assert_check_14,
}


@pytest.mark.parametrize("check", sorted(CHECK_ASSERTIONS))
def test_independent_check_reproduces_its_cited_values(check, check_values):
    assert check in check_values, check
    CHECK_ASSERTIONS[check](check_values[check])


def test_every_cited_check_is_run_and_every_check_is_cited(check_values):
    cited = {match for r in RULES for match in re.findall(r"CHECK-\d+b?", r[5])}
    assert cited <= set(CHECK_ASSERTIONS), sorted(cited - set(CHECK_ASSERTIONS))
    assert set(check_values) == set(CHECK_ASSERTIONS)
    assert cited == set(CHECK_ASSERTIONS), sorted(set(CHECK_ASSERTIONS) - cited)


def test_checks_script_prints_one_line_per_check(tmp_path):
    """The human-readable run a reader reproduces prints every CHECK once, in order."""
    completed = subprocess.run(
        [sys.executable, str(CHECKS_SCRIPT)], cwd=tmp_path, capture_output=True, text=True, encoding="utf-8",
        timeout=600, check=False, env={**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONDONTWRITEBYTECODE": "1"})
    assert completed.returncode == 0, completed.stderr[-3000:]
    printed = [line.split(" ", 1)[0] for line in completed.stdout.splitlines()]
    assert printed == list(CHECK_ASSERTIONS), printed
    assert list(tmp_path.iterdir()) == []
