"""qd_witness_play.py: the witness-reading plays and the fully-witnessed census.

The 2026-07-18/19 night plays with F135+F136 (exact substrate, full-state
closed form + partial trace, reusing qd_letter_gates.py machinery). Three
plays; the third is the census behind the fully-witnessed-worlds corollary
(Corollary 7 of docs/proofs/PROOF_RECORD_LETTER_LAW.md):

PLAY 1  Complete graphs: K_N is a total weave -- every pair holds 1 full Bell
        bit, the letter alternating with N (K3 YY, K4 XX, K5 YY, ...):
        m = N-2 dressers, letter = dresser parity (F136 (ii)).
PLAY 2  The record multiplexer: with S(0)-{1,2,3}, j1=1 adjacent-Bell
        candidate, j2=4 corner-Bell candidate, the single 1-4 bond parity
        ROUTES the testimony: odd -> corner records (YY), adjacent blind;
        even -> adjacent records (XX), corner blind. Adjacent and corner
        Bell witnesses of one S can never coexist (the same bond would need
        both parities).
PLAY 3  The fully-witnessed census: over ALL connected graphs at N=4,5,6
        (38 + 728 + 26704), the worlds where EVERY pair is luminous at
        uniform coupling are EXACTLY the stars and the complete graphs
        (5, 6, 7 winners: the N labeled stars + K_N). This is the from-below
        verification of Corollary 7 (proof in PROOF_RECORD_LETTER_LAW.md):
        uniform => a pair is luminous iff neighborhoods match (Bell) or one
        is the other's leaf (pointer/role-swap); leafless + all-matched
        forces complete, a leaf forces the star. The star is the one world
        showing all three record readings at once (hub-leaf pairs: pointer
        broadcast one way, role-swap the other; leaf-leaf pairs: the Bell
        weave); K_N is fully witnessed with ZERO pointer content anywhere.

All plays are gated (asserts); the script fails loudly on any deviation.
Run: python simulations/qd_witness_play.py   [~10 min, dominated by the N=6 census]
Output: simulations/results/qd_witness/qd_witness_play_out.txt
"""

import numpy as np
from itertools import combinations

exec(open(__file__.replace("qd_witness_play.py", "qd_letter_gates.py"))
     .read().split('print("record-letter GATES')[0])

print("PLAY 1: complete graphs -- the whole-world weave")
for N in range(3, 8):
    bonds = [(a, b, 1.0) for a, b in combinations(range(N), 2)]
    I, corr, td, ev = measure(N, bonds, 0, 1)
    ch = "YY" if abs(corr["YY"]) > 0.5 else ("XX" if abs(corr["XX"]) > 0.5 else "--")
    print(f"  K{N}: m={N - 2}  I={I:.4f}  letter {ch}")
    assert abs(I - 1.0) < 1e-9, f"K{N}: expected 1 full bit, got {I}"
    assert ch == ("YY" if (N - 2) % 2 == 1 else "XX"), f"K{N}: letter {ch} off the dresser parity"

print()
print("PLAY 2: the one-bond record multiplexer")
base = [(0, 1, 1.0), (0, 2, 1.0), (0, 3, 1.0), (1, 2, 1.0), (1, 3, 1.0),
        (4, 2, 1.0), (4, 3, 1.0)]
for r14, tag in [(1.0, "1-4 odd "), (2.0, "1-4 even")]:
    bonds = base + [(1, 4, r14)]
    I1, c1, _, _ = measure(5, bonds, 0, 1)
    I2, c2, _, _ = measure(5, bonds, 0, 4)
    print(f"  {tag}: adjacent j=1 I={I1:.4f}   corner j=4 I={I2:.4f}")
    lit, dark = (I2, I1) if r14 == 1.0 else (I1, I2)
    assert abs(lit - 1.0) < 1e-9 and dark < 1e-9, f"multiplexer at {tag}: {I1}, {I2}"
    letter = c2["YY"] if r14 == 1.0 else c1["XX"]
    assert abs(abs(letter) - 1.0) < 1e-9, f"multiplexer letter at {tag}: {letter}"

print()
print("PLAY 3: the fully-witnessed census (star + K_N and nothing else)")


def census(N):
    edges = list(combinations(range(N), 2))
    winners, n_conn = [], 0
    for mask in range(1 << len(edges)):
        es = [edges[i] for i in range(len(edges)) if mask >> i & 1]
        adj = {i: set() for i in range(N)}
        for a, b in es:
            adj[a].add(b); adj[b].add(a)
        seen, fr = {0}, [0]
        while fr:
            x = fr.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); fr.append(y)
        if len(seen) != N:
            continue
        n_conn += 1
        bonds = [(a, b, 1.0) for a, b in es]
        full = True
        for a, b in combinations(range(N), 2):
            I, corr, td, ev = measure(N, bonds, a, b)
            if I < 0.9999:
                full = False
                break
        if full:
            winners.append(tuple(es))
    return n_conn, winners


N_CONNECTED = {4: 38, 5: 728, 6: 26704}
for N in (4, 5, 6):
    n_conn, winners = census(N)
    kinds = []
    for w in winners:
        degs = sorted(sum(1 for e in w if s in e) for s in range(N))
        kinds.append("STAR" if degs == [1] * (N - 1) + [N - 1]
                     else ("K_N" if degs == [N - 1] * N else "OTHER!"))
    print(f"  N={N}: {n_conn} connected graphs, fully witnessed {len(winners)}: {kinds}")
    assert n_conn == N_CONNECTED[N], f"N={N}: connected count {n_conn} != {N_CONNECTED[N]}"
    assert kinds.count("STAR") == N and kinds.count("K_N") == 1 and "OTHER!" not in kinds, \
        f"N={N}: winners are not exactly the N stars + K_N: {kinds}"

print()
print("witness-play GATES all OK (PLAY 1 letters, PLAY 2 routing, PLAY 3 census)")
