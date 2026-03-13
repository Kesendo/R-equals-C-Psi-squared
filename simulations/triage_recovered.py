"""
TRIAGE: Which recovered files are now mathematically backed?
=============================================================
Maps each of the 44 recovered files to March 13 evidence.
Categorizes: PROVEN, PARTIALLY CONFIRMED, STILL SPECULATIVE, FALLEN
Results -> results_triage.txt
"""
import os, glob

recovered = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\recovered"
outfile = r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations\results_triage.txt"

# What we proved on March 13, 2026
proven_concepts = {
    'mirror_symmetry': {
        'keywords': ['mirror symmetr', 'mirror partner', 'two mirrors', 'facing each other',
                     'we are all mirrors', 'reflection symmetr'],
        'evidence': 'Liouvillian spectrum exactly symmetric around Ng. 100% at N=2-6. Never breaks under any dephasing.',
    },
    'standing_wave': {
        'keywords': ['standing wave', 'stehende welle', 'superposition of two', 'two waves'],
        'evidence': 'c+ and c- are even/odd supermode projections. Two damped sinusoids. Exactly standing wave structure.',
    },
    'between_us': {
        'keywords': ['between us', 'between them', 'between the mirror', 'off-diagonal',
                     'cross-correlation', 'the between'],
        'evidence': 'c+/c- are cross-correlations YZ+ZY and YZ-ZY. Physical content IS between A and B through S.',
    },
    'bidirectional_bridge': {
        'keywords': ['bidirectional', 'two channel', 'two directions', 'two waves meeting'],
        'evidence': 'Two independent information channels: frequency (topology) and decay (noise). Bidirectional confirmed.',
    },
    'quarter_boundary': {
        'keywords': ['1/4', '0.25', 'quarter', 'boundary', 'crossing', 'discriminant'],
        'evidence': 'CPsi=1/4 is Bernoulli variance maximum. z*(1-z*)=CPsi, max at z*=0.5. Not EP, not phase transition.',
    },
    'half_point': {
        'keywords': ['0.5', 'half-fill', 'half full', 'maximum connect', 'fair coin'],
        'evidence': 'z*=0.5 is maximum binary uncertainty. Same 0.5 as IsHalfFull in Stability project.',
    },
    'fixed_point': {
        'keywords': ['fixed point', 'fixpoint', 'convergence', 'iteration converge', 'mandelbrot'],
        'evidence': 'z* = (1-sqrt(1-4CPsi))/2 is the Mandelbrot fixed point. Converges only for CPsi<=1/4.',
    },
    'band_structure': {
        'keywords': ['band structure', 'energy band', 'avoided crossing', 'rate spectrum'],
        'evidence': 'Decay rates form bands at N>=4. Avoided crossings confirmed. Boundary 2g to 2(N-1)g exact.',
    },
    'operator_feedback': {
        'keywords': ['operator feedback', 'state-dependent', 'feedback mechanism'],
        'evidence': 'Documented in Tier 2. State-dependent decoherence preserves delta longer.',
    },
    'pauli_complement': {
        'keywords': ['complement', 'pauli weight', 'k and n-k', 'complementary'],
        'evidence': 'Mirror pairs have complementary Pauli weights: k + (N-k) = N. Proven in Pauli decomposition.',
    },
    'projection': {
        'keywords': ['projection', 'observable filter', 'what you see depends', 'observability'],
        'evidence': 'c+/c- are observable projections, not symmetry sectors. Both parity +1. Confirmed by graph symmetry test.',
    },
    'noise_immunity': {
        'keywords': ['noise immun', 'frequency unchanged', 'robust against', 'topology independent decay'],
        'evidence': 'Frequencies immune to noise (all sweeps). Decay rates topology-independent at N=3.',
    },
    'mirrors_never_break': {
        'keywords': ['never break', 'always symmetric', 'survives', 'robust symmetry', 'unbreakable'],
        'evidence': 'best_sym=100% at every alpha from pure dephasing to pure amplitude damping. Mirrors never break.',
    },
}

# What fell / was disproven
fallen_concepts = {
    'consciousness': ['consciousness is fundamental', 'consciousness required', 'hard problem',
                      'consciousness creates', 'observer creates reality'],
    'gravity': ['gravitational', 'schwarzschild', 'black hole', 'white hole', 'big bang',
                'gravity bridge', 'spacetime'],
    'ftl_signaling': ['faster than light', 'superluminal', 'ftl signal', 'no-signalling violation'],
    'qkd_breaking': ['qkd', 'eavesdrop', 'key distribution'],
    'time_travel': ['time travel', 'backward in time', 'retrocaus'],
}

f = open(outfile, "w")
def log(msg):
    print(msg)
    f.write(msg + "\n")

log("=" * 80)
log("TRIAGE: Recovered Files vs March 13 Evidence")
log("=" * 80)

tier1_proven = []     # Core claims now backed by math
tier2_partial = []    # Mix of proven and unproven claims
tier3_speculative = [] # Still speculative but not disproven
tier4_fallen = []     # Contains disproven claims

for filepath in sorted(glob.glob(os.path.join(recovered, "*.md"))):
    name = os.path.basename(filepath)
    with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
        content = fh.read()
    content_lower = content.lower()
    lines = content.count('\n')
    
    # Count proven concept hits
    proven_hits = {}
    for concept, info in proven_concepts.items():
        count = sum(content_lower.count(kw) for kw in info['keywords'])
        if count > 0:
            proven_hits[concept] = count
    
    # Count fallen concept hits
    fallen_hits = {}
    for concept, keywords in fallen_concepts.items():
        count = sum(content_lower.count(kw) for kw in keywords)
        if count > 0:
            fallen_hits[concept] = count
    
    total_proven = sum(proven_hits.values())
    total_fallen = sum(fallen_hits.values())
    
    # Categorize
    if total_fallen > total_proven and total_fallen > 5:
        tier = "FALLEN"
        tier4_fallen.append((name, lines, proven_hits, fallen_hits))
    elif total_proven > 10 and total_fallen < 3:
        tier = "PROVEN"
        tier1_proven.append((name, lines, proven_hits, fallen_hits))
    elif total_proven > 3:
        tier = "PARTIAL"
        tier2_partial.append((name, lines, proven_hits, fallen_hits))
    else:
        tier = "SPECULATIVE"
        tier3_speculative.append((name, lines, proven_hits, fallen_hits))
    
    log(f"\n{'='*60}")
    log(f"  {name} ({lines} lines) -> {tier}")
    log(f"  Proven hits: {total_proven} | Fallen hits: {total_fallen}")
    if proven_hits:
        top = sorted(proven_hits.items(), key=lambda x: -x[1])[:5]
        log(f"  Proven concepts: {', '.join(f'{k}({v})' for k,v in top)}")
    if fallen_hits:
        top = sorted(fallen_hits.items(), key=lambda x: -x[1])[:3]
        log(f"  Fallen concepts: {', '.join(f'{k}({v})' for k,v in top)}")

log(f"\n{'='*80}")
log("SUMMARY")
log("="*80)

log(f"\n=== TIER 1: PROVEN (restore to active docs) === [{len(tier1_proven)} files]")
for name, lines, ph, fh in tier1_proven:
    log(f"  {name} ({lines} lines)")

log(f"\n=== TIER 2: PARTIALLY CONFIRMED (restore with notes) === [{len(tier2_partial)} files]")
for name, lines, ph, fh in tier2_partial:
    log(f"  {name} ({lines} lines)")

log(f"\n=== TIER 3: STILL SPECULATIVE (keep in recovered/) === [{len(tier3_speculative)} files]")
for name, lines, ph, fh in tier3_speculative:
    log(f"  {name} ({lines} lines)")

log(f"\n=== TIER 4: CONTAINS FALLEN CLAIMS (keep in recovered/, mark fallen) === [{len(tier4_fallen)} files]")
for name, lines, ph, fh in tier4_fallen:
    log(f"  {name} ({lines} lines)")

log(f"\n{'='*80}")
f.close()
print(f"\n>>> Results saved to: {outfile}")
