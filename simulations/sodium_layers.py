"""The sodium nucleus layer by layer: experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md drawn as stacked filters.

Spin 3/2 is two qubits; the quadrupolar jumps are the Pauli strings ZZ, ZX, ZY, XI, YI (Gamma0..Gamma4) at rates
3J0/2, 3J1/2, 3J1/2, 3J2/2, 3J2/2. Each jump is a layer that lights the strings it anticommutes with and charges them
twice its rate; a string lit by every layer and commuting with H = q*ZZ + w*Iz is a colour (F158's far kernel, the
page's mirror at J0 = 0). Before drawing, the script checks the summed rates against the page's rate table.
Writes visualizations/sodium_layers_neon.png and visualizations/sodium_layers_okabe.png (two palettes, a viewer's choice).
"""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

L = "IXYZ"
def anti(a, b):  # two-qubit strings anticommute when an odd number of sites hold two different non-identity letters
    return sum(1 for x, y in zip(a, b) if x != "I" and y != "I" and x != y) % 2 == 1
STRINGS = [a + b for a in L for b in L]
# the sodium jumps (THE_SODIUM_NUCLEUS_TWO_MIRRORS.md): Gamma0..4 = ZZ, ZX, ZY, XI, YI at rates 3J0/2, 3J1/2, 3J1/2, 3J2/2, 3J2/2
JUMPS = [("ZZ", "J0"), ("ZX", "J1"), ("ZY", "J1"), ("XI", "J2"), ("YI", "J2")]
H = ["ZZ", "ZI", "IZ"]   # H_Q = q*ZZ (Gamma0) plus the offset w*Iz, Iz = ZI + IZ/2 in Souza's labelling

PALETTES = {
    "neon": dict(bg="#080b12", panel="#0c1018", text="#7ef9ff", faint="#15303a",
                 layers=["#b14aed", "#08f7fe", "#08f7fe", "#fe53bb", "#fe53bb"], survivor="#ffe700", dark="#05070b"),
    "okabe": dict(bg="#f4f6f8", panel="#ffffff", text="#1b2530", faint="#d5dce3",
                  layers=["#CC79A7", "#0072B2", "#0072B2", "#E69F00", "#E69F00"], survivor="#009E73", dark="#2a2f35"),
}

def rate(s, jumps):
    terms = {}
    for j, k in jumps:
        if anti(s, j):
            terms[k] = terms.get(k, 0) + 1          # 2 * (3J/2) = 3J; the common 3 is in the caption
    sub = {"J0": "J₀", "J1": "J₁", "J2": "J₂"}
    return "+".join((f"{v}" if v > 1 else "") + sub[k] for k, v in sorted(terms.items())) or "0"

def draw(pal, path):
    P = PALETTES[pal]
    rows = [("J₀ = 0: the zero-frequency density off", JUMPS[1:]), ("J₀ > 0: every density on", JUMPS)]
    fig, axes = plt.subplots(2, 6, figsize=(24, 9.2), facecolor=P["bg"])
    for r, (title, jumps) in enumerate(rows):
        steps = [jumps[:k] for k in range(1, len(jumps) + 1)]
        while len(steps) < 5: steps.insert(0, None)
        steps.append("H")
        for c in range(6):
            ax = axes[r, c]; ax.set_facecolor(P["bg"]); ax.set_xlim(-0.6, 4.1); ax.set_ylim(-0.4, 4.6); ax.set_aspect("equal"); ax.axis("off")
            st = steps[c]
            if st is None:
                ax.text(1.75, 2.0, "(no ZZ layer)", ha="center", va="center", color=P["faint"], fontsize=12); continue
            used = jumps if st == "H" else st
            for idx, s in enumerate(STRINGS):
                i, j = idx // 4, idx % 4
                x, y = j, 3 - i
                ax.add_patch(Rectangle((x, y), 0.92, 0.92, facecolor=P["panel"], edgecolor=P["faint"], lw=1))
                lit = [n for n, (jj, _) in enumerate(used) if anti(s, jj)]
                for n in lit:   # one translucent sheet per layer that lights the string, slightly offset
                    k = JUMPS.index(used[n]); off = 0.035 * n
                    ax.add_patch(Rectangle((x + off, y + off), 0.92 - 2 * off, 0.92 - 2 * off, facecolor=P["layers"][k], alpha=0.22, edgecolor="none"))
                every = len(lit) == len(used)
                keep = every and (st != "H" or all(not anti(s, h) for h in H))
                if st == "H" and every and not keep:
                    ax.add_patch(Rectangle((x, y), 0.92, 0.92, facecolor=P["dark"], alpha=0.7, edgecolor="none"))
                if keep:
                    ax.add_patch(Rectangle((x - 0.03, y - 0.03), 0.98, 0.98, facecolor="none", edgecolor=P["survivor"], lw=3.2))
                ax.text(x + 0.46, y + 0.6, s, ha="center", va="center", color=P["text"], fontsize=13, family="monospace", fontweight="bold")
                ax.text(x + 0.46, y + 0.27, rate(s, used), ha="center", va="center", color=P["text"], fontsize=8.6, family="monospace", alpha=0.85)
            if st == "H":
                surv = [s for s in STRINGS if all(anti(s, j) for j, _ in jumps) and all(not anti(s, h) for h in H)]
                head = "+ H = qZZ + ωI_z:  " + (f"colour {surv[0]}" if surv else "black, no colour")
            else:
                j, k = st[-1]; head = f"+ layer {j}  ({k})"
            ax.text(1.75, 4.35, head, ha="center", color=P["text"], fontsize=12.5)
        axes[r, 0].text(-0.55, 4.35 + 0.0, "", color=P["text"])
        fig.text(0.012, 0.93 - r * 0.47, title, color=P["text"], fontsize=15, fontweight="bold")
    fig.text(0.012, 0.015, "Sodium (spin 3/2 as two qubits): each quadrupolar jump is a layer that lights the strings it anticommutes with and charges them its rate (the small numbers, times 3). "
             "Outlined: lit by every layer. Last panel: of those, the ones that also commute with H. Hues are a choice; the outline is the invariant.", color=P["text"], fontsize=10.5)
    fig.savefig(path, dpi=110, facecolor=P["bg"], bbox_inches="tight")

# the page's rate table (THE_SODIUM_NUCLEUS_TWO_MIRRORS.md, "The dissipator is Pauli-string dephasing"), divided by 3,
# in Souza's labelling: Gamma0 = ZZ, Gamma1 = ZX, Gamma2 = ZY, Gamma3 = XI, Gamma4 = YI and their products
PAGE = {"II": "0", "ZZ": "2J₁+2J₂", "IZ": "2J₁", "ZI": "2J₂",
        "YX": "J₁+J₂", "YY": "J₁+J₂", "XX": "J₁+J₂", "XY": "J₁+J₂",
        "ZX": "J₀+J₁+2J₂", "ZY": "J₀+J₁+2J₂", "IY": "J₀+J₁", "IX": "J₀+J₁",
        "XI": "J₀+2J₁+J₂", "YI": "J₀+2J₁+J₂", "YZ": "J₀+J₂", "XZ": "J₀+J₂"}
for s_ in STRINGS:
    assert rate(s_, JUMPS) == PAGE[s_], (s_, rate(s_, JUMPS), PAGE[s_])
survivors = [s_ for s_ in STRINGS if all(anti(s_, j) for j, _ in JUMPS[1:]) and all(not anti(s_, h) for h in H)]
assert survivors == ["ZZ"], survivors
assert not [s_ for s_ in STRINGS if all(anti(s_, j) for j, _ in JUMPS)]

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "visualizations")
for pal in PALETTES:
    draw(pal, os.path.join(out, f"sodium_layers_{pal}.png"))
print("16 rates match the page's table; J0 = 0 leaves the colour ZZ, J0 > 0 none; two images written")
