"""Two spin-3/2 nuclei (Na and Cl) as four qubits: the filter reading of experiments/THE_NUCLEUS_AS_A_KEY.md.

Rows are the Na strings, columns the Cl strings, the shade is the summed decay rate at illustrative densities. A string lit
by every jump of both nuclei and kept by H is a colour. Before drawing, every panel is read two ways: the surviving
words, and the residual of S(rho) = ZZZZ * rho^T against S L S^-1 = -L - 2 sigma, entry by entry (exactly 0 where
the colour survives). The couplings and their verdicts are those of experiments/THE_SODIUM_CHAIN.md and the spin-lock
section of experiments/THE_SODIUM_NUCLEUS_TWO_MIRRORS.md; between Na and Cl the flip-flop is non-secular and drawn
only as the like-nuclei case.
Writes visualizations/nacl_layers_neon.png and visualizations/nacl_layers_okabe.png.
"""
import itertools, os, sys, numpy as np
sys.stdout.reconfigure(encoding="utf-8")   # titles carry superscripts; the default Windows console cannot print them
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
P={'I':np.eye(2),'X':np.array([[0,1],[1,0]]),'Y':np.array([[0,-1j],[1j,0]]),'Z':np.diag([1.,-1])}
def S(w):
    m=np.array([[1.]])
    for c in w: m=np.kron(m,P[c])
    return m
def anti(a,b): return sum(1 for x,y in zip(a,b) if x!='I' and y!='I' and x!=y)%2==1
# spin 3/2 in Souza's labelling |3/2>,|1/2>,|-1/2>,|-3/2> = |00>,|01>,|10>,|11>
Iz=np.diag([1.5,.5,-.5,-1.5]); Ip=np.zeros((4,4)); Ip[0,1]=np.sqrt(3); Ip[1,2]=2; Ip[2,3]=np.sqrt(3)
Ix=(Ip+Ip.T)/2; Iy=(Ip-Ip.T)/(2j); I4=np.eye(4)
na=lambda A: np.kron(A,I4); cl=lambda A: np.kron(I4,A)
JUMPS=["ZZ","ZX","ZY","XI","YI"]
def jumps(J0,J1,J2,where):  # (string on 4 qubits, rate)
    rates=[1.5*J0,1.5*J1,1.5*J1,1.5*J2,1.5*J2]
    return [((j+"II") if where=="na" else ("II"+j), r) for j,r in zip(JUMPS,rates) if r>0]
def L(H,jl):
    d=H.shape[0]; Id=np.eye(d); M=-1j*(np.kron(H,Id)-np.kron(Id,H.T))
    for w,r in jl:
        A=S(w); M+=r*(np.kron(A,A.conj())-np.kron(Id,Id))
    return M
def survivors(H,jl):
    out=[]
    for w in map(''.join,itertools.product('IXYZ',repeat=4)):
        if all(anti(w,j) for j,_ in jl):
            A=S(w)
            if np.allclose(A@H,H@A,atol=0): out.append(w)   # structural: entries of a Pauli conjugation are exact
    return out

W2 = [a + b for a in "IXYZ" for b in "IXYZ"]
def rate2(w, J):   # one nucleus: twice the summed rate of the jumps a two-qubit word anticommutes with
    r = [1.5 * J[0], 1.5 * J[1], 1.5 * J[1], 1.5 * J[2], 1.5 * J[2]]
    return sum(2 * ri for j, ri in zip(JUMPS, r) if ri > 0 and anti(w, j))
JNa, JCl, JClHot = (0, 2, 1), (0, 3, 1), (3.5, 3, 1)   # J0 = 0 is the limit; the last Cl is monotone, J0 > J1 > J2
qNa, qCl, wNa, wCl = 1.3, 0.7, 0.9, -0.4
H0 = qNa * na(np.diag([1, -1, -1, 1.])) + qCl * cl(np.diag([1, -1, -1, 1.])) + wNa * na(Iz) + wCl * cl(Iz)
Hc = H0 + 0.8 * np.kron(Iz, Iz)
PANELS = [
    ("Na⁺ + Cl⁻, coupled by Iz·Sz", Hc, JCl, "colour ZZ|ZZ kept"),
    ("+ flip-flop IxSx + IySy (like nuclei)", Hc + 0.6 * (np.kron(Ix, Ix) + np.kron(Iy, Iy)), JCl, None),
    ("+ rf field on Na", Hc + 0.5 * na(Ix), JCl, None),
    ("Cl⁻ with J₀ > 0", Hc, JClHot, None),
]
PAL = {
    "neon": dict(bg="#080b12", text="#7ef9ff", grid="#15303a", cmap=LinearSegmentedColormap.from_list("n", ["#0c1018", "#0b4f6c", "#08f7fe", "#fe53bb"]), mark="#ffe700", dead="#000000"),
    "okabe": dict(bg="#f4f6f8", text="#1b2530", grid="#d5dce3", cmap=plt.get_cmap("cividis_r"), mark="#D55E00", dead="#1b2530"),
}

def small(ax, J, title, Q, vmax):
    M = np.array([[rate2(a + b, J) for b in "IXYZ"] for a in "IXYZ"])
    ax.imshow(M, cmap=Q["cmap"], vmin=0, vmax=vmax)
    for i, a in enumerate("IXYZ"):
        for j, b in enumerate("IXYZ"):
            ax.text(j, i, a + b, ha="center", va="center", fontsize=10, family="monospace", color=Q["text"] if M[i, j] < 0.55 * vmax else Q["bg"], fontweight="bold")
    lit = [k for k, w in enumerate(W2) if all(anti(w, j) for j, r in zip(JUMPS, [J[0], J[1], J[1], J[2], J[2]]) if r > 0)]
    for k in lit: ax.add_patch(Rectangle((k % 4 - .5, k // 4 - .5), 1, 1, fill=False, ec=Q["mark"], lw=3))
    ax.set_title(title, color=Q["text"], fontsize=12); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_color(Q["grid"])

def big(ax, H, JC, title, Q, vmax):
    M = np.array([[rate2(a, JNa) + rate2(b, JC) for b in W2] for a in W2])
    ax.imshow(M, cmap=Q["cmap"], vmin=0, vmax=vmax, interpolation="nearest")
    for k in range(1, 4):
        ax.axhline(4 * k - .5, color=Q["grid"], lw=.8); ax.axvline(4 * k - .5, color=Q["grid"], lw=.8)
    jl = jumps(*JNa, "na") + jumps(*JC, "cl")
    lit = [(a, b) for a in W2 for b in W2 if all(anti(a + b, j) for j, _ in jl)]
    surv = survivors(H, jl)
    for a, b in lit:
        i, j = W2.index(a), W2.index(b)
        if a + b in surv: ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False, ec=Q["mark"], lw=3.5))
        else:
            ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fc=Q["dead"], ec=Q["mark"], lw=1.5, ls="--"))
    verdict = f"colour {surv[0][:2]}|{surv[0][2:]}" if surv else ("no word lit by every layer" if not lit else "lit, but H turns it: no colour")
    ax.set_title(f"{title}\n{verdict}", color=Q["text"], fontsize=12)
    ax.set_xticks(range(16)); ax.set_xticklabels(W2, rotation=90, fontsize=6.5, family="monospace", color=Q["text"])
    ax.set_yticks(range(16)); ax.set_yticklabels(W2, fontsize=6.5, family="monospace", color=Q["text"])
    ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_color(Q["grid"])

d = 16; Tm = np.zeros((d * d, d * d))
for i in range(d):
    for j in range(d): Tm[j * d + i, i * d + j] = 1
Ssup = np.kron(S("ZZZZ"), np.eye(d)) @ Tm
for title, H, JC, _ in PANELS:
    jl = jumps(*JNa, "na") + jumps(*JC, "cl"); M = L(H, jl); sigma = sum(r for _, r in jl)
    res = np.abs(Ssup @ M @ np.linalg.inv(Ssup) + M + 2 * sigma * np.eye(d * d)).max()
    surv = survivors(H, jl)
    assert (surv == ["ZZZZ"]) == (res == 0.0), (title, surv, res)
    print(f"{title:42s} colour {surv[0] if surv else 'none':5s} residual {res:.3g}")
assert survivors(PANELS[0][1], jumps(*JNa, "na") + jumps(*JCl, "cl")) == ["ZZZZ"]
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "visualizations")

for name, Q in PAL.items():
    fig = plt.figure(figsize=(22, 8.4), facecolor=Q["bg"])
    vmax = max(rate2(a, JNa) + rate2(b, JClHot) for a in W2 for b in W2)
    gs = fig.add_gridspec(2, 5, width_ratios=[0.62, 1, 1, 1, 1], hspace=0.35, wspace=0.25)
    small(fig.add_subplot(gs[0, 0], facecolor=Q["bg"]), JNa, "Na⁺ alone (rows)", Q, vmax)
    small(fig.add_subplot(gs[1, 0], facecolor=Q["bg"]), JCl, "Cl⁻ alone (columns)", Q, vmax)
    for c, (title, H, JC, _) in enumerate(PANELS):
        big(fig.add_subplot(gs[:, c + 1], facecolor=Q["bg"]), H, JC, title, Q, vmax)
    fig.text(0.01, 0.01, "Two spin-3/2 nuclei as four qubits. Row: the Na string, column: the Cl string; shade: the summed decay rate (darker = slower in this palette's scale). "
             "Solid outline: lit by every jump of both nuclei and kept by H, the colour. Dashed dark cell: lit by every jump, but H does not keep it. "
             "Illustrative densities with J₀ = 0 as the limit (a real nucleus has J₀ > 0): Na (J₀, J₁, J₂) = (0, 2, 1), Cl (0, 3, 1), last panel Cl (3.5, 3, 1).", color=Q["text"], fontsize=10)
    fig.savefig(os.path.join(out, f"nacl_layers_{name}.png"), dpi=110, facecolor=Q["bg"], bbox_inches="tight")
print("ok")
