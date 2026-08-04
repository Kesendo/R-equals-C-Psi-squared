"""Detector for the row/column vectorisation-convention mismatch.

A Liouvillian written as -i(H(x)I - I(x)H^T) is the ROW-stack (order='C')
generator; -i(I(x)H - H^T(x)I) is the COLUMN-stack (order='F') one. Pairing
either with the other stacking evolves conj(L) instead of L, i.e. it runs the
Hamiltonian part backwards, and no conjugation-invariant readout can see it.
Eighteen committed scripts carried this (2026-08-04, commits d54c064 + 1da8cff);
docs/CAUGHT_ERRORS.md holds the account.

This is a HEURISTIC pattern matcher, not a proof of absence. Two things it
cannot do, both worth knowing before trusting a zero:

  * it only recognises the two kron spellings above, so a Liouvillian built any
    other way is invisible to it;
  * it pairs a file's generator with any vec order appearing ANYWHERE in the
    same file, so a file that vectorises something other than rho against L can
    read as a false positive.

Its first version reported "0 files" because its regex did not allow `np.`
before the second `kron`. That is why it prints the full inventory and not only
the mismatches: read the residue, not the verdict.

Known deliberate mismatch, expected in the output and NOT a bug:
simulations/framework/diagnostics/polarity_coordinates.py, which mismatches the
stackings on purpose (it carries F113's sign) and never unvecs a state.

Run:  python simulations/vec_convention_sweep.py
"""
import re, subprocess, os, io

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(root)
tracked = subprocess.run(["git", "ls-files", "*.py"], capture_output=True,
                         text=True).stdout.split("\n")

# generator patterns: kron(A, B) - kron(C, D) where one arg is an identity-ish name
ID = r"(?:np\.)?(?:eye|Id|eye_d|I_d|identity)\w*(?:\([^)]*\))?"
HH = r"[A-Za-z_]\w*"
K = r"(?:np\.)?kron"
row_gen = re.compile(K + r"\(\s*(" + HH + r")\s*,\s*(" + ID + r")\s*\)\s*-\s*"
                     + K + r"\(\s*(" + ID + r")\s*,\s*\1\.T\s*\)")
col_gen = re.compile(K + r"\(\s*(" + ID + r")\s*,\s*(" + HH + r")\s*\)\s*-\s*"
                     + K + r"\(\s*\2\.T\s*,\s*(" + ID + r")\s*\)")
f_vec = re.compile(r"order\s*=\s*['\"]F['\"]")
c_vec = re.compile(r"\.flatten\(\s*\)|\.ravel\(\s*\)|\.reshape\(\s*d\s*,\s*d\s*\)"
                   r"|\.reshape\(\s*\(\s*d\s*,\s*d\s*\)\s*\)")

rows = []
for p in tracked:
    p = p.strip()
    if not p or not p.endswith(".py"):
        continue
    try:
        src = io.open(p, encoding="utf-8", errors="replace").read()
    except OSError:
        continue
    gen_row = [i + 1 for i, l in enumerate(src.split("\n")) if row_gen.search(l)]
    gen_col = [i + 1 for i, l in enumerate(src.split("\n")) if col_gen.search(l)]
    if not (gen_row or gen_col):
        continue
    v_f = [i + 1 for i, l in enumerate(src.split("\n")) if f_vec.search(l)]
    v_c = [i + 1 for i, l in enumerate(src.split("\n")) if c_vec.search(l)]
    rows.append((p, gen_row, gen_col, v_f, v_c))

print(f"{len(rows)} tracked files build an identifiable Liouvillian H-term\n")
bad = []
for p, gr, gc, vf, vc in rows:
    gen = "ROW" if gr else "COL"
    uses = []
    if vf: uses.append("F")
    if vc: uses.append("C")
    mismatch = (gen == "ROW" and vf) or (gen == "COL" and vc)
    flag = "  <<< MISMATCH" if mismatch else ""
    if mismatch: bad.append(p)
    print(f"{gen}  gen@{(gr or gc)[:2]}  vec={'+'.join(uses) or 'none':4s}  "
          f"F@{vf[:2]} C@{vc[:2]}  {p}{flag}")
print(f"\n=== {len(bad)} files with a stacking mismatch ===")
for p in bad:
    print("   ", p)
