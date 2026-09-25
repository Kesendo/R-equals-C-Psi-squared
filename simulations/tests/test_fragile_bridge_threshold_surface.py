"""The fragile bridge's threshold: its certificate, its zeros, its sweep, and the pages that quote them.

Everything here is two qubits per chain. The threshold is a second-order exceptional point on the
real gamma axis: certified at J_bridge = 1.0 and 1.9 in the first popcount block of L to go unstable,
and read as one in the symmetry sector of that block that goes first at every coupling of the
producer's 4000-point sweep except the grid point 3/4. At five exact couplings an even and an odd
level of one block already coincide at zero gain and the threshold is zero. The verdict is carried by the tests that compute: they read the committed
artifact and check that its numbers are what they claim to be (the pair's real squared gap changing
sign, one null direction inside the block and in the full L, the singular values attributed block
by block, one Puiseux branch with an O(delta) next order, opposite Krein signatures; one Jordan pair
per block of the first block's orbit under the spin flip and S, with the count recomputed here from
the block label; each zero against its exact form, a linear onset equal to the first-order slope, a
V around each zero; the sweep's runs and maximum, its Jordan rank at every grid coupling with the
exceptions matched against the zeros, and its large-coupling samples against its own last row), and
they recompute the J_bridge = 1 collision and the J_bridge = 3/4 zero from this file's own small
builder (the producer is never imported: its result file is written at the end of a run). What the
pages are held to is numbers and scope: they quote what the producer prints, and they carry "two
qubits per chain" and the certified couplings. No test requires or forbids verdict wording on a
page. Of the producer's own lines the tests require its verdict lines, which it prints only when its
gates pass, and scope sentences such as "These are samples; they establish no limit.", which it
prints unconditionally.
"""
import re
from pathlib import Path

import numpy as np


ROOT = Path(__file__).parents[2]
SIGNATURE_SCRIPT = ROOT / "simulations" / "fragile_bridge_ep_signature.py"
SIGNATURE_RESULT = ROOT / "simulations" / "results" / "fragile_bridge_ep_signature.txt"
PHASE3_SCRIPT = ROOT / "simulations" / "pt_palindrome_breaking.py"
FORMULAS = ROOT / "docs" / "ANALYTICAL_FORMULAS.md"
WHAT_WE_FOUND = ROOT / "docs" / "WHAT_WE_FOUND.md"
GLOSSARY = ROOT / "docs" / "GLOSSARY.md"
FRAGILE = ROOT / "hypotheses" / "FRAGILE_BRIDGE.md"

SUPERSCRIPT = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def _section(text, start, end):
    return text.split(start, 1)[1].split(end, 1)[0]


def _normalized(text):
    return " ".join(text.replace("**", "").split())


def _page_number(printed):
    """A printed 1.8e-14 as the pages write it: 1.8·10⁻¹⁴ (trailing mantissa zeros dropped)."""
    mantissa, exponent = printed.split("e")
    mantissa = mantissa.rstrip("0").rstrip(".") if "." in mantissa else mantissa
    return "{0}·10{1}".format(mantissa, str(int(exponent)).translate(SUPERSCRIPT))


def _artifact():
    return SIGNATURE_RESULT.read_text(encoding="utf-8")


def _certificate_blocks():
    certificate = _section(_artifact(), "1. Threshold certificate", "2. Sector structure")
    blocks = certificate.split("\nJ_bridge = ")[1:]
    return {float(block.split(None, 1)[0]): block for block in blocks}


def _float(pattern, block):
    return float(re.search(pattern, block).group(1))


def _complex_pairs(block):
    """The printed pairs on both sides of gamma*, as (delta, a, b)."""
    pairs = []
    for match in re.finditer(r"delta = ([+-][0-9.e+-]+):\s+\(([+-][0-9.e+-]+) ([+-][0-9.]+)i\)\s+"
                             r"\(([+-][0-9.e+-]+) ([+-][0-9.]+)i\)", block):
        delta = float(match.group(1))
        a = complex(float(match.group(2)), float(match.group(3)))
        b = complex(float(match.group(4)), float(match.group(5)))
        pairs.append((delta, a, b))
    return pairs


def _copies(first_block):
    """Blocks carrying the same spectrum as this one: its orbit under the spin flip (p,q) ->
    (4-p,4-q) and S (p,q) -> (p,4-q), recomputed here from the label alone."""
    if "flip" in first_block:
        return 1
    p, q = (int(x) for x in re.findall(r"\d", first_block))
    return len({(p, q), (4 - p, 4 - q), (p, 4 - q), (4 - p, q)})


def test_certificate_numbers_in_the_artifact_read_as_an_ep2():
    blocks = _certificate_blocks()
    assert sorted(blocks) == [1.0, 1.9]
    for j_bridge, block in blocks.items():
        own_sector = re.search(r"first sector to go unstable: (\(2,2\) flip-(?:even|odd))", block).group(1)
        # The pair's squared gap is real for real gamma (two imaginary eigenvalues below, a mirror
        # pair above, both read off the printed pairs) and changes sign, so its zero is real.
        f_below = _float(r"f\(gamma\*\(1 - 1e-6\)\) = ([+-][0-9.e+-]+)", block)
        f_above = _float(r"f\(gamma\*\(1 \+ 1e-4\)\) = ([+-][0-9.e+-]+)", block)
        assert f_below < 0.0 < f_above, (j_bridge, f_below, f_above)
        pairs = _complex_pairs(block)
        assert len(pairs) == 4
        for delta, a, b in pairs:
            split = abs(a - b)
            if delta < 0:
                assert max(abs(a.real), abs(b.real)) <= 1e-6 * split, (j_bridge, delta, a, b)
            else:
                assert abs(a.imag - b.imag) <= 1e-6 * split, (j_bridge, delta, a, b)
        gap = _float(r"pair gap ([0-9.e+-]+)", block)
        # Algebraic multiplicity 2 in the full L, geometric multiplicity 1 both inside the pair's
        # block and in the full L: at a located Jordan pair the smallest singular value is of the
        # order of the squared pair gap, far below the gap, while the next stays at the scale of the
        # rest of the spectrum; a semisimple pair split by the same gap would put both near the gap.
        assert int(re.search(r"within 1e-5 of lambda\*: (\d+)", block).group(1)) == 2, j_bridge
        sector = [float(x) for x in re.search(
            r"in the sector: smallest singular values of L - lambda\* I: ([0-9.e+-]+), then ([0-9.e+-]+)",
            block).groups()]
        full = [float(x) for x in re.search(
            r"smallest singular values ([0-9.e+-]+)\s+([0-9.e+-]+)\s+([0-9.e+-]+)", block).groups()]
        for sv in (sector, full):
            assert sv[0] <= 1e-3 * gap and sv[1] >= 1e3 * gap, (j_bridge, sv, gap)
        # L is block-diagonal, so the full-L singular values are the blocks' together: the three
        # smallest, attributed block by block, agree with the full-L ones to rounding (the smallest
        # sits at the rounding level in both, the next two agree to the printed digits), and the
        # smallest lives in the colliding pair's own sector.
        attributed = re.search(r"here to ([0-9.e+-]+)\): (.+)", block)
        assert float(attributed.group(1)) <= 1e-14, j_bridge
        entries = re.findall(r"([0-9.]+e[+-]\d+) in (\(\d,\d\)(?: flip-(?:even|odd))?)", attributed.group(2))
        assert len(entries) == 3
        assert float(entries[0][0]) <= 1e-3 * gap, (j_bridge, entries[0])
        assert [float(v) for v, _ in entries[1:]] == full[1:], (j_bridge, entries, full)
        assert entries[0][1] == own_sector, (j_bridge, entries[0], own_sector)
        # One Puiseux branch: the same half-split / sqrt|delta| below and above at 1e-4 and 1e-6.
        both_sides = [float(x) for x in re.findall(r"half-split/sqrt\|delta\| = ([0-9.]+)", block)]
        assert len(both_sides) == 4
        assert max(both_sides) / min(both_sides) < 1.01, (j_bridge, both_sides)
        # The next Puiseux order is O(delta): its deviation shrinks about tenfold per decade.
        above = [float(x) for x in re.findall(r"([0-9.]+) \(1e-0[2-6]\)", block)]
        assert len(above) == 5
        deviation = [value - above[-1] for value in above]
        assert 5.0 < deviation[0] / deviation[1] < 20.0, (j_bridge, above)
        # A Krein collision: opposite signatures under the chain reflection, shrinking toward zero.
        signatures = [(float(a), float(b)) for a, b in re.findall(r"([+-][0-9.]+) / ([+-][0-9.]+) at delta", block)]
        assert len(signatures) == 2
        for plus, minus in signatures:
            assert plus > 0 > minus, (j_bridge, signatures)
        assert abs(signatures[1][0]) < abs(signatures[0][0])
        # The metric is exact: a permutation that commutes with the Hamiltonian part, anticommutes
        # with the gain-loss part, and L is complex symmetric.
        assert re.search(r"\|G LH G - LH\| = 0\.0e\+00, \|G D G \+ D\| = 0\.0e\+00, \|L - L\^T\| = 0\.0e\+00", block)


def _sector_rows(section):
    return re.findall(r"^\s+([0-9.]+)\s+(\(\d,\d\)(?: flip-(?:even|odd))?)\s+(\S+(?: flip-(?:even|odd))?)"
                      r"\s+(\d+)\s+([0-9.]+)\s+([0-9.]+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$", section, re.MULTILINE)


def test_sector_structure_in_the_artifact():
    section = _section(_artifact(), "2. Sector structure", "3. Zeros of the threshold")
    # the symmetry checks are exact
    for check in ("|S L - L S|", "|S^2 - 1|", "|S F - F S|", "|F L F - L|", "|Pi^2 - F|",
                  "|R L(gamma) R - L(-gamma)|", "|G conj(L) G + L|", "|A conj(L) A^-1 + L|", "|A - R|",
                  "|G conj(A) - S|"):
        assert re.search(re.escape(check) + r" = 0\.0e\+00", section), check
    assert section.count("holds") >= 3 and "fails" not in section
    rows = _sector_rows(section)
    assert len(rows) >= 10
    seen = set()
    for j_bridge, first, orbit, copies, gamma_star, omega, within, rounding, positive in rows:
        c = _copies(first)
        assert int(copies) == c, (j_bridge, first, copies)
        assert first.split()[0] in orbit, (j_bridge, first, orbit)
        # one Jordan pair per copy, one quartet per copy
        assert (int(within), int(rounding), int(positive)) == (2 * c, c, 2 * c), (j_bridge, first)
        assert float(omega) > 0.0
        seen.add(c)
    assert seen == {1, 2, 4}
    window = re.search(r"The \(0,1\) orbit goes first for J_bridge between ([0-9.]+) .*?and ([0-9.]+) \(",
                       section, re.DOTALL)
    lo, hi = float(window.group(1)), float(window.group(2))
    inside = [float(j) for j, first, *_ in rows if first == "(0,1)"]
    assert inside and all(lo < j < hi for j in inside)
    assert all(not lo < float(j) < hi for j, first, *_ in rows if first != "(0,1)")
    diagonal = re.search(r"is diagonal, (.+?) \(sites 0 and 1 are chain A\)", section).group(1)
    entries = {int(k): float(v) for v, k in re.findall(r"([+-]\d) gamma on site (\d)", diagonal)}
    assert entries == {0: -2.0, 1: -2.0, 2: 2.0, 3: 2.0}
    crossing = re.search(r"cross at J_bridge = ([0-9.]+), gamma = ([0-9.]+); the lowest other sector there is \S+ at ([0-9.]+)",
                         section)
    j_cross, gamma_cross, other = (float(x) for x in crossing.groups())
    assert other > gamma_cross
    halves = re.findall(r"^\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+\(2,2\) flip-(even|odd)\s+([0-9.]+)\s*$",
                        section, re.MULTILINE)
    assert len(halves) >= 6
    leads = [lead for _, _, _, lead, _ in halves]
    thresholds = [float(value) for _, _, _, _, value in halves]
    switch = leads.index("even")
    assert set(leads[:switch]) == {"odd"} and set(leads[switch:]) == {"even"}
    assert all(b > a for a, b in zip(thresholds[:switch], thresholds[1:switch]))
    assert all(b < a for a, b in zip(thresholds[switch:], thresholds[switch + 1:]))
    assert max(thresholds) < gamma_cross
    below = re.search(r"at J_bridge = 2\.0 the threshold, ([0-9.]+), is ([0-9.]+)% below", section)
    assert abs(float(below.group(2)) - 100.0 * (gamma_cross - float(below.group(1))) / gamma_cross) < 0.01


def _zero_rows():
    section = _section(_artifact(), "3. Zeros of the threshold", "4. Finite-offset diagnostics")
    rows = re.findall(r"^\s+([0-9]\.[0-9]{12}) (x = \S+|real root of [^()]+?)\s{2,}\(.*?"
                      r"([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$", section, re.MULTILINE)
    return section, rows


def _exact_value(form):
    """The coupling J_bridge / J that a printed exact form names (x = J_bridge / J)."""
    if form.startswith("x = "):
        value = form[4:]
        if value == "sqrt(5)/2":
            return np.sqrt(5.0) / 2.0
        num, den = value.split("/")
        return float(num) / float(den)
    poly = form[len("real root of "):].strip()
    coefficients = [0.0, 0.0, 0.0, 0.0]
    for number, power in re.findall(r"([+-]?\d+)(x\^3|x\^2|x)?", poly.replace(" ", "")):
        coefficients[3 - {"x^3": 3, "x^2": 2, "x": 1, "": 0}[power]] = float(number)
    roots = [r.real for r in np.roots(coefficients) if abs(r.imag) < 1e-9]
    assert len(roots) == 1
    return roots[0]


def test_zeros_of_the_threshold_in_the_artifact():
    section, rows = _zero_rows()
    assert len(rows) == 5
    assert "x = J_bridge / J with J = 1" in section
    table = []
    for j0, form, m2, m4, m6 in rows:
        j0 = float(j0)
        # each coupling is where its exact form says
        assert abs(j0 - _exact_value(form.strip())) < 1e-11, (j0, form)
        table.append((j0, float(m2), float(m4), float(m6)))
    # the largest first-order slope per coupling (the s column of the rows that belong to it)
    blocks = re.split(r"^\s+[0-9]\.[0-9]{12} ", section, flags=re.MULTILINE)[1:]
    for (j0, m2, m4, m6), block in zip(table, blocks):
        slopes = [float(s) for s in re.findall(r"[+-][0-9.]+\s+\d\s+([0-9.]+)\s+[0-9.]+", block.split("\n\n")[0])]
        s = max(slopes)
        # the onset is linear: max Re lambda / gamma equals s from 1e-4 down, and the correction at
        # 1e-2 stays well under the 1e-2 relative a correction linear in gamma could give
        assert abs(m4 / s - 1.0) < 2e-6 and abs(m6 / s - 1.0) < 2e-6, (j0, s, m4, m6)
        assert abs(m2 / s - 1.0) < 1e-3, (j0, s, m2)
    control = re.search(r"Control at a generic coupling, J_bridge = 1.0: max Re lambda / gamma = "
                        r"([0-9.e+-]+)\s+([0-9.e+-]+)\s+([0-9.e+-]+)", section)
    for ratio, gamma in zip((float(x) for x in control.groups()), (1e-2, 1e-4, 1e-6)):
        assert ratio * gamma < 1e-12
    positive = re.search(r"Over 61 gammas from 1e-6 to 2 at each of them, max Re lambda stays positive "
                         r"\(smallest max Re lambda / gamma:\s+([0-9.]+),", section)
    assert positive and float(positive.group(1)) > 0.0
    numeric = re.search(r"each refined by brentq: J_bridge = (.+)", section).group(1)
    numeric = sorted(float(x) for x in numeric.split(","))
    assert len(numeric) == 5
    assert all(abs(a - b[0]) < 1e-10 for a, b in zip(numeric, sorted(table)))
    vees = re.findall(r"^\s+([0-9]\.[0-9]{12})\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+"
                      r"([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$", section, re.MULTILINE)
    assert len(vees) == 5
    for j0, far_left, left, right, far_right, slope_left, slope_right, predicted in vees:
        far_left, left, right, far_right = (float(x) for x in (far_left, left, right, far_right))
        # linear rise on both sides: doubling the offset doubles the threshold
        assert 1.9 < far_left / left < 2.1 and 1.9 < far_right / right < 2.1, j0
        mean = 0.5 * (float(slope_left) + float(slope_right))
        assert abs(mean / float(predicted) - 1.0) < 0.02, (j0, mean, predicted)
    take = re.search(r"at J_bridge = ([0-9.]+), gamma = ([0-9.]+)\. It leaves the zero at the V slope ([0-9.]+) "
                     r"and rises on\s+average at ([0-9.]+) per unit", section)
    j_take, g_take, v_slope, mean_slope = (float(x) for x in take.groups())
    assert abs(mean_slope - g_take / (j_take - table[-1][0])) < 1e-4
    assert mean_slope < v_slope


def _sweep_section():
    return _section(_artifact(), "6. The first block over a sweep of couplings", "7. The product")


def test_sweep_of_the_first_block():
    section = _sweep_section()
    runs = re.findall(r"^\s+(\(\d,\d\)(?: flip-(?:even|odd))?)\s+([0-9.]+)\s+([0-9.]+)\s+(\d)\s+(.+)$",
                      section, re.MULTILINE)
    assert len(runs) >= 8
    # contiguous, in order, covering the grid, copies recomputed from the label
    assert float(runs[0][1]) == 0.005 and float(runs[-1][2]) == 20.0
    for (_, _, end, _, _), (_, start, _, _, _) in zip(runs[:-1], runs[1:]):
        assert abs(float(start) - float(end) - 0.005) < 1e-9
    for name, _, _, copies, _ in runs:
        assert int(copies) == _copies(name), name
    four = [(name, a, b) for name, a, b, c, _ in runs if int(c) == 4]
    assert [name for name, _, _ in four] == ["(0,1)"]
    assert "(1,1)" not in {name for name, *_ in runs}
    assert "On this grid (1,1) is never the first block" in section
    largest = re.search(r"The largest threshold on the grid is ([0-9.]+), at J_bridge = ([0-9.]+); the crossing", section)
    crossing = float(re.search(r"\(section 2, gamma = ([0-9.]+)\) lies above every grid value", section).group(1))
    assert float(largest.group(1)) < crossing
    # the run shapes the pages quote
    shapes = {(name, float(a)): shape for name, a, _, _, shape in runs}
    assert shapes[("(2,2) flip-odd", 1.495)].startswith("rises from")
    assert shapes[("(1,2)", 4.46)].startswith("falls from")
    assert "a minimum" in shapes[("(2,2) flip-even", 1.955)]


def _sweep_rank_read():
    section = _sweep_section()
    grid = int(re.search(r"\((\d+) couplings\)", section).group(1))
    count = re.search(r"It is an EP2 \(.*?\) at (\d+) of the (\d+) couplings:", section, re.DOTALL)
    extremes = re.search(r"sigma_1 at most ([0-9.e+-]+); sigma_2 at least ([0-9.e+-]+); the third eigenvalue at "
                         r"least ([0-9.e+-]+) away; the\s+half-split / sqrt\|delta\| below and above equal to "
                         r"([0-9.e+-]+) \(relative\)", section)
    exceptions = re.findall(r"not an EP2 at J_bridge = ([0-9.]+) \((.+?)\): gamma\* = ([0-9.e+-]+), pair gap "
                            r"([0-9.e+-]+), sigma_1 = ([0-9.e+-]+), sigma_2 = ([0-9.e+-]+)(.*)$", section, re.MULTILINE)
    unread = re.findall(r"not read at J_bridge = ([0-9.]+)", section)
    return grid, count, extremes, exceptions, unread


def test_jordan_rank_is_read_at_every_sweep_coupling():
    grid, count, extremes, exceptions, unread = _sweep_rank_read()
    ep2, total = int(count.group(1)), int(count.group(2))
    # every grid coupling is accounted for: read as an EP2, or listed with its reading
    assert total == grid
    assert ep2 + len(exceptions) + len(unread) == total
    assert not unread
    # every exception is a zero of the threshold (section 3), where the pair already coincides at
    # zero gain with two null directions, and every zero that lies on the grid is an exception
    _, zero_rows = _zero_rows()
    zeros = [float(j0) for j0, *_ in zero_rows]
    on_grid = [z for z in zeros if abs(z / 0.005 - round(z / 0.005)) < 1e-9]
    listed = [float(j) for j, *_ in exceptions]
    for j, _, gamma_star, _, _, sigma_2, note in exceptions:
        assert any(abs(float(j) - z) < 1e-9 for z in zeros), j
        assert float(gamma_star) < 1e-12 and float(sigma_2) < 1e-12, (j, gamma_star, sigma_2)
        assert "two null directions at zero gain" in note, j
    assert sorted(listed) == sorted(on_grid)
    # at the EP2 couplings: one null direction at the rounding level with the next six decades above
    # it across the whole sweep, the third eigenvalue far from the pair, and one Puiseux branch
    sigma_1_max, sigma_2_min, third_min, mismatch = (float(x) for x in extremes.groups())
    assert sigma_1_max <= 1e-6 * sigma_2_min, (sigma_1_max, sigma_2_min)
    assert third_min >= 1e6 * sigma_1_max, (third_min, sigma_1_max)
    assert mismatch < 1e-3, mismatch


def _large_coupling_rows():
    section = _artifact().split("7. The product gamma_crit * J_bridge at large couplings", 1)[1]
    return re.findall(r"^\s+(\d+)\s+(\(\d,\d\)(?: flip-(?:even|odd))?)\s+([0-9.]+)\s*$", section, re.MULTILINE)


def test_large_coupling_samples_agree_with_the_sweep():
    rows = _large_coupling_rows()
    assert [int(j) for j, _, _ in rows] == [20, 30, 50, 100, 300, 1000]
    # Computed twice: at J_bridge = 20 the product is the sweep's last threshold times 20, equal to
    # within the rounding of the two printed numbers (half a unit in the fourth and the sixth decimal).
    tail = re.search(r"^\s+\(1,2\)\s+4\.460\s+20\.000\s+2\s+falls from ([0-9.]+) to ([0-9.]+)$", _artifact(), re.MULTILINE)
    assert rows[0][1] == "(1,2)"
    assert abs(float(rows[0][2]) - 20.0 * float(tail.group(2))) <= 5e-5 + 20.0 * 5e-7


def _four_site_liouvillian_parts(j_bridge):
    """Two 2-qubit Heisenberg chains (J = 1) and a Heisenberg bridge: a 4-site chain with bonds
    (1, J_b, 1), Pauli convention, row-major vec; chain A dephases at +gamma, chain B at -gamma."""
    pauli = (np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]], dtype=complex),
             np.array([[1, 0], [0, -1]], dtype=complex))

    def site(op, k):
        out = np.eye(1, dtype=complex)
        for s in range(4):
            out = np.kron(out, op if s == k else np.eye(2))
        return out

    H = sum(J * site(P, b) @ site(P, b + 1) for b, J in enumerate((1.0, j_bridge, 1.0)) for P in pauli)
    identity = np.eye(16)
    LH = -1j * (np.kron(H, identity) - np.kron(identity, H.T))
    z = [np.real(np.diag(site(pauli[2], k))) for k in range(4)]
    D = np.diag(sum(s * (np.kron(z[k], z[k]) - 1.0) for k, s in enumerate((1.0, 1.0, -1.0, -1.0))))
    return LH, D


def test_recomputed_collision_matches_the_artifact_and_is_defective():
    LH, D = _four_site_liouvillian_parts(1.0)
    block = _certificate_blocks()[1.0]
    gamma_artifact = _float(r"collision gamma\* = ([0-9.]+)", block)

    def pair(gamma, anchor):
        w = np.linalg.eigvals(LH + gamma * D)
        order = np.argsort(np.abs(w - anchor))
        return w[order[0]], w[order[1]], w

    g, lam = 0.1873101, complex(0.0, -2.6515)
    for _ in range(40):
        a, b, _ = pair(g, lam)
        lam = 0.5 * (a + b)
        ap, bp, _ = pair(g + 1e-7, lam)
        am, bm, _ = pair(g - 1e-7, lam)
        step = ((a - b) ** 2).real / ((((ap - bp) ** 2).real - ((am - bm) ** 2).real) / 2e-7)
        g -= step
        if abs(step) <= 1e-14 * g:
            break
    a, b, w = pair(g, lam)
    lam = 0.5 * (a + b)
    sv = np.sort(np.linalg.svd(LH + g * D - lam * np.eye(256), compute_uv=False))
    micro_a, micro_b, _ = pair(g * (1 + 1e-6), lam)

    # Located on the collision (by the sqrt-law, within 1e-10 in relative gamma) and at the artifact's
    # gamma*. Both locators find the zero of the same analytic Re f, whose rounding noise moves that
    # zero by its noise over |f'| (8e-16 at J_bridge = 1.0, 2.6e-15 at 1.9, rms), budgeted at 1e-13;
    # the artifact prints gamma* to 12 decimals, which adds half a unit of the last one.
    assert abs(a - b) <= 1e-2 * abs(micro_a - micro_b)
    assert abs(g - gamma_artifact) <= 0.5e-12 + 1e-13
    # One null direction at the rounding level, none at the pair's scale.
    assert sv[0] <= 10 * np.finfo(float).eps * sv[-1]
    assert sv[1] >= 1e3 * abs(a - b)


def test_recomputed_zero_at_three_quarters_matches_the_artifact():
    # An independent route to the first zero: at J_bridge = 3/4 the full L leaves the axis
    # linearly, with the slope the artifact prints; at J_bridge = 1.0 it does not.
    _, rows = _zero_rows()
    printed = {round(float(j0), 6): float(m4) for j0, _, _, m4, _ in rows}
    LH, D = _four_site_liouvillian_parts(0.75)
    slopes = [np.max(np.linalg.eigvals(LH + g * D).real) / g for g in (1e-4, 1e-6)]
    assert abs(slopes[0] / slopes[1] - 1.0) < 1e-6
    assert abs(slopes[0] - printed[0.75]) < 2e-6
    LH, D = _four_site_liouvillian_parts(1.0)
    assert np.max(np.linalg.eigvals(LH + 1e-4 * D).real) < 1e-12


def test_signature_producer_and_artifact_state_what_was_computed():
    for path in (SIGNATURE_SCRIPT, SIGNATURE_RESULT):
        text = path.read_text(encoding="utf-8")
        assert "At J_bridge = 1.0 and 1.9: geometric multiplicity 1 at algebraic multiplicity 2" in text
        assert "In every row the full L holds one Jordan pair per copy" in text
        assert "gamma_crit = 0 there, and any small gain makes the bridge unstable" in text
        assert "sits slightly below the collision" in text
        assert "single-vector K not reported for a degenerate eigenspace" in text
        assert "near-threshold simplicity gate: nearest gap > 1e-10" in text
        assert "gap to the across-axis partner" in text
        assert "These are samples; they establish no limit." in text
        assert "two qubits per chain" in text.lower() or "Two chains of 2 qubits" in text
    source = SIGNATURE_SCRIPT.read_text(encoding="utf-8")
    assert "if gap <= SIMPLE_MODE_GAP_TOL:" in source
    assert "raise RuntimeError(" in source


def test_phase3_producer_fixes_the_centre_by_the_trace():
    source = PHASE3_SCRIPT.read_text(encoding="utf-8")
    assert "np.trace(L) / L.shape[0]" in source
    assert "unique candidate center fixed by the trace" in source


def test_entry_layer_pages_carry_the_scope_and_the_zeros():
    found = WHAT_WE_FOUND.read_text(encoding="utf-8")
    bridge = _normalized(_section(found, "### The stability window is finite", "### The neural palindrome"))
    row = _normalized([line for line in GLOSSARY.read_text(encoding="utf-8").splitlines()
                       if line.startswith("| **Fragile Bridge** |")][0])
    readme = _normalized([line for line in (ROOT / "hypotheses" / "README.md").read_text(encoding="utf-8").splitlines()
                          if line.startswith("- [Fragile Bridge]")][0])
    for text in (bridge, row):
        assert "At J_bridge = 1.0 and 1.9" in text
        for zero in ("3/4", "√5/2", "4/3"):
            assert zero in text, zero
    for text in (bridge, readme):
        assert "two qubits per chain" in text
    assert "sweep of two coupled two-qubit chains" in row


def test_formula_entries_carry_the_certified_number():
    formulas = FORMULAS.read_text(encoding="utf-8")
    f19 = " ".join(_section(formulas, "### F19.", "### F20.").split())
    f40 = " ".join(_section(formulas, "### F40.", "### F41.").split())
    assert "0.187310108345" in f40
    assert "derived symmetry identity, not a third independent trend" in f19


def test_fragile_quotes_what_the_producer_prints():
    artifact = _artifact()
    fragile = " ".join(FRAGILE.read_text(encoding="utf-8").split())
    for value in re.findall(r"collision gamma\* = ([0-9.]+);", artifact):
        assert value in fragile, value
    crossing = re.search(r"cross at J_bridge = ([0-9.]+), gamma = ([0-9.]+);", artifact)
    for value in crossing.groups():
        assert value in fragile, value
    assert "flip-even at 1.0 and flip-odd at 1.9" in fragile
    # the (0,1) window, rounded to four decimals, and its counts at 1.46
    window = re.search(r"between ([0-9.]+) \(where it meets the \(1,2\) orbit.*?and ([0-9.]+) \(where",
                       artifact, re.DOTALL)
    for value in window.groups():
        assert "{0:.4f}".format(float(value)) in fragile, value
    row = re.search(r"^\s+1\.46\s+\(0,1\)\s+\S+\s+4\s+[0-9.]+\s+[0-9.]+\s+8\s+4\s+8\s*$", artifact, re.MULTILINE)
    assert row, "the 1.46 row"
    assert "eight eigenvalues within 10⁻⁵" in fragile and "four quartets" in fragile
    # the threshold frequencies of the table the page shows
    section2 = _section(artifact, "2. Sector structure", "3. Zeros of the threshold")
    omegas = {float(j): omega for j, _, _, _, _, omega, _, _, _ in _sector_rows(section2)}
    for j_bridge in (1.46, 1.5, 1.9, 2.0, 5.0, 10.0):
        assert "| {0} |".format(omegas[j_bridge]) in fragile, (j_bridge, omegas[j_bridge])
    # the zeros: exact forms as the page writes them, in x = J_bridge / J
    forms = re.findall(r"real root of ([0-9x^ +-]+?)\s{2,}", artifact)
    assert len(forms) == 2
    for form in forms:
        page = form.strip().replace("^3", "³").replace("^2", "²").replace(" - ", " − ")
        assert page in fragile, page
    assert "x = J_bridge/J" in fragile
    for zero in ("3/4", "√5/2", "4/3"):
        assert zero in fragile, zero
    # the finer grid's low points, as the zeros section recomputes them
    for j_bridge in ("0.70", "0.80", "1.10", "1.30"):
        value = float(re.search(r"^\s+" + re.escape(j_bridge) + r"\s+([0-9.]+)\s", artifact, re.MULTILINE).group(1))
        assert "{0:.4f}".format(value) in fragile, (j_bridge, value)
    # the takeover beyond the last zero and the two slopes of that rise
    take = re.search(r"at J_bridge = ([0-9.]+), gamma = ([0-9.]+)\. It leaves the zero at the V slope ([0-9.]+) "
                     r"and rises on\s+average at ([0-9.]+) per unit", artifact)
    assert "{0:.4f}".format(float(take.group(1))) in fragile
    assert "{0:.4f}".format(float(take.group(2))) in fragile
    assert "{0:.2f} per unit".format(float(take.group(3))) in fragile
    assert "{0:.3f} per".format(float(take.group(4))) in fragile
    # the sweep: the largest grid value, where it sits, and the tail run
    largest = re.search(r"The largest threshold on the grid is ([0-9.]+), at J_bridge = ([0-9.]+);", artifact)
    assert largest.group(1) in fragile and largest.group(2) in fragile
    tail = re.search(r"^\s+\(1,2\)\s+4\.460\s+20\.000\s+2\s+falls from ([0-9.]+) to ([0-9.]+)$", artifact, re.MULTILINE)
    for value in tail.groups():
        assert "{0:.4f}".format(float(value)) in fragile, value
    below = re.search(r"is ([0-9.]+)% below that crossing value", artifact).group(1)
    assert below + "%" in fragile
    # the sweep's Jordan-rank read: its size, the exception, and the extremes the page quotes
    grid, _, extremes, exceptions, _ = _sweep_rank_read()
    assert "{0} couplings".format(grid) in fragile
    assert [float(j) for j, *_ in exceptions] == [0.75] and "the grid point 3/4" in fragile
    sigma_1_max, sigma_2_min, third_min, _ = extremes.groups()
    for printed in (sigma_1_max, sigma_2_min, third_min):
        assert _page_number(printed) in fragile, (printed, _page_number(printed))
    # the large-coupling samples, in open question 4
    for _, _, product in _large_coupling_rows():
        assert product in fragile, product


def test_bridge_optimum_and_large_coupling_samples_are_quoted():
    fragile = FRAGILE.read_text(encoding="utf-8")
    found = WHAT_WE_FOUND.read_text(encoding="utf-8")
    combined = fragile + found

    assert "J_bridge/J in [1.8, 2.0]" in found
    assert "0.578 at J_bridge=10" in combined
    assert "0.508 at J_bridge=100" in combined


def test_derived_quantities_are_not_counted_as_evidence():
    fragile_normalized = " ".join(FRAGILE.read_text(encoding="utf-8").split())
    signature = SIGNATURE_SCRIPT.read_text(encoding="utf-8") + SIGNATURE_RESULT.read_text(encoding="utf-8")
    formulas = FORMULAS.read_text(encoding="utf-8")
    assert "derived from max Re and supplies no independent EP evidence" in signature
    assert "derived partner gap adds no independent evidence" in fragile_normalized
    assert "derived symmetry identity, not a third independent trend" in formulas
