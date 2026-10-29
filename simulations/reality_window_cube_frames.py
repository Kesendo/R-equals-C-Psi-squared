"""The reality window read on the letter cube, one site, frame by frame.

Second route for the by-hand reading in experiments/REALITY_WINDOW_HARDWARE_PREDICTION.md
section 7 (Read on the cube). One site, letters I, X, Y, Z; the real sector is span{I, X, Z}
(n_Y even), Y the odd letter. A physical one-site generator P that acts when the compiled
circuit has applied the gate product U is seen in the logical frame as Ad_U^dagger o P o Ad_U.
Which carried generators move weight between {I, X, Z} and Y (an odd-Y arrow), exactly, and
at what rate: a Pauli dephasing (a height) in the Clifford frames never, in the frame sx.rz(theta)
through the arrow X <-> Y at sin(2 theta) times the rate; the damping sigma^- (lights plus the
arrow I -> Z along [J, J^dagger]) in the frame sx through I -> Y at its full rate (I -> -Y for the inverse
turn), in the frame h through no odd-Y arrow (the arrow lands on X); a field through the arrow between the two letters perpendicular to its
carried axis, at twice its strength and at no rate.

Exact in sympy; every row compares to 0; two anti-vacuity rows pin the tilted arrow's dependence on the
angle (it vanishes at theta = pi/2 and has modulus 1 at theta = pi/4, so a row that only saw zero could not pass
both). Must print
"cube frame gate: ALL GREEN". Run from the repo root:  python simulations/reality_window_cube_frames.py
"""
import sympy as sp

I2 = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.Matrix([[1, 0], [0, -1]])
P = {'I': I2, 'X': X, 'Y': Y, 'Z': Z}
names = ['I', 'X', 'Y', 'Z']

def ptm(gen):
    """Pauli-transfer matrix of a superoperator gen (a callable on 2x2), rows = output letter."""
    M = sp.zeros(4, 4)
    for j, b in enumerate(names):
        out = gen(P[b])
        for i, a in enumerate(names):
            M[i, j] = sp.simplify((P[a] * out).trace() / 2)
    return M

def dissipator(J, rate=1):
    Jd = J.H
    return lambda rho: rate * (J * rho * Jd - (Jd * J * rho + rho * Jd * J) / 2)

def field(Hm):
    return lambda rho: -sp.I * (Hm * rho - rho * Hm)

def transported(gen, U):
    """The generator seen in the logical frame when the physical gen acts after the gates U."""
    return lambda rho: U.H * gen(U * rho * U.H) * U

th = sp.symbols('theta', real=True)
sx = (I2 - sp.I * X) / sp.sqrt(2)            # exp(-i pi X / 4)
rz = lambda t: sp.diag(sp.exp(-sp.I * t / 2), sp.exp(sp.I * t / 2))
sminus = (X + sp.I * Y) / 2                     # |0><1|, Z|0> = |0>

def odd_y_arrows(M):
    """Entries that carry a real letter into Y or Y into a real letter."""
    iy = names.index('Y')
    into = {names[j]: sp.simplify(M[iy, j]) for j in range(4) if j != iy and sp.simplify(M[iy, j]) != 0}
    outof = {names[i]: sp.simplify(M[i, iy]) for i in range(4) if i != iy and sp.simplify(M[i, iy]) != 0}
    return into, outof

checks = []
def check(label, cond):
    checks.append((label, bool(cond)))
    print(('PASS ' if cond else 'FAIL ') + label)

# 1. Z dephasing is a height in the logical frame: no arrow, in the frame rz(th) (any angle) and in the frame sx.
Dz = dissipator(Z)
for U, lab in [(I2, 'I'), (rz(th), 'rz(theta)'), (sx, 'sx'), (rz(th) * sx, 'rz(theta)*sx (sx first, then rz)')]:
    M = ptm(transported(Dz, U))
    into, outof = odd_y_arrows(M)
    check(f"Z dephasing in frame {lab}: no odd-Y arrow (into={into}, outof={outof})", not into and not outof)
    offdiag = [(names[i], names[j], M[i, j]) for i in range(4) for j in range(4) if i != j and sp.simplify(M[i, j]) != 0]
    check(f"Z dephasing in frame {lab}: pure height (no off-diagonal at all: {offdiag})", not offdiag)

# 2. In the frame sx*rz(th) (rz first, then sx) the dephasing axis is tilted in the XY plane:
#    the arrow X <-> Y appears with rate sin(2 theta) (the dissipator rate is 1 here; D[Z] has heights -2).
M = ptm(transported(Dz, sx * rz(th)))
into, outof = odd_y_arrows(M)
print('Z dephasing in frame sx*rz(theta):'); sp.pprint(sp.simplify(M))
check("tilted dephasing: arrow X->Y has modulus |sin 2theta|", sp.simplify(sp.Abs(M[2, 1]) - sp.Abs(sp.sin(2 * th))) == 0 or sp.simplify(M[2, 1]**2 - sp.sin(2 * th)**2) == 0)
check("tilted dephasing: no arrow out of I or Z into Y", 'I' not in into and 'Z' not in into)
check("tilted dephasing: Z height unchanged (-2)", sp.simplify(M[3, 3] + 2) == 0)
# control that must fail: at theta = pi/2 the axis is a Pauli axis again and the arrow vanishes
check("control: at theta = pi/2 the arrow vanishes", sp.simplify(M[2, 1].subs(th, sp.pi / 2)) == 0)
check("control: at theta = pi/4 the arrow is maximal (modulus 1)", sp.simplify(M[2, 1].subs(th, sp.pi / 4)**2 - 1) == 0)

# 3. sigma^- : logical frame = two lights + arrow I -> Z (real); in frame sx the arrow goes I -> +-Y (odd-Y).
Dm = dissipator(sminus)
M0 = ptm(Dm)
print('sigma^- logical frame:'); sp.pprint(M0)
check("sigma^- logical: arrow I->Z = 1 (rate), no odd-Y arrow", sp.simplify(M0[3, 0] - 1) == 0 and not any(odd_y_arrows(M0)))
Msx = ptm(transported(Dm, sx))
print('sigma^- in frame sx:'); sp.pprint(Msx)
into, outof = odd_y_arrows(Msx)
check(f"sigma^- in frame sx: arrow I->Y of modulus 1 (into={into})", 'I' in into and sp.simplify(into['I']**2 - 1) == 0)
check("sigma^- in frame sx: lights on X and Z (heights -1/2 each), Y height -1", sp.simplify(Msx[1, 1] + sp.Rational(1, 2)) == 0 and sp.simplify(Msx[3, 3] + sp.Rational(1, 2)) == 0 and sp.simplify(Msx[2, 2] + 1) == 0)
Mh = ptm(transported(Dm, (X + Z) / sp.sqrt(2)))
check("sigma^- in frame h (a real gate): no odd-Y arrow", not any(odd_y_arrows(Mh)))

# 4. Coherent letters: a Z field (detuning) is an X<->Y arrow in the logical frame (real-breaking) at rate 2h
#    in any rz frame; in the sx frame it becomes a Y field: X<->Z arrow, real-PRESERVING.
Fz = field(Z)
M = ptm(transported(Fz, rz(th)))
check("Z field in frame rz(theta): odd-Y arrows X<->Y of modulus 2", sp.simplify(M[2, 1]**2 - 4) == 0 and sp.simplify(M[1, 2]**2 - 4) == 0)
M = ptm(transported(Fz, sx))
check("Z field in frame sx: NO odd-Y arrow (it is a Y field, real-preserving)", not any(odd_y_arrows(M)))
M = ptm(transported(Fz, sx * rz(th)))
into, outof = odd_y_arrows(M)
check("Z field in frame sx*rz(theta): a field in the XY plane; odd-Y arrow Z<->Y of modulus 2|sin theta| (its X component), X<->Z of 2|cos theta|, NO X<->Y arrow", sp.simplify(M[2, 3]**2 - 4 * sp.sin(th)**2) == 0 and sp.simplify(M[1, 3]**2 - 4 * sp.cos(th)**2) == 0 and sp.simplify(M[2, 1]) == 0)
# an X field (drive) in a frame rz(phi): tilted field, odd-Y part cos(phi)
ph = sp.symbols('phi', real=True)
M = ptm(transported(field(X), rz(ph)))
check("X field in frame rz(phi): arrow Z<->Y of modulus 2|cos phi| (odd-Y), Z<->X of 2|sin phi|", sp.simplify(M[2, 3]**2 - 4 * sp.cos(ph)**2) == 0 and sp.simplify(M[1, 3]**2 - 4 * sp.sin(ph)**2) == 0)

# 5. The law in one line: a Pauli channel along ANY axis n maps real to real iff n is in the XZ plane or n = Y.
nx, ny, nz = sp.symbols('n_x n_y n_z', real=True)
n = nx * X + ny * Y + nz * Z
M = ptm(dissipator(n))
into, outof = odd_y_arrows(M)
print('Pauli channel along n: odd-Y arrows', into, outof)
check("general axis: the only odd-Y arrows are X<->Y (2 n_x n_y) and Z<->Y (2 n_y n_z); I never moves",
      set(into) <= {'X', 'Z'} and sp.simplify(into.get('X', 0) - 2 * nx * ny) == 0 and sp.simplify(into.get('Z', 0) - 2 * ny * nz) == 0)

fails = [l for l, ok in checks if not ok]
print(f"\n{len(checks)} checks, {len(fails)} failed")
print('cube frame gate: ALL GREEN' if not fails else 'cube frame gate: RED ' + str(fails))
