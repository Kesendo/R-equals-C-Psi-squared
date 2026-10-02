"""Exact checks for the carriers of the section "Blind dephased seats, and palindromes beyond colourings and
site symmetries" in docs/proofs/PROOF_PALINDROME_COMPLEMENT_CONNECTION.md. Run from the repository root:
    python simulations/blind_seat_carriers_gate.py
N = 4: the far element of the pair-channel family in closed form, reconstructed as Pauli strings separately from
its singlet/triplet derivation; [H, W] = 0, {Z0, W} = 0 and W^2 = (a^2 + h^2) I symbolically. N = 5: the 19 x 19
minor that forbids keeping W4 (x) Q when a fifth spin is attached. N = 7: two weighted triangles equivalent as
spin representations, the full 128-dimensional carrier, the trivial weighted automorphism group, the Laplacian.
No producer is imported, no census run, no eigenvalue tolerance used.
"""
import itertools
from functools import lru_cache

import sympy as sp


I2 = sp.eye(2)
PAULI = {
    "I": I2,
    "X": sp.Matrix([[0, 1], [1, 0]]),
    "Y": sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    "Z": sp.diag(1, -1),
}
CHECKS = 0


def check(name, condition):
    global CHECKS
    if not condition:
        raise AssertionError(name)
    CHECKS += 1
    print("PASS", name, flush=True)


def clean(poly):
    return {w: v for w, c in poly.items() if (v := sp.expand(c)) != 0}


def multiply(left, right):
    out = {}
    for u, x in left.items():
        for v, y in right.items():
            word, phase = [], 1
            for p, q in zip(u, v):
                if p == "I":
                    word.append(q)
                elif q == "I":
                    word.append(p)
                elif p == q:
                    word.append("I")
                else:
                    word.append(next(c for c in "XYZ" if c not in (p, q)))
                    phase *= sp.I if p + q in ("XY", "YZ", "ZX") else -sp.I
            label = "".join(word)
            out[label] = out.get(label, 0) + phase * x * y
    return clean(out)


def combine(left, right, sign=1):
    out = dict(left)
    for w, c in right.items():
        out[w] = out.get(w, 0) + sign * c
    return clean(out)


def commutator(left, right):
    return combine(multiply(left, right), multiply(right, left), -1)


def terms(n, bonds, fields):
    out = {}
    for (j, k), c in bonds.items():
        for letter in "XYZ":
            w = ["I"] * n
            w[j] = w[k] = letter
            label = "".join(w)
            out[label] = out.get(label, 0) + c
    for j, letter, c in fields:
        w = ["I"] * n
        w[j] = letter
        label = "".join(w)
        out[label] = out.get(label, 0) + c
    return clean(out)


@lru_cache(None)
def matrix(word):
    return sp.SparseMatrix(sp.kronecker_product(*(PAULI[c] for c in word)))


def dense(poly):
    n = len(next(iter(poly)))
    out = sp.SparseMatrix(2**n, 2**n, {})
    for word, c in poly.items():
        out += c * matrix(word)
    return out


def zero_matrix(value):
    return all(sp.expand(c) == 0 for c in value.todok().values())


def n4_gate():
    a, e, h, g = sp.symbols("a e h g", real=True)
    d = 4*a*a + h*h
    # Q = 2*d*W in site order 0,1,2,3. No denominators or branch divisions.
    Q = {
        "XIII": h*(2*a*a+h*h), "XIXI": -2*a*h*h,
        "XIYY": -4*a*a*h, "XXII": -2*a*h*h,
        "XXXI": -h*(2*a*a+h*h), "XXYY": -2*a*h*h,
        "XYIY": -4*a*a*h, "XYXY": -2*a*h*h,
        "XYYI": -h*(2*a*a-h*h), "XZZI": -h*(2*a*a+h*h),
        "YIIY": a*h*h, "YIXY": -4*a*a*h,
        "YIYI": 2*a*h*h, "YXIY": -4*a*a*h,
        "YXXY": -a*h*h, "YXYI": -4*a*a*h,
        "YYII": 2*a*h*h, "YYXI": -4*a*a*h,
        "YYYY": -a*(8*a*a-h*h), "YZZY": -a*h*h,
    }
    H = terms(4, {(0,1): a, (0,2): -a, (1,2): e,
                  (1,3): -a, (2,3): a}, [(1,"X",-h),(2,"X",h),(3,"Y",g)])
    check("N4 symbolic [H,Q]=0 for free a,e,h,g", not commutator(H, Q))
    check("N4 symbolic {Z0,Q}=0", not combine(multiply({"ZIII":1},Q), multiply(Q,{"ZIII":1})))
    check("N4 symbolic Q^2=4*d^2*(a^2+h^2)*I",
          not combine(multiply(Q,Q), {"IIII":4*d*d*(a*a+h*h)}, -1))
    check("N4 [Y3,Q]=0", not commutator({"IIIY":1}, Q))
    check("N4 [sigma1.dot(sigma2),Q]=0", not commutator(terms(4,{(1,2):1},[]),Q))
    check("N4 h=0 gives Q=-8*a^3*YYYY", clean({w:c.subs(h,0) for w,c in Q.items()}) == {"YYYY":-8*a**3})

    # Reconstruct the displayed dimer-block formula separately from the strings.
    ox = sp.kronecker_product(PAULI["X"],I2)
    oy = sp.kronecker_product(PAULI["Y"],I2)
    oxy = sp.kronecker_product(PAULI["X"],PAULI["Y"])
    oyy = sp.kronecker_product(PAULI["Y"],PAULI["Y"])
    S = h*ox+a*oyy
    An = h**3*ox+a*(h*h-4*a*a)*oyy
    Bn = 2*a*h*h*oxy+4*a*a*h*oy
    Cn = sp.I*(-4*a*a*h*oxy+2*a*h*h*oy)
    En = sp.I*(2*a*h*h*ox+4*a*a*h*oyy)
    triplet = [[An,Bn,Cn],[Bn,-An,En],[-Cn,-En,An]]
    singlet = sp.Matrix([0,1,-1,0])/sp.sqrt(2)
    cartesian = [sp.kronecker_product(PAULI[p],I2)*singlet for p in "XYZ"]
    Wd = sp.kronecker_product(d*S,singlet*singlet.H)
    for i in range(3):
        for j in range(3):
            Wd += sp.kronecker_product(triplet[i][j],cartesian[i]*cartesian[j].H)
    Qreordered = sum((coef*matrix(word[0]+word[3]+word[1]+word[2]) for word,coef in Q.items()),sp.zeros(16))
    check("N4 displayed dimer block equals the independently checked polynomial", zero_matrix(2*Wd-Qreordered))
    channels = [a*(sp.kronecker_product(PAULI[p],I2)-sp.kronecker_product(I2,PAULI[p]))
                -(h*sp.eye(4) if p=="X" else sp.zeros(4)) for p in "XYZ"]
    gram = sum((k*k for k in channels),sp.zeros(4))
    check("N4 singlet map commutes with channel Gram matrix", zero_matrix(gram*S-S*gram))
    span = [sp.kronecker_product(PAULI[p0],PAULI[p3]) for p0,p3 in (("X","I"),("X","Y"),("Y","I"),("Y","Y"))]
    cs = sp.symbols("c0:4")
    trial = sum((c*m for c,m in zip(cs,span)),sp.zeros(4))
    point = {a:sp.Rational(13,10), h:sp.Rational(7,10)}
    eqs = list((gram*trial-trial*gram).subs(point))
    sol = sp.linsolve(eqs, cs)
    free = [c for c in cs if any(c in v.free_symbols for v in next(iter(sol)))]
    check("N4 inside span{X0,X0Y3,Y0,Y0Y3} the Gram commutant is one-dimensional (fixes S up to a factor)", len(free)==1)

    sample = {a:2, e:-3, h:5, g:7}
    Hr = {w:c.subs(sample) for w,c in H.items()}
    Qr = {w:c.subs(sample) for w,c in Q.items()}
    Hm, Qm = dense(Hr), dense(Qr)
    check("N4 independent tensor-matrix commutator", zero_matrix(Hm*Qm-Qm*Hm))
    check("N4 independent tensor-matrix square", zero_matrix(Qm*Qm-4*41**2*29*sp.eye(16)))
    check("N4 the carrier is Hermitian at the sample point", zero_matrix(Qm.H-Qm))
    # Each control changes the guarded physical object or claimed carrier.
    J03 = dense(terms(4,{(0,3):1},[]))
    check("N4 fixed carrier rejects J03=1", not zero_matrix((Hm+J03)*Qm-Qm*(Hm+J03)))
    X2 = matrix("IIXI")
    check("N4 fixed carrier rejects a non-opposite pair field", not zero_matrix((Hm+X2)*Qm-Qm*(Hm+X2)))
    broken = Qm-Qr["XIII"]*matrix("XIII")
    check("N4 omitted carrier term is detected", not zero_matrix(Hm*broken-broken*Hm))
    return Q


def n7_gate():
    dot12 = dense(terms(3,{(0,1):1},[]))
    dot13 = dense(terms(3,{(0,2):1},[]))
    dot23 = dense(terms(3,{(1,2):1},[]))
    E = sp.SparseMatrix(sp.eye(8))
    A = dot12+2*dot13+3*dot23
    B = (11*dot12+22*dot13+9*dot23)/7
    P = (dot12+dot13+dot23+3*E)/6
    D = (-10*dot12+8*dot13+2*dot23)/7
    T, Ti = P+D, P+sp.Rational(7,144)*D
    check("triangle P is a projector and PD=DP=0", zero_matrix(P*P-P) and zero_matrix(P*D) and zero_matrix(D*P))
    check("triangle D^2=(144/7)*(I-P)", zero_matrix(D*D-sp.Rational(144,7)*(E-P)))
    check("triangle T has the displayed two-sided inverse", zero_matrix(T*Ti-E) and zero_matrix(Ti*T-E))
    check("triangle B*T=T*A", zero_matrix(B*T-T*A))
    for axis in "XYZ":
        total = dense(terms(3,{},[(j,axis,1) for j in range(3)]))
        check("triangle T commutes with total "+axis, zero_matrix(total*T-T*total))
    U = P+sp.sqrt(7)*D/12
    check("triangle U is a Hermitian unitary intertwiner", zero_matrix(U.H-U) and zero_matrix(U*U-E) and zero_matrix(B*U-U*A))
    check("triangle detuned B rejects the same T", not zero_matrix((B+dot12/7)*T-T*A))

    # A different (c,h) from the discovering scout; no Liouvillian is formed.
    c, h = 4, 1
    bonds = {(0,j):sp.Integer(c) for j in range(1,7)}
    bonds.update({(1,2):sp.Integer(1),(1,3):sp.Integer(2),(2,3):sp.Integer(3),
                  (4,5):sp.Rational(11,7),(4,6):sp.Rational(22,7),(5,6):sp.Rational(9,7)})
    H = dense(terms(7,bonds,[(j,"X",h) for j in (1,2,3)]+[(j,"Y",h) for j in (4,5,6)]))
    V = sp.SparseMatrix(sp.kronecker_product(sp.eye(16),T))
    Vi = sp.SparseMatrix(sp.kronecker_product(sp.eye(16),Ti))
    swap = sp.SparseMatrix(128,128,{})
    for x in range(128):
        centre, left, right = (x>>6)&1, (x>>3)&7, x&7
        swap[(centre<<6)|(right<<3)|left,x] = 1
    rotation = sp.SparseMatrix(sp.kronecker_product(*([PAULI["X"]+PAULI["Y"]]*7)))
    W = (V*swap*rotation*Vi).applyfunc(sp.expand)
    check("N7 full-space [H,W]=0", zero_matrix(H*W-W*H))
    Z0 = matrix("ZIIIIII")
    check("N7 full-space {Z0,W}=0", zero_matrix(Z0*W+W*Z0))
    check("N7 full-space W^2=128*I", zero_matrix(W*W-128*sp.eye(128)))
    brokenH = H+dense(terms(7,{(4,5):sp.Rational(1,7)},[]))
    check("N7 full-space fixed carrier detects detuned triangle", not zero_matrix(brokenH*W-W*brokenH))
    L = sp.zeros(7)
    for (j,k),w in bonds.items():
        L[j,j] += w
        L[k,k] += w
        L[j,k] -= w
        L[k,j] -= w
    check("N7 positive bonds and exact Laplacian rank 6", all(w>0 for w in bonds.values()) and L.rank()==6)
    automorphisms = []
    for tail in itertools.permutations(range(1,7)):
        perm = (0,)+tail
        if all(bonds.get(tuple(sorted((perm[j],perm[k]))),0)==bonds.get((j,k),0)
               for j in range(7) for k in range(j+1,7)):
            automorphisms.append(perm)
    check("N7 weighted automorphisms fixing 0: identity only", automorphisms==[tuple(range(7))])
    # Connected Heisenberg bonds force a Pauli colouring to be one global letter.
    check("N7 no global Pauli colouring", all(not zero_matrix(H*matrix(p*7)-matrix(p*7)*H) for p in "XY"))
    v = sp.eye(7)[:,0]
    krylov_rank = sp.Matrix.hstack(*(L**k*v for k in range(7))).rank()
    check("N7 watched seat is still blind: Krylov rank 2", krylov_rank==2)
    print("N7 c=4, h=1: Laplacian cofactor =", L[1:,1:].det(), flush=True)


def n5_attachment_gate(Q):
    t = sp.symbols("t", real=True)
    a, h = sp.symbols("a h", real=True)
    q = {w:sp.expand(c.subs({a:1,h:t})) for w,c in Q.items()}
    columns = []
    for j in range(4):
        edge = terms(5,{(j,4):1},[])
        for p in "IXYZ":
            extended = {w+p:c for w,c in q.items()}
            columns.append({w:sp.expand(c/(2*sp.I)) for w,c in commutator(edge,extended).items()})
    # An arbitrary spectator-field commutator lies in these three columns.
    columns.extend({w+p:c for w,c in q.items()} for p in "XYZ")
    words = sorted(set().union(*(col.keys() for col in columns)))
    F = sp.MutableSparseMatrix(len(words),19,{})
    for i,w in enumerate(words):
        for j,col in enumerate(columns):
            if w in col:
                F[i,j] = col[w]
    determinants = []
    selected_words = []
    for label, order in (("forward", list(range(len(words)))), ("reversed", list(reversed(range(len(words)))))):
        pivots = F[order,:].subs(t,3).T.rref()[1]
        check("N5 attachment matrix has 19 pivot columns at t=3 (rows in "+label+" order)", len(pivots)==19)
        rows = [order[i] for i in pivots]
        determinants.append(sp.factor(F[rows,:].det(method="domain-ge")))
        selected_words.append([words[i] for i in rows])
    common = sp.Poly(sp.gcd(*determinants),t).monic().as_expr()
    print("N5 first symbolic minor =", determinants[0], flush=True)
    print("N5 second symbolic minor =", determinants[1], flush=True)
    print("N5 second minor rows =", ",".join(selected_words[1]), flush=True)
    print("N5 gcd of symbolic minors =", sp.factor(common), flush=True)
    # Count real roots after removing the only permitted exceptional point t=0.
    poly = sp.Poly(common,t)
    while poly.eval(0)==0:
        poly = poly.exquo(sp.Poly(t,t))
    check("N5 minors have no common nonzero real root", poly.count_roots(-sp.oo,sp.oo)==0)
    check("N5 displayed minor = 256*t^34*(t^2+2)^2", sp.expand(determinants[1]-256*t**34*(t*t+2)**2)==0)
    check("N5 obstruction excludes h=0: attachment rank drops", F.subs(t,0).rank()<19)
    coloured = terms(5,{(0,1):2,(0,2):-2,(1,2):3,(1,3):-2,(2,3):2,(0,4):1},[(3,"Y",7)])
    check("N5 h=0 admits a connected attachment with carrier Y^5", not commutator(coloured,{"YYYYY":1}))


if __name__ == "__main__":
    Q = n4_gate()
    n5_attachment_gate(Q)
    n7_gate()
    print(f"ALL PASS: {CHECKS} exact checks")
