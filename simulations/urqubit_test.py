"""
Urqubit-Test: Erzeugt die partielle Spur über einen Paritätssektor
den Z-Dephasing-Dissipator?

Die algebraische Frage:
- Der Pauli-Raum eines Qubits ist 4-dimensional: {I, X, Y, Z}
- Π² = X (bit-flip), teilt in +1 Parität {I,X} und -1 Parität {Y,Z}  
- Z-Dephasing teilt in immune {I,Z} und decaying {X,Y}
- DIESE SPLITS SIND VERSCHIEDEN. Kreuzstruktur.

Tensorprodukt-Zerlegung: Pauli-Index = (a, b)
  σ_{a,b} = i^{ab} X^a Z^b  (Weyl-Form)
  (0,0) = I    (0,1) = Z
  (1,0) = X    (1,1) = iY → Y (bis auf Phase)

  Bit a: immune (0) vs decaying (1) unter Z-dephasing
  Bit b: Π²-even (0) vs Π²-odd (1)

Frage: Was passiert wenn man über Bit b spurt?
"""

import numpy as np
from itertools import product

# Pauli matrices
I = np.array([[1, 0], [0, 1]], dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

paulis = {'I': I, 'X': X, 'Y': Y, 'Z': Z}
pauli_list = [I, X, Y, Z]
pauli_names = ['I', 'X', 'Y', 'Z']

print("=" * 60)
print("TEIL 1: Die Kreuzstruktur des Pauli-Raums")
print("=" * 60)

# Die zwei Splits
print("\nSplit 1 (Z-Dephasing): immune {I,Z} vs decaying {X,Y}")
print("Split 2 (Π²-Parität): even {I,X} vs odd {Y,Z}")
print()
print("Kreuzstruktur als 2×2 Tabelle:")
print()
print("              Π²-even (b=0)    Π²-odd (b=1)")
print("immune (a=0):     I                Z")
print("decaying(a=1):    X                Y")
print()

# Verify: Z-dephasing eigenvalues
print("Z-Dephasing-Eigenvalues (σ_z · σ · σ_z - σ):")
for name, sigma in paulis.items():
    result = Z @ sigma @ Z - sigma
    eigenval = np.trace(result @ sigma.conj().T) / np.trace(sigma @ sigma.conj().T)
    print(f"  D({name}) = {eigenval.real:.1f} · {name}")

print()
print("=" * 60)
print("TEIL 2: Die Tensorprodukt-Struktur")
print("=" * 60)

# Index as (a, b): a = dephasing-bit, b = parity-bit
# σ_{00} = I, σ_{10} = X, σ_{01} = Z, σ_{11} = Y (up to phase)
print()
print("Pauli-Operatoren als (a,b) mit a=dephasing, b=parität:")
mapping = {(0,0): 'I', (1,0): 'X', (0,1): 'Z', (1,1): 'Y'}
for (a,b), name in sorted(mapping.items()):
    deph = "immune" if a == 0 else "decaying"
    par = "even" if b == 0 else "odd"
    print(f"  ({a},{b}) = {name:2s}  [{deph:8s}, {par:4s}]")

print()
print("Der Dissipator in dieser Struktur:")
print("  D = -2γ · |a=1⟩⟨a=1| ⊗ I_b")
print("  → wirkt NUR auf Faktor a (dephasing-bit)")
print("  → Faktor b (parität) ist unsichtbar für den Dissipator")

print()
print("=" * 60)
print("TEIL 3: Partielle Spur über Faktor b (Parität)")
print("=" * 60)
print()
print("Wenn wir über den Paritäts-Faktor (b) spuren:")
print("  I und Z → beide a=0 (immune) → kollabieren zu 'immune'")
print("  X und Y → beide a=1 (decaying) → kollabieren zu 'decaying'")
print()
print("Ergebnis: ein 2-dimensionaler Raum {immune, decaying}")
print("mit Dynamik: immune bleibt, decaying zerfällt mit Rate 2γ")
print()
print("Das IST der Dissipator, projiziert auf den Dephasing-Faktor.")
print("Aber es ERZEUGT ihn nicht. γ muss von außen gegeben werden.")

print()
print("=" * 60)
print("TEIL 4: Die andere Richtung - Spur über Faktor a (Dephasing)")  
print("=" * 60)
print()
print("Wenn wir über den Dephasing-Faktor (a) spuren:")
print("  I und X → beide b=0 (Π²-even) → 'even'")
print("  Z und Y → beide b=1 (Π²-odd) → 'odd'")
print()
print("Ergebnis: ein 2-dimensionaler Raum {even, odd}")
print("mit Dynamik: even enthält Σ(immune+decaying), odd ebenso")
print("Der Dissipator wird herausgemittelt (I und X mitteln zu γ)")

print()
print("=" * 60)
print("TEIL 5: Der tiefere Test - N=2, gekoppeltes System")
print("=" * 60)
print()

def build_liouvillian_n2(J, gamma):
    """Build 16x16 Liouvillian for 2-qubit Heisenberg + Z-dephasing."""
    dim = 4  # 2^N
    dim2 = dim * dim  # 4^N = 16
    
    # Pauli operators for 2 sites
    def kron(A, B):
        return np.kron(A, B)
    
    # Hamiltonian: H = J(XX + YY + ZZ)
    H = J * (kron(X,X) + kron(Y,Y) + kron(Z,Z))
    
    # Liouvillian: L(rho) = -i[H,rho] + D(rho)
    # In superoperator form
    L = np.zeros((dim2, dim2), dtype=complex)
    
    # Hamiltonian part: -i(H⊗I - I⊗H^T)
    L += -1j * (np.kron(H, np.eye(dim)) - np.kron(np.eye(dim), H.T))
    
    # Dissipator: γ_k (σ_z^k ρ σ_z^k - ρ) for each site
    for site in range(2):
        if site == 0:
            Lk = kron(Z, I)
        else:
            Lk = kron(I, Z)
        
        # D(rho) = gamma * (Lk rho Lk† - rho)  [since Lk=Lk†=σ_z]
        D = gamma * (np.kron(Lk, Lk.conj()) - np.kron(np.eye(dim), np.eye(dim)))
        L += D
    
    return L

# Build at gamma=0 and gamma>0
J = 1.0
gamma = 0.05
L0 = build_liouvillian_n2(J, 0.0)
Lg = build_liouvillian_n2(J, gamma)

evals_0 = np.linalg.eigvals(L0)
evals_g = np.linalg.eigvals(Lg)

# Sort by real part
evals_0 = np.sort(evals_0)[::-1]
evals_g = np.sort(evals_g)[::-1]

print("Eigenvalues bei Σγ = 0 (der Spiegel, Null):")
print("  λ ↔ -λ Paarung (Π·L·Π⁻¹ = -L):")
evals_0_sorted = sorted(evals_0, key=lambda x: (round(x.real, 10), round(x.imag, 10)))
for ev in evals_0_sorted:
    partner = -ev
    # Find closest match
    dists = [abs(ev2 - partner) for ev2 in evals_0_sorted]
    min_dist = min(dists)
    print(f"  λ = {ev.real:+8.4f} {ev.imag:+8.4f}i  →  -λ Fehler: {min_dist:.2e}")

print()
print(f"\nEigenvalues bei Σγ = {2*gamma} (γ={gamma} pro Site):")
print(f"  λ ↔ -(λ+2Σγ) Paarung, Σγ = {2*gamma}:")
sum_gamma = 2 * gamma
for ev in sorted(evals_g, key=lambda x: x.real):
    partner = -(ev + 2*sum_gamma)
    dists = [abs(ev2 - partner) for ev2 in evals_g]
    min_dist = min(dists)
    if min_dist < 1e-10:
        sym = "✓"
    else:
        sym = f"✗ ({min_dist:.2e})"
    print(f"  λ = {ev.real:+8.4f} {ev.imag:+8.4f}i  Palindrom: {sym}")

print()
print("=" * 60)
print("TEIL 6: Die eigentliche Frage - Selbstkonsistenz")
print("=" * 60)
print()
print("Tom's Hypothese: Der Urqubit sitzt auf beiden Seiten.")
print("Null ist die Spiegelgrenze. Ihr Rauschen ist unser Rauschen.")
print()
print("Algebraische Übersetzung:")
print("  Bei Σγ = 0: Π·L·Π⁻¹ = -L")
print("  Die Eigenvalues sind rein imaginär: ±iω")
print("  Keine Seite zerfällt. Beide sind eins.")
print()
print("  Bei Σγ > 0: Π·L·Π⁻¹ = -L - 2Σγ·I")
print("  Die Eigenvalues haben Realteile: Zerfall.")
print("  EINE Seite 'lebt', die andere ist 'Rauschen'.")
print()
print("Die Kreuzstruktur zeigt:")
print("  Faktor a (dephasing) und Faktor b (parität) sind UNABHÄNGIG.")
print("  Der Dissipator wirkt nur auf a.")
print("  Die Parität lebt in b.")
print("  Sie sehen sich nicht direkt.")
print()
print("ABER: Π selbst mischt a und b!")
print("  Π: (a,b) → (1-a, ...) mit Phasen")
print("  Π tauscht immune↔decaying UND dreht die Parität.")
print("  Π IST die Brücke zwischen den beiden Faktoren.")

# Verify Pi mixes the factors
print()
print("Π-Wirkung auf (a,b):")
# Pi per site: I→X (+1), X→I (+1), Y→iZ (+i), Z→iY (+i)
# In (a,b): (0,0)→(1,0), (1,0)→(0,0), (1,1)→(0,1)·i, (0,1)→(1,1)·i
pi_map = {
    (0,0): ((1,0), 1),    # I → X
    (1,0): ((0,0), 1),    # X → I
    (0,1): ((1,1), 1j),   # Z → iY
    (1,1): ((0,1), 1j),   # Y → iZ
}
for (a,b), ((a2,b2), phase) in sorted(pi_map.items()):
    name_in = mapping[(a,b)]
    name_out = mapping[(a2,b2)]
    print(f"  Π: ({a},{b})={name_in} → ({a2},{b2})={name_out}  (×{phase})")
    da = a2 - a
    db = b2 - b
    print(f"       Δa={da:+d}, Δb={db:+d}  → Π mischt {'BEIDE Faktoren' if da != 0 and db != 0 else 'nur a' if da != 0 else 'nur b' if db != 0 else 'keinen'}")

print()
print("=" * 60)
print("TEIL 7: Ergebnis")
print("=" * 60)
print()
print("Die partielle Spur über einen Paritätssektor ERZEUGT NICHT")
print("den Z-Dephasing-Dissipator. Die Sektoren stehen senkrecht")
print("aufeinander: Parität ({I,X} vs {Y,Z}) ist ORTHOGONAL zu")
print("Dephasing ({I,Z} vs {X,Y}).")
print()
print("ABER die Kreuzstruktur offenbart etwas Tieferes:")
print()
print("Der Pauli-Raum IST ein Tensorprodukt C² ⊗ C²:")
print("  Faktor a = Dephasing-Empfindlichkeit (immune/decaying)")
print("  Faktor b = Π²-Parität (even/odd)")
print()
print("Und Π VERSCHRÄNKT diese beiden Faktoren.")
print("Π ist kein Produkt-Operator auf C² ⊗ C².")
print("Π erzeugt Verschränkung zwischen a und b.")
print()
print("Das heißt: die Dephasing-Struktur und die Paritäts-Struktur")
print("sind nicht unabhängig. Sie werden durch Π miteinander verwoben.")
print("Die eine KANN ohne die andere nicht existieren.")
print()

# Final check: is Pi entangling?
# Pi maps product states to... let's check
print("Ist Π verschränkend auf C²⊗C²?")
print("  |a=0,b=0⟩ = I → Π → |a=1,b=0⟩ = X : Produktzustand → Produktzustand")
print("  |a=0,b=1⟩ = Z → Π → i·|a=1,b=1⟩ = iY : Produktzustand → Produktzustand")
print()
print("  Π ist NICHT verschränkend in dieser Basis.")
print("  Π ist eine kontrollierte Operation: flip a, phase auf b.")
print()
print("  ABER: Π² = X^N ist NICHT trivial (es negiert b=1).")
print("  Die Gruppe ⟨Π⟩ hat Ordnung 4 (pro Site).")
print("  Π ist eine Z₄-Symmetrie, nicht Z₂.")
print()
print("=" * 60)
print("FAZIT FÜR HOMEWORK")  
print("=" * 60)
print()
print("Die direkte Rechnung zeigt: partielle Spur → Dissipator")
print("funktioniert NICHT in der naiven Form.")
print()
print("Was sie ZEIGT:")
print("1. Der Pauli-Raum hat eine echte C²⊗C² Tensorstruktur")
print("2. Dephasing lebt in Faktor a, Parität in Faktor b")
print("3. Π mischt die Faktoren (a→1-a, b rotiert)")
print("4. γ kann nicht aus der Struktur allein abgeleitet werden")
print()
print("Die OFFENE Frage (für Homework):")
print("Gibt es eine selbstkonsistente Konstruktion auf dem")
print("VERDOPPELTEN Raum (C²⊗C²)_uns ⊗ (C²⊗C²)_spiegel,")
print("bei der die partielle Spur über 'spiegel' den Dissipator")
print("auf 'uns' erzeugt, UND UMGEKEHRT?")
print()
print("D.h.: nicht EIN Qubit das sein eigenes Rauschen erzeugt,")
print("sondern der Urqubit ALS VERDOPPLUNG, dessen zwei Hälften")
print("das Rauschen der jeweils anderen sind?")
