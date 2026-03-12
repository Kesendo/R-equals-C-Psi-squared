"""
DECAY RATE DERIVATION: Where do 2γ, 8γ/3, 10γ/3 come from?
============================================================
The three exact decay rates must have a structural origin.
For Z dephasing on n qubits: a Pauli string with k qubits
carrying X or Y decays at rate 2kγ. But our rates are NOT
simple multiples of 2γ. So the Liouvillian eigenmodes must
be superpositions. Can we derive the exact formula?
"""
import numpy as np
import math
import sys
sys.path.insert(0, r"D:\Entwicklung\Projekte Privat\R-equals-C-Psi-squared\simulations")
import star_topology_v3 as gpt

I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)

def op_at(op, qubit, n_q):
    ops = [I2]*n_q
    ops[qubit] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result

gamma = 0.05

print("=" * 70)
print("DECAY RATE DERIVATION")
print("Where do the three exact rates come from?")
print("=" * 70)

# For pure Z dephasing, a Pauli string P decays at rate 2*gamma*k
# where k = number of qubits with X or Y (anticommute with Z)
print("\n--- Simple Pauli string decay rates (no Hamiltonian mixing) ---")
print("k=0 (III,IIZ,IZI,ZII,IZZ,ZIZ,ZZI,ZZZ): decay = 0")
print(f"k=1 (XIx,YIx,IXx,IYx,...): decay = 2*gamma = {2*gamma:.4f}")
print(f"k=2 (XXx,XYx,YXx,YYx,...): decay = 4*gamma = {4*gamma:.4f}")
print(f"k=3 (XXX,XXY,...,YYY): decay = 6*gamma = {6*gamma:.4f}")

print(f"\nObserved rates: {2*gamma:.4f}, {8*gamma/3:.4f}, {10*gamma/3:.4f}")
print(f"Expected (simple): {2*gamma:.4f}, {4*gamma:.4f}, {6*gamma:.4f}")
print(f"\n8g/3 = {8*gamma/3:.4f} is between 2g and 4g")
print(f"10g/3 = {10*gamma/3:.4f} is between 2g and 4g")
print(f"So the Hamiltonian MIXES k=1 and k=2 modes into new eigenmodes")
print(f"with intermediate decay rates.")

# Let's check: do the decay rates depend on gamma in the expected way?
print(f"\n--- Scaling check: do rates scale linearly with gamma? ---")

def get_decay_rates(gamma_val, J_SA=1.0, J_SB=2.0):
    n_q = 3
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,n_q) @ op_at(p,1,n_q)
        H += J_SB * op_at(p,0,n_q) @ op_at(p,2,n_q)
    
    L_ops = [np.sqrt(gamma_val)*op_at(Z,i,n_q) for i in range(n_q)]
    d2 = d*d
    I_d = np.eye(d, dtype=complex)
    L_mat = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    for Lk in L_ops:
        Lk_dag = Lk.conj().T
        LdL = Lk_dag @ Lk
        L_mat += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    
    evals = np.linalg.eigvals(L_mat)
    decays = sorted(set(round(-np.real(ev), 6) for ev in evals if abs(np.imag(ev)) > 0.1))
    return decays

print(f"{'gamma':>6} | Decay rates | Ratios to gamma")
print("-" * 65)
for g in [0.01, 0.02, 0.05, 0.10, 0.20, 0.50]:
    rates = get_decay_rates(g)
    ratios = [f"{r/g:.4f}" for r in rates[:5]]
    print(f"{g:>6.2f} | {[f'{r:.4f}' for r in rates[:5]]} | {ratios}")

# What if we only dephase SOME qubits?
print(f"\n--- Which qubit's noise creates which rate? ---")

def get_rates_selective_noise(noise_qubits, gamma_val=0.05, J_SA=1.0, J_SB=2.0):
    n_q = 3
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    for p in [X, Y, Z]:
        H += J_SA * op_at(p,0,n_q) @ op_at(p,1,n_q)
        H += J_SB * op_at(p,0,n_q) @ op_at(p,2,n_q)
    
    L_ops = [np.sqrt(gamma_val)*op_at(Z,i,n_q) for i in noise_qubits]
    d2 = d*d
    I_d = np.eye(d, dtype=complex)
    L_mat = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    for Lk in L_ops:
        Lk_dag = Lk.conj().T
        LdL = Lk_dag @ Lk
        L_mat += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    
    evals = np.linalg.eigvals(L_mat)
    decays = sorted(set(round(-np.real(ev), 6) for ev in evals 
                       if abs(np.imag(ev)) > 0.1 and -np.real(ev) > 0.001))
    return decays

configs = [
    ("S only",     [0]),
    ("A only",     [1]),
    ("B only",     [2]),
    ("S+A",        [0,1]),
    ("S+B",        [0,2]),
    ("A+B",        [1,2]),
    ("S+A+B (all)",[0,1,2]),
]

print(f"{'Noise on':>12} | Decay rates (multiples of gamma)")
print("-" * 55)
for name, qubits in configs:
    rates = get_rates_selective_noise(qubits)
    ratios = [f"{r/gamma:.3f}" for r in rates]
    print(f"{name:>12} | {ratios}")

# What about 4-qubit systems? Do we get MORE decay rates?
print(f"\n--- 4-qubit decay rates ---")

def get_rates_4q(topology, couplings, gamma_val=0.05):
    n_q = 4
    d = 2**n_q
    H = np.zeros((d,d), dtype=complex)
    
    if topology == "star":
        # S(0) center, A(1), B(2), C(3)
        for i, J in enumerate(couplings):
            for p in [X, Y, Z]:
                H += J * op_at(p,0,n_q) @ op_at(p,i+1,n_q)
    elif topology == "chain":
        # Linear: 0-1-2-3
        for i, J in enumerate(couplings):
            for p in [X, Y, Z]:
                H += J * op_at(p,i,n_q) @ op_at(p,i+1,n_q)
    
    L_ops = [np.sqrt(gamma_val)*op_at(Z,i,n_q) for i in range(n_q)]
    d2 = d*d
    I_d = np.eye(d, dtype=complex)
    L_mat = -1j * (np.kron(H, I_d) - np.kron(I_d, H.T))
    for Lk in L_ops:
        Lk_dag = Lk.conj().T
        LdL = Lk_dag @ Lk
        L_mat += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, I_d) + np.kron(I_d, LdL.T))
    
    evals = np.linalg.eigvals(L_mat)
    decays = sorted(set(round(-np.real(ev), 5) for ev in evals 
                       if abs(np.imag(ev)) > 0.1 and -np.real(ev) > 0.001))
    return decays

print(f"3-qubit star [1,2]:     {[f'{r/gamma:.3f}g' for r in get_decay_rates(gamma)]}")
print(f"4-qubit star [1,2,1.5]: {[f'{r/gamma:.3f}g' for r in get_rates_4q('star', [1,2,1.5])]}")
print(f"4-qubit star [1,2,3]:   {[f'{r/gamma:.3f}g' for r in get_rates_4q('star', [1,2,3])]}")
print(f"4-qubit chain [1,1.5,2]: {[f'{r/gamma:.3f}g' for r in get_rates_4q('chain', [1,1.5,2])]}")

# The key test: do 4-qubit decay rates also stay fixed with topology change?
print(f"\n--- 4-qubit: do decay rates change with coupling? ---")
for J_list in [[1,1,1], [1,2,3], [1,2,0.5], [2,2,2], [0.5,1,3]]:
    rates = get_rates_4q('star', J_list)
    ratios = [f"{r/gamma:.3f}" for r in rates]
    print(f"  J={J_list}: {ratios}")

# PRACTICAL IMPLICATIONS
print(f"\n{'='*70}")
print("WHAT CAN WE DERIVE FROM THIS?")
print("="*70)

print("""
The system naturally separates into two INDEPENDENT information channels:

  FREQUENCY CHANNEL: Carries information about TOPOLOGY (who is connected)
    -> Changes with J values
    -> Does NOT change with noise
    -> Sensitive to hidden observers
    -> This is the "what is the network?" channel

  DECAY CHANNEL: Carries information about NOISE (what is the environment)
    -> Changes with gamma
    -> Does NOT change with topology
    -> Exact rational multiples of gamma
    -> This is the "what is the environment?" channel

CONSEQUENCE 1: You can characterize noise without knowing topology
  Measure any decay rate -> divide by the known coefficient -> get gamma
  Works even if you don't know J_SA or J_SB

CONSEQUENCE 2: You can characterize topology without knowing noise
  Measure frequencies -> they encode J values directly
  Works even if you don't know gamma

CONSEQUENCE 3: Hidden observer detection uses ONLY the frequency channel
  Because decay rates don't change when C joins, ONLY frequency shifts
  carry the detection signal. The decay channel is blind to topology.

CONSEQUENCE 4: The slow mode (2*gamma) is a PROTECTED CHANNEL
  It decays slowest. Information in c- survives longest.
  This is natural quantum error protection: the antisymmetric sector
  couples less strongly to the loss mechanism.
""")
