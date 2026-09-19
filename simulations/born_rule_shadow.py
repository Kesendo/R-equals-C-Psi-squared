"""
Finite Born-shadow decomposition (captured-result producer)
============================================================
This script evaluates a finite selected readout and decomposes its density
operator into named slow/fast spectral-coordinate pieces. Standard Born
probabilities are assumed, not derived. The diagonal probabilities decompose
linearly; the purity also contains a cross term. Photography and shutter
language is an interpretation, not a measurement mechanism.

The historical search used the purity target 0.27. Historical CΨ labels do not
name this observable. A trajectory that misses the target is reported as a
terminal no-hit sample rather than as a crossing.

The tracked result is a digest-locked historical capture and is not regenerated
by Task 6. This source is import-inert and writes only through main(output_path).

Output: simulations/results/born_rule_shadow.txt
"""

import argparse
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

DEFAULT_OUTPUT = Path(__file__).parent / "results" / "born_rule_shadow.txt"
J = 1.0
GAMMA = 0.05

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Ym = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)

def kron_chain(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r

def build_liouvillian(N, gammas):
    d = 2**N
    Id = np.eye(d, dtype=complex)
    H = np.zeros((d, d), dtype=complex)
    for i in range(N - 1):
        for P in [Xm, Ym, Zm]:
            ops = [I2] * N; ops[i] = P; ops[i + 1] = P
            H += J * kron_chain(ops)
    L = -1j * (np.kron(H, Id) - np.kron(Id, H.T))
    for k in range(N):
        ops = [I2] * N; ops[k] = Zm
        Lk = np.sqrt(gammas[k]) * kron_chain(ops)
        LdL = Lk.conj().T @ Lk
        L += np.kron(Lk, Lk.conj()) - 0.5 * (np.kron(LdL, Id) + np.kron(Id, LdL.T))
    return L

def purity_readout(rho):
    """The scalar actually evaluated by this historical calculation: Tr(rho^2)."""
    return np.real(np.trace(rho @ rho))

def propagate_to(rho0_vec, eigvals, R, L_left, d, t_target):
    """Propagate to a specific time."""
    coeffs = L_left.conj().T @ rho0_vec
    exp_l = np.exp(eigvals * t_target)
    rv = R @ (coeffs * exp_l)
    rho = rv.reshape(d, d)
    rho = (rho + rho.conj().T) / 2
    rho /= np.trace(rho).real
    return rho

def sample_purity_target_or_terminal(
        rho0_vec, eigvals, R, L_left, d, target_purity=0.27, t_max=50):
    """Return (reached, time, rho); false means the terminal sample, not a hit."""
    coeffs = L_left.conj().T @ rho0_vec
    for t in np.linspace(0.01, t_max, 10000):
        exp_l = np.exp(eigvals * t)
        rv = R @ (coeffs * exp_l)
        rho = rv.reshape(d, d)
        rho = (rho + rho.conj().T) / 2
        rho /= np.trace(rho).real
        p = purity_readout(rho)
        if p <= target_purity:
            return True, t, rho
    return False, t_max, rho


def decompose_spectral_coordinates(rho_past, rho_future):
    """Return the linear diagonal and quadratic purity decomposition."""
    I_past = np.real(np.diag(rho_past))
    I_future = np.real(np.diag(rho_future))
    I_cross = I_past + I_future
    purity_past = np.real(np.trace(rho_past @ rho_past))
    purity_future = np.real(np.trace(rho_future @ rho_future))
    purity_cross = 2 * np.real(np.trace(rho_past @ rho_future))
    purity_total = purity_past + purity_future + purity_cross
    return I_cross, purity_past, purity_future, purity_cross, purity_total


def diagonal_reconstruction_error(diagonal_sum, probabilities):
    """Compare the normalized decomposed diagonal with the live probabilities."""
    normalized = diagonal_sum / np.sum(diagonal_sum)
    return np.max(np.abs(normalized - probabilities))

def main(output_path=DEFAULT_OUTPUT):
    destination = Path(output_path)
    if sys.platform == "win32":
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    out = []
    def log(msg=""):
        print(msg)
        out.append(msg)

    log("=" * 75)
    log("FINITE BORN-SHADOW DECOMPOSITION")
    log("=" * 75)
    log()

    # ─────────────────────────────────────────────
    # N=2 Analysis
    # ─────────────────────────────────────────────

    N = 2
    d = 2**N
    d2 = d * d
    sigma_gamma = N * GAMMA
    gammas = [GAMMA] * N

    # Build Liouvillian and eigendecompose
    L_mat = build_liouvillian(N, gammas)
    eigvals, R = np.linalg.eig(L_mat)
    _, Lf = np.linalg.eig(L_mat.T)
    ov = Lf.conj().T @ R
    for j in range(d2):
        Lf[:, j] /= ov[j, j]

    # Classify eigenmodes as "past" (slow absorption, |Re| < Σγ) and "future" (fast, |Re| ≥ Σγ)
    past_mask = (-eigvals.real) < sigma_gamma  # Re closer to 0
    future_mask = (-eigvals.real) >= sigma_gamma  # Re closer to -2Σγ

    log(f"N={N}: {d2} eigenvalues, Σγ = {sigma_gamma}")
    log(f"  Past modes (|Re| < Σγ): {np.sum(past_mask)}")
    log(f"  Future modes (|Re| ≥ Σγ): {np.sum(future_mask)}")
    log()

    # ─── Test with multiple initial states ───

    states = {
        "Bell+": None,  # will construct below
        "|01>": None,
        "|++>": None,
    }

    # Bell+ = (|00> + |11>)/sqrt(2)
    psi_bell = np.zeros(d, dtype=complex)
    psi_bell[0] = 1/np.sqrt(2)  # |00>
    psi_bell[3] = 1/np.sqrt(2)  # |11>
    states["Bell+"] = np.outer(psi_bell, psi_bell.conj())

    # |01>
    psi_01 = np.zeros(d, dtype=complex)
    psi_01[1] = 1  # |01>
    states["|01>"] = np.outer(psi_01, psi_01.conj())

    # |++> = |+> tensor |+>
    plus = np.array([1, 1], dtype=complex) / np.sqrt(2)
    psi_pp = np.kron(plus, plus)
    states["|++>"] = np.outer(psi_pp, psi_pp.conj())

    basis_labels = ["|00>", "|01>", "|10>", "|11>"]

    for state_name, rho0 in states.items():
        rho0_vec = rho0.ravel()
        coeffs = Lf.conj().T @ rho0_vec

        log("=" * 75)
        log(f"INITIAL STATE: {state_name}")
        log("=" * 75)

        reached, t_readout, rho_readout = sample_purity_target_or_terminal(
            rho0_vec, eigvals, R, Lf, d)
        status = "purity target reached" if reached else "terminal no-hit"
        purity = purity_readout(rho_readout)
        log(f"  Readout status: {status}")
        log(f"  Sample time: t = {t_readout:.4f}")
        log(f"  Purity at sample (target 0.27): {purity:.6f}")
        log()

        # Standard Born probabilities at this selected sample.
        P = np.real(np.diag(rho_readout))
        log(f"  Standard Born probabilities P(i) at the selected sample:")
        for i, label in enumerate(basis_labels):
            log(f"    P({label}) = {P[i]:.6f}")
        log()

        # ─── Decompose into past and future contributions ───

        exp_l = np.exp(eigvals * t_readout)

        # Past contribution: modes with slow absorption
        rho_past_vec = R @ (coeffs * exp_l * past_mask)
        rho_past = rho_past_vec.reshape(d, d)
        rho_past = (rho_past + rho_past.conj().T) / 2

        # Future contribution: modes with fast absorption
        rho_future_vec = R @ (coeffs * exp_l * future_mask)
        rho_future = rho_future_vec.reshape(d, d)
        rho_future = (rho_future + rho_future.conj().T) / 2

        # Verify: rho_past + rho_future = rho_readout (up to normalization)
        rho_sum = rho_past + rho_future
        rho_sum /= np.trace(rho_sum).real
        recon_err = np.linalg.norm(rho_sum - rho_readout, 'fro')
        log(f"  Reconstruction: ||rho_past + rho_future - rho_readout|| = {recon_err:.2e}")

        # Diagonal elements and purity are computed by the same live helper
        # exercised by the independent source-copy gate.
        I_past = np.real(np.diag(rho_past))
        I_future = np.real(np.diag(rho_future))
        I_cross, purity_past, purity_future, purity_cross, purity_total = \
            decompose_spectral_coordinates(rho_past, rho_future)
        diag_recon_err = diagonal_reconstruction_error(I_cross, P)
        log(f"  Diagonal reconstruction error: {diag_recon_err:.2e}")

        # Interference term: P(i) - (I_past + I_future) after normalization
        # In the density matrix formalism, there is no explicit interference term
        # because rho = rho_past + rho_future (linear superposition of operators).
        # But we can look at how the cross-terms in rho_past * rho_future contribute
        # to this calculation's selected purity readout.

        # For the Born rule P(i) = <i|rho|i>, the decomposition is exactly linear:
        # P(i) = <i|rho_past|i> + <i|rho_future|i>
        # No interference term in the probabilities themselves!

        # However, there IS interference in the COHERENCES (off-diagonal elements):
        # rho[i,j] = rho_past[i,j] + rho_future[i,j]
        # The purity Tr(rho^2) contains cross-terms Tr(rho_past * rho_future)

        log(f"\n  Purity decomposition at the selected sample:")
        log(f"    Tr(rho_past²)   = {purity_past:.6f}  ({purity_past/purity_total*100:.1f}%)")
        log(f"    Tr(rho_future²) = {purity_future:.6f}  ({purity_future/purity_total*100:.1f}%)")
        log(f"    2Tr(rho_p·rho_f)= {purity_cross:.6f}  ({purity_cross/purity_total*100:.1f}%) <- interference")
        log(f"    Total            = {purity_total:.6f}")
        log()

        # P(i) decomposition (diagonal: no interference possible)
        log(f"  P(i) decomposition (past + future, NO interference in diagonal):")
        log(f"  {'basis':>6s} {'P(i)':>10s} {'past':>10s} {'future':>10s} {'past%':>8s}")
        for i, label in enumerate(basis_labels):
            frac = I_past[i] / P[i] * 100 if P[i] > 1e-10 else 0
            log(f"  {label:>6s} {P[i]:10.6f} {I_past[i]:10.6f} {I_future[i]:10.6f} {frac:7.1f}%")
        log()

        # C_i analysis: effective "mirror quality" per basis state
        # C_i = P(i) / (1/d) = d * P(i) (how much more probable than uniform)
        log(f"  Effective mirror quality C_i = d × P(i):")
        for i, label in enumerate(basis_labels):
            C_i = d * P[i]
            log(f"    C({label}) = {C_i:.4f}")
        log()

        # Coherence analysis: where IS the interference?
        off_diag_past = np.sum(np.abs(rho_past)**2) - np.sum(np.abs(np.diag(rho_past))**2)
        off_diag_future = np.sum(np.abs(rho_future)**2) - np.sum(np.abs(np.diag(rho_future))**2)
        off_diag_total = np.sum(np.abs(rho_readout)**2) - np.sum(np.abs(np.diag(rho_readout))**2)

        log(f"  Coherence (off-diagonal) contribution:")
        log(f"    Past coherence:   {off_diag_past:.6f}")
        log(f"    Future coherence: {off_diag_future:.6f}")
        log(f"    Total coherence:  {off_diag_total:.6f}")
        log()

    # ─────────────────────────────────────────────
    # Summary
    # ─────────────────────────────────────────────

    log("=" * 75)
    log("SUMMARY")
    log("=" * 75)
    log()
    log("1. The diagonal probabilities decompose linearly into the named slow")
    log("   and fast spectral-coordinate pieces. Standard Born probabilities are")
    log("   assumed, not derived; linearity is the finite algebraic statement:")
    log("   <i|rho_past + rho_future|i> = <i|rho_past|i> + <i|rho_future|i>.")
    log()
    log("2. The purity Tr(rho²) contains the cross-term")
    log("   2Tr(rho_past · rho_future). This changes the finite selected readout")
    log("   time but supplies no universal measurement or irreversibility law.")
    log()
    log("3. Photography and shutter language is an interpretation.")
    log("   It is retained only as an invitation for reading the finite decomposition,")
    log("   not as a derived Born rule or a physical measurement mechanism.")
    log()
    log("4. The search target is purity 0.27, not CΨ = 1/4.")
    log("   Historical CΨ labels do not name this observable. Bell+ and |01>")
    log("   terminate without a hit at t=50; only |++> reaches the purity target.")
    log()
    log("5. P(i) IS STATE-DEPENDENT. Different initial states produce")
    log("   different probability distributions at the selected readout.")
    log("   This is a finite state-dependent comparison, not a universal law.")

    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    print(f"\n>>> Results saved to: {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()
    main(args.output)
