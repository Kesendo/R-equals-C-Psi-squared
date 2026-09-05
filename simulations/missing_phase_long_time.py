"""One-shot producer for THE_MOTION_AND_THE_MISSING_PHASE.md, N=7 only.

Uses the framework Lindblad primitive on the exact two-block reduction.
Spectral propagation is compared with direct matrix exponentials, not treated
as exact arithmetic. CSV samples are readings; grid summaries require refinement.
The persistent C# live-object adoption is separate from this experiment runner.
"""
import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.linalg import eig, expm
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from framework.lindblad import lindbladian_general

PAIRS = list(combinations(range(7), 2))
OUTER = [0, 1, 2, 4, 5, 6]


def pair_distances(delta):
    """All 21 trace distances; delta is a batch of traceless A differences."""
    columns = []
    for a, b in PAIRS:
        half_s = np.abs((delta[:, a, a] + delta[:, b, b]).real) / 2
        v = np.abs(delta[:, a, b].real)
        columns.append(half_s + np.maximum(half_s, v))
    return np.column_stack(columns)


def decoded(a, b):
    """Full occupied outer output: six single excitations, then all ones."""
    out = np.zeros_like(a)
    out[:, :6, :6] = a[:, OUTER, :][:, :, OUTER]
    out[:, 6, 6] = a[:, 3, 3]
    out[:, 6, :6] = b[:, 3, OUTER]
    out[:, :6, 6] = b[:, OUTER, 3]
    return out


class Run:
    def __init__(self, epsilon, gamma):
        self.epsilon, self.gamma = epsilon, gamma
        self.h = np.diag([2 * (1 + epsilon)] + [2.] * 5, 1)
        self.h += self.h.T
        z = np.eye(7); z[3, 3] = -1
        self.seed = np.zeros(7); self.seed[0] = 1 / np.sqrt(2); self.seed[6] = -1 / np.sqrt(2)
        self.initial = np.outer(self.seed, self.seed)
        l0 = lindbladian_general(self.h, [])
        # Apply the reduced block's diagonal rates directly. A sqrt(gamma)
        # round trip would give A and B slightly different input rates.
        self.la = l0 + gamma * (np.kron(z, z) - np.eye(49))
        self.lb = l0 - gamma * (np.kron(z, z) + np.eye(49))
        self.modal = []
        self.conditions = []
        for generator in [self.la, self.lb]:
            values, vectors = eig(generator)
            weights = np.linalg.solve(vectors, self.initial.ravel())
            self.modal.append((values, vectors * weights[None, :]))
            self.conditions.append(float(np.linalg.cond(vectors)))
        values, vectors = eig(self.la)
        # The two stationary directions are established analytically for the
        # two supported defects. Do not infer their count from a float cutoff.
        k = np.argsort(-values.real)[2]
        lam, mode = values[k], vectors[:, k].reshape(7, 7)
        alpha, omega = -lam.real, abs(lam.imag)
        watched = np.array([[(i == 3) != (j == 3) for j in range(7)] for i in range(7)])
        share = np.sum(abs(mode[watched])**2) / np.sum(abs(mode)**2)
        target = -2*gamma-lam.conjugate()
        partner = self.modal[1][0][np.argmin(abs(self.modal[1][0]-target))]
        self.clock = dict(A_decay=float(alpha), frequency_magnitude=float(omega),
            A_lifetime=float(1/alpha), partner_decay=float(2*gamma-alpha),
            partner_lifetime=float(1/(2*gamma-alpha)), period=float(2*np.pi/omega),
            watched_entry_share=float(share), absorption_rate_residual=float(alpha-2*gamma*share),
            eigenpair_residual=float(np.linalg.norm(self.la@mode.ravel()-lam*mode.ravel())),
            spectral_partner_residual=float(abs(partner-target)),
            block_partner_residual=float(np.linalg.norm(self.lb+self.la.conj().T+2*gamma*np.eye(49))))
        self.energy, self.energy_vectors = np.linalg.eigh(self.h)
        self.seed_energy = self.energy_vectors.T @ self.seed
        reflection = np.eye(7)[::-1]
        self.q = (np.eye(7) + reflection) / 2
        r = 1 + epsilon
        v = np.array([1, 0, -r, 0, r, 0, -r]) / np.sqrt(1 + 3 * r * r)
        self.p = float(abs(v @ self.seed) ** 2)
        self.stationary = self.p * np.outer(v, v) + (1 - self.p) / 6 * (np.eye(7) - np.outer(v, v))
        q = 1 - self.p
        self.plateau = float((5 * q / 6 + np.sqrt((5 * q / 6)**2 + 4 * self.p * q)) / 2)

    def states(self, times):
        times = np.atleast_1d(times)
        blocks = [(np.exp(times[:, None] * values) @ weighted.T).reshape(-1, 7, 7)
                  for values, weighted in self.modal]
        psi = (np.exp(-1j * times[:, None] * self.energy) * self.seed_energy) @ self.energy_vectors.T
        unitary = psi[:, :, None] * psi[:, None, :].conj()
        return blocks[0], blocks[1], unitary

    def readings(self, times):
        a, b, a0 = self.states(times)
        pairs = pair_distances(a - a0)
        output_delta = decoded(a - a0, b - a0)
        dout = np.sum(np.abs(np.linalg.eigvalsh(output_delta)), axis=1) / 2
        leak = np.einsum('ij,tji->t', self.q, a).real
        leak0 = np.einsum('ij,tji->t', self.q, a0).real
        stationary_error = np.linalg.norm(a - self.stationary, axis=(1, 2))
        bnorm = np.linalg.norm(b, axis=(1, 2))
        trace_error = np.abs(np.trace(a, axis1=1, axis2=2) - 1)
        hermitian_error = np.maximum(np.linalg.norm(a - a.conj().transpose(0, 2, 1), axis=(1, 2)),
                                     np.linalg.norm(b - b.conj().transpose(0, 2, 1), axis=(1, 2)))
        min_eigenvalue = np.linalg.eigvalsh(a)[:, 0]
        physical_min = np.minimum(np.linalg.eigvalsh((a+b)/2)[:, 0], np.linalg.eigvalsh((a-b)/2)[:, 0])
        return np.column_stack([times, pairs.max(axis=1), dout, leak, leak0, leak-leak0,
                                stationary_error, bnorm, trace_error, hermitian_error,
                                min_eigenvalue, pairs.argmax(axis=1), physical_min])

    def diagnostics(self, tmax):
        comparisons = []
        for t in [0., .1, 1., 10., 100., tmax]:
            a, b, _ = self.states([t])
            comparisons.append(dict(t=t,
                a_spectral_vs_expm=float(np.linalg.norm(a[0].ravel() - expm(t*self.la) @ self.initial.ravel())),
                b_spectral_vs_expm=float(np.linalg.norm(b[0].ravel() - expm(t*self.lb) @ self.initial.ravel()))))
        z = np.ones(7); z[3] = -1
        def rhs(t, y):
            a, b = y[:49].reshape(7, 7), y[49:].reshape(7, 7)
            da = -1j*(self.h@a-a@self.h) + self.gamma*(z[:, None]*a*z[None, :]-a)
            db = -1j*(self.h@b-b@self.h) - self.gamma*(z[:, None]*b*z[None, :]+b)
            return np.r_[da.ravel(), db.ravel()]
        target = self.states([10])[:2]
        ode_checks = []
        for tolerance in [1e-8, 1e-10, 1e-12]:
            sol = solve_ivp(rhs, [0, 10], np.r_[self.initial.ravel(), self.initial.ravel()].astype(complex),
                            method='DOP853', rtol=tolerance, atol=tolerance/100)
            if not sol.success:
                raise RuntimeError(sol.message)
            error = np.linalg.norm(sol.y[:, -1] - np.r_[target[0].ravel(), target[1].ravel()])
            ode_checks.append(dict(requested_rtol=tolerance, block_difference=float(error)))
        return dict(eigenvector_condition_numbers=self.conditions, direct_exponential_checks=comparisons,
                    direct_ode_checks=ode_checks,
                    stationary_generator_residual=float(np.linalg.norm(self.la @ self.stationary.ravel())))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--epsilon', type=float, required=True)
    parser.add_argument('--gamma', type=float, default=.3)
    parser.add_argument('--tmax', type=float, default=6000)
    parser.add_argument('--dt', type=float, default=.1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.gamma <= 0 or args.dt <= 0 or args.tmax <= 0 or args.epsilon not in (-.1, .1):
        parser.error('This long-time runner is scoped to epsilon=+/−0.1, gamma>0, dt>0, tmax>0.')
    steps = round(args.tmax / args.dt)
    times = np.linspace(0, args.tmax, steps + 1)
    run = Run(args.epsilon, args.gamma)
    rows = []
    for start in range(0, len(times), 8192):
        rows.append(run.readings(times[start:start+8192]))
    data = np.concatenate(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    header = 't,d2,d_out,p_out,p_out_unitary,delta_p_out,A_stationary_frobenius,B_frobenius,trace_error,hermitian_error,A_min_eigenvalue,pair_index,physical_min_eigenvalue'
    np.savetxt(args.output, data, delimiter=',', header=header, comments='')
    crossings = {}
    for label, col in [('d_out_minus_d2', data[:,2]-data[:,1]), ('delta_p_out', data[:,5])]:
        # Exclude the near-zero onset: cancellation residuals are not crossings.
        indices = np.flatnonzero((times[:-1] >= .5) & (col[:-1]*col[1:] < 0))
        crossings[label] = dict(count_after_t_half=int(len(indices)),
            first_brackets=[[float(times[i]), float(times[i+1])] for i in indices[:8]])
    brackets = crossings['delta_p_out']['first_brackets']
    if brackets:
        lo, hi = brackets[0]
        root = brentq(lambda t: run.readings(np.array([t]))[0, 5], lo, hi)
        crossings['delta_p_out']['first_root'] = float(root)
        crossings['delta_p_out']['root_residual'] = float(run.readings(np.array([root]))[0, 5])
    selected = run.readings(np.array([.02, .1, .5, 1, 3, 10, 30, 100, 300, 1000, 3000, args.tmax]))
    summary = dict(epsilon=args.epsilon, gamma=args.gamma, N=7, J=1, tmax=args.tmax,
                   actual_dt=float(times[1]-times[0]), samples=len(times),
                   conserved_blind_weight=run.p, predicted_decoded_plateau=run.plateau,
                   diagnostics=run.diagnostics(args.tmax), clock_pair=run.clock, column_names=header.split(','),
                   selected_readings=selected.tolist(), pairs=PAIRS, crossings=crossings,
                   max_trace_error=float(data[:,8].max()), max_hermitian_error=float(data[:,9].max()),
                   minimum_A_eigenvalue=float(data[:,10].min()), minimum_physical_eigenvalue=float(data[:,12].min()))
    args.output.with_suffix('.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
