"""Virtual ternary readout at the N4 equal-end EP. No hardware model.

BLAS/OpenMP threads are pinned to one. Time is gamma*t; gamma=1, q=2,
XY hopping=2*q*weight. Masks use site zero as the least significant bit.
"""
import json
import os
from pathlib import Path

os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import expm_multiply
from route_b_other_n_unfolding import parts
from route_b_artifact_provenance import (
    canonical_text_sha256 as _canonical_text_sha256,
    artifact_provenance as _artifact_provenance,
    verify_artifact_provenance as _verify_artifact_provenance,
)

ROOT = Path(__file__).resolve().parents[1]
EP = np.sqrt(2)-2
TIMES = np.linspace(0, 2, 401)
SOURCE = ROOT/'simulations/results/route_b_n4_self_fold.json'


def generator(epsilon, q=2):
    # Independent full Hilbert-space construction, row-major density vectors.
    h = np.zeros((16, 16))
    for site, weight in enumerate([1+epsilon/2, 1, 1+epsilon/2]):
        for mask in range(16):
            if (mask >> site & 1) != (mask >> (site+1) & 1):
                h[mask ^ (3 << site), mask] = 2*q*weight
    decay = np.array([-2*(a ^ b).bit_count()
                      for a in range(16) for b in range(16)])
    return -1j*(np.kron(h, np.eye(16))-np.kron(np.eye(16), h.T))+np.diag(decay)


def verify_generator_contract(builder=generator):
    """Check the q=2 end-bias convention against the independent (1,2) pencil.

    parts() builds the coherence block directly, without this Hilbert-space
    generator. At q=2 its entries use exact binary scalings of the integer
    pencil, so equality is entrywise, with no eigensolver tolerance. The two
    displaced inputs expose a rescaling about EP that the EP alone cannot see.
    """
    d, c, left, right, _ = parts(4)
    indices = [16*a+b for a in (1, 2, 4, 8) for b in (3, 5, 6, 9, 10, 12)]
    for shift in (0., -.05, .05):
        epsilon = EP+shift
        expected = d+2*(c+epsilon*(left+right)/2)
        actual = builder(epsilon, q=2)[np.ix_(indices, indices)]
        if not np.array_equal(actual, expected):
            residual = float(np.max(abs(actual-expected)))
            raise ValueError(f'generator contract failed at epsilon={epsilon}: '
                             f'max entry residual={residual}')


def artifact_provenance(script):
    """Canonical source and transitive-import manifest for the N4 artifact."""
    return _artifact_provenance(script, SOURCE)


def verify_artifact_provenance(artifact, script):
    """Reject absent/stale source, script or dependency hashes before writing."""
    _verify_artifact_provenance(artifact, script, SOURCE)


def run(epsilon, q=2):
    psi = np.zeros(16); psi[[1, 3]] = 1/np.sqrt(2)
    initial = np.outer(psi, psi).ravel().astype(complex)
    g = generator(epsilon, q)
    states = expm_multiply(g, initial, start=0, stop=2, num=len(TIMES)).reshape(-1, 16, 16)
    signal = 2*states[:, 1, 3].real
    acceptance = (states[:, 1, 1]+states[:, 3, 3]).real
    probabilities = np.array([(acceptance+signal)/2, (acceptance-signal)/2, 1-acceptance]).T
    assert np.min(probabilities) > -2e-12
    assert np.max(abs(np.trace(states, axis1=1, axis2=2)-1)) < 2e-12
    assert np.max(abs(states-states.conj().transpose(0, 2, 1))) < 2e-12
    assert min(la.eigvalsh(r)[0] for r in states[::20]) > -2e-12
    assert la.norm(la.expm(g*.375)@initial-states[75].ravel()) < 2e-12
    return states, signal, acceptance, probabilities


def main():
    verify_generator_contract(builder=generator)
    d, c, left, right, _ = parts(4)
    block = d+2*(c+EP*(left+right)/2)
    pairs = [(a, b) for a in range(16) if a.bit_count() == 1
             for b in range(16) if b.bit_count() == 2]
    indices = [16*a+b for a, b in pairs]
    states, signal, acceptance, probs = run(EP)
    initial = states[0].ravel()[indices]
    block_states = expm_multiply(block, initial, start=0, stop=2, num=len(TIMES))
    assert np.max(abs(block_states-states.reshape(-1, 256)[:, indices])) < 2e-12
    # Riesz projector avoids diagonalizing the defective pair. Two quadratures
    # independently check contour convergence and the nilpotent component.
    lam = -4+2j
    def projector(count):
        result = np.zeros_like(block)
        for theta in 2*np.pi*(np.arange(count)+.5)/count:
            z = .2*np.exp(1j*theta)
            result += z*la.solve((lam+z)*np.eye(24)-block, np.eye(24))/count
        return result
    p = projector(128)
    assert la.norm(p-projector(256)) < 1e-10
    nil = (block-lam*np.eye(24))@p
    assert abs(np.trace(p)-2) < 1e-10
    assert la.norm(p@p-p) < 1e-10
    assert la.norm(nil@nil) < 1e-10
    assert la.norm(nil) > .1
    read = pairs.index((1, 3))
    a, b = (p@initial)[read], (nil@initial)[read]
    for time in [.1, .375, 1.0]:
        assert la.norm(la.expm(block*time)@p-np.exp(lam*time)*(p+time*nil)) < 2e-11
    ep_component = 2*np.real(np.exp(lam*TIMES)*(a+TIMES*b))
    # Hamiltonian-off control knows the actual one-bit coherence decay.
    _, zero_signal, _, _ = run(EP, q=0)
    assert np.max(abs(zero_signal-np.exp(-2*TIMES))) < 2e-12
    assert np.max(abs(zero_signal-np.exp(-4*TIMES))) > .2  # wrong-rate control
    rows = []
    curves = []
    for shift in [-.05, .05]:
        _, other, acc, _ = run(EP+shift)
        difference = other-signal
        # Independent datasets: M shots EACH, unconditioned ternary variance.
        variance_sum = acceptance-signal**2+acc-other**2
        snr_per_sqrt_m = abs(difference)/np.sqrt(np.maximum(variance_sum, 1e-30))
        snr_per_sqrt_m[0] = 0
        best = int(np.argmax(snr_per_sqrt_m))
        rows.append(dict(epsilon_shift=shift, best_grid_time=float(TIMES[best]),
                         signal_ep=float(signal[best]), signal_other=float(other[best]),
                         difference=float(difference[best]),
                         acceptance_ep=float(acceptance[best]),
                         shots_per_setting_for_mean_snr_5=float(25/snr_per_sqrt_m[best]**2)))
        curves.append((shift, other))
    output = dict(model=dict(n=4, q=2, gamma=1, epsilon_ep=float(EP),
                             initial_masks=[1, 3], time_grid=[0, 2, 401]),
                  jordan_readout=dict(a=[a.real, a.imag], b=[b.real, b.imag],
                                      nilpotent_norm=float(la.norm(nil))),
                  comparisons=rows,
                  caveat='Grid optimum, ideal independent ternary shots; mean SNR is not a test power or EP certificate.',
                  **artifact_provenance(__file__))
    path = ROOT/'simulations/results/route_b_n4_virtual_readout.json'
    verify_artifact_provenance(output, __file__)
    path.write_text(json.dumps(output, indent=2)+'\n', encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    axes[0].plot(TIMES, signal, label='Full signal at EP', color='black')
    for shift, other in curves:
        axes[0].plot(TIMES, other, label=f'End bias {shift:+.2f}', linestyle='--')
        axes[1].plot(TIMES, other-signal, label=f'End bias {shift:+.2f}')
    axes[0].plot(TIMES, ep_component, label='Selected Jordan pair contribution', linestyle=':')
    axes[0].set_ylabel('Unconditioned conditional-X expectation')
    axes[1].set_ylabel('Expectation difference from EP')
    for ax in axes:
        ax.set_xlabel('Dimensionless time gamma t'); ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.suptitle('Four-qubit virtual readout: fixed q=2, local Z dephasing')
    fig.savefig(ROOT/'visualizations/route_b_n4_virtual_readout.png', dpi=180)
    print(json.dumps(output, indent=2))
    print('ALL PASS')


if __name__ == '__main__':
    main()
