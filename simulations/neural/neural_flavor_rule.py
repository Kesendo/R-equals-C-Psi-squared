"""Neural inheritance test for the flavor-resolved relaxation rule.

This is not a claim that neural networks are quantum systems.  It tests the
abstract rule in a Wilson-Cowan-style linearization:

    local non-commuting E/I interaction + leak/damping + graph topology
        -> flavor-resolved slow modes with different lifetimes.

Each graph node carries an excitatory/inhibitory pair (E_i, I_i).  The
linearized generator is

    dE/dt = (-gamma_E I + c A) E - h I
    dI/dt =  h E + (-gamma_I I - c A) I

where A is the graph adjacency.  The signs are chosen as a minimal stable
E/I push-pull model: E recurrence is positive, I recurrence is suppressive,
and local E/I cross-coupling creates the non-commuting axis.

Modes are classified by their one-node E vs I weight.  The reported ratio is

    T(E-dominant) / T(I-dominant) = rate(I-dominant) / rate(E-dominant).

Ratio > 1 means E-dominant activity is the longer-lived channel.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def adjacency(N, topology):
    A = np.zeros((N, N), dtype=float)
    if topology == "empty":
        return A
    if topology == "chain":
        for i in range(N - 1):
            A[i, i + 1] = A[i + 1, i] = 1.0
    elif topology == "ring":
        for i in range(N):
            A[i, (i + 1) % N] = A[(i + 1) % N, i] = 1.0
    elif topology == "star":
        for i in range(1, N):
            A[0, i] = A[i, 0] = 1.0
    elif topology == "complete":
        A[:] = 1.0
        np.fill_diagonal(A, 0.0)
    else:
        raise ValueError(topology)
    # Normalize by spectral radius so c has comparable meaning per topology.
    radius = max(abs(np.linalg.eigvals(A))) if N else 1.0
    if radius > 0:
        A = A / radius
    return A


def build_jacobian(N, topology, h, coupling, gamma_e, gamma_i):
    A = adjacency(N, topology)
    I = np.eye(N)
    top_left = -gamma_e * I + coupling * A
    top_right = -h * I
    bottom_left = h * I
    bottom_right = -gamma_i * I - coupling * A
    return np.block([[top_left, top_right], [bottom_left, bottom_right]])


def classify_mode(v, N, tol=1e-12):
    e_w = float(np.sum(np.abs(v[:N]) ** 2))
    i_w = float(np.sum(np.abs(v[N:]) ** 2))
    total = e_w + i_w
    if total < tol:
        return "dark", e_w, i_w
    i_frac = i_w / total
    if i_frac < 0.35:
        return "E-dom", e_w, i_w
    if i_frac > 0.65:
        return "I-dom", e_w, i_w
    return "mixed", e_w, i_w


def slow_flavor_rates(J, N, n_modes):
    vals, vecs = np.linalg.eig(J)
    stable = vals.real < -1e-10
    order = np.argsort(vals.real)[::-1]
    slow_e = None
    slow_i = None
    inspected = []
    for k in order[: min(n_modes, len(order))]:
        if not stable[k]:
            inspected.append((k, vals[k], "unstable", 0.0, 0.0))
            continue
        flavor, e_w, i_w = classify_mode(vecs[:, k], N)
        inspected.append((k, vals[k], flavor, e_w, i_w))
        if flavor == "E-dom" and slow_e is None:
            slow_e = (k, -float(vals[k].real), float(vals[k].imag))
        elif flavor == "I-dom" and slow_i is None:
            slow_i = (k, -float(vals[k].real), float(vals[k].imag))
        if slow_e is not None and slow_i is not None:
            break
    return slow_e, slow_i, inspected


def parse_int_list(raw):
    return [int(x.strip()) for x in raw.split(",") if x.strip()]


def run(args):
    Ns = parse_int_list(args.N)
    topologies = [x.strip() for x in args.topology.split(",") if x.strip()]
    lines = []

    def emit(s=""):
        print(s, flush=True)
        lines.append(s)

    emit("=" * 108)
    emit("Neural flavor-rule inheritance test")
    emit(
        f"N={Ns}, topology={topologies}, h={args.h}, coupling={args.coupling}, "
        f"gamma_E={args.gamma_e}, gamma_I={args.gamma_i}"
    )
    emit("=" * 108)
    emit("")
    emit(
        f"{'N':>3}  {'topology':>8}  {'dim':>5}  {'k_E':>5}  {'rate_E':>12}  "
        f"{'k_I':>5}  {'rate_I':>12}  {'T ratio':>10}  {'max Re':>12}  status"
    )
    emit("  " + "-" * 96)

    for N in Ns:
        for topo in topologies:
            J = build_jacobian(N, topo, args.h, args.coupling, args.gamma_e, args.gamma_i)
            vals = np.linalg.eigvals(J)
            max_re = float(np.max(vals.real))
            slow_e, slow_i, inspected = slow_flavor_rates(J, N, args.n_modes)
            if max_re >= -1e-10:
                status = "unstable"
            elif slow_e is None or slow_i is None:
                status = "classifier-failed"
            else:
                status = "ok"

            if status != "ok":
                emit(
                    f"{N:>3}  {topo:>8}  {2*N:>5}  {'-':>5}  {'-':>12}  "
                    f"{'-':>5}  {'-':>12}  {'-':>10}  {max_re:>12.6e}  {status}"
                )
                if args.verbose:
                    for k, val, flavor, e_w, i_w in inspected:
                        total = e_w + i_w
                        i_frac = i_w / total if total else float("nan")
                        emit(
                            f"       mode {k:>3}: lambda={val.real:+.6e}{val.imag:+.6e}i "
                            f"flavor={flavor:>8} I-frac={i_frac:.4f}"
                        )
                continue

            k_e, rate_e, im_e = slow_e
            k_i, rate_i, im_i = slow_i
            ratio = rate_i / rate_e
            emit(
                f"{N:>3}  {topo:>8}  {2*N:>5}  {k_e:>5}  {rate_e:>12.6e}  "
                f"{k_i:>5}  {rate_i:>12.6e}  {ratio:>10.6f}  {max_re:>12.6e}  {status}"
            )
            if args.verbose:
                emit(f"       slow E imag={im_e:+.6e}, slow I imag={im_i:+.6e}")

    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("\n".join(lines), encoding="utf-8")
        emit(f"\nWritten to {out}")


def main():
    parser = argparse.ArgumentParser(description="Test flavor-rule inheritance in neural E/I networks.")
    parser.add_argument("--N", default="4,5,6,8")
    parser.add_argument("--topology", default="empty,chain,ring,star,complete")
    parser.add_argument("--h", type=float, default=0.5)
    parser.add_argument("--coupling", type=float, default=0.25)
    parser.add_argument("--gamma-e", type=float, default=1.0)
    parser.add_argument("--gamma-i", type=float, default=2.0)
    parser.add_argument("--n-modes", type=int, default=64)
    parser.add_argument("--output", default="simulations/results/neural_flavor_rule.txt")
    parser.add_argument("--verbose", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
