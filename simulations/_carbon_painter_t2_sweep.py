"""Sweep flavor-resolved T2 anisotropy across larger N and parameters.

This is the scalable companion to _carbon_painter_t2_anisotropy.py.  The
detailed script propagates time traces, which is useful at N=4 but quickly
becomes expensive.  Here we keep the same Hückel + h_y ΣY + Holstein-Z
Liouvillian and read the T2 anisotropy directly from the slow flavor-resolved
eigenmodes:

    T2(non-field) / T2(field-only) = rate(field-only) / rate(non-field).

The single-site Pauli classifier is intentionally lightweight: it projects
each eigenoperator onto one-site X/Y/Z observables and classifies a mode by
dominant field-axis-only or non-field-axis weight.  Modes with negligible one-site weight or
mixed weight are skipped.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np


PAULI = {
    "I": np.eye(2, dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "Z": np.array([[1, 0], [0, -1]], dtype=complex),
}


def pauli_op(letters):
    op = PAULI[letters[0]]
    for letter in letters[1:]:
        op = np.kron(op, PAULI[letter])
    return op


def site_op(N, site, letter):
    letters = ["I"] * N
    letters[site] = letter
    return pauli_op(letters)


def two_site_op(N, a, b, la, lb):
    letters = ["I"] * N
    letters[a] = la
    letters[b] = lb
    return pauli_op(letters)


def bonds_for(N, topology):
    if topology == "ring":
        return [(a, (a + 1) % N) for a in range(N)]
    if topology == "chain":
        return [(a, a + 1) for a in range(N - 1)]
    if topology == "star":
        return [(0, a) for a in range(1, N)]
    if topology == "complete":
        return [(a, b) for a in range(N) for b in range(a + 1, N)]
    if topology == "empty":
        return []
    raise ValueError(f"unknown topology: {topology}")


def hueckel_H(N, topology):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for a, b in bonds_for(N, topology):
        H += two_site_op(N, a, b, "X", "X")
        H += two_site_op(N, a, b, "Y", "Y")
    return H


def zeeman_total(N, axis):
    d = 2**N
    H = np.zeros((d, d), dtype=complex)
    for site in range(N):
        H += site_op(N, site, axis)
    return H


def commutator_superop_vec(H):
    d = H.shape[0]
    I = np.eye(d, dtype=complex)
    return -1j * (np.kron(H, I) - np.kron(I, H.T))


def dissipator_superop_vec(c):
    d = c.shape[0]
    I = np.eye(d, dtype=complex)
    c_dag_c = c.conj().T @ c
    return np.kron(c, c.conj()) - 0.5 * (
        np.kron(c_dag_c, I) + np.kron(I, c_dag_c.T)
    )


def lindbladian_vec(H, c_list, gammas):
    L = commutator_superop_vec(H)
    for c, gamma in zip(c_list, gammas):
        L += gamma * dissipator_superop_vec(c)
    return L


def build_problem(N, topology, field_axis):
    H0 = hueckel_H(N, topology)
    Hfield = zeeman_total(N, field_axis)
    z_jumps = [site_op(N, site, "Z") for site in range(N)]
    site_ops = {
        (site, letter): site_op(N, site, letter)
        for site in range(N)
        for letter in ("X", "Y", "Z")
    }
    return H0, Hfield, z_jumps, site_ops


def build_L(H0, Hfield, z_jumps, h_field, gamma):
    return lindbladian_vec(H0 + h_field * Hfield, z_jumps, [gamma] * len(z_jumps))


def classify_mode(rho, N, site_ops, field_axis):
    d = 2**N
    field = 0.0
    other = 0.0
    for site in range(N):
        for letter in ("X", "Y", "Z"):
            coeff = np.trace(site_ops[(site, letter)].conj().T @ rho) / d
            if letter == field_axis:
                field += abs(coeff) ** 2
            else:
                other += abs(coeff) ** 2
    total = other + field
    if total < 1e-14:
        return "dark", other, field
    field_frac = field / total
    if field_frac < 0.01:
        return "non-field", other, field
    if field_frac > 0.99:
        return "field-only", other, field
    return "mixed", other, field


def slow_flavor_rates(L, N, site_ops, field_axis, n_modes, zero_tol=1e-9):
    d = 2**N
    eigvals, eigvecs = np.linalg.eig(L)
    order = np.argsort(eigvals.real)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    slow_non_field = None
    slow_field = None
    inspected = []

    for k in range(min(n_modes, len(eigvals))):
        eig = eigvals[k]
        if abs(eig.real) < zero_tol:
            continue
        rho = eigvecs[:, k].reshape(d, d, order="F")
        flavor, other, field = classify_mode(rho, N, site_ops, field_axis)
        inspected.append((k, eig.real, eig.imag, flavor, other, field))
        if flavor == "non-field" and slow_non_field is None:
            slow_non_field = (k, abs(eig.real), eig.imag)
        elif flavor == "field-only" and slow_field is None:
            slow_field = (k, abs(eig.real), eig.imag)
        if slow_non_field is not None and slow_field is not None:
            break

    return slow_non_field, slow_field, inspected


def parse_float_list(raw):
    return [float(part.strip()) for part in raw.split(",") if part.strip()]


def parse_int_list(raw):
    return [int(part.strip()) for part in raw.split(",") if part.strip()]


def run_sweep(args):
    Ns = parse_int_list(args.N)
    h_values = parse_float_list(args.h_y)
    gamma_values = parse_float_list(args.gamma)
    rows = []
    out_lines = []

    def emit(line=""):
        print(line, flush=True)
        out_lines.append(line)

    emit("=" * 118)
    emit("Carbon Painter spectral T2 anisotropy sweep")
    emit(
        "T2(non-field)/T2(field-only) = rate(field-only)/rate(non-field); "
        "time-domain propagation is intentionally skipped."
    )
    emit(
        f"N={Ns}, topology={args.topology}, field_axis={args.field_axis}, "
        f"h_field={h_values}, gamma={gamma_values}, n_modes={args.n_modes}"
    )
    emit("=" * 118)
    emit("")
    emit(
        f"{'N':>3}  {'topo':>5}  {'field':>5}  {'h':>8}  {'gamma':>8}  {'dim':>7}  "
        f"{'k_non':>6}  {'rate_non':>13}  {'k_fld':>5}  {'rate_fld':>13}  "
        f"{'T2 ratio':>10}  {'rate/g ratios':>15}  {'sec':>7}  status"
    )
    emit("  " + "-" * 114)

    for N in Ns:
        H0, Hfield, z_jumps, site_ops = build_problem(N, args.topology, args.field_axis)
        dim = 4**N
        for gamma in gamma_values:
            for h_field in h_values:
                t0 = time.time()
                L = build_L(H0, Hfield, z_jumps, h_field, gamma)
                slow_non_field, slow_field, inspected = slow_flavor_rates(
                    L, N, site_ops, args.field_axis, args.n_modes
                )
                elapsed = time.time() - t0
                if slow_non_field is None or slow_field is None:
                    status = "classifier-failed"
                    rows.append((N, args.topology, args.field_axis, h_field, gamma, None, None, None, status))
                    emit(
                        f"{N:>3}  {args.topology:>5}  {args.field_axis:>5}  {h_field:>8.4g}  {gamma:>8.4g}  "
                        f"{dim:>7}  {'-':>6}  {'-':>13}  {'-':>5}  {'-':>13}  "
                        f"{'-':>10}  {'-':>15}  {elapsed:>7.1f}  {status}"
                    )
                    continue

                k_non_field, rate_non_field, _ = slow_non_field
                k_field, rate_field, _ = slow_field
                t2_ratio = rate_field / rate_non_field
                scaled = f"{rate_non_field/gamma:.6g}/{rate_field/gamma:.6g}"
                status = "ok"
                rows.append((N, args.topology, args.field_axis, h_field, gamma, rate_non_field, rate_field, t2_ratio, status))
                emit(
                    f"{N:>3}  {args.topology:>5}  {args.field_axis:>5}  {h_field:>8.4g}  {gamma:>8.4g}  "
                    f"{dim:>7}  {k_non_field:>6}  {rate_non_field:>13.6e}  "
                    f"{k_field:>5}  {rate_field:>13.6e}  {t2_ratio:>10.6f}  "
                    f"{scaled:>15}  {elapsed:>7.1f}  {status}"
                )

                if args.verbose:
                    for k, re, im, flavor, other, field in inspected:
                        total = other + field
                        field_frac = field / total if total else float("nan")
                        emit(
                            f"       mode {k:>3}: Re={re:+.6e} Im={im:+.6e} "
                            f"flavor={flavor:>10} field_frac={field_frac:.6f}"
                        )
        emit("")

    emit("=" * 118)
    emit("Compact ratio table")
    emit("=" * 118)
    for N, topo, field_axis, h_field, gamma, rate_non_field, rate_field, ratio, status in rows:
        if ratio is None:
            emit(f"N={N} {topo} field={field_axis} h={h_field:g} gamma={gamma:g}: {status}")
        else:
            emit(
                f"N={N} {topo} field={field_axis} h={h_field:g} gamma={gamma:g}: "
                f"rate_non_field={rate_non_field:.10e}, rate_field={rate_field:.10e}, "
                f"T2_ratio={ratio:.10f}"
            )

    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("\n".join(out_lines), encoding="utf-8")
        emit(f"\nWritten to {out_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Sweep larger-N spectral T2 anisotropy for the Carbon Painter model."
    )
    parser.add_argument("--N", default="4,5", help="comma-separated N values")
    parser.add_argument(
        "--topology",
        default="ring",
        choices=("ring", "chain", "star", "complete", "empty"),
    )
    parser.add_argument("--h-y", default="0.25,0.5,0.75", help="comma-separated field-strength values")
    parser.add_argument("--field-axis", default="Y", choices=("X", "Y", "Z"))
    parser.add_argument("--gamma", default="0.5,1.0,2.0", help="comma-separated gamma values")
    parser.add_argument(
        "--n-modes",
        type=int,
        default=80,
        help="number of slowest eigenmodes inspected for flavor classification",
    )
    parser.add_argument("--output", default="simulations/results/carbon_painter_t2_sweep.txt")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    run_sweep(args)


if __name__ == "__main__":
    main()
