#!/usr/bin/env python3
"""Sampled leading-angle categories of neural Jacobian spectra.

At each grid point, sort eigenvalues by descending real part and report
the first two rates (-Re) and angles atan2(|Im|, |Re|). Categorize the first
angle by whether it exceeds 2 degrees. A sampled leading-angle category
change marks different categories at adjacent samples; no mode tracking
or exceptional-point detection is performed. Ties follow the eigensolver
and sort order, and negative reported rates indicate growth.

An unchanged category on the sampled grid establishes no absence of a
crossing between samples or of a mode-identity change. Coupling, external
drive P and the leak-time ratio are the three displayed axes. The sigmoid
builder uses 500 updates without a fixed-point residual check.
"""
import numpy as np

TAU_E, TAU_I = 5.0, 10.0


def make_balanced_network(N, density=0.3, seed=42):
    rng = np.random.RandomState(seed)
    signs = np.ones(N)
    signs[rng.choice(N, N // 2, replace=False)] = -1
    mask = rng.random((N, N)) < density
    np.fill_diagonal(mask, False)
    weights = rng.exponential(0.3, (N, N))
    W = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if mask[i, j]:
                W[i, j] = signs[j] * weights[i, j]
    mx = np.max(np.abs(W))
    if mx > 0:
        W /= mx
    return W, signs


def build_linear_jacobian(W, signs, tau_E, tau_I, alpha):
    n = len(signs)
    J = np.zeros((n, n))
    for i in range(n):
        tau_i = tau_E if signs[i] > 0 else tau_I
        J[i, i] = -1.0 / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] / tau_i
    return J


def sigmoid(x, a, theta):
    return 1.0 / (1.0 + np.exp(np.clip(-a * (x - theta), -500, 500)))


def dsigmoid(x, a, theta):
    s = sigmoid(x, a, theta)
    return a * s * (1.0 - s)


def build_jacobian_with_sigmoid(W, signs, tau_E, tau_I, alpha, P):
    n = len(signs)
    x = np.ones(n) * 0.3
    for _ in range(500):
        inputs = alpha * W @ x + P
        for i in range(n):
            a_i = 1.3 if signs[i] > 0 else 2.0
            th_i = 4.0 if signs[i] > 0 else 3.7
            x[i] = sigmoid(inputs[i], a_i, th_i)
    inputs = alpha * W @ x + P
    J = np.zeros((n, n))
    for i in range(n):
        a_i = 1.3 if signs[i] > 0 else 2.0
        th_i = 4.0 if signs[i] > 0 else 3.7
        tau_i = tau_E if signs[i] > 0 else tau_I
        dS = dsigmoid(inputs[i], a_i, th_i)
        J[i, i] = (-1.0 + alpha * W[i, i] * dS) / tau_i
        for j in range(n):
            if i != j:
                J[i, j] = alpha * W[i, j] * dS / tau_i
    return J


def slowest_reading(J):
    """Read rates and angles of the first two roots in descending Re order.

    The reading does not match mode identities across parameter samples.
    """
    ev = np.linalg.eigvals(J)
    order = np.argsort(ev.real)[::-1]            # largest real part first, including growth
    lam0 = ev[order[0]]
    lam1 = ev[order[1]]
    def th(l):
        return float(np.degrees(np.arctan2(abs(l.imag), abs(l.real))))
    return (-lam0.real, th(lam0), -lam1.real, th(lam1))


def sweep(label, jac_of, axis_vals, axis_name):
    print(f"\n  {label}  (leading = largest real part)")
    print(f"  {axis_name:>8s}  {'rate':>7s}  {'theta':>7s}  {'2nd rate':>8s}  {'2nd th':>7s}   angle category")
    prev = None
    for v in axis_vals:
        r0, t0, r1, t1 = slowest_reading(jac_of(v))
        char = "above 2 deg" if t0 > 2.0 else "at/below 2 deg"
        switch = ""
        if prev is not None and (prev > 2.0) != (t0 > 2.0):
            switch = "   <-- sampled leading-angle category change"
        prev = t0
        print(f"  {v:8.3f}  {r0:7.4f}  {t0:7.2f}  {r1:8.4f}  {t1:7.2f}   {char}{switch}")


def main():
    W, signs = make_balanced_network(50, density=0.3, seed=42)
    print("=" * 72)
    print("Neural leading-angle census: threshold = 2 deg, on the sampled grid")
    print("  Adjacent angle categories only; no mode tracking or EP detection.")
    print(f"  Wilson-Cowan balanced E-I net, N=50, tau_E={TAU_E}, tau_I={TAU_I}")
    print("=" * 72)

    sweep("axis 1: coupling alpha (linear Jacobian, no drive)",
          lambda a: build_linear_jacobian(W, signs, TAU_E, TAU_I, a),
          [0.05, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0], "alpha")

    sweep("axis 2: external drive P (sigmoid Jacobian), alpha=0.3",
          lambda P: build_jacobian_with_sigmoid(W, signs, TAU_E, TAU_I, 0.3, P),
          [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0], "P")

    sweep("axis 3: tau_I/tau_E ratio (linear Jacobian, alpha=0.5)",
          lambda r: build_linear_jacobian(W, signs, TAU_E, TAU_E * r, 0.5),
          [1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0], "tauI/E")

    print("\n  Each marker is a sampled leading-angle category change across the 2 deg threshold.")
    print("  No marker means no category change between these sampled endpoints.")
    print("  Crossings between samples and mode-identity changes are not tested.")


if __name__ == "__main__":
    main()
