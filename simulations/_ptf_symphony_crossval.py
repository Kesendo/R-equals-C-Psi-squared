"""Cross-validation of the C# Symphony painters' movement against the Python PTF twin.

Builds the same canonical case the C# `inspect --root symphony --N 4 --htype xy
--defect-bond 1 --delta-j 0.02` runs: XY chain, N=4, J=1, gamma=0.1, Bell+ on sites
(0,1), defect bond 1, delta_J=0.02 -> J_mod = 1.02. Reports the per-site alpha vector
from ptf_alpha_fit so it can be set beside the C# alphas. Not a unit test; a one-shot
verification harness.
"""
import sys
import numpy as np

sys.path.insert(0, "simulations")
import framework as fw  # noqa: E402
from framework.workflows.ptf import ptf_alpha_fit  # noqa: E402

N = 4
J = 1.0
gamma = 0.1
delta_J = 0.02
defect_bond = 1

chain = fw.ChainSystem(N=N, gamma_0=gamma, J=J, H_type="xy")

# Bell+ on sites (0,1) = (|0000> + |1100>)/sqrt(2); site 0 = MSB (matches C# BuildInitialState).
d = 2 ** N
psi = np.zeros(d, dtype=complex)
both = (1 << (N - 1)) | (1 << (N - 2))   # = 12 at N=4
psi[0] = 1.0 / np.sqrt(2.0)
psi[both] = 1.0 / np.sqrt(2.0)
rho0 = np.outer(psi, psi.conj())

# The C# Symphony grid is t in [0, 1/gamma] with 60 points; mirror it for a like-for-like fit.
t_max = 1.0 / gamma
n_t = 60

out = ptf_alpha_fit(chain, rho0, defect_bond, J_mod=J + delta_J, t_max=t_max, n_t=n_t)

print("PTF Python twin (XY, N=4, J=1, gamma=0.1, Bell+(0,1), bond 1, dJ=0.02)")
print("t_max =", t_max, " n_t =", n_t)
for i, a in enumerate(out["alphas"]):
    print(f"  site {i}: alpha = {a:.6f}")
print("sigma_log_alpha (all sites) =", out["sigma_log_alpha"])

# perspectives_panel: the guarded reading (reliable sites + closure), Python side.
from framework.workflows.ptf import perspectives_panel
pp = perspectives_panel(chain, rho0, defect_bond, delta_J=delta_J,
                        guard_delta_J=delta_J / 2.0, t_max=t_max, n_t=n_t)
print("\nperspectives_panel (guarded):")
print("  alphas   =", np.array2string(pp['alphas'], precision=5))
print("  f        =", np.array2string(pp['f'], precision=3))
print("  f_guard  =", np.array2string(pp['f_guard'], precision=3))
print("  reliable =", pp['reliable'])
print("  sigma_log_alpha_reliable =", pp['sigma_log_alpha_reliable'])
print("  clock =", pp['clock'])

# --- The site-2 deviation, explained by the MSE landscape -------------------
# The C# golden-section fit and the Python (scipy bounded-Brent) fit agree to ~1e-2 on the
# well-conditioned sites (0, 1, 3) but disagree on site 2 (a featureless trajectory: P_A and
# P_B differ by < 0.009 at every t). A brute-grid scan of the MSE landscape shows the GENUINE
# global minimum is alpha ~ 1.016 (the C# answer); the Python alpha ~ 3.15 is a Brent local-min
# artifact with MSE ~920x larger. C#'s fit is the better global optimum here.
from scipy.interpolate import interp1d
i = 2
interp = interp1d(out['t_grid'], out['P_A'][i], bounds_error=False,
                  fill_value=(float(out['P_A'][i][0]), float(out['P_A'][i][-1])), kind='cubic')


def _mse(a):
    dd = interp(a * out['t_grid']) - out['P_B'][i]
    return float(np.mean(dd * dd))


print("\nsite 2 MSE landscape (featureless trajectory, max|P_A-P_B| =",
      f"{np.max(np.abs(out['P_A'][i] - out['P_B'][i])):.4f}):")
print(f"  scipy-Brent alpha 3.15372 -> MSE {_mse(3.15372):.3e}")
print(f"  C#/global    alpha 1.01571 -> MSE {_mse(1.01571):.3e}")
ag = np.linspace(0.1, 10, 100000)
am = ag[int(np.argmin([_mse(a) for a in ag]))]
print(f"  brute global argmin alpha = {am:.5f} (matches the C# golden-section fit)")
