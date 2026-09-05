# Neural Quantum-Translation Rebuild Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rebuild the Neural bridge from repository-owned quantum objects, with gates that distinguish exact transfer, scoped numerical sighting, and failed/open translation.

**Architecture:** Put the neural matrix identity and its separating diagnostics in one import-safe Python module and one small `MirrorWorld` typed object. The documents consume those owners in dependency order: theorem and formula registry first, then the five Neural pages, then every sibling carrying the same claims. Event records remain records; current-facing guides and producer text carry only current truth.

**Tech Stack:** Python 3, NumPy, SciPy, pytest, C#/.NET 10, xUnit, Markdown, repository Claim/F registries.

---

### Task 1: Canonical Python matrix identity and separating residuals

**Files:**
- Create: `simulations/neural/neural_palindrome.py`
- Create: `simulations/neural/tests/test_neural_palindrome.py`

- [ ] **Step 1: Write failing tests for the exact identity, the old false-positive instrument, and complex multiset pairing**

```python
import pathlib
import sys

import numpy as np

NEURAL = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(NEURAL))

from neural_palindrome import (
    fitted_offdiagonal_residual,
    is_involution,
    scalar_center_residual,
    spectral_pairing_error,
)


def test_exact_complex_pair_has_zero_scalar_residual_and_pairing_error():
    J = np.array([[-0.2, -0.1], [0.1, -0.1]])
    perm = np.array([1, 0])
    assert is_involution(perm)
    assert scalar_center_residual(J, perm, s=0.15) == 0.0
    assert spectral_pairing_error(np.linalg.eigvals(J), s=0.15) < 1e-13
    assert np.max(np.abs(np.linalg.eigvals(J).imag)) > 0.08


def test_fixed_seat_is_rejected_by_scalar_residual_but_not_old_fit():
    J = np.diag([-0.2, -0.1, -0.2])
    perm = np.array([1, 0, 2])
    assert fitted_offdiagonal_residual(J, perm) == 0.0
    assert scalar_center_residual(J, perm, s=0.15) > 0.09


def test_pairing_error_preserves_multiplicity():
    values = np.array([-0.2 + 0j, -0.1 + 0j, -0.1 + 0j])
    assert spectral_pairing_error(values, s=0.15) > 0.09
```

- [ ] **Step 2: Run the tests and verify that collection fails because the module does not exist**

Run: `python -m pytest simulations/neural/tests/test_neural_palindrome.py -q`

Expected: FAIL during collection with `ModuleNotFoundError: No module named 'neural_palindrome'`.

- [ ] **Step 3: Implement the minimal canonical matrix helpers**

```python
from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment


def permutation_matrix(perm: np.ndarray) -> np.ndarray:
    perm = np.asarray(perm, dtype=int)
    n = len(perm)
    if sorted(perm.tolist()) != list(range(n)):
        raise ValueError("perm must contain each seat exactly once")
    Q = np.zeros((n, n), dtype=float)
    Q[np.arange(n), perm] = 1.0
    return Q


def is_involution(perm: np.ndarray) -> bool:
    perm = np.asarray(perm, dtype=int)
    return sorted(perm.tolist()) == list(range(len(perm))) and np.array_equal(perm[perm], np.arange(len(perm)))


def _normalised_norm(residual: np.ndarray, J: np.ndarray) -> float:
    scale = np.linalg.norm(J)
    return float(np.linalg.norm(residual) / scale) if scale else float(np.linalg.norm(residual))


def scalar_center_residual(J: np.ndarray, perm: np.ndarray, s: float) -> float:
    Q = permutation_matrix(perm)
    return _normalised_norm(Q @ J @ Q.T + J + 2.0 * s * np.eye(len(J)), J)


def fitted_offdiagonal_residual(J: np.ndarray, perm: np.ndarray) -> float:
    Q = permutation_matrix(perm)
    QJQ = Q @ J @ Q.T
    fitted = -(np.diag(QJQ) + np.diag(J)) / 2.0
    R = QJQ + J + 2.0 * np.diag(fitted)
    return _normalised_norm(R - np.diag(np.diag(R)), J)


def spectral_pairing_error(values: np.ndarray, s: float) -> float:
    values = np.asarray(values, dtype=complex)
    targets = -values - 2.0 * s
    cost = np.abs(values[:, None] - targets[None, :])
    rows, cols = linear_sum_assignment(cost)
    return float(cost[rows, cols].max(initial=0.0))
```

- [ ] **Step 4: Run the focused tests and verify all three pass**

Run: `python -m pytest simulations/neural/tests/test_neural_palindrome.py -q`

Expected: `3 passed`.

- [ ] **Step 5: Commit the canonical identity helpers**

```bash
git add simulations/neural/neural_palindrome.py simulations/neural/tests/test_neural_palindrome.py
git commit -m "test(neural): separate scalar palindrome from fitted residual"
```

### Task 2: Exact-network construction, ensemble census, and mode-transport gate

**Files:**
- Modify: `simulations/neural/neural_palindrome.py`
- Modify: `simulations/neural/tests/test_neural_palindrome.py`
- Create: `simulations/neural/neural_translation_gate.py`

- [ ] **Step 1: Add failing tests for the deterministic ensemble and `Q`-transported partner eigenspaces**

```python
from neural_palindrome import exact_ensemble_census, make_exact_network, partner_subspace_error


def test_exact_ensemble_reproduces_the_committed_counts():
    assert exact_ensemble_census() == [
        (0.5, 200, 24, 0),
        (1.5, 200, 110, 1),
        (3.0, 200, 149, 15),
        (5.0, 200, 159, 24),
        (10.0, 200, 167, 45),
    ]


def test_q_transports_each_isolated_eigenspace_to_its_partner():
    J, perm, s = make_exact_network(n=10, tau_e=5.0, tau_i=10.0, alpha=0.5, seed=42)
    assert scalar_center_residual(J, perm, s) < 1e-13
    assert partner_subspace_error(J, perm, s, cluster_tol=1e-8) < 1e-8


def test_q_transport_is_basis_invariant_at_a_degeneracy():
    J = -0.15 * np.eye(4)
    perm = np.array([1, 0, 3, 2])
    assert partner_subspace_error(J, perm, s=0.15, cluster_tol=1e-10) < 1e-12
```

- [ ] **Step 2: Run the focused tests and verify the new imports fail**

Run: `python -m pytest simulations/neural/tests/test_neural_palindrome.py -q`

Expected: FAIL with an import error naming `exact_ensemble_census`.

- [ ] **Step 3: Add the exact builder and deterministic census**

Move the existing construction from `veffect_exact.py` into import-safe functions, preserving `RandomState`, iteration order, normalization, and Dale-column signs. Add:

```python
def build_exact_weights(n, n_exc, tau_e, tau_i, density=0.3, seed=42):
    rng = np.random.RandomState(seed)
    signs = np.ones(n)
    signs[rng.choice(n, n - n_exc, replace=False)] = -1
    e_idx = np.where(signs > 0)[0]
    i_idx = np.where(signs < 0)[0]
    perm = np.arange(n)
    for e, i in zip(e_idx, i_idx):
        perm[e], perm[i] = i, e
    W = np.zeros((n, n))
    mask = rng.random((n, n)) < density
    np.fill_diagonal(mask, False)
    for i in range(n):
        for j in range(n):
            if i == j or i >= j or not mask[i, j]:
                continue
            qi, qj = perm[i], perm[j]
            base = rng.exponential(0.3)
            W[i, j] = signs[j] * base
            row_tau = tau_e if signs[i] > 0 else tau_i
            partner_tau = tau_e if signs[qi] > 0 else tau_i
            W[qi, qj] = -(partner_tau / row_tau) * W[i, j]
    scale = np.max(np.abs(W))
    if scale:
        W /= scale
    return W, signs, perm


def build_linear_jacobian(W, signs, tau_e, tau_i, alpha):
    tau = np.where(signs > 0, tau_e, tau_i)
    J = alpha * W / tau[:, None]
    J[np.diag_indices_from(J)] = -1.0 / tau
    return J


def make_exact_network(n, tau_e, tau_i, alpha, seed, density=0.3):
    W, signs, perm = build_exact_weights(n, n // 2, tau_e, tau_i, density, seed)
    return build_linear_jacobian(W, signs, tau_e, tau_i, alpha), perm, 0.5 * (1.0 / tau_e + 1.0 / tau_i)


def exact_ensemble_census():
    rows = []
    for alpha in (0.5, 1.5, 3.0, 5.0, 10.0):
        complex_count = unstable_count = exact_count = 0
        for seed in range(200):
            J, perm, s = make_exact_network(10, 5.0, 10.0, alpha, seed)
            if scalar_center_residual(J, perm, s) < 1e-13:
                exact_count += 1
            eig = np.linalg.eigvals(J)
            complex_count += int(np.max(np.abs(eig.imag)) > 1e-8)
            unstable_count += int(np.max(eig.real) > 0.0)
        rows.append((alpha, exact_count, complex_count, unstable_count))
    return rows
```

- [ ] **Step 4: Implement basis-invariant eigenspace transport**

Cluster eigenvalues by complex distance, obtain each invariant subspace with an
ordered complex Schur decomposition, transport it with `Q`, compare it to the
partner invariant subspace using principal angles, and return the largest sine.
The function compares subspaces, not individual eigenvectors, so a degenerate
basis rotation cannot change the verdict.

```python
from scipy.linalg import schur, subspace_angles


def cluster_eigenvalues(values, tol):
    remaining = set(range(len(values)))
    clusters = []
    while remaining:
        first = min(remaining)
        members = {i for i in remaining if abs(values[i] - values[first]) <= tol}
        remaining -= members
        clusters.append((complex(np.mean(values[list(members)])), len(members)))
    return clusters


def schur_subspace(J, value, multiplicity, tol):
    _, vectors, selected = schur(
        np.asarray(J, dtype=complex),
        output="complex",
        sort=lambda candidate: abs(candidate - value) <= tol,
    )
    if selected != multiplicity:
        raise AssertionError(f"expected {multiplicity} Schur vectors, found {selected}")
    return vectors[:, :selected]


def partner_subspace_error(J, perm, s, cluster_tol=1e-8):
    Q = permutation_matrix(perm)
    values = np.linalg.eigvals(J)
    clusters = cluster_eigenvalues(values, cluster_tol)
    worst = 0.0
    for value, multiplicity in clusters:
        partner = -value - 2.0 * s
        source = schur_subspace(J, value, multiplicity, cluster_tol)
        target = schur_subspace(J, partner, multiplicity, cluster_tol)
        angles = subspace_angles(Q @ source, target)
        worst = max(worst, float(np.sin(angles).max(initial=0.0)))
    return worst
```

- [ ] **Step 5: Add a self-validating command that prints named gates and exits nonzero on failure**

`neural_translation_gate.py` must run the two-seat complex counterexample, the odd fixed-seat rejection, multiset-multiplicity rejection, five ensemble rows, and exact mode transport. End with exactly `ALL NEURAL TRANSLATION GATES PASS` after assertions.

- [ ] **Step 6: Run the tests and command**

Run: `python -m pytest simulations/neural/tests/test_neural_palindrome.py -q`

Expected: all tests pass.

Run: `python simulations/neural/neural_translation_gate.py`

Expected: the five census rows match `24/110/149/159/167` complex and `0/1/15/24/45` unstable; final line `ALL NEURAL TRANSLATION GATES PASS`.

- [ ] **Step 7: Commit the from-below gate suite**

```bash
git add simulations/neural/neural_palindrome.py simulations/neural/neural_translation_gate.py simulations/neural/tests/test_neural_palindrome.py
git commit -m "feat(neural): gate exact pairing and mode transport"
```

### Task 3: Repair existing producers and gate the drive sweep

**Files:**
- Modify: `simulations/neural/veffect_exact.py`
- Modify: `simulations/neural/veffect_and_heat.py`
- Modify: `simulations/neural/cpsi_two_perspectives.py`
- Modify: `simulations/neural/celegans_trichotomy.py`
- Modify: `simulations/neural/validation_checks.py`
- Modify: `simulations/neural/tests/test_neural_palindrome.py`

- [ ] **Step 1: Add a failing fixed-point convergence test**

```python
from neural_palindrome import solve_fixed_point


def test_drive_fixed_point_reports_a_small_equation_residual():
    W, signs = make_balanced_dale_network(n=50, n_exc=25, density=0.3, seed=42)
    _, residual = solve_fixed_point(W, signs, alpha=0.3, drive=4.0, tol=1e-12, max_iter=5000)
    assert residual < 1e-10
```

- [ ] **Step 2: Run the test and verify `solve_fixed_point` is missing**

Run: `python -m pytest simulations/neural/tests/test_neural_palindrome.py::test_drive_fixed_point_reports_a_small_equation_residual -q`

Expected: FAIL with an import error.

- [ ] **Step 3: Implement a declared fixed-point solver**

Add the random Dale builder from `veffect_and_heat.py`, then return both the
iterate and `max(abs(x - sigmoid(alpha*W@x + P)))`; stop only when that residual
is below `tol`, otherwise raise at `max_iter`:

```python
def make_balanced_dale_network(n, n_exc, density=0.3, seed=42):
    rng = np.random.RandomState(seed)
    signs = np.ones(n)
    signs[rng.choice(n, n - n_exc, replace=False)] = -1
    mask = rng.random((n, n)) < density
    np.fill_diagonal(mask, False)
    weights = rng.exponential(0.3, (n, n))
    W = mask * weights * signs[None, :]
    scale = np.max(np.abs(W))
    if scale:
        W /= scale
    return W, signs


def _sigmoid_by_type(inputs, signs):
    a = np.where(signs > 0, 1.3, 2.0)
    theta = np.where(signs > 0, 4.0, 3.7)
    argument = np.clip(-a * (inputs - theta), -500, 500)
    return 1.0 / (1.0 + np.exp(argument))


def solve_fixed_point(W, signs, alpha, drive, tol=1e-12, max_iter=5000):
    x = np.full(len(signs), 0.3)
    for _ in range(max_iter):
        candidate = _sigmoid_by_type(alpha * W @ x + drive, signs)
        residual = float(np.max(np.abs(candidate - x)))
        x = candidate
        if residual < tol:
            equation_residual = float(np.max(np.abs(x - _sigmoid_by_type(alpha * W @ x + drive, signs))))
            return x, equation_residual
    raise RuntimeError(f"fixed point did not converge below {tol:g}; last residual={residual:g}")
```

- [ ] **Step 4: Make existing scripts import-safe and consume the canonical helpers**

Wrap executable sections in `main()`/`if __name__ == "__main__":`. Replace local fitted residuals where a theorem verdict is printed with `scalar_center_residual`; retain `fitted_offdiagonal_residual` only under that explicit name. Use one `frequency_tolerance` for both activity and correlation counts in each table.

- [ ] **Step 5: Replace stale producer narratives with current measurements**

Use these exact meanings in module docstrings and console headings:

```text
veffect_exact.py: coupling and drive sweeps on constructed networks; the coupled odd-seat system is not an exact palindrome and no V-effect mechanism is claimed.
veffect_and_heat.py: synthetic external-drive sweep; P is not temperature or metabolic energy; frequency counts are measurements of this model.
cpsi_two_perspectives.py: historical greedy real-part matcher, retained as an instrument control; it is not a test of F36 or Qv transport.
celegans_trichotomy.py: equality with the degree-preserved null does not establish degree-distribution causality because the withdrawn metric is insensitive to the rewiring.
```

- [ ] **Step 6: Run every affected cheap producer**

Run:

```powershell
python simulations/neural/veffect_exact.py
python simulations/neural/veffect_and_heat.py
$env:PYTHONIOENCODING='utf-8'; python simulations/neural/cpsi_two_perspectives.py
python simulations/neural/celegans_trichotomy.py
python simulations/neural/validation_checks.py
```

Expected: each exits 0; no output says exact palindrome implies silence, coupling releases symmetry-held modes, P is heat/metabolism, or the degree-preserving result establishes causality.

- [ ] **Step 7: Commit the producer repair**

```bash
git add simulations/neural/veffect_exact.py simulations/neural/veffect_and_heat.py simulations/neural/cpsi_two_perspectives.py simulations/neural/celegans_trichotomy.py simulations/neural/validation_checks.py simulations/neural/tests/test_neural_palindrome.py
git commit -m "fix(neural): make producers report the measured objects"
```

### Task 4: Carry F36/F37 into the sober typed base

**Files:**
- Create: `compute/MirrorWorld/NeuralPalindrome.cs`
- Create: `compute/MirrorWorld.Tests/NeuralPalindromeTests.cs`
- Modify: `compute/MirrorWorld/Formulas.cs`
- Modify: `compute/MirrorWorld/Program.cs`
- Modify: `compute/MirrorWorld/README.md`

- [ ] **Step 1: Write failing xUnit tests for involution, scalar residual, fixed-seat rejection, and F37**

```csharp
namespace MirrorWorld.Tests;

public class NeuralPalindromeTests
{
    [Fact]
    public void ExactComplexPair_MeetsF36EntryWise()
    {
        double[,] j = { { -0.2, -0.1 }, { 0.1, -0.1 } };
        Assert.Equal(0.0, NeuralPalindrome.MaxResidual(j, new[] { 1, 0 }, 0.15));
    }

    [Fact]
    public void FixedSeat_FailsTheScalarCentreCondition()
    {
        double[,] j = { { -0.2, 0, 0 }, { 0, -0.1, 0 }, { 0, 0, -0.2 } };
        Assert.Equal(0.1, NeuralPalindrome.MaxResidual(j, new[] { 1, 0, 2 }, 0.15), 12);
    }

    [Fact]
    public void F37PairSum_IsMinusTwiceTheF36Centre()
    {
        double s = NeuralPalindrome.Centre(5.0, 10.0);
        Assert.Equal(-2.0 * s, Formulas.F37_NeuralPairSum(5.0, 10.0), 12);
    }
}
```

- [ ] **Step 2: Run the focused tests and verify the type is missing**

Run: `dotnet test compute/MirrorWorld.Tests --filter FullyQualifiedName~NeuralPalindromeTests`

Expected: build FAIL naming `NeuralPalindrome`.

- [ ] **Step 3: Implement the typed identity without an eigensolver**

```csharp
namespace MirrorWorld;

public static class NeuralPalindrome
{
    public static double Centre(double tauE, double tauI)
    {
        if (tauE <= 0 || tauI <= 0) throw new ArgumentOutOfRangeException();
        return 0.5 * (1.0 / tauE + 1.0 / tauI);
    }

    public static bool IsInvolution(int[] permutation) =>
        permutation.Order().SequenceEqual(Enumerable.Range(0, permutation.Length)) &&
        Enumerable.Range(0, permutation.Length).All(i => permutation[permutation[i]] == i);

    public static double MaxResidual(double[,] j, int[] permutation, double s)
    {
        int n = j.GetLength(0);
        if (j.GetLength(1) != n || permutation.Length != n || !IsInvolution(permutation))
            throw new ArgumentException("F36 requires a square J and an involutive seat permutation Q");
        double max = 0;
        for (int i = 0; i < n; i++)
        for (int k = 0; k < n; k++)
        {
            double residual = j[permutation[i], permutation[k]] + j[i, k] + (i == k ? 2 * s : 0);
            max = Math.Max(max, Math.Abs(residual));
        }
        return max;
    }
}
```

- [ ] **Step 4: Narrow the F36 comment and expose a `neural` run mode**

The run mode prints the exact two-seat example, its complex eigenvalues from the closed quadratic, the odd fixed-seat residual, and `F37_NeuralPairSum(5,10)=-0.3`. It must say `conditional matrix identity`, not `brain symmetry`.

- [ ] **Step 5: Run typed tests and the sober witness**

Run: `dotnet test compute/MirrorWorld.Tests --filter FullyQualifiedName~NeuralPalindromeTests`

Expected: all focused tests PASS.

Run: `dotnet run --project compute/MirrorWorld -- neural`

Expected: exact example residual `0`, complex pair `-0.15 +/- 0.0866i`, fixed-seat residual `0.1`, pair sum `-0.3`.

- [ ] **Step 6: Commit the typed owner**

```bash
git add compute/MirrorWorld/NeuralPalindrome.cs compute/MirrorWorld.Tests/NeuralPalindromeTests.cs compute/MirrorWorld/Formulas.cs compute/MirrorWorld/Program.cs compute/MirrorWorld/README.md
git commit -m "feat(mirrorworld): carry the conditional neural palindrome"
```

### Task 5: Repair the theorem and formula owners

**Files:**
- Modify: `docs/neural/proofs/PROOF_PALINDROME_NEURAL.md`
- Modify: `docs/ANALYTICAL_FORMULAS.md`
- Modify: `compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs`

- [ ] **Step 1: Rewrite the theorem statement around the exact operator object**

State:

```text
Let J = D + W_eff be a finite complex neural Jacobian, let Q be an involutive
permutation matrix, and let s be a scalar. The identity
    Q J Q + J + 2 s I = 0
holds iff its diagonal and off-diagonal entries hold separately. It implies
multiset invariance under mu -> -mu - 2s. If W_eff = alpha T W and alpha is
nonzero, the off-diagonal entry condition may equivalently be divided through
by alpha and written on W; at alpha=0 the W condition is not necessary.
```

Make clear that complex pairs and instability are allowed, Dale fixes only signs on nonzero support, and `Q^2=I` is a theorem hypothesis.

- [ ] **Step 2: Replace the C. elegans “Verification” authority with a scoped record**

Keep the historical March table only as an experiment record linked to the withdrawal. The current result is the support-count obstruction `253` non-empty excitatory rows versus `18` inhibitory rows; no qualifying swap exists. Do not call unreconstructable 1.6-percent rows informative evidence.

- [ ] **Step 3: Narrow F36/F37 in the registry**

F36 must say `Tier 1 derived algebra + constructed numerical gate`, valid only for an involutive `Q`, scalar `s`, and the full matrix identity. F37 remains the exact conditional pair sum. Remove the stale sentence that F36 replaces a connectome quality assessment; that assessment is withdrawn rather than replaced by the theorem.

- [ ] **Step 4: Close only the typed half-carriage portion of the OpenArc**

Update the existing `f_registry_meets_the_typed_layer` item to record that F36/F37 now have `MirrorWorld.NeuralPalindrome`, `Formulas.F37_NeuralPairSum`, focused tests, and the `neural` run mode. Do not close unrelated missing F-formulas in that umbrella arc.

- [ ] **Step 5: Run owner verification**

Run:

```powershell
python simulations/neural/neural_translation_gate.py
dotnet test compute/MirrorWorld.Tests --filter FullyQualifiedName~NeuralPalindromeTests
rg -n "arbitrary permutation|Dale.s Law.*sufficient|purely real|automatically satisfied" docs/neural/proofs/PROOF_PALINDROME_NEURAL.md docs/ANALYTICAL_FORMULAS.md
```

Expected: gates PASS; phrase search returns no live overclaim.

- [ ] **Step 6: Commit the canonical owners**

```bash
git add docs/neural/proofs/PROOF_PALINDROME_NEURAL.md docs/ANALYTICAL_FORMULAS.md compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs
git commit -m "docs(neural): state the conditional palindrome exactly"
```

### Task 6: Rebuild the five Neural documents from the translation map

**Files:**
- Modify: `docs/neural/README.md`
- Modify: `docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md`
- Modify: `docs/neural/V_EFFECT_NEURAL.md`
- Modify: `docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md`

- [ ] **Step 1: Rewrite the README as the current front door**

Use this order: `What transfers exactly`, `What was tested`, `What failed`, `What remains open`, `How to rerun`. Say that no biological brain is known here to satisfy F36; C. elegans fails the tested support condition. Define centre as `-s` and reflection shift as `-2s`. Remove the biological thermal-window and Yerkes-Dodson mechanism. Describe `V_EFFECT_NEURAL.md` as a coupling/drive census with an open mechanism, not a 2-times law.

- [ ] **Step 2: Rebuild the algebra document around four authority layers**

The four layers are: exact conditional theorem; constructed exact networks; C. elegans support-count null; gated mode transport. Delete the 96/97/74-percent greedy-matcher result from current findings. Replace it with the exact statement `Jv=mu v => J(Qv)=(-mu-2s)(Qv)` and the basis-invariant subspace gate. Keep withdrawn connectome tables only in a clearly bounded event-record section.

- [ ] **Step 3: Rebuild the V-effect document around surviving measurements**

Lead with: coupling changes activity/correlation frequency counts non-monotonically in the specified synthetic model; no palindrome mechanism is established. Separate the odd-seat coupled sweep from the exact single-network ensemble. Call `P` external drive. Scope the `124` correlation-frequency maximum to `N=50`, seed 42, declared tolerance, and a converged fixed point. Remove `life operates`, `metabolic energy`, `oscillation requires breaking`, `strong coupling overwhelms the palindrome`, and the quantum-V mechanism transfer.

- [ ] **Step 4: Keep the mechanism proof as the falsification owner and link the new gates**

Preserve the two-seat complex counterexample, scalar-versus-fitted residual, odd fixed seat, correlation-count threshold issue, and drive null. Replace ad-hoc reproducibility wording with links to `neural_translation_gate.py` and its exact assertions. Do not append a replacement causal story.

- [ ] **Step 5: Normalize living-document dates and scripts**

Each living Neural page gets one line: `last refreshed 2026-09-05 (the change history lives in git)`. Document both POSIX and PowerShell UTF-8 commands where needed. Do not retain struck-through, withdrawn, or superseded current prose; git owns the editing history.

- [ ] **Step 6: Run the Neural gate and local link checker**

Run the gate, then use a read-only Markdown link parser over `docs/neural/**/*.md` and require zero missing local targets.

- [ ] **Step 7: Commit the Neural subtree rebuild**

```bash
git add docs/neural/README.md docs/neural/ALGEBRAIC_PALINDROME_NEURAL.md docs/neural/V_EFFECT_NEURAL.md docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md
git commit -m "docs(neural): rebuild the bridge from gated objects"
```

### Task 7: Carry the corrected trail through every living sibling

**Files:**
- Modify: `experiments/NEURAL_CLOCK_TWO_HANDS.md`
- Modify: `docs/WHAT_WE_FOUND.md`
- Modify: `docs/READING_GUIDE.md`
- Modify: `hypotheses/UNIVERSAL_PALINDROME_CONDITION.md`
- Modify: `hypotheses/THE_PATTERN_RECOGNIZES_ITSELF.md`
- Modify: `hypotheses/COMPLEXITY_THRESHOLD.md`
- Modify: `hypotheses/EVOLUTION_AS_CRYSTALLIZATION.md`
- Modify: `review/EMERGING_QUESTIONS.md`
- Modify: `review/OPEN_QUESTIONS_INDEX.md`
- Modify: `simulations/neural/README.md`
- Modify as found by the residual sweep: other living Markdown/C#/Python files carrying the same claims

- [ ] **Step 1: Repair the Clock sibling without erasing its event**

Keep `trace(J)/N` as the mean real eigenvalue and clock-rate reading. Delete claims that all graph effects live in imaginary parts, that exact palindromes are silent, that the mediator awakens a quantum V-effect, or that worm subcircuits possess the palindrome. Point mode pairing to F36 and state that trace alone cannot witness it.

- [ ] **Step 2: Repair overview and guide claims**

In `WHAT_WE_FOUND.md` and `READING_GUIDE.md`, replace broad transfer lists with the exact conditional identity, constructed gate, connectome null, and open neural mechanism. The quarter threshold, quantum V-effect, Hopf/fold identity, and greedy character-swap percentage do not transfer.

- [ ] **Step 3: Repair hypotheses without promoting or burying them**

Remove current-facing assertions that V-effect is inevitable under palindromic symmetry, that the data show a life mechanism, or that a thermal window is biological. Preserve genuinely open questions as questions with explicit discriminating gates. Delete old editing-history paragraphs; do not turn false claims into discounted-tier truths.

- [ ] **Step 4: Repair operational documentation and open-question indexes**

Make `simulations/neural/README.md` describe what each producer currently measures and list `neural_translation_gate.py` as the canonical entry point. Update question indexes to point at the surviving open mechanism rather than withdrawn conclusions.

- [ ] **Step 5: Run a whole-repository tell-tale sweep**

Run:

```powershell
rg -n -i "neural noble gas|perfectly dead|perfect balance is silence|symmetry.*releases.*oscillat|oscillation requires.*palindrom|life operates|metabolic energy|biological thermal window|96%|97%|character.swap fidelity|neural fold.*Hopf|V.Effect.*inevitable|2.? decay law" . -g '*.md' -g '*.py' -g '*.cs'
```

Classify every remaining match as a current claim, an explicit falsification/caught-error record, or an unrelated quantum statement. Repair all current Neural matches; leave event records intact.

- [ ] **Step 6: Commit the propagated current truth**

Stage only files actually changed by this task, inspect `git diff --cached --check`, then:

```bash
git commit -m "docs(neural): carry the corrected translation repo-wide"
```

### Task 8: Append the error ledger and private review index

**Files:**
- Modify: `docs/CAUGHT_ERRORS.md`
- Modify: `docs/superpowers/REVIEW_LOG.md`

- [ ] **Step 1: Append one new caught-error entry**

Record, without editing older rows: stale biological front door; synthetic drive promoted to thermal/life; F36 scoped to Dale alone; theorem omitted involution/zero-gain scope; 96-percent matcher tested the wrong object; current producer text repeated the refuted V-effect mechanism. Name the from-below corrections and the separating gates.

- [ ] **Step 2: Append the private review-log row**

Use one row with date, the five Neural files plus siblings swept, physics/math/operational/future-us lenses, genuine-break verdict, and anchors `F36/F37 + neural_translation_gate.py + MirrorWorld NeuralPalindrome`.

- [ ] **Step 3: Verify append-only behavior**

Run: `git diff -- docs/CAUGHT_ERRORS.md docs/superpowers/REVIEW_LOG.md`

Expected: only additions at the end of each file; no earlier line changes.

- [ ] **Step 4: Commit the review record**

```bash
git add docs/CAUGHT_ERRORS.md
git add -f docs/superpowers/REVIEW_LOG.md
git commit -m "docs: record the neural translation repair"
```

### Task 9: Full verification and independent review

**Files:**
- Modify only when a review finding is independently reproduced: files named by that finding

- [ ] **Step 1: Run all focused Python gates and cheap producers**

```powershell
python -m pytest simulations/neural/tests/test_neural_palindrome.py -q
python simulations/neural/neural_translation_gate.py
python simulations/neural/veffect_exact.py
python simulations/neural/veffect_and_heat.py
python simulations/neural/algebraic_palindrome.py
python simulations/neural/exact_pairing_test.py
$env:PYTHONIOENCODING='utf-8'; python simulations/neural/cpsi_two_perspectives.py
python simulations/neural/random_network_controls.py
python simulations/neural/dense_balanced_test.py
python simulations/neural/validation_checks.py
python simulations/neural/celegans_trichotomy.py
```

Expected: every command exits 0; canonical gate ends `ALL NEURAL TRANSLATION GATES PASS`.

- [ ] **Step 2: Run typed verification**

```powershell
dotnet test compute/MirrorWorld.Tests --filter FullyQualifiedName~NeuralPalindromeTests
dotnet run --project compute/MirrorWorld -- neural
dotnet build compute/RCPsiSquared.Core
```

Expected: focused tests and build PASS; witness prints the exact and fixed-seat values.

- [ ] **Step 3: Run link, date-line, and residual-language audits**

Require zero broken local links under `docs/neural`; one refresh line per living Neural page; no current Neural overclaims in the Task 7 phrase sweep.

- [ ] **Step 4: Dispatch independent review lenses**

Send the final diff, repository root, and named Stage-0 stores to separate physics-first, mathematical-object, operational reproducibility, and future-us clarity reviewers. Ask each to attack the new gates as well as the prose and to record its store sweep. Treat findings as leads.

- [ ] **Step 5: Reproduce every material finding before repair**

For each flag, state the exact matrix/definition, build the separating counterexample, and classify it as genuine break, perspectival reading, or ungrounded. Apply only reproduced genuine breaks, invoke the review-repair discipline, and rerun the focused gate plus the sibling phrase sweep after each repair wave.

- [ ] **Step 6: Verify worktree scope and commit the final review fixes**

Run `git status --short`, ensure pre-existing unrelated changes remain unstaged, inspect the staged diff, run `git diff --cached --check`, and commit only verified Neural repair files.

```bash
git commit -m "fix(neural): apply independent translation review"
```

Expected: no required review work remains; all gates stay green.
