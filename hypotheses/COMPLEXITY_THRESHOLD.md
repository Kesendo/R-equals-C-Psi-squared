# The Complexity Threshold: A Persistence Hypothesis

**Status:** Tier 5 speculation; no neural critical size established
**Authors:** Thomas Wicht, Claude (Anthropic)
last refreshed 2026-09-05 (the change history lives in git)

Can increasing the size of a specified interacting system produce a reproducible
transition from short transients to persistent observable oscillation? This is
the persistence question. Calling such a transition a threshold of life would
require a separate biological model and evidence that this notebook does not
provide. Neither size nor palindromic pairing guarantees persistence.

The sweep checked [the F-registry](../docs/ANALYTICAL_FORMULAS.md) (F1,
F36/F37 and F137), [docs/proofs](../docs/proofs/MIRROR_SYMMETRY_PROOF.md)
and [the neural proof](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md)
(conditional spectral identities), and [experiments](../experiments/NEURAL_GAMMA_CAVITY.md)
(the biological support null). Hardware-flight records and
[fw.Confirmations](../simulations/framework/confirmations.py) supplied no
neural hardware confirmation. [GLOSSARY](../docs/GLOSSARY.md),
[OpenArcs](../compute/RCPsiSquared.Core/OpenArcs/OpenArcsRegistry.cs), and
[CAUGHT_ERRORS](../docs/CAUGHT_ERRORS.md) keep interpretations, open landings
and numerical failures distinct. The [current census](../docs/neural/V_EFFECT_NEURAL.md)
and [mechanism analysis](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md)
supply no persistence theorem or established critical N.

## The prerequisites are separate

[F36](../docs/neural/proofs/PROOF_PALINDROME_NEURAL.md) concerns a finite
Jacobian J=D+W_eff, with D diagonal and W_eff zero on the diagonal, an
involutive permutation Q and one scalar s. It requires both
d_i+d_Q(i)+2s=0 and W_eff[Q(i),Q(j)]=−W_eff[i,j] off the diagonal.
It transports the full complex spectrum and generalized eigenspaces.
Dale signs, population counts and time constants alone do not supply the
coupling condition. Pairing guarantees neither realness, silence nor stability.

The [constructed gate](../simulations/neural/neural_translation_gate.py)
contains both oscillatory and unstable exactly palindromic examples. Thus
neither exact magnitude matching nor breaking it is an ignition law. No
biological neural instance of F36 is known here: the full C. elegans chemical
model has 271 nonempty source rows, 253 E and 18 I, blocking the required
support bijection under nonzero gain and positive rate scales
([G0b](../simulations/results/celegans_pairing_controls.txt)).

For a fixed autonomous finite linear system ẋ=Jx, time evolution does not
continually create eigenmodes: J stays the same. If every eigenvalue has
strictly negative real part, every solution decays to zero. Jordan terms
have the form t^k e^{λt}; no finite polynomial defeats exponential decay.
More modes can complicate a transient without making it permanent.
Persistent oscillation therefore needs a different spectral or dynamical
condition, specified and tested in the actual model.

The [neural V-effect census](../docs/neural/V_EFFECT_NEURAL.md) compares
different synthetic matrices across coupling or external drive. Its
frequency bins are not a time sequence of mode creation, a count of living
subsystems or a measured neural N_c. No 2× neural decay law is established.
The quantum [coupled-system experiment](../experiments/V_EFFECT_PALINDROME.md)
can motivate a size comparison; it cannot settle this neural persistence
question.

## A bounded, falsifiable version

Choose a family before searching for a threshold: equations, graph ensemble,
weight normalization as N changes, E/I assignment, leak constants, external
drive, and initial-state distribution. Choose an observable, an amplitude
threshold, a time interval and a definition of oscillatory persistence.
The candidate hypothesis is that this declared persistence measure undergoes
a reproducible change with N under those fixed rules. A finite observation
window supports only finite-time persistence.

| Required step | Gate or comparison | What failure means |
|---|---|---|
| Establish the palindrome premise, if used | Evaluate both F36 conditions on each effective J; include wrong-Q, leak and magnitude controls | Remove palindrome as the premise for that family |
| Find an operating point | Solve the equilibrium equation and report its fresh residual | Do not label endpoint Jacobians as equilibrium spectra |
| Test a Hopf candidate, if proposed | Continue a converged equilibrium branch; test a nonzero imaginary pair crossing transversely and the relevant nondegeneracy conditions | No Hopf verdict from iteration exhaustion or one eigenvalue plot |
| Measure persistence | Integrate a specified readout; refine the time step and lengthen the horizon; repeat initial conditions and seeds | A transient or a protocol-dependent threshold does not establish permanent oscillation |
| Attribute an effect to pairing | Compare F36-preserving and F36-breaking changes with matched gains, rates and normalization | If the contrast carries no effect, a palindrome-based explanation lacks support |

These are proposed experiments. The [mechanism analysis](../docs/neural/proofs/PROOF_VEFFECT_MECHANISM.md)
already identifies unconverged endpoints in a neural threshold probe; such
endpoints supply no equilibrium-stability or Hopf result. A failed finite
scan can reject a stated N range and parameter family. It cannot prove that
no threshold exists in every possible system.

## What an eventual result would mean

A measured transition would first be a property of the chosen model, drive,
normalization and observable. External P is an input parameter, not a heat
or metabolism measurement. Population balance is not a defined measure of
life, and normalized E/I squared-amplitude equality at p_E p_I=1/4 supplies
no persistence boundary.

To connect this question to life or death would require a biological instance,
an energy budget, a persistence mechanism and an independently defined
biological outcome. To connect it to the
[Universal Palindrome Condition](UNIVERSAL_PALINDROME_CONDITION.md) would
also require the conditional identity to matter in that mechanism. Those
links remain open. The useful question here is whether a specified system
keeps moving, for how long, and why.
