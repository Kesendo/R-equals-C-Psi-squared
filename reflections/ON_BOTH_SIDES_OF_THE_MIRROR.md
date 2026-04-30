# On Both Sides of the Mirror

**Status:** Reflection. After F80 entered the framework as a callable primitive, while reading what its 2i factor actually says geometrically. Revised 2026-04-30 after a third-party algebra check found that the through-line claim holds for the spectrum, not for the operator: M and Π·M·Π⁻¹ are Frobenius-orthogonal partners for the Π²-odd cases F80 covers, sharing only their eigenvalues.
**Date:** 2026-04-30
**Authors:** Thomas Wicht, Claude (Opus 4.7)
**Context:** Tom asked a question that the algebra had already been answering for months without anyone hearing it. If the mirror projects everything 90° onto itself so that it does not forget, must not what enters on one side leave on the other? From our view minus zero, from theirs plus zero.

---

The palindromic equation Π·L·Π⁻¹ + L + 2Σγ·I = M has been read for months as a single equation with a single sign convention. The Π conjugates L; the result plus L plus the dissipation shift gives the residual. M is what the mirror cannot reflect away. That reading is correct, but it conceals a structure that becomes visible the moment one steps to the other side and reads the same equation backward.

Conjugate M itself by Π. The first term Π·L·Π⁻¹ becomes Π²·L·Π⁻². For Π²-even L (which includes truly H = XX+YY+ZZ and Π²-even non-truly H = YZ+ZY), Π² commutes with L and Π²·L·Π⁻² = L. The full chain then collapses: Π·M·Π⁻¹ = L + Π·L·Π⁻¹ + 2Σγ·I = M. M is its own Π-conjugate; the residual is literally the through-line operator that both sides of the mirror read identically.

For Π²-odd L (the case F80 covers, where H is XY+YX or any Π²-odd non-truly chain bilinear), Π²·L·Π⁻² ≠ L. Direct computation gives Π·M·Π⁻¹ - M = Π²·L·Π⁻² - L; numerically, ‖Π·M·Π⁻¹ - M‖ = √2·‖M‖, which is the Frobenius distance between two equal-norm operators that are perpendicular in operator space. M and Π·M·Π⁻¹ are different matrices. They are not the same operator read from two sides; they are two operators of equal magnitude and identical spectrum sitting at 90° to each other in the operator sphere.

The eigenvalues, however, are common to both. Spec(Π·M·Π⁻¹) = Spec(M) by unitary invariance of the spectrum. What both readings share is the list of allowed values; what differs is which eigenvectors carry which values. The through-line is the spectrum, not the matrix.

[F80](../docs/ANALYTICAL_FORMULAS.md) says what that spectrum is:

```
Spec(M)_{nontrivial} = 2i · Spec(H_non-truly)
```

The factor i is the 90° rotation in the complex plane. H's eigenvalues are real numbers; they are the energies our side measures. M's eigenvalues are purely imaginary; they are the frequencies and decay rates that govern time evolution. The mirror is the rotation that maps energy into time. From one side an eigenvalue is a number on the real axis. The 90° turn places the same eigenvalue on the imaginary axis. Same magnitude, different axis. Reading-direction is what changes.

The factor 2 comes from H's particle-hole pair structure: the JW-mapped Majorana bilinear is chirally symmetric, so for every +λ there is a partner -λ in Spec(H). The 2i applied to the pair gives ±2iλ. What looks like "doubling because there are two copies of L in M" is shorthand for the deeper fact that H itself is already mirrored within its own spectrum, and the i factor lifts that mirroring into M's imaginary spectrum.

Truly Hamiltonians have M = 0 from both sides. Our side reads zero as "no negative residue." The mirror's side reads zero as "no positive residue." The two zeros are the same point in operator space, labeled by two opposite conventions. Minus zero from us, plus zero from them; both correct, both the absence of the same thing. Here Π·M·Π⁻¹ = M holds trivially because both are zero.

For non-truly Hamiltonians the spectrum opens. Chain N=3 with H = J(XY+YX) has Spec(M) = {-5.66i, 0, +5.66i} with multiplicities 16, 32, 16, summing to 4³ = 64. The +5.66i mode is what we call "outward"; the -5.66i is what the mirror calls "outward." They are paired by H's particle-hole symmetry, lifted by the 2i factor onto opposite sides of the imaginary axis. The middle zero, with multiplicity 32, is not a "truly subspace inside soft" (XY+YX has no algebraic truly component). It is the free-fermion zero-mode degeneracy of the JW spectrum: at N=3, mode k=2 has dispersion ε(2) = 2cos(π·2/4) = 0, plus particle-hole-paired fillings of modes 1 and 3 that cancel. These zeros are not the mirror still closing perfectly; they are a property of the chiral H spectrum that survives the 2i lift unchanged.

This is the geometric content of [Zero Is the Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md) carried one level higher and corrected. At γ = 0 with truly H, Π·L·Π⁻¹ = -L exactly: same operator, opposite sign, mirror closes. At γ = 0 with non-truly H, the closure breaks; the residual M ≠ 0 captures the breaking. For Π²-even non-truly H (like YZ+ZY) the residual is its own Π-conjugate; the defect is mirror-symmetric. For Π²-odd non-truly H (the F80 cases) the residual and its Π-conjugate are distinct orthogonal partners with shared spectrum. The mirror does not vanish when symmetry breaks; it relocates. Sometimes onto the operator (Π²-even cases). Sometimes onto the spectrum alone (Π²-odd cases). Always onto something.

What the framework remembers, it remembers by sharing eigenvalues across the mirror. The 2i factor is the channel: 90° on one side, 90° on the other, summing into an imaginary spectrum on which both readings agree. Reading the framework is reading from one of the two sides. Reading it from the other gives the same eigenvalue list, possibly differently distributed across operator-space directions. There is no canonical orientation; only the choice of which side to call ours.

---

*"Müsste dann nicht das was rein kommt auf der andere seite wieder rauskommen?" Tom, 2026-04-30.*
*"Spec(Π·M·Π⁻¹) = Spec(M) always; Π·M·Π⁻¹ = M only for Π²-even L." The corrected through-line statement.*
*"Spec(M) = 2i · Spec(H_non-truly)." F80, the 90° rotation that maps energy into time.*
*"Two mirrors, perfectly aligned, zero distance apart." From [Zero Is The Mirror](../hypotheses/ZERO_IS_THE_MIRROR.md), March 2026.*
