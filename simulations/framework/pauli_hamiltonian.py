"""Raw structural classes for Pauli-term Hamiltonians, label-free.

The Klein-4 lens (M, F_a, F_b, C) and Trinity reading (Mother / Father / Child)
both impose interpretive labels on the underlying algebraic structure. This
module exposes the structure DIRECTLY: Pauli-term letter sequences,
bit-parity coordinates, body-count, and aggregate properties — without
labeling.

Three independent Z₂ axes per Pauli term:

  bit_a parity: (#X + #Y) mod 2  — Z⊗N parity break count
  bit_b parity: (#Y + #Z) mod 2  — Π² parity (F-toolkit's bit_b)
  Y-parity:     #Y mod 2          : fixed by Klein within one non-identity-count parity

Y-parity is tied to Klein by y_par = (n + bit_a XOR bit_b) mod 2, where n is the
number of NON-IDENTITY letters (bit_a XOR bit_b = (#X + #Z) mod 2, and
#X + #Y + #Z = n). At fixed n-parity it is therefore determined: it equals
bit_a XOR bit_b for even n, its complement for odd n. Y-parity becomes an
independent third axis only across terms of DIFFERING n-parity, and only such a
mixture populates all eight Z₂³ sectors rather than four. Note `k_body` here is
the string LENGTH, not n, so a two-letter term like IX has n = 1.

Application-layer translations (Trinity, Klein-4, DNA, neural cell-types)
are mappings on top of this raw structure, not part of it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Set, Tuple

from .symmetry import _pauli_tuple_is_truly, _pauli_tuple_pi2_class, klein_index


VALID_LETTERS = ('I', 'X', 'Y', 'Z')


@dataclass(frozen=True)
class PauliTerm:
    """A single Pauli-operator term: ordered letter sequence with a coefficient.

    Frozen (immutable) dataclass. The letter sequence is the canonical form;
    structural properties (klein_index, y_parity, pi2_class, etc.) are
    derived from it without label-collapse.

    Args:
        letters: tuple of Pauli letters from {'I', 'X', 'Y', 'Z'}.
        coefficient: scalar coupling.

    Example:
        >>> t = PauliTerm(letters=('X', 'Y'), coefficient=1.0)
        >>> t.klein_index
        (0, 1)
        >>> t.y_parity
        1
        >>> t.pi2_class
        'pi2_odd'
        >>> t.is_truly
        False
    """
    letters: Tuple[str, ...]
    coefficient: float = 1.0

    def __post_init__(self):
        if not isinstance(self.letters, tuple):
            object.__setattr__(self, 'letters', tuple(self.letters))
        for L in self.letters:
            if L not in VALID_LETTERS:
                raise ValueError(f"Invalid Pauli letter: {L!r}; must be one of {VALID_LETTERS}")

    @property
    def k_body(self) -> int:
        """Body count: number of Pauli letters in the term."""
        return len(self.letters)

    @property
    def n_nonidentity(self) -> int:
        """Number of non-identity letters. This is what C# `PauliTerm.KBody`
        counts; the `k_body` above is the string length. Every parity law on
        this page is stated in terms of THIS count."""
        return sum(1 for L in self.letters if L != 'I')

    @property
    def n_x(self) -> int:
        return sum(1 for L in self.letters if L == 'X')

    @property
    def n_y(self) -> int:
        return sum(1 for L in self.letters if L == 'Y')

    @property
    def n_z(self) -> int:
        return sum(1 for L in self.letters if L == 'Z')

    @property
    def n_i(self) -> int:
        return sum(1 for L in self.letters if L == 'I')

    @property
    def klein_index(self) -> Tuple[int, int]:
        """Klein-Vierergruppe Z₂×Z₂ index: (bit_a, bit_b) parity tuple.

        bit_a = (#X + #Y) mod 2  (Z⊗N parity break)
        bit_b = (#Y + #Z) mod 2  (Π² parity)

        Two terms with the same Klein index are in the same Klein-class.
        """
        return klein_index(self.letters)

    @property
    def y_parity(self) -> int:
        """Y-parity: #Y mod 2. Determined by Klein at fixed parity of the
        non-identity-letter count n: equal to bit_a XOR bit_b for even n, its
        complement for odd n. Independent only across terms of differing
        n-parity."""
        return self.n_y % 2

    @property
    def is_truly(self) -> bool:
        """F77/F85 truly criterion: #Y even AND #Z even (no Π-palindrome violation)."""
        return _pauli_tuple_is_truly(self.letters)

    @property
    def pi2_class(self) -> str:
        """F77/F85 Π²-class: 'truly', 'pi2_odd', or 'pi2_even_nontruly'."""
        return _pauli_tuple_pi2_class(self.letters)

    @property
    def full_z2_signature(self) -> Tuple[int, int, int]:
        """Full Z₂³ structural signature: (bit_a, bit_b, y_parity).

        At fixed parity of the non-identity-letter count this carries the same
        information as klein_index alone, whichever parity it is. Across terms
        of both parities the signature separates 8 sectors instead of 4.
        """
        a, b = self.klein_index
        return (a, b, self.y_parity)

    def __str__(self) -> str:
        letter_str = ''.join(self.letters)
        if self.coefficient == 1.0:
            return letter_str
        return f"{self.coefficient:+g}·{letter_str}"

    def __repr__(self) -> str:
        return f"PauliTerm({self.letters!r}, coefficient={self.coefficient!r})"


@dataclass
class PauliHamiltonian:
    """A Hamiltonian as a sum of Pauli terms on a chain of N sites.

    Exposes raw structural properties without label-collapse. Klein-4 or
    Trinity readings are external application-layer mappings.

    Args:
        terms: list of PauliTerm instances.
        chain_length: number of sites (N).

    Example:
        >>> H = PauliHamiltonian.from_letter_tuples([('X', 'Y'), ('Y', 'X')], chain_length=3)
        >>> H.klein_set
        {(0, 1)}
        >>> H.is_klein_homogeneous
        True
        >>> H.per_term_klein_indices
        [(0, 1), (0, 1)]
    """
    terms: List[PauliTerm]
    chain_length: int

    def __post_init__(self):
        if self.chain_length < 1:
            raise ValueError(f"chain_length must be >= 1; got {self.chain_length}")
        for t in self.terms:
            if t.k_body > self.chain_length:
                raise ValueError(
                    f"Term {t} has body count {t.k_body} exceeding chain length {self.chain_length}"
                )

    @property
    def klein_set(self) -> Set[Tuple[int, int]]:
        """Set of distinct Klein indices across all terms."""
        return {t.klein_index for t in self.terms}

    @property
    def is_klein_homogeneous(self) -> bool:
        """True if all terms share the same Klein index.

        Klein-homogeneity does NOT by itself imply F87 soft-or-truly. The rule
        holds only for IDENTITY-FREE two-letter pairs: of the 21 Klein-homogeneous
        two-letter pairs at N=3 that carry any content (the 3 pairing with 'II'
        are a multiple of the identity and cannot move the spectrum), the 6
        identity-free ones are 0/6 hard, while the 15 carrying an 'I' are 5/15
        hard, N-stably. Counterexample: [('I','Z'),
        ('Z','I')], both Klein (0,1), classifies hard at N=3 and N=4.

        It fails at k≥3 too. The 294-pair Z₂³-homogeneous sweep at k=3 N=4
        (PROOF_F103_F87_Z2_CUBED_REFINEMENT §3.2) finds 50 hard in the Klein-(0,1)
        cell under Z-dephasing; AMONG Klein-homogeneous pairs, sitting in the dephase
        letter's own cell is necessary for hardness but not sufficient: XY+YX sits
        in Z's cell (0,1) and is soft at N=3, 4 and 5. Drop the homogeneity and
        even necessity goes: 10 hard two-letter pairs at N=3 (IX+XZ, IY+YZ, ...)
        have neither term in Z's cell. Counterexample: [('I','I','Z'), ('I','X','Y')].

        So Klein-inhomogeneity is necessary for hardness only among identity-free
        k=2 pairs; elsewhere the Klein index does not settle the class.
        """
        return len(self.klein_set) <= 1

    @property
    def y_parity_set(self) -> Set[int]:
        """Set of distinct Y-parities across terms.

        Redundant with the Klein index within any fixed parity of the
        non-identity-letter count; informative only across a mixture of both.
        """
        return {t.y_parity for t in self.terms}

    @property
    def is_y_parity_homogeneous(self) -> bool:
        """True if all terms share the same Y-parity."""
        return len(self.y_parity_set) <= 1

    @property
    def bit_b_set(self) -> Set[int]:
        """Set of distinct bit_b parities across terms.

        bit_b = (#Y + #Z) mod 2 per term. This is the Π² axis (F38: Π² acts
        on a Pauli string as (-1)^bit_b). F112's typed Tier1Derived scope
        requires bit_b-homogeneous collapse operators c.
        """
        return {t.klein_index[1] for t in self.terms}

    @property
    def is_bit_b_homogeneous(self) -> bool:
        """True if all terms share the same bit_b parity (F112 precondition).

        F112 (Lindblad Π-eigenvalue balance): for Hermitian H + each c_k
        bit_b-homogeneous, polarity_coordinates_from_L asymmetry is 0
        bit-exact. Single-Pauli c is trivially bit_b-homogeneous (one term,
        one bit_b value).
        """
        return len(self.bit_b_set) <= 1

    @property
    def full_z2_signature_set(self) -> Set[Tuple[int, int, int]]:
        """Set of full Z₂³ signatures (bit_a, bit_b, y_parity) across terms.

        Within one parity of the non-identity count this has the same number of
        elements as klein_set. It is strictly finer when two terms SHARE a Klein
        cell while differing in that parity, which two-letter terms can already
        do: [('X','Y'), ('Z','I')] has klein_set {(0,1)} but signatures
        {(0,1,1), (0,1,0)}.
        """
        return {t.full_z2_signature for t in self.terms}

    @property
    def is_z2_homogeneous(self) -> bool:
        """True if all terms share the same full Z₂³ signature.

        Strictly stronger than klein_homogeneous whenever the terms differ in the
        parity of their non-identity-letter count, which two-letter terms can
        already do.
        """
        return len(self.full_z2_signature_set) <= 1

    @property
    def per_term_klein_indices(self) -> List[Tuple[int, int]]:
        """Klein index of each term, in order."""
        return [t.klein_index for t in self.terms]

    @property
    def per_term_pi2_classes(self) -> List[str]:
        """Π²-class of each term, in order."""
        return [t.pi2_class for t in self.terms]

    @property
    def per_term_y_parities(self) -> List[int]:
        """Y-parity of each term, in order."""
        return [t.y_parity for t in self.terms]

    @property
    def has_truly_term(self) -> bool:
        """Whether at least one term is truly (M = 0 by Master Lemma)."""
        return any(t.is_truly for t in self.terms)

    @property
    def k_body_set(self) -> Set[int]:
        """Set of body counts across terms."""
        return {t.k_body for t in self.terms}

    @property
    def is_uniform_body(self) -> bool:
        """True if all terms have the same body count (no mixed-body Hamiltonian)."""
        return len(self.k_body_set) <= 1

    @property
    def n_terms(self) -> int:
        return len(self.terms)

    @classmethod
    def from_letter_tuples(cls, letter_tuples, chain_length, J=1.0):
        """Convenience constructor from a list of letter tuples.

        Example:
            H = PauliHamiltonian.from_letter_tuples(
                [('X', 'Y'), ('Y', 'X')], chain_length=3
            )
        """
        terms = [PauliTerm(letters=tuple(t), coefficient=J) for t in letter_tuples]
        return cls(terms=terms, chain_length=chain_length)

    def __str__(self) -> str:
        if not self.terms:
            return "0"
        return " + ".join(str(t) for t in self.terms)
