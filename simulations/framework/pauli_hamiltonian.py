"""Raw structural classes for Pauli-term Hamiltonians, label-free.

The Klein-4 lens (M, F_a, F_b, C) and Trinity reading (Mother / Father / Child)
both impose interpretive labels on the underlying algebraic structure. This
module exposes the structure DIRECTLY: Pauli-term letter sequences,
bit-parity coordinates, body-count, and aggregate properties — without
labeling.

Three independent Z₂ axes per Pauli term:

  bit_a parity: (#X + #Y) mod 2  — Z⊗N parity break count
  bit_b parity: (#Y + #Z) mod 2  — Π² parity (F-toolkit's bit_b)
  Y-parity:     #Y mod 2          — independent at k ≥ 3

At k = 2 (bilinears), Y-parity = (bit_a XOR bit_b), redundant with Klein.
At k ≥ 3, Y-parity carries independent information; the full per-letter
symmetry alphabet is Z₂³ (8 sectors), not Z₂² (4 Klein slots).

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
        """Y-parity: #Y mod 2. At k≥3 this is an independent third Z₂ axis;
        at k=2 it equals bit_a XOR bit_b (redundant with Klein index)."""
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

        At k=2 this carries the same information as klein_index alone.
        At k≥3, y_parity becomes an independent axis and the full signature
        labels 8 sectors instead of 4.
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

        Empirical structural fact (verified at k=2 full enumeration; k=3 sample):
        Klein-homogeneous Hamiltonians are always F77 soft or F77 truly,
        never F77 hard. Klein-inhomogeneity is necessary (not sufficient)
        for F77 hardness.
        """
        return len(self.klein_set) <= 1

    @property
    def y_parity_set(self) -> Set[int]:
        """Set of distinct Y-parities across terms.

        At k≥3 this is an independent classification axis; at k=2 it is
        determined by Klein index (redundant).
        """
        return {t.y_parity for t in self.terms}

    @property
    def is_y_parity_homogeneous(self) -> bool:
        """True if all terms share the same Y-parity."""
        return len(self.y_parity_set) <= 1

    @property
    def full_z2_signature_set(self) -> Set[Tuple[int, int, int]]:
        """Set of full Z₂³ signatures (bit_a, bit_b, y_parity) across terms.

        For k=2 Hamiltonians this has the same number of elements as klein_set
        (Y-parity redundant). For k≥3 Hamiltonians, may be finer.
        """
        return {t.full_z2_signature for t in self.terms}

    @property
    def is_z2_homogeneous(self) -> bool:
        """True if all terms share the same full Z₂³ signature.

        Stronger than klein_homogeneous at k≥3.
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
