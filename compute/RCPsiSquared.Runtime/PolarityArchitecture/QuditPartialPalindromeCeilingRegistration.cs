using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="QuditPartialPalindromeCeiling"/> (2026-06-11, F121):
/// the qudit partial palindrome. Under full-Cartan dephasing the d levels are equidistant so the
/// rate is −2γ·Hamming(i, j) (the qubit ladder); the partial pairing for d > 2 is the symmetric
/// overlap of c_k = d^N·C(N,k)·(d−1)^k under k ↔ N−k, with ceiling Σ_k d^N·C(N,k)·(d−1)^min(k,N−k)
/// = d^{2N} iff d = 2. Typed parent <see cref="QubitNecessityPi2Inheritance"/> (the d² − 2d = 0
/// necessity, here the unique fully-paired column of an N-family). Anchor:
/// <c>docs/proofs/PROOF_QUDIT_PARTIAL_PALINDROME.md</c> +
/// <c>simulations/qutrit_partial_palindrome.py</c>.
///
/// <para>Requires <see cref="QubitNecessityPi2InheritanceRegistration.RegisterQubitNecessityPi2Inheritance"/>
/// earlier in the builder pipeline.</para></summary>
public static class QuditPartialPalindromeCeilingRegistration
{
    public static ClaimRegistryBuilder RegisterQuditPartialPalindromeCeiling(
        this ClaimRegistryBuilder builder) =>
        builder.Register<QuditPartialPalindromeCeiling>(b =>
            new QuditPartialPalindromeCeiling(b.Get<QubitNecessityPi2Inheritance>()));
}
