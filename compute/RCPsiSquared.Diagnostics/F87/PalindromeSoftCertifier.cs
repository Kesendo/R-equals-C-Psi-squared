using System;
using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.Pauli;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>A Liouvillian-free, one-sided SOFT-certifier (PROOF_F103 §7.12). The true soft criterion is
/// the bipartiteness of the BASIS-STATE hopping graph (a 2^N, letter-dependent object); this class tries
/// scalable SUFFICIENT conditions ("soft colourings") and certifies soft if any holds. It never claims
/// hard: NotCertified means no scalable strategy applies (the chain-scope hard proxy stays in
/// <see cref="PalindromeMaskClassifier"/>). A certificate is correct for any N and any topology.</summary>
public static class PalindromeSoftCertifier
{
    /// <summary>Which scalable soft-colouring certified the Hamiltonian (None = not certified).</summary>
    public enum SoftStrategy { None, LinearSiteColoring, ExcitationPairing }

    /// <summary>Result of <see cref="Certify"/>: whether soft is certified, and by which strategy.</summary>
    public readonly record struct SoftCertificate(bool Certified, SoftStrategy Strategy);

    /// <summary>True iff the summed Hamiltonian is a pure pairing (every basis-edge Δn = ±2), detected
    /// by a σ± decomposition: the mixed (hopping) pieces must cancel. N-independent.</summary>
    public static bool IsPurePairing(IReadOnlyList<PauliTerm> terms) => throw new NotImplementedException();

    /// <summary>The excitation-number strategy: certify soft iff the Hamiltonian is a pure pairing
    /// (then ⌊n/2⌋ mod 2 two-colours the basis-state graph, soft on any topology).</summary>
    public static bool CertifyByExcitationPairing(IReadOnlyList<PauliTerm> terms) => throw new NotImplementedException();

    /// <summary>The linear site-colouring strategy: certify soft iff the chain flip-mask set is bipartite
    /// (the chiral K). Reuses <see cref="PalindromeMaskClassifier"/>.</summary>
    public static bool CertifyByLinearSiteColoring(IReadOnlyList<PauliTerm> terms, int n) => throw new NotImplementedException();

    /// <summary>Try the stronger excitation strategy first, then the linear one; return the certificate.</summary>
    public static SoftCertificate Certify(IReadOnlyList<PauliTerm> terms, int n) => throw new NotImplementedException();
}
