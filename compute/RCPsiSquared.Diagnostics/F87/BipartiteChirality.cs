using System;
using System.Collections.Generic;
using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Pauli;
using CMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.F87;

/// <summary>The F103 §7 bipartite-chirality criterion for the F87 diagonal cell.
///
/// <para>A diagonal-cell pair's Hamiltonian H, read in the dephasing letter's eigenbasis
/// (computational basis for Z; rotated by the single-qubit gate that sends the letter to Z
/// for X and Y), is a real hopping matrix. Its graph G_H is bipartite iff a chiral sign
/// operator K = diag(±1) with K·H·K = −H exists (the 2-colouring), and K, being diagonal in
/// the dephasing basis, commutes with the dephasing dissipator. That second mirror, composed
/// with Π, restores the palindrome Spec(L) = Spec(−L − 2σ): the pair is <b>soft</b>. A
/// non-bipartite G_H (odd hopping cycle, or a lifted diagonal from a pure-D template) admits
/// no such K and the pair is <b>hard</b>.</para>
///
/// <para><b>bipartite ⟹ soft is derived</b> (the chiral K construction, see
/// <see cref="Core.Symmetry.ChiralKClaim"/>); the converse <b>non-bipartite ⟹ hard</b> is
/// now derived (PROOF_F103 §7.5, 2026-06-04) modulo the first-order-block premise (the K3 triangle
/// obstructs the chiral functional supplying the gain channel's −N reflection mode; any
/// palindromizer forces a spectral palindrome); that premise is itself closed (§7.6, degenerate PT +
/// analyticity), so the converse is fully derived modulo standard perturbation theory. Tier1Derived-eligible,
/// kept Tier1Candidate pending the formal promotion pass.
/// This primitive returns both the criterion's verdict (<see cref="BipartiteChiralityResult.PredictedClass"/>)
/// and the actual F87 verdict (<see cref="PauliPairTrichotomy.Classify"/>) so a witness can
/// check that they agree.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F103_F87_Z2_CUBED_REFINEMENT.md</c> §7 +
/// <c>simulations/f87_42_8_bipartite_fullcell.py</c> +
/// <c>simulations/f87_bipartite_chiral_witness.py</c>.</para></summary>
public static class BipartiteChirality
{
    /// <summary>Classify a k-body diagonal-cell pair by the bipartite criterion and compare to
    /// the actual F87 trichotomy. <paramref name="termTemplates"/> are k-body Pauli templates
    /// (placed across the chain by the sliding-window builder, as in
    /// <see cref="PauliPairTrichotomy.Classify(ChainSystem, IReadOnlyList{PauliTerm}, double, double, PauliLetter)"/>).</summary>
    public static BipartiteChiralityResult Classify(
        ChainSystem chain, IReadOnlyList<PauliTerm> termTemplates,
        PauliLetter dephaseLetter = PauliLetter.Z, double tol = 1e-9)
    {
        if (chain is null) throw new ArgumentNullException(nameof(chain));
        if (termTemplates is null || termTemplates.Count == 0)
            throw new ArgumentException("need at least one term template", nameof(termTemplates));

        var hRotated = RotateToDephaseBasis(termTemplates.ChainKBody(chain.N), chain.N, dephaseLetter);
        var (isBipartite, chiralK) = BipartiteChiralK(hRotated, tol);

        bool kVerified = false;
        if (chiralK is not null)
            kVerified = (chiralK * hRotated * chiralK + hRotated).FrobeniusNorm() < 1e-7;

        var predicted = isBipartite ? TrichotomyClass.Soft : TrichotomyClass.Hard;
        var actual = PauliPairTrichotomy.Classify(chain, termTemplates, dephaseLetter: dephaseLetter);
        return new BipartiteChiralityResult(isBipartite, kVerified, predicted, actual);
    }

    /// <summary>Rotate H by U^⊗N where U sends the dephasing Pauli letter to Z (Hadamard for
    /// X, the Y→Z rotation for Y, identity for Z), so the chiral K becomes diagonal in the
    /// computational basis there.</summary>
    private static CMatrix RotateToDephaseBasis(CMatrix h, int n, PauliLetter dephaseLetter)
    {
        if (dephaseLetter == PauliLetter.Z) return h;
        var u1 = SingleSiteRotation(dephaseLetter);
        var u = u1;
        for (int s = 1; s < n; s++) u = u.KroneckerProduct(u1);
        return u * h * u.ConjugateTranspose();
    }

    private static CMatrix SingleSiteRotation(PauliLetter dephaseLetter)
    {
        double r = 1.0 / Math.Sqrt(2.0);
        Complex i = Complex.ImaginaryOne;
        return dephaseLetter switch
        {
            // Hadamard: H·X·H = Z
            PauliLetter.X => CMatrix.Build.DenseOfArray(new Complex[,] { { r, r }, { r, -r } }),
            // U·Y·U† = Z
            PauliLetter.Y => CMatrix.Build.DenseOfArray(new Complex[,] { { r, -r * i }, { r, r * i } }),
            _ => CMatrix.Build.DenseIdentity(2),
        };
    }

    /// <summary>Two-colour H's hopping graph (basis states as nodes, nonzero off-diagonal
    /// entries as edges). Returns (false, null) if the diagonal is nonzero (a lifted diagonal
    /// admits no K with K·H·K = −H) or an odd cycle is found; otherwise (true, K) with
    /// K = diag(±1) the 2-colouring.</summary>
    private static (bool isBipartite, CMatrix? chiralK) BipartiteChiralK(CMatrix h, double tol)
    {
        int d = h.RowCount;
        for (int a = 0; a < d; a++)
            if (h[a, a].Magnitude > tol) return (false, null);

        var color = new int[d];
        for (int a = 0; a < d; a++) color[a] = -1;
        var queue = new Queue<int>();
        for (int start = 0; start < d; start++)
        {
            if (color[start] != -1) continue;
            color[start] = 0;
            queue.Enqueue(start);
            while (queue.Count > 0)
            {
                int u = queue.Dequeue();
                for (int v = 0; v < d; v++)
                {
                    if (v == u || h[u, v].Magnitude <= tol) continue;
                    if (color[v] == -1) { color[v] = color[u] ^ 1; queue.Enqueue(v); }
                    else if (color[v] == color[u]) return (false, null);
                }
            }
        }

        var k = CMatrix.Build.Dense(d, d);
        for (int a = 0; a < d; a++) k[a, a] = color[a] == 0 ? Complex.One : -Complex.One;
        return (true, k);
    }
}

/// <summary>Result of <see cref="BipartiteChirality.Classify"/>: the bipartite criterion's
/// prediction next to the actual F87 trichotomy verdict.</summary>
public sealed class BipartiteChiralityResult
{
    /// <summary>Whether H's hopping graph (in the dephasing letter's eigenbasis) is bipartite,
    /// i.e. admits a chiral K = diag(±1) with K·H·K = −H.</summary>
    public bool IsBipartite { get; }

    /// <summary>For a bipartite pair, whether the extracted 2-colouring K satisfies K·H·K = −H
    /// to 1e-7 (a sanity check on the chiral construction). False when not bipartite.</summary>
    public bool ChiralKVerified { get; }

    /// <summary>The bipartite criterion's verdict: Soft if bipartite, else Hard. (The criterion
    /// distinguishes soft from hard among non-truly diagonal-cell pairs; truly is M=0, a
    /// separate bucket the criterion does not address.)</summary>
    public TrichotomyClass PredictedClass { get; }

    /// <summary>The actual F87 verdict from <see cref="PauliPairTrichotomy.Classify"/>.</summary>
    public TrichotomyClass ActualClass { get; }

    /// <summary>True iff the bipartite criterion's verdict equals the actual F87 verdict.</summary>
    public bool Agrees => PredictedClass == ActualClass;

    public BipartiteChiralityResult(
        bool isBipartite, bool chiralKVerified,
        TrichotomyClass predictedClass, TrichotomyClass actualClass)
    {
        IsBipartite = isBipartite;
        ChiralKVerified = chiralKVerified;
        PredictedClass = predictedClass;
        ActualClass = actualClass;
    }
}
