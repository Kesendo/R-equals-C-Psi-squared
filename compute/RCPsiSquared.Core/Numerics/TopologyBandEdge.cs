using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;

namespace RCPsiSquared.Core.Numerics;

/// <summary>The XY single-excitation band edge of a hopping graph: the largest |E_k| (frequency) of the
/// vac↔single-excitation coherences. The single-excitation block of H = (J/2)·Σ_bond(XX+YY) is exactly
/// J × the graph adjacency matrix, so the band edge = J × the adjacency spectral radius ρ. Closed forms:
/// chain (path P_N) 2cos(π/(N+1)), star (K_{1,N−1}) √(N−1), ring (cycle C_N) 2.
///
/// <para>This is the Im/frequency (L_H) side only. The RATE side — that this coherence sits at Re=−2γ for
/// any topology — is the Absorption Theorem (n_XY=1 ⟹ Re=−2γ, PROOF_ABSORPTION_THEOREM.md §4.3/§4.5), not
/// recomputed here. XY-specific: under Heisenberg the ZZ term adds a diagonal and the band edge ≠ J·ρ.</para></summary>
public static class TopologyBandEdge
{
    /// <summary>The adjacency spectral radius ρ = max|eigenvalue| of the N×N unit adjacency matrix built
    /// from the bonds (bond couplings are ignored — the J factor is applied in <see cref="BandEdge"/>).
    /// The adjacency is real symmetric, so its eigenvalues are real.</summary>
    public static double SpectralRadius(int n, IReadOnlyList<Bond> bonds)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "N must be >= 1.");
        var a = Matrix<double>.Build.Dense(n, n);
        foreach (var b in bonds)
        {
            a[b.Site1, b.Site2] = 1.0;
            a[b.Site2, b.Site1] = 1.0;
        }
        double rho = 0.0;
        foreach (var ev in a.Evd().EigenValues)
            rho = Math.Max(rho, Math.Abs(ev.Real));
        return rho;
    }

    /// <summary>The band edge J·ρ — the largest single-excitation frequency of the XY hopping graph.</summary>
    public static double BandEdge(int n, IReadOnlyList<Bond> bonds, double j) => j * SpectralRadius(n, bonds);
}
