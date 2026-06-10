using MathNet.Numerics.LinearAlgebra.Factorization;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Ptf;

/// <summary>The slow subspace's per-site light distribution as a computed, basis-free object:
/// the s*-boundary surface made measurable.
///
/// <para><b>The question this answers</b>: <see cref="Foundation.PostEpFlowField"/> decides
/// sterile zone vs birth canal by a two-point rate probe (rate at Q=1000 minus rate at Q=1.5).
/// The mechanism beneath that probe is the Absorption Theorem
/// (<see cref="RCPsiSquared.Core.Symmetry.AbsorptionTheoremClaim"/>,
/// <c>docs/proofs/PROOF_ABSORPTION_THEOREM.md</c>): the slowest rate is
/// 2·Σ_l γ_l·⟨light at site l⟩ of the slow mode, exactly. So the zone question IS whether the
/// slow subspace's per-site light distribution is Q-frozen: sterile ⟺ this distribution is
/// Q-invariant; birth canal ⟺ the Hamiltonian redistributes it with Q and the rate follows.
/// The boundary in γ-profile space is a SURFACE, not a number: s* crossings range 0.11–0.77
/// across endpoint pairs (<c>simulations/birth_canal_boundary_pathdependence.py</c>, 2026-06-04;
/// the once-celebrated 0.709 is one hand-picked line's crossing of it). This primitive computes
/// the distribution itself, so the surface can be probed pointwise.</para>
///
/// <para><b>The object</b>: take the slowest non-kernel rate cluster (all eigenvalues with
/// Re λ within <paramref name="clusterTolerance"/> of the slowest non-kernel Re λ; they share
/// one rate but generally differ in frequency). Its biorthogonal spectral projector
/// P = Σ_k |M_k⟩⟨W_k| (left covectors = rows of R⁻¹, the construction of
/// <c>F8PartnerLightComplementarityTests</c>) is basis-free: it does not depend on the choice
/// of eigenbasis within the cluster, by the same spectral argument
/// <see cref="StationaryManifold"/> makes for the kernel count. The per-site light is read
/// through the ORTHOGONAL projector Π_V onto range(P), the slow invariant subspace:
///
/// <code>
///   light_l = Tr(Π_V · Δ_l) / g,   Δ_l(x) = 1 iff bit l differs between bra and ket
///                                  of the coherence index x = a·d + b (row-major vec),
///   g = Degeneracy = dim range(P).
/// </code></para>
///
/// <para><b>Why the orthogonal projector and not the raw diagonal of P</b> (pinned in Python at
/// N=5, 2026-06-10): the diagonal of the biorthogonal P fails the absorption cross-check in the
/// birth canal (residual 1.7·10⁻¹ on the flat-bulk-edge profile) and can go negative; the
/// per-eigenvector |v|² average (the former <c>ReadAssembly</c> carrier) is absorption-exact but
/// basis-dependent under exact degeneracy, the gauge caveat this primitive retires. The Π_V trace
/// is basis-free (Π_V is the unique orthogonal projector onto the L-invariant subspace range(P)),
/// stays in [0, 1] per site (Π_V is PSD with a {0, 1} spectrum), and satisfies the absorption
/// cross-check to machine precision: in an orthonormal basis adapted to the invariant subspace V,
/// L is block upper-triangular (L·V ⊆ V), so Tr(Π_V·L) = Σ_cluster λ_k; the Hamiltonian part is
/// anti-Hermitian and the dissipator is the real diagonal −2·Σ_l γ_l·Δ_l, hence
/// Re Tr(Π_V·L) = Tr(Π_V·D) gives rate = 2·Σ_l γ_l·light_l exactly. Where Degeneracy = 1 the
/// projector value coincides with the per-eigenvector Rayleigh light, so nothing the old reading
/// got right is lost.</para>
///
/// <para><b>Pinned behavior (N=5 dimensionless XY chain, γ-profiles summing to N)</b>: uniform
/// [1,1,1,1,1] is sterile, distribution exactly [0.2]⁵ at Q=1.5 and Q=1000 (drift ~10⁻¹⁶);
/// peaked-V [0.25,0.75,3.0,0.75,0.25] is sterile with light [¼,¼,0,¼,¼] (the light avoids the
/// strong center, rate 1); flat-bulk-edge [0.25,1.5,1.5,1.5,0.25] is in the canal: the
/// distribution drifts from [0.3503, 0.0218, 0.2557, …] at Q=1.5 to [⅓, 0, ⅓, 0, ⅓] at Q=1000
/// (max per-site drift 7.8·10⁻²) and the rate follows it, 1.2483 → 4/3.</para>
///
/// <para><b>Honest scope</b>: the absorption identity behind the cross-check needs a Hermitian
/// Hamiltonian and pure Z-dephasing (the dissipator diagonal in the coherence basis); the
/// block-triangularity argument then holds for ANY such Liouvillian. The pinned surface story
/// (sterile/canal verdicts) is established for the dimensionless XY-chain family of
/// <see cref="Foundation.PostEpFlowField"/> at N ≤ 6. Near an exceptional point of the cluster
/// the eigendecomposition that locates range(P) is ill-conditioned, the same caveat
/// <see cref="Foundation.PostEpFlowField"/> documents for its spectral propagation. Not a Claim:
/// a compute primitive of the same standing as <see cref="SlowModeMixing"/> and
/// <see cref="StationaryManifold"/>.</para></summary>
public static class SlowLightDistribution
{
    /// <summary>Modes with |λ| below this are kernel (steady-state) modes and never part of the
    /// slow cluster. Matches <see cref="Foundation.PostEpFlowField"/>'s kernel tolerance.</summary>
    public const double DefaultKernelTolerance = 1e-7;

    /// <summary>Eigenvalues with Re λ within this of the slowest non-kernel Re λ belong to the
    /// slow cluster. Matches the degeneracy tolerance the former averaged carrier used.</summary>
    public const double DefaultClusterTolerance = 1e-6;

    /// <summary>The reading: <paramref name="Rate"/> = −Re λ of the slowest non-kernel cluster;
    /// <paramref name="PerSiteLight"/> = the basis-free distribution light_l = Tr(Π_V·Δ_l)/g,
    /// each in [0, 1]; <paramref name="TotalLight"/> = Σ_l light_l = the drain depth ⟨n_XY⟩;
    /// <paramref name="AbsorptionRate"/> = 2·Σ_l γ_l·light_l (equal to the rate to machine
    /// precision, the Absorption cross-check); <paramref name="Degeneracy"/> = g, the number of
    /// modes sharing the slowest rate (= dim of the slow invariant subspace).</summary>
    public sealed record Reading(
        double Rate,
        IReadOnlyList<double> PerSiteLight,
        double TotalLight,
        double AbsorptionRate,
        int Degeneracy)
    {
        /// <summary>|Rate − AbsorptionRate|, the absorption cross-check residual. Machine
        /// precision at moderate Q; at very large Q it tracks the eigensolver's own accuracy
        /// on the rate (‖L‖ grows with Q), not a failure of the identity.</summary>
        public double AbsorptionResidual => Math.Abs(Rate - AbsorptionRate);
    }

    /// <summary>Compute the slow subspace's light distribution of <paramref name="liouvillian"/>
    /// (a dense 4^n × 4^n operator-space matrix, row-major vec as in
    /// <see cref="RCPsiSquared.Core.Lindblad.PauliDephasingDissipator"/>: vec(ρ)[a·d+b] = ρ[a,b],
    /// site l ↔ bit n−1−l). <paramref name="gammaProfile"/> (length n) is the per-site dephasing
    /// weight vector used for the absorption cross-check; pass the same profile the Liouvillian
    /// was built with. Throws if L has no non-kernel mode.</summary>
    public static Reading Compute(ComplexMatrix liouvillian, int n, IReadOnlyList<double> gammaProfile,
        double kernelTolerance = DefaultKernelTolerance, double clusterTolerance = DefaultClusterTolerance)
    {
        if (n < 1) throw new ArgumentOutOfRangeException(nameof(n), n, "n must be >= 1");
        int d = 1 << n;
        if (liouvillian.RowCount != d * d || liouvillian.ColumnCount != d * d)
            throw new ArgumentException($"liouvillian must be {d * d}×{d * d} for n={n}, got {liouvillian.RowCount}×{liouvillian.ColumnCount}", nameof(liouvillian));
        if (gammaProfile is null) throw new ArgumentNullException(nameof(gammaProfile));
        if (gammaProfile.Count != n)
            throw new ArgumentException($"gammaProfile must have length n={n}, got {gammaProfile.Count}", nameof(gammaProfile));

        var evd = liouvillian.Evd();
        var vals = evd.EigenValues;

        double maxRe = double.NegativeInfinity;
        for (int k = 0; k < vals.Count; k++)
            if (vals[k].Magnitude > kernelTolerance && vals[k].Real > maxRe) maxRe = vals[k].Real;
        if (double.IsNegativeInfinity(maxRe))
            throw new InvalidOperationException("L has no non-kernel mode to read.");

        var cluster = new List<int>();
        for (int k = 0; k < vals.Count; k++)
            if (vals[k].Magnitude > kernelTolerance && Math.Abs(vals[k].Real - maxRe) <= clusterTolerance)
                cluster.Add(k);
        int g = cluster.Count;

        // Orthonormal basis of the slow invariant subspace V = range(P_slow): thin QR of the
        // cluster's right eigenvectors. Any basis of V gives the same Π_V = U·U†; only the
        // span enters, which is exactly what makes the reading basis-free.
        var m = ComplexMatrix.Build.Dense(d * d, g, (row, j) => evd.EigenVectors[row, cluster[j]]);
        var u = m.QR(QRMethod.Thin).Q;

        // diag(Π_V)[x] = Σ_a |U[x,a]|², then light_l = Σ_x Δ_l(x)·diag(Π_V)[x] / g.
        var perSite = new double[n];
        for (int x = 0; x < d * d; x++)
        {
            double w = 0.0;
            for (int a = 0; a < g; a++)
            {
                var c = u[x, a];
                w += c.Real * c.Real + c.Imaginary * c.Imaginary;
            }
            if (w == 0.0) continue;
            int diff = (x / d) ^ (x % d);
            for (int l = 0; l < n; l++)
                if (((diff >> (n - 1 - l)) & 1) != 0) perSite[l] += w;
        }

        double total = 0.0, absorption = 0.0;
        for (int l = 0; l < n; l++)
        {
            perSite[l] /= g;
            total += perSite[l];
            absorption += 2.0 * gammaProfile[l] * perSite[l];
        }
        return new Reading(-maxRe, perSite, total, absorption, g);
    }
}
