using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.F86.JordanWigner;

/// <summary>F86 Item 1' Direction (b'') JW track, T9: foundational basis-transformation
/// primitive from the computational basis of the c=2 (n=1, n+1=2) coherence block to the
/// JW-mode basis indexed by triples (k, k₁, k₂) with k ∈ [1, N], 1 ≤ k₁ &lt; k₂ ≤ N.
///
/// <para>The transformation comes from textbook Wick contraction with sine-mode JW
/// operators η_k = Σ_p ψ_k(p) c_p (T1 <see cref="XyJordanWignerModes"/>):</para>
/// <list type="bullet">
///   <item>Popcount-1 part: ⟨p | k⟩ = ψ_k(p) on the single-excitation site p.</item>
///   <item>Popcount-2 part: ⟨p, q | k₁, k₂⟩ = ψ_{k₁}(p)·ψ_{k₂}(q) − ψ_{k₁}(q)·ψ_{k₂}(p) for
///         p &lt; q (Slater determinant; <b>no 1/√2 prefactor</b> — adding one breaks
///         unitarity of U).</item>
///   <item>Coherence-block matrix element:
///         <c>U[(p, q) flat, (k, k₁, k₂) α] = ψ_k(p) · [ψ_{k₁}(p)·ψ_{k₂}(q) − ψ_{k₁}(q)·ψ_{k₂}(p)]</c>.</item>
/// </list>
///
/// <para><b>Tier outcome: <see cref="Tier.Tier1Derived"/>.</b> Wick contraction is exact
/// textbook XY-JW + free-fermion algebra. The witnesses verify the algebra at runtime:</para>
/// <list type="bullet">
///   <item><see cref="OrthonormalityResidual"/> = <c>‖U · U^† − I‖_F</c> (Wick-contracted
///         Slater determinants are orthonormal under sine-mode orthonormality).</item>
///   <item><see cref="MhTotalDiagonalityResidual"/> = <c>‖(U⁻¹ · MhTotal · U)_off‖_F</c>
///         (free-fermion XY is diagonal under JW).</item>
///   <item><see cref="MhTotalEigenvalueMatchResidual"/> = <c>max_α |MhTotal_JW[α, α] −
///         (−i·(ε_k − ε_{k₁} − ε_{k₂}))|</c> (textbook free-fermion dispersion identity).</item>
/// </list>
/// All three witnesses are verified below <see cref="Tolerance"/> = 1e-10 at construction.
///
/// <para>This primitive is the foundation for the open Tier1-promotion path of F86 Item 1'
/// Direction (b''): the analytical eigenstructure of L(Q) in JW basis (D + Q·γ·MhTotal with
/// MhTotal diagonal, D 20%-off-diagonal in JW basis), which gives the EP-relevant
/// dispersion-degenerate mode pairs needed for the closed-form Q_peak / HWHM derivation.</para>
///
/// <para>Anchor: <c>docs/proofs/PROOF_F86_QPEAK.md</c> Item 1' Direction (b'') (JW track) +
/// textbook XY Jordan-Wigner sine-mode basis.</para>
/// </summary>
public sealed class JwBlockBasis : Claim
{
    /// <summary>Tolerance for the three orthonormality / diagonality / eigenvalue-match
    /// witnesses verified at construction: <c>1e-10</c>.</summary>
    public const double Tolerance = 1e-10;

    public CoherenceBlock Block { get; }

    /// <summary>Composition: T1 sine-mode basis + dispersion. Built once at <see cref="Build"/>
    /// time and reused for every triple.</summary>
    public XyJordanWignerModes Modes { get; }

    /// <summary>Triples (k, k₁, k₂) enumerating the JW-mode basis of the coherence block.
    /// Length Mtot = N · C(N, 2). Ordered by k (outer), then k₁ &lt; k₂ (inner).</summary>
    public IReadOnlyList<JwTriple> Triples { get; }

    /// <summary>Mtot × Mtot unitary basis-transformation matrix. Columns are JW basis
    /// vectors expressed in the computational basis: column α is the JW state (k, k₁, k₂)
    /// of triple α, row index is the flat (p, q) index in
    /// <see cref="BlockBasis.FlatIndex"/> order.</summary>
    public ComplexMatrix U { get; }

    /// <summary><c>U^†</c>. Witness <see cref="OrthonormalityResidual"/> verifies
    /// <c>U · U^† = I</c>, so <c>U^† = U⁻¹</c>.</summary>
    public ComplexMatrix Uinv { get; }

    /// <summary><c>‖U · U^† − I‖_F</c>. Verified below <see cref="Tolerance"/> at
    /// construction.</summary>
    public double OrthonormalityResidual { get; }

    /// <summary><c>‖(U^† · MhTotal · U)_off-diagonal‖_F</c>. Verified below
    /// <see cref="Tolerance"/> at construction (free-fermion XY is diagonal under JW).</summary>
    public double MhTotalDiagonalityResidual { get; }

    /// <summary><c>max_α |MhTotal_JW[α, α] − (−i·(ε_k − ε_{k₁} − ε_{k₂}))|</c>.
    /// Verified below <see cref="Tolerance"/> at construction (free-fermion dispersion
    /// identity).</summary>
    public double MhTotalEigenvalueMatchResidual { get; }

    /// <summary>Public factory: validates c=2, builds the T1 sine-mode basis, enumerates
    /// triples (k, k₁, k₂), assembles U from Wick-contracted Slater overlaps, and verifies
    /// the three runtime witnesses. Throws if any witness exceeds <see cref="Tolerance"/>.
    /// </summary>
    public static JwBlockBasis Build(CoherenceBlock block)
    {
        if (block.C != 2)
            throw new ArgumentException(
                $"JwBlockBasis applies only to the c=2 stratum (n=1, n+1=2); got c={block.C} (N={block.N}, n={block.LowerPopcount}).",
                nameof(block));

        const double J = 1.0;
        int N = block.N;
        var modes = XyJordanWignerModes.Build(N, J);
        BlockBasis basis = block.Basis;
        int Mtot = basis.MTotal;

        // 1) Enumerate triples (k, k₁, k₂): for k = 1..N, for k₁ = 1..N, for k₂ = k₁+1..N.
        //    Length = N · C(N, 2) = N · N · (N - 1) / 2.
        var triples = new List<JwTriple>(N * N * (N - 1) / 2);
        for (int k = 1; k <= N; k++)
            for (int k1 = 1; k1 <= N; k1++)
                for (int k2 = k1 + 1; k2 <= N; k2++)
                    triples.Add(new JwTriple(k, k1, k2));

        if (triples.Count != Mtot)
            throw new InvalidOperationException(
                $"Triple-count mismatch: enumerated {triples.Count} triples but coherence block has Mtot = {Mtot}.");

        // 2) Build U: column α = JW basis vector |k⟩⟨k₁,k₂| in computational basis.
        //    No 1/√2 prefactor on the Slater — Wick-contracted bra-ket already absorbs the
        //    antisymmetry through the single (b<c) representative. Adding 1/√2 breaks
        //    orthonormality; this is the documented trap from the 2026-05-08 exploration.
        var uRaw = new Complex[Mtot, Mtot];
        for (int pIdx = 0; pIdx < basis.Mp; pIdx++)
        {
            long pState = basis.StatesP[pIdx];
            int pSite = SingleExcitationSite(pState, N);
            for (int qIdx = 0; qIdx < basis.Mq; qIdx++)
            {
                long qState = basis.StatesQ[qIdx];
                var (qSiteB, qSiteC) = TwoExcitationSites(qState, N);
                int flat = basis.FlatIndex(pState, qState);
                for (int alpha = 0; alpha < triples.Count; alpha++)
                {
                    var (k, k1, k2) = (triples[alpha].K, triples[alpha].K1, triples[alpha].K2);
                    double psik_p = modes.SineMode(k, pSite);
                    double psik1_b = modes.SineMode(k1, qSiteB);
                    double psik2_c = modes.SineMode(k2, qSiteC);
                    double psik1_c = modes.SineMode(k1, qSiteC);
                    double psik2_b = modes.SineMode(k2, qSiteB);
                    double slater = psik1_b * psik2_c - psik1_c * psik2_b;
                    uRaw[flat, alpha] = new Complex(psik_p * slater, 0.0);
                }
            }
        }

        var U = Matrix<Complex>.Build.DenseOfArray(uRaw);
        var Uinv = U.ConjugateTranspose();

        // 3) Orthonormality witness: ‖U · U^† − I‖_F.
        var identity = Matrix<Complex>.Build.DenseIdentity(Mtot);
        double orthoResidual = (U * Uinv - identity).FrobeniusNorm();
        if (orthoResidual > Tolerance)
            throw new InvalidOperationException(
                $"JwBlockBasis (N={N}): orthonormality residual {orthoResidual:G3} exceeds tolerance {Tolerance:G3}.");

        // 4) MhTotal diagonality witness: ‖(U^† · MhTotal · U)_off‖_F.
        var mhTotalJw = Uinv * block.Decomposition.MhTotal * U;
        double diagSqSum = 0.0;
        for (int i = 0; i < Mtot; i++)
            for (int j = 0; j < Mtot; j++)
            {
                if (i == j) continue;
                Complex z = mhTotalJw[i, j];
                diagSqSum += z.Real * z.Real + z.Imaginary * z.Imaginary;
            }
        double diagResidual = Math.Sqrt(diagSqSum);
        if (diagResidual > Tolerance)
            throw new InvalidOperationException(
                $"JwBlockBasis (N={N}): MhTotal diagonality residual {diagResidual:G3} exceeds tolerance {Tolerance:G3}.");

        // 5) Eigenvalue-match witness: max_α |MhTotal_JW[α, α] − (−i·(ε_k − ε_{k₁} − ε_{k₂}))|.
        double eigMaxDev = 0.0;
        for (int alpha = 0; alpha < triples.Count; alpha++)
        {
            var (k, k1, k2) = (triples[alpha].K, triples[alpha].K1, triples[alpha].K2);
            double expectedImag = -(modes.Dispersion[k - 1] - modes.Dispersion[k1 - 1] - modes.Dispersion[k2 - 1]);
            Complex expected = new Complex(0.0, expectedImag);
            Complex actual = mhTotalJw[alpha, alpha];
            double dev = (actual - expected).Magnitude;
            if (dev > eigMaxDev) eigMaxDev = dev;
        }
        if (eigMaxDev > Tolerance)
            throw new InvalidOperationException(
                $"JwBlockBasis (N={N}): MhTotal eigenvalue-match residual {eigMaxDev:G3} exceeds tolerance {Tolerance:G3}.");

        return new JwBlockBasis(block, modes, triples, U, Uinv, orthoResidual, diagResidual, eigMaxDev);
    }

    private JwBlockBasis(
        CoherenceBlock block,
        XyJordanWignerModes modes,
        IReadOnlyList<JwTriple> triples,
        ComplexMatrix u,
        ComplexMatrix uinv,
        double orthoResidual,
        double diagResidual,
        double eigMaxDev)
        : base("c=2 JW-mode basis transformation U from computational basis to (k, k₁, k₂) triples",
               Tier.Tier1Derived,
               "docs/proofs/PROOF_F86_QPEAK.md Item 1' Direction (b'') (JW track) + textbook XY Jordan-Wigner sine-mode basis")
    {
        Block = block;
        Modes = modes;
        Triples = triples;
        U = u;
        Uinv = uinv;
        OrthonormalityResidual = orthoResidual;
        MhTotalDiagonalityResidual = diagResidual;
        MhTotalEigenvalueMatchResidual = eigMaxDev;
    }

    /// <summary>Site index (0-indexed, big-endian: site 0 = MSB) of the single excitation in
    /// a popcount-1 bitmask state. Throws if popcount(state) ≠ 1.</summary>
    private static int SingleExcitationSite(long state, int N)
    {
        int site = -1;
        for (int s = 0; s < N; s++)
        {
            if ((state & (1L << (N - 1 - s))) != 0)
            {
                if (site != -1)
                    throw new ArgumentException(
                        $"State 0x{state:X} has popcount > 1; expected single excitation.");
                site = s;
            }
        }
        if (site == -1)
            throw new ArgumentException(
                $"State 0x{state:X} has popcount = 0; expected single excitation.");
        return site;
    }

    /// <summary>Site indices (b, c) with b &lt; c (0-indexed, big-endian: site 0 = MSB) of the
    /// two excitations in a popcount-2 bitmask state. Throws if popcount(state) ≠ 2.</summary>
    private static (int b, int c) TwoExcitationSites(long state, int N)
    {
        int b = -1, c = -1;
        for (int s = 0; s < N; s++)
        {
            if ((state & (1L << (N - 1 - s))) != 0)
            {
                if (b == -1) b = s;
                else if (c == -1) c = s;
                else throw new ArgumentException(
                    $"State 0x{state:X} has popcount > 2; expected two excitations.");
            }
        }
        if (c == -1)
            throw new ArgumentException(
                $"State 0x{state:X} has popcount < 2; expected two excitations.");
        return (b, c);
    }

    public override string DisplayName =>
        $"JW-mode basis (N={Block.N}, Mtot={Triples.Count})";

    public override string Summary =>
        $"U: computational ↔ (k, k₁, k₂); ortho={OrthonormalityResidual:G3}, " +
        $"diag={MhTotalDiagonalityResidual:G3}, eig={MhTotalEigenvalueMatchResidual:G3} " +
        $"(tol {Tolerance:G3}) ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("N", summary: Block.N.ToString());
            yield return new InspectableNode("Mtot", summary: Triples.Count.ToString());
            yield return Modes;
            yield return InspectableNode.RealScalar("orthonormality residual", OrthonormalityResidual, "G3");
            yield return InspectableNode.RealScalar("MhTotal diagonality residual", MhTotalDiagonalityResidual, "G3");
            yield return InspectableNode.RealScalar("MhTotal eigenvalue-match residual", MhTotalEigenvalueMatchResidual, "G3");
            yield return InspectableNode.RealScalar("tolerance", Tolerance, "G3");
            yield return InspectableNode.Group("triples (k, k₁, k₂)", Triples, Triples.Count);
        }
    }
}

/// <summary>JW-mode triple (k, k₁, k₂) with k ∈ [1, N] the popcount-1 mode index and
/// 1 ≤ k₁ &lt; k₂ ≤ N the popcount-2 mode index pair. Visible in the inspection tree under
/// <see cref="JwBlockBasis.Triples"/>.</summary>
public sealed record JwTriple(int K, int K1, int K2) : IInspectable
{
    public string DisplayName => $"(k={K}, k₁={K1}, k₂={K2})";
    public string Summary => DisplayName;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            yield return new InspectableNode("k", summary: K.ToString());
            yield return new InspectableNode("k₁", summary: K1.ToString());
            yield return new InspectableNode("k₂", summary: K2.ToString());
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;
}
