using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.States;
using RCPsiSquared.Core.Symmetry;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Receivers;

/// <summary>F71 receiver-engineering signature for a state ψ on a chain. See
/// experiments/J_BLIND_RECEIVER_CLASSES.md.</summary>
public sealed record ReceiverSignature(
    int N,
    int? F71Eigenvalue,
    int NSym,
    int NAsym,
    bool BondBlockBalanced,
    string Prediction);

/// <summary>State-bearing wrapper for a quantum state ψ on a chain, with cached F71 classification
/// and density matrix.
///
/// Selecting <see cref="BondingMode"/> on both ends IS the per-side handshake (no exchange step
/// needed) — see HANDSHAKE_ALGEBRA.md.
/// </summary>
public sealed class Receiver
{
    public ChainSystem? Chain { get; }
    public ComplexVector Psi { get; }
    public int N { get; }

    private int? _f71Class;
    private bool _f71Computed;
    private ComplexMatrix? _rho;

    public Receiver(ComplexVector psi, ChainSystem? chain = null, double normTolerance = 1e-6)
    {
        N = (int)Math.Round(Math.Log2(psi.Count));
        if ((1 << N) != psi.Count) throw new ArgumentException($"psi length {psi.Count} is not a power of 2");
        if (chain is not null && chain.N != N) throw new ArgumentException($"chain.N ({chain.N}) does not match psi N ({N})");
        double norm = Math.Sqrt(psi.ConjugateDotProduct(psi).Real);
        if (Math.Abs(norm - 1.0) > normTolerance)
            throw new ArgumentException($"psi must be normalized (‖psi‖=1); got ‖psi‖={norm:G6}.");
        Psi = psi;
        Chain = chain;
    }

    public static Receiver FromUnnormalized(ComplexVector psi, ChainSystem? chain = null)
    {
        double norm = Math.Sqrt(psi.ConjugateDotProduct(psi).Real);
        if (norm == 0) throw new ArgumentException("psi is the zero vector; cannot normalize.");
        return new Receiver(psi / norm, chain);
    }

    /// <summary>F65 single-excitation bonding-mode receiver |ψ_k⟩ (or pair state with vacuum).</summary>
    public static Receiver BondingMode(ChainSystem chain, int k, bool withVacuum = false)
    {
        var psi = withVacuum ? Core.States.BondingMode.PairState(chain.N, k) : Core.States.BondingMode.Build(chain.N, k);
        return new Receiver(psi, chain);
    }

    /// <summary>F71 eigenstate class of ψ (+1 / −1 / null for mixed). Cached.</summary>
    public int? F71Eigenvalue
    {
        get
        {
            if (!_f71Computed)
            {
                _f71Class = ChainMirror.EigenstateClass(Psi);
                _f71Computed = true;
            }
            return _f71Class;
        }
    }

    /// <summary>Density matrix ρ = |ψ⟩⟨ψ|. Cached.</summary>
    public ComplexMatrix Rho => _rho ??= DensityMatrix.FromStateVector(Psi);

    /// <summary>F71-based receiver-engineering favorability forecast.
    /// Tier 2: verified at N=5 (balanced, F71-eigenstate optimal) and N=6 (unbalanced, F71-breaking wins).
    /// See experiments/J_BLIND_RECEIVER_CLASSES.md.</summary>
    public ReceiverSignature Signature()
    {
        var (nSym, nAsym) = ChainMirror.BondMirrorCounts(N);
        bool balanced = nSym == nAsym;
        var eig = F71Eigenvalue;

        string prediction;
        if (eig is null) prediction = "no-prediction (F71-mixed state — no block structure)";
        else if (balanced) prediction = "capacity-optimal (balanced split, F71-eigenstate)";
        else prediction = "capacity-suboptimal (unbalanced split, F71-eigenstate)";

        return new ReceiverSignature(N, eig, nSym, nAsym, balanced, prediction);
    }
}
