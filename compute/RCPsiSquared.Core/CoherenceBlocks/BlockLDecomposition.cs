using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.CoherenceBlocks;

/// <summary>Block-restricted Liouvillian split L = D + Σ_b J_b · M_H_per_bond[b].
///
/// D is the J-independent dephasing diagonal: D[i, i] = −2γ₀ · HD(p, q) on the (p, q)
/// basis pair at flat index i. M_H_per_bond[b] is the J-coefficient contribution to L_H
/// from bond b.
///
/// Convention: H_b = (J/2)·(X_b X_{b+1} + Y_b Y_{b+1}). Big-endian: site 0 = MSB.
/// Each bond flips a single adjacent-bit pair where the two bits are opposite, contributing
/// ±i to the corresponding off-diagonal entries.
/// </summary>
public sealed class BlockLDecomposition
{
    public CoherenceBlock Block { get; }
    public BlockBasis Basis => Block.Basis;
    public ComplexMatrix D { get; }
    public IReadOnlyList<ComplexMatrix> MhPerBond { get; }
    public int NumBonds => MhPerBond.Count;

    private ComplexMatrix? _mhTotal;

    /// <summary>Σ_b M_H_per_bond[b], lazy and cached. Used for uniform-J assembly and
    /// inter-channel SVD without rebuilding from scratch.</summary>
    public ComplexMatrix MhTotal
    {
        get
        {
            if (_mhTotal is null)
            {
                var sum = MhPerBond[0].Clone();
                for (int b = 1; b < MhPerBond.Count; b++) sum = sum + MhPerBond[b];
                _mhTotal = sum;
            }
            return _mhTotal;
        }
    }

    private BlockLDecomposition(CoherenceBlock block, ComplexMatrix d, IReadOnlyList<ComplexMatrix> mh)
    {
        Block = block;
        D = d;
        MhPerBond = mh;
    }

    public static BlockLDecomposition Build(CoherenceBlock block)
    {
        int N = block.N;
        double gamma0 = block.GammaZero;
        BlockBasis basis = block.Basis;
        int Mtot = basis.MTotal;
        int numBonds = N - 1;

        // Build dense first — sparse insertions on MathNet.SparseMatrix are O(nnz) per write.
        // Each per-bond M_H is ~2·Mtot nonzeros; building dense and (optionally) sparsifying
        // later is faster for the small N range used here.
        var dRaw = new Complex[Mtot, Mtot];
        var mhRaw = new Complex[numBonds][,];
        for (int b = 0; b < numBonds; b++) mhRaw[b] = new Complex[Mtot, Mtot];

        var pFlips = new Dictionary<int, (int Bond, int Flipped)[]>(basis.Mp);
        foreach (int p in basis.StatesP) pFlips[p] = BondFlipTargets(p, N).ToArray();
        var qFlips = new Dictionary<int, (int Bond, int Flipped)[]>(basis.Mq);
        foreach (int q in basis.StatesQ) qFlips[q] = BondFlipTargets(q, N).ToArray();

        foreach (int p in basis.StatesP)
        {
            foreach (int q in basis.StatesQ)
            {
                int i = basis.FlatIndex(p, q);
                int hd = BitOperations.PopCount((uint)(p ^ q));
                dRaw[i, i] = new Complex(-2.0 * gamma0 * hd, 0.0);

                foreach (var (bond, pFlipped) in pFlips[p])
                {
                    int j = basis.FlatIndex(pFlipped, q);
                    mhRaw[bond][j, i] += new Complex(0.0, -1.0);
                }
                foreach (var (bond, qFlipped) in qFlips[q])
                {
                    int j = basis.FlatIndex(p, qFlipped);
                    mhRaw[bond][j, i] += new Complex(0.0, +1.0);
                }
            }
        }

        var d = Matrix<Complex>.Build.DenseOfArray(dRaw);
        var mh = new ComplexMatrix[numBonds];
        for (int b = 0; b < numBonds; b++) mh[b] = Matrix<Complex>.Build.DenseOfArray(mhRaw[b]);
        return new BlockLDecomposition(block, d, mh);
    }

    /// <summary>Yield (bond, flipped_state) for every bond where the state's two adjacent bits
    /// are opposite (i.e. an XX+YY swap acts non-trivially). Big-endian: site 0 = MSB.</summary>
    private static IEnumerable<(int Bond, int Flipped)> BondFlipTargets(int state, int N)
    {
        for (int bond = 0; bond < N - 1; bond++)
        {
            int bitHi = (state >> (N - 1 - bond)) & 1;
            int bitLo = (state >> (N - 2 - bond)) & 1;
            if (bitHi != bitLo)
            {
                int mask = (1 << (N - 1 - bond)) | (1 << (N - 2 - bond));
                yield return (bond, state ^ mask);
            }
        }
    }

    /// <summary>L = D + Σ_b J_b · M_H_per_bond[b], with one J per bond.</summary>
    public ComplexMatrix AssembleAt(IReadOnlyList<double> bondCouplings)
    {
        if (bondCouplings.Count != NumBonds)
            throw new ArgumentException($"expected {NumBonds} bond couplings, got {bondCouplings.Count}.");
        var L = D.Clone();
        for (int b = 0; b < NumBonds; b++)
            L = L + (Complex)bondCouplings[b] * MhPerBond[b];
        return L;
    }

    /// <summary>L = D + J · MhTotal for uniform J across all bonds.</summary>
    public ComplexMatrix AssembleUniform(double J) => D + (Complex)J * MhTotal;

    /// <summary>Build the uniform-J Liouvillian L = D + J · Σ_b M_h_per_bond[b] DIRECTLY,
    /// without ever materialising the per-bond matrices. Memory cost: 1 × Mtot² complex
    /// (the L matrix itself) instead of (NumBonds + 2) × Mtot² for the full
    /// <see cref="BlockLDecomposition"/>. Use when only a single uniform-J L is needed
    /// (no per-bond Hellmann-Feynman): e.g. <see cref="F86.Item1Derivation.C2FullBlockSigmaAnatomy"/>
    /// at high N where the per-bond storage exceeds available RAM.
    /// <para>Numerically equivalent to <c>block.Decomposition.AssembleUniform(J)</c>,
    /// but skips the BlockLDecomposition allocation.</para></summary>
    public static ComplexMatrix BuildUniformLAt(CoherenceBlock block, double J)
    {
        int N = block.N;
        double gamma0 = block.GammaZero;
        BlockBasis basis = block.Basis;
        int Mtot = basis.MTotal;

        var lRaw = new Complex[Mtot, Mtot];

        var pFlips = new Dictionary<int, (int Bond, int Flipped)[]>(basis.Mp);
        foreach (int p in basis.StatesP) pFlips[p] = BondFlipTargets(p, N).ToArray();
        var qFlips = new Dictionary<int, (int Bond, int Flipped)[]>(basis.Mq);
        foreach (int q in basis.StatesQ) qFlips[q] = BondFlipTargets(q, N).ToArray();

        Complex hopP = new Complex(0.0, -J);
        Complex hopQ = new Complex(0.0, +J);

        foreach (int p in basis.StatesP)
        {
            foreach (int q in basis.StatesQ)
            {
                int i = basis.FlatIndex(p, q);
                int hd = BitOperations.PopCount((uint)(p ^ q));
                lRaw[i, i] = new Complex(-2.0 * gamma0 * hd, 0.0);

                foreach (var (_, pFlipped) in pFlips[p])
                    lRaw[basis.FlatIndex(pFlipped, q), i] += hopP;
                foreach (var (_, qFlipped) in qFlips[q])
                    lRaw[basis.FlatIndex(p, qFlipped), i] += hopQ;
            }
        }

        return Matrix<Complex>.Build.DenseOfArray(lRaw);
    }
}
