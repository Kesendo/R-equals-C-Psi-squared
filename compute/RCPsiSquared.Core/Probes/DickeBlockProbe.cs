using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.CoherenceBlocks;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Probes;

/// <summary>Dicke-state coherence probe ρ₀ = |S_n⟩⟨S_{n+1}| / 2 expressed as a flat coefficient
/// vector in the (n, n+1) coherence-block basis. Canonical initial state for the
/// **F86** K_CC_pr Q-scan (see docs/ANALYTICAL_FORMULAS.md, F86 entry, and
/// docs/proofs/PROOF_F86_QPEAK.md Statement 2).
///
/// <para>|S_n⟩ = (1/√C(N, n)) · Σ_{popcount(p)=n} |p⟩, so |S_n⟩⟨S_{n+1}| has every basis-pair
/// coefficient equal to 1/√(C(N,n)·C(N,n+1)). The factor 1/2 accounts for the Hermitian
/// combination convention used in K_CC_pr observables. The probe lives entirely in the
/// channel-uniform 2D subspace span{|c_1⟩, |c_3⟩} for c=2 chains — orthogonal to the SVD
/// EP-partner modes (PROOF_F86_QPEAK structural finding 2026-05-02 night).</para>
/// </summary>
public static class DickeBlockProbe
{
    public static ComplexVector Build(CoherenceBlock block)
    {
        BlockBasis basis = block.Basis;
        double coeff = 1.0 / (2.0 * Math.Sqrt((double)basis.Mp * basis.Mq));
        var v = MathNet.Numerics.LinearAlgebra.Vector<Complex>.Build.Dense(basis.MTotal);
        for (int i = 0; i < basis.MTotal; i++)
            v[i] = new Complex(coeff, 0.0);
        return v;
    }
}
