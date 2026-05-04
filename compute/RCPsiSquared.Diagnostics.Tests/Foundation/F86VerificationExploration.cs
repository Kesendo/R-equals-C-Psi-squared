using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Core.Decomposition;
using RCPsiSquared.Core.States;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>F86 verification via the trio's state-level diagnostics. Two tasks (per
/// 2026-05-04 conversation):
///
/// <list type="number">
///   <item><b>MemoryAxisRho on F86-EP-regime test states.</b> The K_CC_pr observable
///         lives in the (n, n+1) popcount-coherence block. A representative coherence
///         state should have non-zero Π²-odd memory content (the framework's memory
///         channel for Π²-odd Hamiltonian dynamics, F81's M_anti). This verifies that
///         hardware-prepared test states sit in the right sector for F86 measurement.</item>
///   <item><b>BlochAxisReading + MemoryAxisRho on c=2 4-mode-basis-derived density
///         matrices.</b> The {|c_1⟩, |c_3⟩, |u_0⟩, |v_0⟩} basis (FourModeBasis at c=2)
///         spans the F86 minimal effective subspace. Embed each basis vector as a
///         coherence-block operator + I/d perturbation to form a positive density
///         matrix; the per-qubit Bloch and Frobenius readings give independent
///         state-level fingerprints of the 4-mode subspace structure.</item>
/// </list>
///
/// <para>Neither test PROVES F86's open elements (4-mode cross-coupling closed form,
/// σ_0(c, N → ∞), multi-k extension repair). They confirm the framework's state-level
/// Π²-odd reading is consistent with F86's premises.</para>
/// </summary>
public class F86VerificationExploration
{
    private readonly ITestOutputHelper _output;
    public F86VerificationExploration(ITestOutputHelper output) => _output = output;

    [Fact]
    public void DumpF86CoherenceState_HasPi2OddMemoryContent_AtN5C2()
    {
        const int N = 5;
        int d = 1 << N;
        var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);

        _output.WriteLine($"=== F86 coherence states at N=5, c=2 (popcounts 1↔2) ===");
        _output.WriteLine("(K_CC_pr observable measures coherence between adjacent popcount sectors)");
        _output.WriteLine("");

        // HD = 1: |00001⟩ ↔ |00011⟩ (single bit flip at site 3)
        DumpCoherenceState(N, d, chain, "|00001⟩+|00011⟩ (HD=1)", 0b00001, 0b00011);
        // HD = 1, different position: |00010⟩ ↔ |00110⟩
        DumpCoherenceState(N, d, chain, "|00010⟩+|00110⟩ (HD=1)", 0b00010, 0b00110);
        // HD = 3: |00001⟩ ↔ |11100⟩ (popcounts 1, 3 — wrong popcount diff)
        // Need HD=3 with popcount diff 1: |00001⟩ ↔ |01110⟩ (popcount 1, 3 — diff 2, also wrong)
        // For popcount diff 1 + HD=3: e.g. |00001⟩ ↔ |11000⟩? popcount 1, 2, HD=3 ✓
        DumpCoherenceState(N, d, chain, "|00001⟩+|11000⟩ (HD=3)", 0b00001, 0b11000);
        DumpCoherenceState(N, d, chain, "|00010⟩+|01100⟩ (HD=3)", 0b00010, 0b01100);
    }

    [Fact]
    public void DumpFourModeBasisVectors_AsTestStates_AtN5C2()
    {
        const int N = 5;
        const int n = 1; // c=min(1, 3)+1 = 2
        int d = 1 << N;

        var chain = new ChainSystem(N, J: 1.0, GammaZero: 0.05);
        var block = new CoherenceBlock(N, n, gammaZero: chain.GammaZero);
        Assert.Equal(2, block.C); // chromaticity 2

        var fmb = FourModeBasis.Build(block);

        _output.WriteLine($"=== 4-mode basis as test states, N={N}, n={n}, c={block.C} ===");
        _output.WriteLine($"orthonormality residual: {fmb.OffOrthonormalityResidual:E2}");
        _output.WriteLine("(each basis column embedded as a coherence operator + I/d, normalised)");
        _output.WriteLine("");

        var labels = new[] { "|c_1⟩  channel-uniform HD=1", "|c_3⟩  channel-uniform HD=3", "|u_0⟩  SVD-top right", "|v_0⟩  SVD-top left" };
        for (int col = 0; col < 4; col++)
        {
            var v = fmb.BasisMatrix.Column(col);
            var rho = EmbedAsDensityMatrix(v, block, d);

            var memReading = MemoryAxisRho.Decompose(rho, chain);
            var blochReading = BlochAxisReading.Compute(rho, N);

            _output.WriteLine($"[{col}] {labels[col]}");
            _output.WriteLine($"    static={memReading.StaticFraction:F4}  memory={memReading.MemoryFraction:F4}  Π²-odd/mem={memReading.Pi2OddFractionWithinMemory:F4}");
            _output.WriteLine($"    per-qubit |r|: {string.Join("  ", blochReading.Qubits.Select(q => $"q{q.QubitIndex}={q.RMagnitude:F3}"))}");
            _output.WriteLine($"    dominant axes: {string.Join(" ", blochReading.Qubits.Select(q => q.DominantAxis == 'I' ? "I" : $"{q.DominantAxis}{(q.DominantSign>0?"+":"-")}"))}");
        }
    }

    private void DumpCoherenceState(int N, int d, ChainSystem chain, string name, int pBits, int qBits)
    {
        var psi = ComplexVector.Build.Dense(d);
        psi[pBits] = 1.0 / Math.Sqrt(2);
        psi[qBits] = 1.0 / Math.Sqrt(2);
        var rho = DensityMatrix.FromStateVector(psi);

        var memReading = MemoryAxisRho.Decompose(rho, chain);
        var blochReading = BlochAxisReading.Compute(rho, N);

        _output.WriteLine($"  {name}");
        _output.WriteLine($"    static={memReading.StaticFraction:F4}  memory={memReading.MemoryFraction:F4}  Π²-odd/mem={memReading.Pi2OddFractionWithinMemory:F4}");
        _output.WriteLine($"    per-qubit |r|: {string.Join("  ", blochReading.Qubits.Select(q => $"q{q.QubitIndex}={q.RMagnitude:F3}"))}");
        _output.WriteLine("");
    }

    /// <summary>Embed a coherence-block vector v ∈ ℂ^(Mp·Mq) as a Hermitian operator
    /// A = Σ v[i,j]·(|p_i⟩⟨q_j| + h.c.) on d×d, then form ρ = (I/d + α·A) / Tr with α
    /// chosen small enough for positivity and renormalised.</summary>
    private static ComplexMatrix EmbedAsDensityMatrix(ComplexVector v, CoherenceBlock block, int d)
    {
        var A = ComplexMatrix.Build.Dense(d, d);
        for (int i = 0; i < block.Basis.Mp; i++)
        {
            int p = block.Basis.StatesP[i];
            for (int j = 0; j < block.Basis.Mq; j++)
            {
                int q = block.Basis.StatesQ[j];
                Complex coeff = v[i * block.Basis.Mq + j];
                A[p, q] += coeff;
                A[q, p] += Complex.Conjugate(coeff);
            }
        }
        // Normalise A to operator-norm 1, then add to I/d with small alpha for positivity
        double opNorm = A.L2Norm();
        if (opNorm > 0) A = A / opNorm;
        double alpha = 0.4 / d; // ρ = I/d + α·A, eigenvalue ≥ 1/d − α ≥ 1/d − 0.4/d = 0.6/d > 0
        var rho = ComplexMatrix.Build.DiagonalIdentity(d) / d + alpha * A;
        // Renormalise trace to 1 (alpha·A has zero trace, but numerical drift)
        double tr = rho.Trace().Real;
        if (Math.Abs(tr - 1.0) > 1e-10) rho = rho / tr;
        return rho;
    }
}
