using System.Numerics;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.States;
using RCPsiSquared.Diagnostics.DZero;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>F80-lift dynamics exploration: propagate a Y-state under a Π²-odd
/// Hamiltonian + Z-dephasing and watch the trio's fractions evolve. Tests whether
/// the Π²-odd content visible in the static MemoryAxisRho reading manifests as
/// time-evolution dynamics (the F80 ±2iλ lift made observable).
///
/// <para>The trio's dynamic question (open thread from project_y_phase_pi2_odd_lens.md):
/// at t=0 a Y-state has ~57% Π²-odd-fraction-within-memory at N=3. As ρ(t) evolves
/// under L (XY-like Π²-odd Hamiltonian + Z-dephasing), does the Π²-odd fraction
/// stay constant (inheritance preserves the structural identity through dynamics),
/// or does it decay differently from Π²-even?</para>
/// </summary>
public class F80LiftDynamicsExploration
{
    private readonly ITestOutputHelper _output;
    public F80LiftDynamicsExploration(ITestOutputHelper output) => _output = output;

    [Fact]
    public void DumpYStateUnderPi2OddH_DynamicsTrajectory_AtN3()
    {
        const int N = 3;
        const double J = 1.0;
        const double gamma = 0.05;

        _output.WriteLine($"=== F80-lift dynamics: Y-state under H = J·Σ X_b Y_{{b+1}}  +  γ Z-dephasing ===");
        _output.WriteLine($"N = {N}, J = {J}, γ_Z per site = {gamma}");
        _output.WriteLine($"Π²-odd 2-body H (F79 universal class): expect Spec(M) = ±2i·Spec(H_non-truly)");
        _output.WriteLine("");

        // ChainSystem skeleton (only used for bonds + chain.N + dephase letter; we override L)
        var chain = new ChainSystem(N, J: J, GammaZero: gamma);

        // Π²-odd 2-body Hamiltonian on chain bonds
        var H = PauliHamiltonian.Bilinear(N, chain.Bonds, new[]
        {
            (PauliLetter.X, PauliLetter.Y, (Complex)J),
        }).ToMatrix();

        var gammas = Enumerable.Repeat(gamma, N).ToArray();
        var L = PauliDephasingDissipator.BuildZ(H, gammas);

        // Stationary modes for THIS L (kernel structure may differ from default chain L)
        var sm = StationaryModes.Compute(chain, Loverride: L);
        _output.WriteLine($"kernel dim of custom L: {sm.KernelDimension}");
        _output.WriteLine("");

        // Eigendecomp of L for time evolution
        var evd = L.Evd();
        var lambdas = evd.EigenValues;
        var R = evd.EigenVectors;
        var Rinv = R.Inverse();

        // Initial states to compare
        var initialStates = new (string Name, ComplexVector Psi)[]
        {
            ("|+i,+i,+i⟩ (Y)",  YBasisProduct(N, new[] { +1, +1, +1 })),
            ("|+++⟩ (X)",        PolarityState.Uniform(N, sign: +1)),
        };

        double[] times = { 0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0 };

        foreach (var (name, psi0) in initialStates)
        {
            _output.WriteLine($"--- initial: {name} ---");
            _output.WriteLine($"{"t",6} {"static",10} {"memory",10} {"Π²-odd/mem",12}");
            _output.WriteLine(new string('-', 42));

            var rho0 = DensityMatrix.FromStateVector(psi0);
            var vecRho0 = VecF(rho0);

            // c[k] = (R^{-1} · vec(ρ_0))[k]; trajectory = R · diag(exp(λ_k · t)) · c
            var c = Rinv * vecRho0;

            int d = 1 << N;
            foreach (double t in times)
            {
                var expLambdaTimesC = ComplexVector.Build.Dense(c.Count);
                for (int k = 0; k < c.Count; k++)
                    expLambdaTimesC[k] = c[k] * Complex.Exp(lambdas[k] * t);
                var vecRhoT = R * expLambdaTimesC;
                var rhoT = UnvecF(vecRhoT, d);
                rhoT = (rhoT + rhoT.ConjugateTranspose()) / 2.0; // symmetrize against numerical drift

                var result = MemoryAxisRho.Decompose(rhoT, chain, sm);
                _output.WriteLine($"{t,6:F2} {result.StaticFraction,10:F4} {result.MemoryFraction,10:F4} {result.Pi2OddFractionWithinMemory,12:F4}");
            }
            _output.WriteLine("");
        }
    }

    // === Helpers ===

    private static ComplexVector YBasisProduct(int N, IReadOnlyList<int> signs)
    {
        int d = 1 << N;
        var vec = ComplexVector.Build.Dense(d);
        double norm = 1.0 / Math.Sqrt(d);
        for (int idx = 0; idx < d; idx++)
        {
            Complex amp = Complex.One;
            for (int k = 0; k < N; k++)
            {
                int bit = (idx >> (N - 1 - k)) & 1;
                if (bit == 1)
                    amp *= signs[k] == +1 ? Complex.ImaginaryOne : -Complex.ImaginaryOne;
            }
            vec[idx] = amp * norm;
        }
        return vec;
    }

    private static ComplexVector VecF(ComplexMatrix m)
    {
        int d = m.RowCount;
        var v = ComplexVector.Build.Dense(d * d);
        for (int j = 0; j < d; j++)
            for (int i = 0; i < d; i++)
                v[j * d + i] = m[i, j];
        return v;
    }

    private static ComplexMatrix UnvecF(ComplexVector v, int d)
    {
        var m = ComplexMatrix.Build.Dense(d, d);
        for (int j = 0; j < d; j++)
            for (int i = 0; i < d; i++)
                m[i, j] = v[j * d + i];
        return m;
    }
}
