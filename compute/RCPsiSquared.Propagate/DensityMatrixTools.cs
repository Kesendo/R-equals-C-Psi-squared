using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Propagate;

/// <summary>
/// Quantum information metrics on density matrices.
/// </summary>
public static class DensityMatrixTools
{
    /// <summary>
    /// Partial trace: keep specified qubits, trace out the rest.
    /// Uses precomputed index map for O(d^2) performance.
    /// </summary>
    public static Matrix<Complex> PartialTrace(Matrix<Complex> rho, int nQubits, int[] keep)
    {
        int dim = 1 << nQubits;
        int nKeep = keep.Length;
        int dKeep = 1 << nKeep;
        var traced = Enumerable.Range(0, nQubits).Except(keep).ToArray();
        int nTraced = traced.Length;

        // Precompute index mapping: for each full index, store (kept_index, traced_bits)
        var keptIdx = new int[dim];
        var tracedBits = new int[dim];
        for (int i = 0; i < dim; i++)
        {
            int ki = 0;
            for (int m = 0; m < nKeep; m++)
                ki |= ((i >> (nQubits - 1 - keep[m])) & 1) << (nKeep - 1 - m);
            keptIdx[i] = ki;

            int tb = 0;
            for (int m = 0; m < nTraced; m++)
                tb |= ((i >> (nQubits - 1 - traced[m])) & 1) << (nTraced - 1 - m);
            tracedBits[i] = tb;
        }

        var result = DenseMatrix.Create(dKeep, dKeep, Complex.Zero);
        var rhoVals = ((DenseMatrix)rho).Values;

        for (int i = 0; i < dim; i++)
        {
            int ki = keptIdx[i];
            int tbi = tracedBits[i];
            for (int j = 0; j < dim; j++)
            {
                if (tracedBits[j] != tbi) continue;
                result[ki, keptIdx[j]] += rhoVals[j * dim + i]; // column-major
            }
        }
        return result;
    }

    /// <summary>
    /// Von Neumann entropy S = -Tr(rho log2 rho).
    /// </summary>
    public static double VonNeumannEntropy(Matrix<Complex> rho)
    {
        var evd = rho.Evd();
        double S = 0;
        foreach (var ev in evd.EigenValues)
        {
            double p = ev.Real;
            if (p > 1e-15)
                S -= p * Math.Log2(p);
        }
        return S;
    }

    /// <summary>
    /// Mutual information I(A:B) = S(A) + S(B) - S(AB).
    /// Traces from full state to AB first, then A and B from AB.
    /// </summary>
    public static double MutualInformation(Matrix<Complex> rhoFull, int nQubits,
        int[] keepA, int[] keepB)
    {
        var keepAB = keepA.Union(keepB).OrderBy(x => x).ToArray();

        // Trace full state to AB subsystem first (smaller matrix)
        var rhoAB = PartialTrace(rhoFull, nQubits, keepAB);

        // Map keepA/keepB indices into the AB subsystem's local indices
        int nAB = keepAB.Length;
        var localA = keepA.Select(q => Array.IndexOf(keepAB, q)).ToArray();
        var localB = keepB.Select(q => Array.IndexOf(keepAB, q)).ToArray();

        var rhoA = PartialTrace(rhoAB, nAB, localA);
        var rhoB = PartialTrace(rhoAB, nAB, localB);

        return VonNeumannEntropy(rhoA) + VonNeumannEntropy(rhoB)
             - VonNeumannEntropy(rhoAB);
    }

    /// <summary>
    /// Purity Tr(rho^2) = sum_ij |rho_ij|^2 (Frobenius norm squared). O(d^2).
    /// </summary>
    public static double Purity(Matrix<Complex> rho)
    {
        var vals = ((DenseMatrix)rho).Values;
        double sum = 0;
        for (int i = 0; i < vals.Length; i++)
        {
            var v = vals[i];
            sum += v.Real * v.Real + v.Imaginary * v.Imaginary;
        }
        return sum;
    }

    /// <summary>
    /// L1 coherence: sum of absolute values of off-diagonal elements.
    /// </summary>
    public static double L1Coherence(Matrix<Complex> rho)
    {
        int d = rho.RowCount;
        var vals = ((DenseMatrix)rho).Values;
        double l1 = 0;
        for (int j = 0; j < d; j++)
            for (int i = 0; i < d; i++)
            {
                if (i != j) l1 += vals[j * d + i].Magnitude;
            }
        return l1;
    }

    /// <summary>
    /// CΨ = Purity x normalized L1 coherence.
    /// </summary>
    public static (double CPsi, double C, double Psi) ComputeCPsi(Matrix<Complex> rho)
    {
        double C = Purity(rho);
        int d = rho.RowCount;
        double Psi = d > 1 ? L1Coherence(rho) / (d - 1) : 0;
        return (C * Psi, C, Psi);
    }

    public static double Theta(double cPsi)
    {
        if (cPsi > 0.25)
            return Math.Atan(Math.Sqrt(4 * cPsi - 1)) * 180 / Math.PI;
        return 0;
    }

    private static readonly Matrix<Complex> _SY2 = PauliOps.Y.KroneckerProduct(PauliOps.Y);

    /// <summary>
    /// Wootters concurrence for a 2-qubit (4x4) density matrix.
    /// </summary>
    public static double Concurrence2Q(Matrix<Complex> rho)
    {
        if (rho.RowCount != 4) throw new ArgumentException("Must be 4x4");

        var rhoTilde = _SY2 * rho.Conjugate() * _SY2;
        var R = rho * rhoTilde;

        var evd = R.Evd();
        var lambdas = evd.EigenValues
            .Select(ev => Math.Sqrt(Math.Max(ev.Real, 0)))
            .OrderByDescending(x => x)
            .ToArray();

        return Math.Max(0, lambdas[0] - lambdas[1] - lambdas[2] - lambdas[3]);
    }

    /// <summary>
    /// Expectation value Tr(rho * op) computed without full matrix multiply.
    /// </summary>
    public static Complex ExpectationValue(Matrix<Complex> rho, Matrix<Complex> op)
    {
        int d = rho.RowCount;
        var rv = ((DenseMatrix)rho).Values;
        var ov = ((DenseMatrix)op).Values;
        Complex trace = Complex.Zero;
        // Tr(A*B) = sum_ij A_ij * B_ji = sum_ij rv[j*d+i] * ov[i*d+j]
        // With column-major: A[i,j] = rv[j*d+i], B[j,i] = ov[i*d+j]
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                trace += rv[j * d + i] * ov[i * d + j];
        return trace;
    }

    /// <summary>
    /// Build a pure state density matrix |psi><psi|.
    /// </summary>
    public static Matrix<Complex> PureState(Complex[] psi)
    {
        int d = psi.Length;
        var rho = DenseMatrix.Create(d, d, Complex.Zero);
        var vals = ((DenseMatrix)rho).Values;
        // Column-major: rho[i,j] = vals[j*d+i] = psi[i] * conj(psi[j])
        for (int j = 0; j < d; j++)
        {
            var pjc = Complex.Conjugate(psi[j]);
            for (int i = 0; i < d; i++)
                vals[j * d + i] = psi[i] * pjc;
        }
        return rho;
    }

    /// <summary>
    /// Bell state fidelity Tr(rho * |b><b|) = b† rho b for a 4x4 reduced state.
    /// </summary>
    public static double BellFidelity(Matrix<Complex> rho4x4, Complex[] bellVec)
    {
        if (rho4x4.RowCount != 4 || bellVec.Length != 4)
            throw new ArgumentException("BellFidelity requires 4x4 rho and length-4 vec");
        Complex acc = Complex.Zero;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                acc += Complex.Conjugate(bellVec[i]) * rho4x4[i, j] * bellVec[j];
        return acc.Real;
    }

    /// <summary>
    /// Phase angle of the (0,3) element, normalized by pi. In [-1, 1].
    /// Returns 0 if the magnitude is below tolerance.
    /// </summary>
    public static double Ph03(Matrix<Complex> rho4x4, double tol = 1e-12)
    {
        var z = rho4x4[0, 3];
        if (z.Magnitude < tol) return 0.0;
        return Math.Atan2(z.Imaginary, z.Real) / Math.PI;
    }

    /// <summary>
    /// Extract the 9 cockpit features from a 4x4 reduced density matrix.
    /// Order: PhiPlus, PhiMinus, PsiPlus, PsiMinus, Purity, SvN, Concurrence, PsiNorm, ph03.
    /// This order must match simulations/cockpit_universality.py feat_names.
    /// </summary>
    public static double[] ExtractCockpitFeatures(Matrix<Complex> rhoPair)
    {
        if (rhoPair.RowCount != 4)
            throw new ArgumentException("ExtractCockpitFeatures requires a 4x4 reduced state");

        double phiP  = BellFidelity(rhoPair, BellStates.PhiPlus);
        double phiM  = BellFidelity(rhoPair, BellStates.PhiMinus);
        double psiP  = BellFidelity(rhoPair, BellStates.PsiPlus);
        double psiM  = BellFidelity(rhoPair, BellStates.PsiMinus);
        double pur   = Purity(rhoPair);
        double svn   = VonNeumannEntropy(rhoPair);
        double conc  = Concurrence2Q(rhoPair);
        double psiN  = L1Coherence(rhoPair) / 3.0; // d - 1 = 3
        double ph03v = Ph03(rhoPair);

        return new[] { phiP, phiM, psiP, psiM, pur, svn, conc, psiN, ph03v };
    }
}

/// <summary>
/// The four Bell state vectors as Complex[].
/// </summary>
public static class BellStates
{
    private static readonly double S = 1.0 / Math.Sqrt(2);
    public static readonly Complex[] PhiPlus  = { S, 0, 0,  S };
    public static readonly Complex[] PhiMinus = { S, 0, 0, -S };
    public static readonly Complex[] PsiPlus  = { 0, S,  S, 0 };
    public static readonly Complex[] PsiMinus = { 0, S, -S, 0 };
}
