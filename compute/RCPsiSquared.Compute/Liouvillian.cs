using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Complex;

namespace RCPsiSquared.Compute;

/// <summary>
/// Constructs the Lindblad Liouvillian superoperator.
/// L(rho) = -i[H,rho] + sum_k (L_k rho L_k^dag - 1/2 {L_k^dag L_k, rho})
/// In vectorized form: L_mat |rho>> where |rho>> = vec(rho).
/// </summary>
public static class Liouvillian
{
    /// <summary>
    /// Build the full (dense) Liouvillian superoperator.
    /// Size: d^2 x d^2 where d = 2^nQubits.
    /// </summary>
    public static Matrix<Complex> Build(int nQubits, Bond[] bonds, double[] gammaPerQubit)
    {
        int d = 1 << nQubits;
        int d2 = d * d;

        var H = Topology.BuildHamiltonian(nQubits, bonds);
        var Id = DenseMatrix.CreateIdentity(d).ToComplex();

        // Hamiltonian part: -i(H⊗I - I⊗H^T)
        var L = new Complex(0, -1) * (H.KroneckerProduct(Id) - Id.KroneckerProduct(H.Transpose()));

        // Dissipator: sum_k [ L_k⊗L_k* - 1/2(L_k^dag*L_k⊗I + I⊗L_k^T*L_k*) ]
        for (int k = 0; k < nQubits; k++)
        {
            var Lk = Math.Sqrt(gammaPerQubit[k]) * PauliOps.At(PauliOps.Z, k, nQubits);
            var LkDag = Lk.ConjugateTranspose();
            var LdL = LkDag * Lk;

            L += Lk.KroneckerProduct(Lk.Conjugate())
               - 0.5 * (LdL.KroneckerProduct(Id) + Id.KroneckerProduct(LdL.Transpose()));
        }

        return L;
    }

    /// <summary>
    /// Extract oscillatory decay rates from Liouvillian eigenvalues.
    /// Returns sorted list of decay rates (positive, in units of gamma).
    /// </summary>
    public static List<double> GetOscillatoryRates(Matrix<Complex> L, double threshold = 0.05)
    {
        var evals = L.Evd().EigenValues;
        var rates = new List<double>();

        foreach (var ev in evals)
        {
            if (Math.Abs(ev.Imaginary) > threshold)
            {
                double rate = -ev.Real;
                if (rate > 0.0001)
                    rates.Add(Math.Round(rate, 6));
            }
        }

        rates.Sort();
        return rates;
    }
}
