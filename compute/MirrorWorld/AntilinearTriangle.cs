using System.Numerics;

namespace MirrorWorld;

// The antilinear triangle (adopted 2026-07-04 from F119, docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md):
// the transpose theta(A) = A^T, the entrywise conjugation conj(A) = A-bar, and the adjoint
// dagger(A) = A-dagger satisfy dagger = theta o conj and form, with the identity, a Klein
// four-group graded by two +-1 characters: ell (linearity) and m (multiplicativity) --
// theta = (linear, reversing), conj = (antilinear, preserving), dagger = (antilinear, reversing).
// On a Pauli string: theta(sigma) = conj(sigma) = (-1)^nY sigma, dagger(sigma) = sigma.
//
// The engine is the TRANSPORT LAW
//
//   mu o L_H o mu = ell(mu)*m(mu) * L_{mu(H)},   L_H = -i[H, .],  any H,
//
// one sign from the -i (antilinearity), one from the commutator's order (reversal); the dephasing
// dissipator is fixed by all three. The fixed-point collapse ties the vertices: each vertex's
// fixed-point set is where the other two agree (H = H-dagger iff H^T = H-bar).
//
// In the Pauli basis the triangle docks onto the mirror group: theta = D, dagger = the antilinear
// unit K, conj = D o K, and the closure <R, D, K> is the ANTILINEAR DOUBLE, D4 x Z2 -- order 16,
// exactly eight antiunitary members. (conj commutes with R; dagger does not -- the side flips.)
public sealed class AntilinearTriangle : GameObject
{
    public int N { get; }

    public AntilinearTriangle(World world, int n) : base(world) => N = n;

    // left: what the triangle itself produces.
    public override IReadOnlyList<string> Own => new[] { "vertices", "transport", "double" };

    // the three vertices, as members of the (doubled) mirror group:
    public static MirrorGroup.Member Theta => MirrorGroup.D;                       // linear, reversing
    public static readonly MirrorGroup.Member Conj =
        MirrorGroup.Compose(MirrorGroup.D, MirrorGroup.K);                         // antilinear, preserving
    public static MirrorGroup.Member Dagger => MirrorGroup.K;                      // antilinear, reversing

    // ell: +1 linear, -1 antilinear.
    public static int Ell(MirrorGroup.Member vertex) => vertex.Antilinear ? -1 : +1;

    // m: +1 multiplicative (preserving), -1 antimultiplicative (reversing). The four V4 vertices only.
    public static int Em(MirrorGroup.Member vertex)
    {
        if (vertex.Equals(MirrorGroup.Identity) || vertex.Equals(Conj)) return +1;
        if (vertex.Equals(Theta) || vertex.Equals(Dagger)) return -1;
        throw new ArgumentException("m is graded on the triangle's four vertices only.");
    }

    // ---- the transport law on a deterministic NON-Hermitian H: for each vertex mu,
    // mu(L_H(mu(rho))) = ell*m * L_{mu(H)}(rho), machine-exact (integer entries). ----
    public double TransportWorstResidual()
    {
        int d = 1 << N;
        var h = Dense.TestMatrix(d, 3, 5);
        var rho = Dense.TestMatrix(d, 7, 2);
        var vertices = new (Func<Complex[,], Complex[,]> Mu, Complex[,] MuH, int Sign)[]
        {
            (Dense.Transpose, Dense.Transpose(h), Ell(Theta) * Em(Theta)),         // theta -> -L_{H^T}
            (Dense.Conjugate, Dense.Conjugate(h), Ell(Conj) * Em(Conj)),           // conj  -> -L_{H-bar}
            (Dense.Dagger, Dense.Dagger(h), Ell(Dagger) * Em(Dagger)),             // dagger -> +L_{H-dagger}
        };
        double worst = 0;
        foreach (var (mu, muH, sign) in vertices)
        {
            var lhs = mu(Dense.CommutatorFlow(h, mu(rho)));
            var rhs = Dense.Scale(sign, Dense.CommutatorFlow(muH, rho));
            worst = Math.Max(worst, Dense.MaxAbsDiff(lhs, rhs));
        }
        return worst;
    }

    // ---- the dephasing dissipator is fixed by all three vertices, verified on the full Pauli basis. ----
    public double DissipatorFixedWorstResidual(double[] gammas)
    {
        var zs = Enumerable.Range(0, N).Select(l => Dense.SiteZ(N, l)).ToArray();
        var vertices = new Func<Complex[,], Complex[,]>[] { Dense.Transpose, Dense.Conjugate, Dense.Dagger };
        double worst = 0;
        foreach (var letters in Dense.AllStrings(N))
        {
            var sigma = Dense.PauliString(letters);
            var ldiss = Dense.Dephasing(zs, gammas, sigma);
            foreach (var mu in vertices)
                worst = Math.Max(worst, Dense.MaxAbsDiff(mu(Dense.Dephasing(zs, gammas, mu(sigma))), ldiss));
        }
        return worst;
    }

    // ---- the fixed-point collapse: H^T = H-bar iff H = H-dagger. A Hermitian H agrees exactly;
    // the deterministic non-Hermitian H splits at O(1). ----
    public (double Hermitian, double NonHermitian) FixedPointCollapse()
    {
        int d = 1 << N;
        var m = Dense.TestMatrix(d, 3, 5);
        var hermitian = new Complex[d, d];
        var mDagger = Dense.Dagger(m);
        for (int i = 0; i < d; i++)
            for (int j = 0; j < d; j++)
                hermitian[i, j] = m[i, j] + mDagger[i, j];
        return (Dense.MaxAbsDiff(Dense.Transpose(hermitian), Dense.Conjugate(hermitian)),
                Dense.MaxAbsDiff(Dense.Transpose(m), Dense.Conjugate(m)));
    }
}
