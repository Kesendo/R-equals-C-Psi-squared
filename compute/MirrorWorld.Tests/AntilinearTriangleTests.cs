using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the antilinear triangle (F119, adopted 2026-07-04 from
// docs/proofs/PROOF_ANTILINEAR_TRIANGLE.md via the registry entry): the transpose theta, the
// entrywise conjugation conj, and the adjoint dagger = theta o conj form, with the identity, a
// Klein four-group graded by two characters -- ell (linearity) and m (multiplicativity) -- and the
// engine is the transport law mu o L_H o mu = ell*m * L_{mu(H)} for every H, Hermitian or not.
// In the Pauli basis the triangle docks onto F118's mirror group: theta = D, dagger = the
// antilinear unit K, conj = D o K, and the closure <R, D, K> is the antilinear double
// D4 x Z2 -- order 16, exactly eight antiunitary members.
public class AntilinearTriangleTests
{
    static readonly World W = new();

    // the V4 table: dagger = theta o conj = conj o theta, and every non-identity vertex squares
    // to the identity. The transport sign ell*m is -1, -1, +1 for theta, conj, dagger.
    [Fact]
    public void The_Triangle_Is_A_Klein_FourGroup_With_Product_Characters()
    {
        Assert.True(MirrorGroup.Compose(AntilinearTriangle.Theta, AntilinearTriangle.Conj)
            .Equals(AntilinearTriangle.Dagger));
        Assert.True(MirrorGroup.Compose(AntilinearTriangle.Conj, AntilinearTriangle.Theta)
            .Equals(AntilinearTriangle.Dagger));
        foreach (var v in new[] { AntilinearTriangle.Theta, AntilinearTriangle.Conj, AntilinearTriangle.Dagger })
            Assert.True(MirrorGroup.Compose(v, v).Equals(MirrorGroup.Identity));
        Assert.Equal((+1, -1), (AntilinearTriangle.Ell(AntilinearTriangle.Theta), AntilinearTriangle.Em(AntilinearTriangle.Theta)));
        Assert.Equal((-1, +1), (AntilinearTriangle.Ell(AntilinearTriangle.Conj), AntilinearTriangle.Em(AntilinearTriangle.Conj)));
        Assert.Equal((-1, -1), (AntilinearTriangle.Ell(AntilinearTriangle.Dagger), AntilinearTriangle.Em(AntilinearTriangle.Dagger)));
    }

    // on a Pauli string the odd move reads nY: theta(sigma) = conj(sigma) = (-1)^nY sigma,
    // dagger(sigma) = sigma (every Pauli string is Hermitian). All strings at N=3.
    [Fact]
    public void On_A_Pauli_String_The_Odd_Move_Reads_nY()
    {
        foreach (var m in PauliMode.Enumerate(W, 3, 0.5))
        {
            Complex expected = m.Ny % 2 == 1 ? -Complex.One : Complex.One;
            var (lt, pt) = AntilinearTriangle.Theta.Apply(m.Letters);
            var (lc, pc) = AntilinearTriangle.Conj.Apply(m.Letters);
            var (ld, pd) = AntilinearTriangle.Dagger.Apply(m.Letters);
            Assert.Equal(m.Letters, lt); Assert.Equal(expected, pt);
            Assert.Equal(m.Letters, lc); Assert.Equal(expected, pc);
            Assert.Equal(m.Letters, ld); Assert.Equal(Complex.One, pd);
        }
    }

    // the engine: mu o L_H o mu = ell(mu)*m(mu) * L_{mu(H)} for a NON-Hermitian H -- theta pays
    // -L_{H^T}, conj pays -L_{H-bar}, dagger's two signs cancel into +L_{H-dagger}. Machine-exact.
    [Fact]
    public void The_Transport_Law_Carries_The_Product_Character()
    {
        foreach (int n in new[] { 2, 3 })
            Assert.True(new AntilinearTriangle(W, n).TransportWorstResidual() < 1e-12,
                $"the transport law broke at N={n}");
    }

    // the dephasing dissipator is fixed by all three vertices (site-dependent gamma, full basis).
    [Fact]
    public void The_Dephasing_Dissipator_Is_Fixed_By_All_Three()
    {
        var triangle = new AntilinearTriangle(W, 3);
        Assert.True(triangle.DissipatorFixedWorstResidual(new[] { 0.3, 0.5, 0.7 }) < 1e-12);
    }

    // the fixed-point collapse: each vertex's fixed-point set is where the other two agree --
    // H = H-dagger iff H^T = H-bar. A Hermitian H agrees exactly; a non-Hermitian H splits at O(1).
    [Fact]
    public void The_Fixed_Point_Collapse_Ties_The_Three_Vertices()
    {
        var (hermitian, nonHermitian) = new AntilinearTriangle(W, 2).FixedPointCollapse();
        Assert.True(hermitian < 1e-12, $"H^T = H-bar must hold for Hermitian H: {hermitian:E2}");
        Assert.True(nonHermitian > 0.1, $"H^T = H-bar must fail at O(1) for non-Hermitian H: {nonHermitian:E2}");
    }

    // the antilinear double: <R, D, K> closes at order 16 with exactly eight antilinear members --
    // D4 x Z2, the eight antiunitary companions of the mirror group.
    [Fact]
    public void The_Antilinear_Double_Closes_At_Sixteen()
    {
        var members = MirrorGroup.Closure(MirrorGroup.R, MirrorGroup.D, MirrorGroup.K);
        Assert.Equal(16, members.Count);
        Assert.Equal(8, members.Count(m => m.Antilinear));
    }

    // the docking: theta IS D, dagger IS the antilinear unit K, conj = D o K -- and the subtlety
    // that keeps the double honest: conj commutes with R (conj(A F) = A-bar F), dagger does NOT
    // ((A F)-dagger = F A-dagger, the side flips).
    [Fact]
    public void The_Triangle_Docks_Onto_The_Mirror_Group()
    {
        Assert.True(AntilinearTriangle.Theta.Equals(MirrorGroup.D));
        Assert.True(AntilinearTriangle.Dagger.Equals(MirrorGroup.K));
        Assert.True(AntilinearTriangle.Conj.Equals(MirrorGroup.Compose(MirrorGroup.D, MirrorGroup.K)));
        Assert.True(MirrorGroup.Compose(AntilinearTriangle.Conj, MirrorGroup.R)
            .Equals(MirrorGroup.Compose(MirrorGroup.R, AntilinearTriangle.Conj)));
        Assert.False(MirrorGroup.Compose(AntilinearTriangle.Dagger, MirrorGroup.R)
            .Equals(MirrorGroup.Compose(MirrorGroup.R, AntilinearTriangle.Dagger)));
    }
}
