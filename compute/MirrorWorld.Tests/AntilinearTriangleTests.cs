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

    // ---- one level up: the generator's transpose (PROOF_ANTILINEAR_TRIANGLE §7) ----

    // the transposed generator is the generator of the transposed system: for a NON-Hermitian H, Pauli-string
    // jumps (Z on every site, an all-X string, a string carrying Y) and, from N = 2, the phased cycle (normal,
    // neither symmetric nor unitary) at integer rates, L(E_ij)[k,l] = L'(E_kl)[i,j] with L' built from H^T
    // and the transposed jumps. Integer entries, so the comparison is exact.
    [Fact]
    public void The_Generator_Transposes_To_The_Transposed_System()
    {
        foreach (int n in new[] { 1, 2, 3 })
        {
            double residual = new AntilinearTriangle(W, n).GeneratorTransposeWorstResidual();
            Assert.True(residual == 0.0, $"the transpose law broke at N={n}: {residual:E2}");
        }
    }

    // the jumps must be TRANSPOSED on the way back: carried back unchanged, conjugated, or daggered, the
    // phased cycle breaks the law at O(1). (On Pauli strings the four maps agree up to a sign that squares
    // away, so only a jump like the cycle tells them apart.)
    [Fact]
    public void The_Jumps_Must_Be_Transposed_Not_Kept_Conjugated_Or_Daggered()
    {
        foreach (int n in new[] { 2, 3 })
        {
            var (kept, conjugated, daggered) = new AntilinearTriangle(W, n).WrongJumpMapResiduals();
            Assert.True(kept > 0.1 && conjugated > 0.1 && daggered > 0.1,
                $"N={n}: kept {kept:E2}, conjugated {conjugated:E2}, daggered {daggered:E2}");
        }
    }

    // Lindblad is a Lindbladian: the trace is preserved exactly with sigma-minus and the phased cycle among
    // the jumps (a missing 1/2 or a missing dagger would not be), and with H = 0 and Z jumps it is the
    // world's dephasing dissipator, entry for entry.
    [Fact]
    public void Lindblad_Preserves_The_Trace_And_Is_The_Dephasing_At_H_Zero()
    {
        foreach (int n in new[] { 1, 2, 3 })
        {
            var (trace, dephasing) = new AntilinearTriangle(W, n).LindbladFormResidual();
            Assert.True(trace == 0.0, $"N={n}: Tr L(rho) = {trace:E2}");
            Assert.True(dephasing == 0.0, $"N={n}: Lindblad vs Dephasing {dephasing:E2}");
        }
    }

    // the law asks L to be unital; amplitude damping (sigma-minus alone) is not, and breaks it at O(1),
    // even with the jump transposed on the way back.
    [Fact]
    public void A_NonUnital_Jump_Breaks_The_Transpose_Law()
    {
        double residual = new AntilinearTriangle(W, 2).NonUnitalJumpTransposeResidual();
        Assert.True(residual > 0.1, $"sigma-minus must break the transpose law: {residual:E2}");
    }

    // unitality, not normality: sigma-minus and sigma-plus (neither normal) at equal rates keep L unital and
    // the law exact; at unequal rates it breaks.
    [Fact]
    public void Unitality_Not_Normality_Is_The_Condition()
    {
        foreach (int n in new[] { 1, 2 })
        {
            var (equal, unequal) = new AntilinearTriangle(W, n).UnitalPairTransposeResidual();
            Assert.True(equal == 0.0, $"N={n}: equal rates must keep the law: {equal:E2}");
            Assert.True(unequal > 0.1, $"N={n}: unequal rates must break it: {unequal:E2}");
        }
    }

    // a real H (here generic: a dense integer symmetric matrix, not number-conserving)
    // under Z dephasing gives a complex-symmetric generator, M^T = M exactly. A flux through a ring (the
    // XY - YX term on the wrap bond, H Hermitian but not real) does not. The rates carry no teeth here: Z
    // dephasing is diagonal in the matrix-unit basis, so no rate profile can break the symmetry; H can.
    [Fact]
    public void A_Real_Hamiltonian_Makes_The_Generator_Complex_Symmetric()
    {
        var rates = new[] { 1.0, 2.0, 3.0, 5.0 };
        var real = AntilinearTriangle.RealSymmetricTestMatrix(16, 3, 5);
        var zs = AntilinearTriangle.SiteZs(4);
        Func<Complex[,], Complex[,]> l = r => AntilinearTriangle.Lindblad(real, zs, rates, r);
        double residual = AntilinearTriangle.TransposeResidual(l, l, 16);
        Assert.True(residual == 0.0, $"a real H must give M^T = M: {residual:E2}");

        var ring = AntilinearTriangle.FluxRing(3, +1);
        var zs3 = AntilinearTriangle.SiteZs(3);
        Func<Complex[,], Complex[,]> lr = r => AntilinearTriangle.Lindblad(ring, zs3, rates[..3], r);
        Assert.True(AntilinearTriangle.TransposeResidual(lr, lr, 8) > 0.1, "a flux ring must break M^T = M");
    }

    // reciprocity, order by order in t: Tr(B L^k(A)) = Tr(A L^k(B)) for k = 0..4 on a real H, for real
    // symmetric A and B (preparation and readout trade places exactly; the unequal rates cannot break it,
    // Z dephasing being diagonal in the matrix-unit basis). The condition on A is load-bearing: a Hermitian
    // A that is not real (a |+i>-type coherence) breaks it.
    [Fact]
    public void Preparation_And_Readout_Trade_Places_On_A_Real_Hamiltonian()
    {
        var rates = new[] { 1.0, 2.0, 3.0 };
        var h = AntilinearTriangle.RealSymmetricTestMatrix(8, 2, 3);
        var zs = AntilinearTriangle.SiteZs(3);
        Func<Complex[,], Complex[,]> l = r => AntilinearTriangle.Lindblad(h, zs, rates, r);
        var a = AntilinearTriangle.RealSymmetricTestMatrix(8, 1, 4);
        var b = AntilinearTriangle.RealSymmetricTestMatrix(8, 4, 1);
        double residual = AntilinearTriangle.ReciprocityResidual(l, l, a, b, 4);
        Assert.True(residual == 0.0, $"prepare A, read B must equal prepare B, read A on a real H: {residual:E2}");

        var hermitianNotReal = AntilinearTriangle.HermitianTestMatrix(8, 3, 2);
        Assert.True(AntilinearTriangle.ReciprocityResidual(l, l, hermitianNotReal, b, 4) > 0.1,
            "a Hermitian but non-real payload must break the swap");
    }

    // reciprocity swaps OPERATORS, not positions: on a generic real H, "prepare site 0 over a quiet rest,
    // read site N-1 ignoring the rest" against its positional mirror differs already in a linear readout,
    // and against its operator swap (prepare what was read, read what was prepared) agrees exactly.
    [Fact]
    public void Reciprocity_Swaps_Operators_Not_Positions()
    {
        foreach (int n in new[] { 2, 3 })
        {
            var (positional, swap) = new AntilinearTriangle(W, n).BackgroundSwapResiduals(4);
            Assert.True(swap == 0.0, $"N={n}: the operator swap must agree exactly: {swap:E2}");
            Assert.True(positional > 0.1, $"N={n}: the positional mirror must differ: {positional:E2}");
        }
    }

    // ...and when positions DO trade: on a number-conserving real chain (uneven Heisenberg and XXZ, unequal
    // rates) one excitation sent from end site to end site over the rest in |0>, or a |+>, is its own
    // operator swap, so the positional reversal agrees exactly (to order t^10; the excitation's signal first
    // appears at t^8, four hops each way). Read at an interior site it differs. A Bell pair is not its own
    // swap (the far pair's |00> population collects the preparation's own return probability): at N = 4 it
    // agrees only through the spin flip X^N, and an uneven longitudinal field breaks that; at N = 5 it
    // differs (first at t^6); on the XX chain it agrees through the free-fermion closure of the one-body
    // correlations, field or not; on the XX star the Jordan-Wigner strings break the closure.
    [Fact]
    public void Inside_One_Excitation_Block_Positions_Do_Trade()
    {
        var triangle = new AntilinearTriangle(W, 5);
        foreach (int delta in new[] { 1, 2 })
        {
            var (single, interior, plus, bell) = triangle.NumberConservingPositionalResiduals(delta, 10);
            Assert.True(single == 0.0, $"delta={delta}: one excitation must travel end to end the same both ways: {single:E2}");
            Assert.True(interior > 0.1, $"delta={delta}: an interior readout must differ under the positional mirror: {interior:E2}");
            Assert.True(plus == 0.0, $"delta={delta}: a |+> must travel the same both ways: {plus:E2}");
            Assert.True(bell > 0.1, $"delta={delta}: a Bell pair over a spectator must differ between the ends: {bell:E2}");
        }
        var (_, _, _, bellXxChain) = triangle.NumberConservingPositionalResiduals(0, 10);
        Assert.True(bellXxChain == 0.0, $"on the XX chain the Bell pair agrees through t^10: {bellXxChain:E2}");
        var (_, _, _, bellXxField) = triangle.NumberConservingPositionalResiduals(0, 10, field: true);
        Assert.True(bellXxField == 0.0, $"on the XX chain with a field the Bell pair still agrees: {bellXxField:E2}");

        var four = new AntilinearTriangle(W, 4);
        var (_, _, _, bellFlip) = four.NumberConservingPositionalResiduals(1, 10);
        Assert.True(bellFlip == 0.0, $"N=4: the spin flip carries |1100> to |0011>, the Bell pair agrees: {bellFlip:E2}");
        var (_, _, _, bellFieldBroken) = four.NumberConservingPositionalResiduals(1, 10, field: true);
        Assert.True(bellFieldBroken > 0.1, $"N=4: a longitudinal field breaks the spin flip and the Bell pair differs: {bellFieldBroken:E2}");
        var (singleStar, _, _, bellXxStar) = triangle.NumberConservingPositionalResiduals(0, 10, star: true);
        Assert.True(singleStar == 0.0, $"on the star one excitation still travels leaf to leaf the same both ways: {singleStar:E2}");
        Assert.True(bellXxStar > 0.1, $"on the XX star the Bell pair must differ: {bellXxStar:E2}");
    }

    // the silence reading: a Bell pair sent from one pair to another over the rest in |0>, the receiving pair
    // read with the rest projected on |0...0>, gives the same pair state entry for entry whichever pair sends,
    // so every function of it agrees, CPsi included; read by habit (the rest traced out) the two directions
    // differ. Order by order through t^12 on the uneven Heisenberg and XX chains and on a star: the |11>
    // transfer first appears at t^12 on the chain (six hops on each side, pure hopping at that order, so
    // there the anisotropy and the rates enter only through the coherence) and at t^8 on the star, the
    // vacuum-to-|11> coherence at t^6 and t^4. Every entry is a dyadic number whose odd part stays far
    // below 2^53 at these orders, so the comparison is exact.
    [Fact]
    public void Read_With_The_Silence_The_Channel_Is_The_Same_From_Both_Ends()
    {
        var triangle = new AntilinearTriangle(W, 5);
        foreach (var (delta, star) in new[] { (1, false), (0, false), (1, true) })
        {
            var (silence, habit) = triangle.SilenceReadingResiduals(delta, 12, star);
            Assert.True(silence == 0.0, $"delta={delta}, star={star}: the silence-read pair must agree: {silence:E2}");
            Assert.True(habit > 0.1, $"delta={delta}, star={star}: the habit-read pair must differ: {habit:E2}");
        }
    }

    // with a flux the swap lands on the ring with the flux reversed (H^T = H-bar): exact against the
    // reversed ring, broken against the same ring.
    [Fact]
    public void On_A_Flux_Ring_The_Way_Back_Reverses_The_Flux()
    {
        var rates = new[] { 1.0, 2.0, 3.0 };
        var zs = AntilinearTriangle.SiteZs(3);
        var forward = AntilinearTriangle.FluxRing(3, +1);
        var reversed = AntilinearTriangle.FluxRing(3, -1);
        Assert.True(AntilinearTriangle.TransposeEqualsResidual(forward, reversed) == 0.0, "H^T must be the reversed ring");
        Func<Complex[,], Complex[,]> lf = r => AntilinearTriangle.Lindblad(forward, zs, rates, r);
        Func<Complex[,], Complex[,]> lb = r => AntilinearTriangle.Lindblad(reversed, zs, rates, r);
        var a = AntilinearTriangle.RealSymmetricTestMatrix(8, 1, 4);
        var b = AntilinearTriangle.RealSymmetricTestMatrix(8, 4, 1);
        double back = AntilinearTriangle.ReciprocityResidual(lf, lb, a, b, 4);
        Assert.True(back == 0.0, $"B->A must equal A->B with the flux reversed: {back:E2}");
        Assert.True(AntilinearTriangle.ReciprocityResidual(lf, lf, a, b, 4) > 0.1,
            "on the same flux ring the two directions must differ");
    }
}
