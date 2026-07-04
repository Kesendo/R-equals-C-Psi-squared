using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the golden/metallic router (F116, adopted 2026-07-04 from
// docs/proofs/PROOF_CEILING_GOLDEN_ROUTER.md sections 1-3 and 8 via the registry entry): the two
// Z-middle ceiling cases are palindromized by a LOCAL period-4 per-site router W = (x)_l q_{l mod 4}
// with the [a,a,b,b] frame on the golden locus (a = phi X + Y, b = X - phi Y), and the whole
// construction transports to the metallic family r(c) = (c + sqrt(c^2+4))/2 on the soft line
// t2 = t3. Verified here in MirrorWorld's genre: the window lemma on the 64-dim window space
// (coefficient algebra), and the dense two-sided end-to-end W L W^-1 = -L - 2 sigma on the full
// Pauli basis (previously Python-only) -- no eigensolver anywhere.
public class RouterTests
{
    static readonly World W = new();
    static readonly double Phi = (1.0 + Math.Sqrt(5.0)) / 2.0;

    // F116's metallic mean: r(c) = (c + sqrt(c^2+4))/2, the positive root of r^2 = c r + 1.
    // Golden at c=1, the silver ratio 1+sqrt2 at c=2 (the third silver sighting after F99's 45
    // degrees), bronze at c=3, the 45-degree frame r=1 at c=0, and r(-c) = 1/r(c).
    [Fact]
    public void The_Metallic_Mean_Closes_The_Family()
    {
        Assert.Equal(Phi, Formulas.F116_MetallicMean(1.0), 12);
        Assert.Equal(1.0 + Math.Sqrt(2.0), Formulas.F116_MetallicMean(2.0), 12);
        Assert.Equal((3.0 + Math.Sqrt(13.0)) / 2.0, Formulas.F116_MetallicMean(3.0), 12);
        Assert.Equal(1.0, Formulas.F116_MetallicMean(0.0), 12);
        foreach (double c in new[] { 0.5, 1.0, 2.0, Math.PI })
        {
            double r = Formulas.F116_MetallicMean(c);
            Assert.Equal(c * r + 1.0, r * r, 10);                            // the defining quadratic
            Assert.Equal(1.0 / r, Formulas.F116_MetallicMean(-c), 10);       // the reciprocal twin
        }
    }

    // each site map is a scalar times a unitary: q_l^2 = -(1+r^2) I (golden: -(2+phi) I), and it
    // is class-swapping ({I,Z} <-> {X,Y}; the diagonal blocks vanish) -- the dissipator leg.
    [Fact]
    public void Each_Site_Map_Squares_To_The_Metallic_Scalar_And_Swaps_The_Classes()
    {
        foreach (double r in new[] { Phi, 1.0 + Math.Sqrt(2.0) })
            foreach (var q in Router.SiteMaps(r))
            {
                var sq = new Complex[4, 4];
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        for (int k = 0; k < 4; k++)
                            sq[i, j] += q[i, k] * q[k, j];
                for (int i = 0; i < 4; i++)
                    for (int j = 0; j < 4; j++)
                        Assert.True((sq[i, j] - (i == j ? -(1.0 + r * r) : 0.0)).Magnitude < 1e-12,
                            $"q^2 must be -(1+r^2) I; broke at ({i},{j})");
                foreach (var (i, j) in new[] { (0, 0), (0, 3), (3, 0), (3, 3), (1, 1), (1, 2), (2, 1), (2, 2) })
                    Assert.True(q[i, j].Magnitude < 1e-12, "the class-diagonal blocks must vanish");
            }
    }

    // the window lemma: the window-summed anticommutator {Q_k, S} vanishes at ALL FOUR offsets --
    // for the golden case, the X<->Y sibling (via the mirrored maps), and metallic c on the line.
    [Fact]
    public void The_Window_Lemma_Holds_At_All_Four_Offsets()
    {
        Assert.True(Router.WindowAnticommutatorNorm(1.0, Phi, sibling: false) < 1e-9);
        Assert.True(Router.WindowAnticommutatorNorm(1.0, Phi, sibling: true) < 1e-9);
        foreach (double c in new[] { 0.0, 0.5, 2.0, 3.0, -1.0 })
            Assert.True(Router.WindowAnticommutatorNorm(c, Formulas.F116_MetallicMean(c), sibling: false) < 1e-9,
                $"the window lemma broke at c={c}");
    }

    // the metallic locus gates the frame: off the locus (r shifted by 0.25) the window-summed
    // anticommutator is O(1), not zero -- the router exists exactly on r(c).
    [Fact]
    public void The_Metallic_Locus_Gates_The_Frame()
    {
        foreach (double c in new[] { 1.0, 2.0 })
        {
            double r = Formulas.F116_MetallicMean(c);
            Assert.True(Router.WindowAnticommutatorNorm(c, r, sibling: false) < 1e-9);
            Assert.True(Router.WindowAnticommutatorNorm(c, r + 0.25, sibling: false) > 0.1,
                $"off the locus the lemma must fail at O(1), c={c}");
        }
    }

    // the crown (previously Python-only): the dense two-sided end-to-end. With P and Q the
    // period-4 product unitaries of the proof's section 2, P H P^-1 = Q H Q^-1 = -H (the chiral
    // pair), and W L W^-1 = -L - 2 sigma holds on the FULL Pauli basis with site-dependent rates:
    // the ceiling cases are soft, and locally so.
    [Fact]
    public void The_Router_Conjugates_The_Full_Lindbladian()
    {
        var router = new Router(W);
        foreach (var (n, c) in new[] { (4, 1.0), (5, 1.0), (4, 2.0) })
        {
            var gammas = Enumerable.Range(0, n).Select(l => 0.3 + 0.15 * l).ToArray();
            var (chiralP, chiralQ, conjugation) = router.DenseResiduals(n, c, gammas);
            Assert.True(chiralP < 1e-9, $"P H P^-1 = -H broke at N={n}, c={c}: {chiralP:E2}");
            Assert.True(chiralQ < 1e-9, $"Q H Q^-1 = -H broke at N={n}, c={c}: {chiralQ:E2}");
            Assert.True(conjugation < 1e-9, $"W L W^-1 = -L - 2 sigma broke at N={n}, c={c}: {conjugation:E2}");
        }
    }
}
