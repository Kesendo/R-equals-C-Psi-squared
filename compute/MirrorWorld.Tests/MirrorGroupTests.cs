using System.Numerics;
using MirrorWorld;

namespace MirrorWorldTests;

// From-below guard for the mirror group (F118, adopted 2026-07-04 from
// docs/proofs/PROOF_PI_FACTORS_AS_R_TIMES_D.md via the registry entry): the canonical palindromizer
// factors, Pi_Z = R o D, and the two generators close into the dihedral D4 -- eight signed
// permutations of the Pauli basis, every phase in {1, -1, i, -i}, so every pin compares EXACTLY.
// The palindrome splits along the generators (D carries the Hamiltonian sign, R carries the
// constant), and the polarity cube's three axes are characters of two conjugations and the transpose.
public class MirrorGroupTests
{
    static readonly World W = new();

    // |<R, D>| = 8: the dihedral D4, closed by brute composition from the two generators.
    [Fact]
    public void R_And_D_Close_Into_Eight()
    {
        var eight = MirrorGroup.Closure(MirrorGroup.R, MirrorGroup.D);
        Assert.Equal(8, eight.Count);
        foreach (var named in new[] { MirrorGroup.Identity, MirrorGroup.PiZ, MirrorGroup.F, MirrorGroup.PiY,
                                      MirrorGroup.D, MirrorGroup.FD, MirrorGroup.R, MirrorGroup.FR })
            Assert.Contains(eight, m => m.Equals(named));
    }

    // Pi_Z = R o D (D applied first), and per site the April rule with no extra phase:
    // I -> X, X -> I, Y -> iZ, Z -> iY (the hard-won i falls out of Y^T = -Y meeting YX = -iZ).
    [Fact]
    public void PiZ_Is_R_After_D_And_Walks_The_April_Rule()
    {
        Assert.True(MirrorGroup.Compose(MirrorGroup.R, MirrorGroup.D).Equals(MirrorGroup.PiZ));
        Assert.Equal(('X', Complex.One), MirrorGroup.PiZ.ApplySite('I'));
        Assert.Equal(('I', Complex.One), MirrorGroup.PiZ.ApplySite('X'));
        Assert.Equal(('Z', Complex.ImaginaryOne), MirrorGroup.PiZ.ApplySite('Y'));
        Assert.Equal(('Y', Complex.ImaginaryOne), MirrorGroup.PiZ.ApplySite('Z'));
    }

    // the dihedral relations: Pi_Z has order 4 through the center F = Pi_Z^2, and the inversion
    // s r s = r^-1 in disguise: D o Pi_Z o D = Pi_Y = Pi_Z^3 = Pi_Z^-1.
    [Fact]
    public void The_Dihedral_Relations_Hold()
    {
        var r2 = MirrorGroup.Compose(MirrorGroup.PiZ, MirrorGroup.PiZ);
        Assert.True(r2.Equals(MirrorGroup.F));
        Assert.False(r2.Equals(MirrorGroup.Identity));
        Assert.True(MirrorGroup.Compose(r2, r2).Equals(MirrorGroup.Identity));
        var srs = MirrorGroup.Compose(MirrorGroup.D, MirrorGroup.Compose(MirrorGroup.PiZ, MirrorGroup.D));
        Assert.True(srs.Equals(MirrorGroup.PiY));
        Assert.True(MirrorGroup.Compose(MirrorGroup.PiZ, MirrorGroup.PiY).Equals(MirrorGroup.Identity));
    }

    // the diagonal mirrors are literally the Klein cell signs: D = diag((-1)^nY), FD = diag((-1)^nZ),
    // letters untouched -- the quartering PauliMode already carries, now as group members.
    [Fact]
    public void The_Diagonal_Mirrors_Are_The_Klein_Cell_Signs()
    {
        foreach (var m in PauliMode.Enumerate(W, 3, 0.5))
        {
            var (lettersD, phD) = MirrorGroup.D.Apply(m.Letters);
            var (lettersFD, phFD) = MirrorGroup.FD.Apply(m.Letters);
            Assert.Equal(m.Letters, lettersD);
            Assert.Equal(m.Letters, lettersFD);
            Assert.Equal(m.Ny % 2 == 1 ? -Complex.One : Complex.One, phD);
            Assert.Equal(m.Nz % 2 == 1 ? -Complex.One : Complex.One, phFD);
        }
    }

    // the spine's involution set {1, F, R, FR} is a Klein four-subgroup, and F is the center:
    // it commutes with all eight (the charge conjugation).
    [Fact]
    public void The_Spine_Involutions_Are_A_Klein_Subgroup_With_Central_F()
    {
        var klein = new[] { MirrorGroup.Identity, MirrorGroup.F, MirrorGroup.R, MirrorGroup.FR };
        foreach (var a in klein)
        {
            Assert.True(MirrorGroup.Compose(a, a).Equals(MirrorGroup.Identity));
            foreach (var b in klein)
                Assert.Contains(klein, c => c.Equals(MirrorGroup.Compose(a, b)));
        }
        foreach (var g in MirrorGroup.Closure(MirrorGroup.R, MirrorGroup.D))
            Assert.True(MirrorGroup.Compose(MirrorGroup.F, g).Equals(MirrorGroup.Compose(g, MirrorGroup.F)),
                $"F must be central; it fails against {g.Name}");
    }

    // all eight Pauli-basis forms, compared against their operator-level definitions
    // (rho^T*F, F*rho*F, rho*F, ...) on every Pauli string: signed permutations compare exactly.
    [Fact]
    public void The_Eight_Forms_Match_Their_Operators_Exactly()
    {
        foreach (int n in new[] { 2, 3 })
            Assert.True(new MirrorGroup(W, n).EightFormsWorstResidual() < 1e-12,
                $"a Pauli-basis form broke at N={n}");
    }

    // the action on a full deterministic rho: Pi_Z(rho) = rho^T * F, and the wrong-sided F * rho^T
    // is rejected at O(1) -- the side matters.
    [Fact]
    public void PiZ_Acts_As_Transpose_Times_F_And_The_Wrong_Side_Is_Rejected()
    {
        foreach (int n in new[] { 2, 3 })
        {
            var (right, wrong) = new MirrorGroup(W, n).PiZOnRho();
            Assert.True(right < 1e-12, $"rho^T * F broke at N={n}: {right:E2}");
            Assert.True(wrong > 0.1, $"F * rho^T must be rejected at O(1) at N={n}: {wrong:E2}");
        }
    }

    // the palindrome splits along the generators (XXZ Delta=0.7, site-dependent gamma, verified on
    // the full Pauli basis): D flips L_H and fixes the dissipator; R fixes L_H and reflects the
    // dissipator, carrying the entire constant -2*sigma.
    [Fact]
    public void The_Palindrome_Splits_Along_The_Generators()
    {
        var group = new MirrorGroup(W, 3);
        var (dFlipsH, dFixesDiss, rFixesH, rReflectsDiss) =
            group.PalindromeSplitResiduals(j: 1.0, delta: 0.7, gammas: new[] { 0.3, 0.5, 0.7 });
        Assert.True(dFlipsH < 1e-12, $"D L_H D = -L_H broke: {dFlipsH:E2}");
        Assert.True(dFixesDiss < 1e-12, $"D L_diss D = L_diss broke: {dFixesDiss:E2}");
        Assert.True(rFixesH < 1e-12, $"R L_H R = L_H broke: {rFixesH:E2}");
        Assert.True(rReflectsDiss < 1e-12, $"R L_diss R = -L_diss - 2 sigma broke: {rReflectsDiss:E2}");
    }

    // the cube of characters: bit_a is the character of Ad_{Z^N}, bit_b the character of
    // Ad_{X^N} = F, y_par the character of the transpose D -- and the two mixed rows
    // theta o Ad_{Z^N} = (-1)^nX, theta o Ad_{X^N} = (-1)^nZ. All 64 strings at N=3, exact.
    [Fact]
    public void The_Cube_Of_Characters_Reads_The_Three_Axes()
    {
        var thetaAdZ = MirrorGroup.Compose(MirrorGroup.D, MirrorGroup.AdZ);
        var thetaAdX = MirrorGroup.Compose(MirrorGroup.D, MirrorGroup.F);
        foreach (var m in PauliMode.Enumerate(W, 3, 0.5))
        {
            Assert.Equal(Sign((m.Nx + m.Ny) % 2), MirrorGroup.AdZ.Apply(m.Letters).Phase);
            Assert.Equal(Sign((m.Ny + m.Nz) % 2), MirrorGroup.F.Apply(m.Letters).Phase);
            Assert.Equal(Sign(m.Ny % 2), MirrorGroup.D.Apply(m.Letters).Phase);
            Assert.Equal(Sign(m.Nx % 2), thetaAdZ.Apply(m.Letters).Phase);
            Assert.Equal(Sign(m.Nz % 2), thetaAdX.Apply(m.Letters).Phase);
        }
    }

    // the truly criterion (nY even AND nZ even) is the joint-fixed cell of the diagonal mirror pair:
    // sigma is truly iff both D and FD fix it. 64/64 strings at N=3.
    [Fact]
    public void The_Truly_Cell_Is_The_Joint_Fixed_Cell_Of_The_Diagonal_Mirrors()
    {
        foreach (var m in PauliMode.Enumerate(W, 3, 0.5))
        {
            bool truly = m.Ny % 2 == 0 && m.Nz % 2 == 0;
            bool jointFixed = MirrorGroup.D.Apply(m.Letters).Phase == Complex.One
                           && MirrorGroup.FD.Apply(m.Letters).Phase == Complex.One;
            Assert.Equal(truly, jointFixed);
        }
    }

    static Complex Sign(int parity) => parity == 1 ? -Complex.One : Complex.One;
}
