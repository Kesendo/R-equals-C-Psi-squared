using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.CoherenceBlocks;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The spectator W(ρ) = Σ_l c_l†ρc_l (JW strings included, <see cref="SpectatorIntertwiner.BuildW"/>)
/// on the XY RING, where H is not one quadratic form in the fermions. With hopping matrices h_ket, h_bra held
/// fixed on the two sides, pushing W through the two-sided commutator leaves Σ_{l,m} (h_ket − h_bra)_{ml}·c_m†ρc_l.
/// On the ring the wrap bond's Jordan-Wigner sign follows the parity of the particle number it acts on, so the
/// ket and the bra see one matrix when p ≡ q (mod 2), the coherence's height k_Z = n_X + n_Y being even; and W,
/// adding one fermion to each side, flips that shared sign. Hence, for p ≡ q: W carries the ring's (p,q) block
/// onto the (p+1,q+1) block of the ring with its wrap bond negated (the π-flux ring), and W² carries it back
/// onto the ring's own (p+2,q+2) block. For p ≢ q (the Δn = 1 blocks, the floor's, among them) neither holds,
/// with one exception, also gated: at even N the blocks (0, N−1) and (N−1, 0), whose empty starting side and
/// full arriving side impose no hopping matrix, go under W onto the ring itself.
///
/// <para>Exact route (no tolerance): integer bonds, integer per-site γ, W entries ±1, so every entry of
/// L₂W − WL₁ is an integer computed exactly in double; the zero residuals are compared with == 0.0 and the
/// nonzero ones are asserted nonzero. Profiles are non-uniform and non-palindromic so that no reflection
/// symmetry can hide a convention error; rows with a negative wrap coupling start from the π-flux ring; the open chain is the control (Theorem B of
/// <c>docs/proofs/PROOF_CODIM1_BY_ADDITIVITY.md</c>, W exact on every block).</para>
///
/// <para>The adjoint step on the diagonal, Ψ = W† from the rung p of one spin ring to the rung p − 1 of the
/// other, the descent whose kernel is F170's lowest-weight space, is checked live in
/// <see cref="RingFrozenBandWitness"/>
/// (<c>inspect --root ringband</c>).</para>
///
/// <para>Conventions as <see cref="SpectatorIntertwinerGateTests"/>: site l = bit N−1−l, block flat index
/// pIdx·Mq + qIdx over <see cref="BlockBasis.PopcountStates"/>, Liouville flat = ket·2^N + bra,
/// H = Σ_b J_b (σ⁺_iσ⁻_j + h.c.) built here in the spin basis (no strings), so one excitation hops with J_b.</para></summary>
public class SpectatorRingParityGateTests
{
    private readonly ITestOutputHelper _out;
    public SpectatorRingParityGateTests(ITestOutputHelper o) => _out = o;

    private static ComplexMatrix SpinXY(int n, IEnumerable<(int I, int J, int Coupling)> bonds)
    {
        int d = 1 << n;
        var h = Matrix<Complex>.Build.Dense(d, d);
        var list = bonds.ToList();
        for (long s = 0; s < d; s++)
            foreach (var (i, j, c) in list)
                foreach (var (x, y) in new[] { (i, j), (j, i) })
                    if (SpectatorIntertwiner.Occupied(n, s, y) && !SpectatorIntertwiner.Occupied(n, s, x))
                    {
                        long t = s ^ SpectatorIntertwiner.SiteMask(n, x) ^ SpectatorIntertwiner.SiteMask(n, y);
                        h[(int)t, (int)s] += c;
                    }
        return h;
    }

    private static ComplexMatrix Block(int n, ComplexMatrix h, IReadOnlyList<double> gamma, int p, int q)
    {
        int d = 1 << n;
        var sp = BlockBasis.PopcountStates(n, p);
        var sq = BlockBasis.PopcountStates(n, q);
        var flat = new int[sp.Count * sq.Count];
        for (int i = 0; i < sp.Count; i++)
            for (int j = 0; j < sq.Count; j++)
                flat[i * sq.Count + j] = (int)(sp[i] * d + sq[j]);
        return PerBlockLiouvillianBuilder.BuildBlockZ(h, gamma, flat);
    }

    private static double MaxAbs(ComplexMatrix m)
    {
        double r = 0;
        for (int i = 0; i < m.RowCount; i++)
            for (int j = 0; j < m.ColumnCount; j++)
                r = Math.Max(r, m[i, j].Magnitude);
        return r;
    }

    [Theory(DisplayName = "Ring parity: for p ≡ q, W carries the ring onto the π-flux ring and W² onto itself, exactly; for p ≢ q neither")]
    [InlineData(4, new[] { 2, 3, 1, 4 }, new[] { 1, 3, 2, 2 })]
    [InlineData(5, new[] { 1, 2, 4, 3, 5 }, new[] { 1, 2, 3, 1, 3 })]
    [InlineData(6, new[] { 5, 1, 2, 3, 1, 4 }, new[] { 3, 2, 1, 1, 2, 1 })]
    [InlineData(7, new[] { 2, 4, 1, 3, 5, 1, 3 }, new[] { 2, 1, 3, 1, 2, 3, 1 })]
    // A negative wrap coupling makes the base the π-flux ring, so these rows gate the step from the π-flux ring
    // onto the ring as a source in its own right (the even rungs of experiments/XY_FROZEN_BAND.md start there).
    [InlineData(5, new[] { 1, 2, 4, 3, -5 }, new[] { 1, 2, 3, 1, 3 })]
    [InlineData(6, new[] { 5, 1, 2, 3, 1, -4 }, new[] { 3, 2, 1, 1, 2, 1 })]
    public void Ring_W_And_WSquared_ByHeightParity(int n, int[] j, int[] g)
    {
        var gamma = g.Select(v => (double)v).ToList();
        var chainBonds = Enumerable.Range(0, n - 1).Select(l => (l, l + 1, j[l])).ToList();
        var hChain = SpinXY(n, chainBonds);
        var hRing = SpinXY(n, chainBonds.Append((n - 1, 0, j[n - 1])));
        var hPiRing = SpinXY(n, chainBonds.Append((n - 1, 0, -j[n - 1])));

        int checkedEven = 0, checkedOdd = 0;
        for (int p = 0; p + 1 <= n; p++)
            for (int q = 0; q + 1 <= n; q++)
            {
                var w = SpectatorIntertwiner.BuildW(n, p, q);
                double chain = MaxAbs(Block(n, hChain, gamma, p + 1, q + 1) * w - w * Block(n, hChain, gamma, p, q));
                Assert.True(chain == 0.0, $"control: chain W residual {chain:R} at ({p},{q}), N={n}");

                var lRing = Block(n, hRing, gamma, p, q);
                double toPi = MaxAbs(Block(n, hPiRing, gamma, p + 1, q + 1) * w - w * lRing);
                bool even = (p - q) % 2 == 0;
                string line = $"N={n} ({p},{q}) {(even ? "even" : "odd ")}: chain {chain:R}, W ring→π-ring {toPi:R}";

                if (p + 2 <= n && q + 2 <= n)
                {
                    var w2 = SpectatorIntertwiner.BuildW(n, p + 1, q + 1) * w;
                    double sq = MaxAbs(Block(n, hRing, gamma, p + 2, q + 2) * w2 - w2 * lRing);
                    line += $", W² ring→ring {sq:R}";
                    if (even) Assert.True(sq == 0.0, $"W² ring residual {sq:R} at even ({p},{q}), N={n}");
                    else Assert.True(sq > 0.0, $"W² ring residual vanished at odd ({p},{q}), N={n}: the gate cannot fail");
                }
                if (even) { Assert.True(toPi == 0.0, $"W ring→π-ring residual {toPi:R} at even ({p},{q}), N={n}"); checkedEven++; }
                else
                {
                    Assert.True(toPi > 0.0, $"W ring→π-ring residual vanished at odd ({p},{q}), N={n}");
                    double toRing = MaxAbs(Block(n, hRing, gamma, p + 1, q + 1) * w - w * lRing);
                    bool edgePair = n % 2 == 0 && Math.Min(p, q) == 0 && Math.Max(p, q) == n - 1;
                    line += $", W ring→ring {toRing:R}";
                    if (edgePair) Assert.True(toRing == 0.0, $"edge pair ({p},{q}) at even N={n}: W ring→ring residual {toRing:R}");
                    else Assert.True(toRing > 0.0, $"W ring→ring residual vanished at odd ({p},{q}), N={n}: an exception beyond the edge pair");
                    checkedOdd++;
                }
                _out.WriteLine(line);
            }
        Assert.True(checkedEven > 0 && checkedOdd > 0);
    }
}
