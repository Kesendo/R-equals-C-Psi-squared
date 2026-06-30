using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using RCPsiSquared.Core.ChainSystems;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.States;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;
using ComplexVector = MathNet.Numerics.LinearAlgebra.Vector<System.Numerics.Complex>;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>GROUNDING SCOUT (play-session investigation, 2026-06-30): how does the decoder's α-dictionary
/// degeneracy scale with N? Measured directly across N=3..6 via the SAME Symphony/PaintersMovement
/// pipeline DefectDecoder.Calibrate uses (one source of α semantics), bypassing the decoder's N≤5 cap
/// (Symphony itself allows N≤6). Two DISTINCT degeneracy channels emerge:
///
/// <para>(1) ANTI-COLLINEARITY, the min SIGNED cos of the α-dictionary, the sign+location confusion the
/// decoder's de-loss addresses ("weakening bond a reads like strengthening bond b", cos → −1). It follows
/// the PARITY of N: strong at ODD N, weak at EVEN N. Measured: N=3 −0.976, N=4 −0.541, N=5 −0.965,
/// N=6 −0.378 (HIGH-LOW-HIGH-LOW). This FALSIFIES the interior-bond rule (which, like the deviation
/// dictionary, predicts N=5 weak); the α-dictionary at N=5 is as anti-collinear as N=3. The realizing
/// pairs are NOT mirror pairs: N=3 the single (0,1); N=5 the distance-2 (0,2)/(1,3), the (1,3) is exactly
/// the demo's "bond 3 weakened ≈ bond 1 strengthened".</para>
///
/// <para>(2) COLLINEARITY, the max POSITIVE cos: a SECOND, distinct degeneracy where adjacent EDGE bonds
/// become near-identical (same sign, pure location confusion). Absent at N≤4, +0.846 at N=5 (pairs
/// (0,1),(2,3)), +0.982 at N=6 (pairs (0,1),(3,4)). This GROWS with N rather than alternating.</para>
///
/// <para>Mechanism for the anti-collinearity parity is OPEN; the empirical pattern holds on 4 points
/// (N=3..6, Symphony's ceiling, N=7's 16384² eig is hours). The decoder-relevant quantity is the SIGNED
/// minimum cos (anti-collinearity), NOT worst |cos| (which at N=6 catches the COLLINEAR pair +0.982 and
/// conflates the two channels).</para></summary>
public class DictionaryParityInvestigationTests
{
    private readonly ITestOutputHelper _out;
    public DictionaryParityInvestigationTests(ITestOutputHelper output) => _out = output;

    // The canonical PTF protocol the painters require (same as DefectDecoderTests / DefectDecoder.Calibrate).
    const double J = 1.0, Gamma = 0.05, DeltaJCal = 0.02;

    /// <summary>Build both per-bond dictionaries at N through the SAME pipeline DefectDecoder.Calibrate
    /// uses (Symphony → PaintersMovement.F / .DeviationResponse). Not routed through DefectDecoder.Calibrate
    /// so it can reach N=6 (Symphony.MaxN), past the decoder's N≤5 cap.</summary>
    static (double[][] Alpha, double[][] Dev) BuildDictionaries(int n)
    {
        int bonds = n - 1;
        var alpha = new double[bonds][];
        var dev = new double[bonds][];
        for (int b = 0; b < bonds; b++)
        {
            var s = new Symphony(n: n, j: J, gamma: Gamma, hType: HamiltonianType.XY,
                                 initialState: InitialStateKind.BondingMode, defectBond: b, deltaJ: DeltaJCal);
            var pm = ((IInspectable)s).Children.OfType<PaintersMovement>().Single();
            Assert.True(pm.HasLenses, $"N={n} bond {b}: painters declined ({pm.DeclineReason})");
            alpha[b] = pm.F.ToArray();
            dev[b] = pm.DeviationResponse.ToArray();
        }
        return (alpha, dev);
    }

    static double Cos(double[] a, double[] b)
    {
        double dot = 0, na = 0, nb = 0;
        for (int i = 0; i < a.Length; i++) { dot += a[i] * b[i]; na += a[i] * a[i]; nb += b[i] * b[i]; }
        double d = Math.Sqrt(na * nb);
        return d > 0 ? dot / d : 0;
    }

    /// <summary>The two-channel worst-pair statistics: the most anti-collinear pair (min SIGNED cos, the
    /// sign-confusion channel) and the most collinear pair (max POSITIVE cos, the location-confusion channel).</summary>
    static (double MinSigned, (int, int) MinPair, double MaxPos, (int, int) MaxPosPair) Worst(double[][] dict)
    {
        double minSigned = double.PositiveInfinity, maxPos = double.NegativeInfinity;
        (int, int) minPair = (-1, -1), maxPosPair = (-1, -1);
        for (int i = 0; i < dict.Length; i++)
            for (int k = i + 1; k < dict.Length; k++)
            {
                double c = Cos(dict[i], dict[k]);
                if (c < minSigned) { minSigned = c; minPair = (i, k); }
                if (c > maxPos) { maxPos = c; maxPosPair = (i, k); }
            }
        return (minSigned, minPair, maxPos, maxPosPair);
    }

    void Dump(int n, double[][] alpha, double[][] dev)
    {
        var a = Worst(alpha);
        var d = Worst(dev);
        _out.WriteLine($"N={n} ({n - 1} bonds): " +
                       $"ALPHA antiCos={a.MinSigned:F4} {a.MinPair}, collinearCos={a.MaxPos:F4} {a.MaxPosPair}  ||  " +
                       $"DEV antiCos={d.MinSigned:F4} {d.MinPair}, collinearCos={d.MaxPos:F4} {d.MaxPosPair}");
        for (int i = 0; i < alpha.Length; i++)
        {
            var sb = new StringBuilder($"   α cos[{i}]: ");
            for (int k = 0; k < alpha.Length; k++) sb.Append($"{Cos(alpha[i], alpha[k]),7:F3} ");
            _out.WriteLine(sb.ToString());
        }
    }

    /// <summary>The discriminator (N=3,4,5, fast, ≤1024² eig). The α-dictionary's anti-collinearity (min
    /// signed cos) follows N-parity: strong at ODD N (3,5), weak at EVEN N (4). N=5 is the discriminator:
    /// parity predicts strong anti (≤ −0.9); the interior-bond rule predicts weak (like the deviation
    /// dictionary, which resolves N=5). Strong anti at N=5 confirms parity, falsifies interior-bond.</summary>
    [Fact]
    public void AlphaAntiCollinearity_FollowsParity_N3to5()
    {
        var anti = new Dictionary<int, double>();
        for (int n = 3; n <= 5; n++)
        {
            var (alpha, dev) = BuildDictionaries(n);
            anti[n] = Worst(alpha).MinSigned;
            Dump(n, alpha, dev);
        }

        Assert.True(anti[3] < -0.9, $"N=3 (odd) α dictionary should be strongly anti-collinear; got {anti[3]:F4}");
        Assert.True(anti[4] > -0.7, $"N=4 (even) α dictionary should be weakly anti-collinear; got {anti[4]:F4}");
        Assert.True(anti[5] < -0.9,
            $"PARITY predicts N=5 (odd) α dictionary strongly anti-collinear (the discriminator vs the " +
            $"interior-bond rule, which predicts > −0.7 like the deviation dictionary); got {anti[5]:F4}");
    }

    /// <summary>The 4th point at Symphony's ceiling N=6 (5 bonds, 4096² eig,~30 min). Confirms BOTH
    /// channels: the anti-collinearity stays WEAK (even N, parity holds, min signed cos ≈ −0.38), while the
    /// SECOND channel, collinear adjacent-edge bonds (max positive cos), has GROWN to ≈ +0.98 (a pure
    /// location degeneracy that the parity-alternating sign channel does not see). ~30 min (4096² eig ×5
    /// bonds), tagged SLOW_DICTPARITY so it skips by default; run via
    /// <c>--filter "Category=SLOW_DICTPARITY"</c>. Observed: α antiCos=−0.378, collinearCos=+0.982.</summary>
    [Fact]
    [Trait("Category", "SLOW_DICTPARITY")]
    public void N6_BothChannels_AntiWeak_CollinearGrown()
    {
        var (alpha, dev) = BuildDictionaries(6);
        Dump(6, alpha, dev);
        var a = Worst(alpha);
        Assert.True(a.MinSigned > -0.7,
            $"N=6 (even) α anti-collinearity should stay weak (parity); got {a.MinSigned:F4}");
        Assert.True(a.MaxPos > 0.95,
            $"N=6 collinear edge-bond channel should have grown to near +1; got {a.MaxPos:F4}");
    }

    /// <summary>BARE-vs-PAINTED reconciliation (the key to landing without tripping). The BARE mode-basis
    /// location dictionary M[b,k] = ⟨ψ_k|V_b|ψ_1⟩ (the object <see cref="DefectReadingEquivarianceClaim"/>
    /// types: pure single-excitation algebra, NO propagation, Q-INDEPENDENT) reaches high N cheaply (no
    /// Liouvillian eig, N=7,8,9 give the discriminators the painted path cannot). Question: does the
    /// PARITY the painted α-dictionary shows (strong anti at odd N) live in the BARE dictionary too, or is
    /// it a Q=20 PAINTING effect? The claim's note says bare mirror-pair cos = −0.33 at N=5 (weak) vs
    /// painted −0.97 (strong), which predicts the bare object does NOT show the strong odd-N parity. This
    /// prints the bare worst-anti and mirror-pair-anti per N to settle it.</summary>
    [Fact]
    public void BareModeDictionary_AntiCollinearity_AcrossN_3to9()
    {
        _out.WriteLine("BARE mode-basis dictionary M[b,k]=<psi_k|V_b|psi_1> (Q-independent, no propagation):");
        for (int n = 3; n <= 9; n++)
        {
            var M = BuildBareModeDictionary(n);
            var w = Worst(M);
            int bonds = n - 1;
            double worstMirror = double.PositiveInfinity; (int, int) mp = (-1, -1);
            for (int b = 0; b < bonds; b++)
            {
                int bm = bonds - 1 - b;
                if (bm > b) { double c = Cos(M[b], M[bm]); if (c < worstMirror) { worstMirror = c; mp = (b, bm); } }
            }
            _out.WriteLine($"N={n} ({bonds} bonds): worst antiCos={w.MinSigned:F4} {w.MinPair}, " +
                           $"mirror-pair worst antiCos={worstMirror:F4} {mp}, collinearCos={w.MaxPos:F4} {w.MaxPosPair}");
        }
    }

    /// <summary>The bare mode-basis location dictionary M[b][kIdx] = ⟨ψ_{kIdx+2}|V_b|ψ_1⟩, built exactly as
    /// <see cref="DefectReadingEquivarianceClaim"/> does (bonding-mode carrier ψ_1, modes k=2..N, bond term
    /// V_b = ½(X_bX_{b+1}+Y_bY_{b+1})). Pure algebra, no Liouvillian/propagation.</summary>
    static double[][] BuildBareModeDictionary(int n)
    {
        ComplexMatrix BondTerm(int l) =>
            0.5 * (PauliString.SiteOp(n, l, PauliLetter.X) * PauliString.SiteOp(n, l + 1, PauliLetter.X)
                 + PauliString.SiteOp(n, l, PauliLetter.Y) * PauliString.SiteOp(n, l + 1, PauliLetter.Y));
        double Element(ComplexVector bra, ComplexMatrix Vb, ComplexVector ket) =>
            bra.Conjugate().DotProduct(Vb * ket).Real;

        var psi1 = BondingMode.Build(n, 1);
        int bonds = n - 1, locModes = n - 1;   // b = 0..N−2, k = 2..N
        var M = new double[bonds][];
        for (int b = 0; b < bonds; b++)
        {
            var Vb = BondTerm(b);
            M[b] = new double[locModes];
            for (int kIdx = 0; kIdx < locModes; kIdx++)
                M[b][kIdx] = Element(BondingMode.Build(n, kIdx + 2), Vb, psi1);
        }
        return M;
    }
}
