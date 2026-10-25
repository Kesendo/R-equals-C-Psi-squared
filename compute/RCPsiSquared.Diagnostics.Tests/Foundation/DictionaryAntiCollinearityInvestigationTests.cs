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

/// <summary>The decoder's α-dictionary through the full 4^N pipeline at N = 3..6 (Symphony's ceiling), the
/// same pipeline DefectDecoder.Calibrate uses. Two channels: (1) ANTI-COLLINEARITY, the min SIGNED cos, the
/// sign+location confusion the de-loss addresses ("weakening bond a reads like strengthening bond b"); over
/// all sites it reads N=3 −0.976, N=4 −0.541, N=5 −0.965, N=6 −0.378, and over the sites the painters' own
/// reliability rule accepts N = 6 reads −0.962, so the odd/even alternation of the all-site reading is
/// carried by N = 6's two rejected sites. The single-excitation block dictionary
/// (SingleExcitationPainterDictionary, handshaken with this pipeline at N = 3..5) reads the worst signed cos
/// ≤ −0.91 at every N = 3, 5..12 over reliable sites and −0.541 at N = 4, the seesaw's anti-collinearity
/// (cos = 1 − 2w on mirror pairs; DecoderAntiCollinearitySeesawClaim). (2) COLLINEARITY, the max POSITIVE
/// cos: adjacent EDGE bonds becoming near-identical (same sign, pure location confusion), absent at N ≤ 4,
/// +0.846 at N=5 (pairs (0,1),(2,3)), +0.982 at N=6 (pairs (0,1),(3,4)).</summary>
public class DictionaryAntiCollinearityInvestigationTests
{
    private readonly ITestOutputHelper _out;
    public DictionaryAntiCollinearityInvestigationTests(ITestOutputHelper output) => _out = output;

    // The canonical PTF protocol the painters require (same as DefectDecoderTests / DefectDecoder.Calibrate).
    const double J = 1.0, Gamma = 0.05, DeltaJCal = 0.02;

    /// <summary>Build both per-bond dictionaries at N through the SAME pipeline DefectDecoder.Calibrate
    /// uses (Symphony → PaintersMovement.F / .DeviationResponse). Not routed through DefectDecoder.Calibrate
    /// so it can reach N=6 (Symphony.MaxN), past the decoder's N≤5 cap.</summary>
    static (double[][] Alpha, double[][] Dev, bool[][] Reliable) BuildDictionaries(int n)
    {
        int bonds = n - 1;
        var alpha = new double[bonds][];
        var dev = new double[bonds][];
        var reliable = new bool[bonds][];
        for (int b = 0; b < bonds; b++)
        {
            var s = new Symphony(n: n, j: J, gamma: Gamma, hType: HamiltonianType.XY,
                                 initialState: InitialStateKind.BondingMode, defectBond: b, deltaJ: DeltaJCal);
            var pm = ((IInspectable)s).Children.OfType<PaintersMovement>().Single();
            Assert.True(pm.HasLenses, $"N={n} bond {b}: painters declined ({pm.DeclineReason})");
            alpha[b] = pm.F.ToArray();
            dev[b] = pm.DeviationResponse.ToArray();
            reliable[b] = pm.Reliable.ToArray();
        }
        return (alpha, dev, reliable);
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

    /// <summary>N = 3, 4, 5 through the full pipeline (fast, ≤ 1024² eig): strong anti-collinearity at N = 3
    /// and 5, weak at N = 4, where one of the three bonds is the mirror-even centre.</summary>
    [Fact]
    public void AlphaAntiCollinearity_StrongAtThreeAndFive_WeakAtFour()
    {
        var anti = new Dictionary<int, double>();
        for (int n = 3; n <= 5; n++)
        {
            var (alpha, dev, _) = BuildDictionaries(n);
            anti[n] = Worst(alpha).MinSigned;
            Dump(n, alpha, dev);
        }
        Assert.True(anti[3] < -0.9, $"N=3 α dictionary should be strongly anti-collinear; got {anti[3]:F4}");
        Assert.True(anti[4] > -0.7, $"N=4 α dictionary should be weakly anti-collinear; got {anti[4]:F4}");
        Assert.True(anti[5] < -0.9, $"N=5 α dictionary should be strongly anti-collinear; got {anti[5]:F4}");
    }

    /// <summary>N = 6 at Symphony's ceiling (5 bonds, 4096² eig, about 30 min; tagged SLOW_DICTANTICOLLINEAR so it
    /// skips by default; run via <c>--filter "Category=SLOW_DICTANTICOLLINEAR"</c>). Over all sites the α
    /// anti-collinearity reads weak (−0.378); over the sites the painters' reliability rule accepts for every
    /// bond (four of six) the same dictionary reads strong (−0.962), so the weak reading is carried by the two
    /// rejected sites. The collinear edge-bond channel has grown to ≈ +0.98.</summary>
    [Fact]
    [Trait("Category", "SLOW_DICTANTICOLLINEAR")]
    public void N6_AllSitesWeak_ReliableSitesStrong_CollinearGrown()
    {
        var (alpha, dev, reliable) = BuildDictionaries(6);
        Dump(6, alpha, dev);
        var a = Worst(alpha);
        Assert.True(a.MinSigned > -0.7, $"N=6 over all sites: {a.MinSigned:F4}");
        Assert.True(a.MaxPos > 0.95, $"N=6 collinear edge-bond channel: {a.MaxPos:F4}");
        var keep = Enumerable.Range(0, 6).Where(i => reliable.All(r => r[i])).ToArray();
        Assert.Equal(4, keep.Length);
        var masked = alpha.Select(f => keep.Select(i => f[i]).ToArray()).ToArray();
        Assert.True(Worst(masked).MinSigned < -0.95, $"N=6 over reliable sites: {Worst(masked).MinSigned:F4}");
    }

    /// <summary>BARE-vs-PAINTED reconciliation (the key to landing without tripping). The BARE mode-basis
    /// location dictionary M[b,k] = ⟨ψ_k|V_b|ψ_1⟩ (the object <see cref="DefectReadingEquivarianceClaim"/>
    /// types: pure single-excitation algebra, NO propagation, Q-INDEPENDENT) reaches high N cheaply (no
    /// Liouvillian eig, N=7,8,9 give the discriminators the painted path cannot). Question: does the
    /// strong anti-collinearity of the painted α-dictionary live in the BARE dictionary too, or is it a
    /// Q=20 PAINTING effect? The bare mirror-pair cos = −0.33 at N=5 (weak) against the painted mirror pairs'
    /// −0.78 and −0.89 says painting. This prints the bare worst-anti and mirror-pair-anti per N.</summary>
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
