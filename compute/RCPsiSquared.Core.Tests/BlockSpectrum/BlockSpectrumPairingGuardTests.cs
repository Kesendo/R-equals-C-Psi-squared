using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.BlockSpectrum;
using RCPsiSquared.Core.Lindblad;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.SymmetryFamily;
using Xunit;
using Xunit.Abstractions;
using ComplexMatrix = MathNet.Numerics.LinearAlgebra.Matrix<System.Numerics.Complex>;

namespace RCPsiSquared.Core.Tests.BlockSpectrum;

/// <summary>When the sector pairing of <see cref="LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(ComplexMatrix, IReadOnlyList{double}, int)"/>
/// is exact, and what the engine does when it is not. The pairing fills three sectors of
/// every orbit of the F1 mirror Π on joint-popcount labels from the fourth: the Π²-image by
/// copying, the Π- and Π³-images by the reflection λ ↦ −2σ − λ. Both rules hold whenever H
/// commutes with X^⊗N (ρ ↦ ρ·X^⊗N carries L to −L† − 2σ, and ρ ↦ ρ† maps a sector to its
/// transposed label with conjugated spectrum), and, with every site dephased at a positive
/// rate, wherever the spectrum is palindromic (the engine's class summary); the engine checks
/// [H, X^⊗N] = 0 exactly (<see cref="LiouvillianBlockSpectrum.CommutesWithXN"/>) and solves
/// every sector when it fails. Among the nine rows below, a Z field on one site or on every
/// site, an XY − YX bond, and XZX + YZY or ZZZ added to the chain fail the check; of these
/// only the XY − YX bond pairs.</summary>
public class BlockSpectrumPairingGuardTests
{
    private readonly ITestOutputHelper _out;

    public BlockSpectrumPairingGuardTests(ITestOutputHelper output) => _out = output;

    private static ComplexMatrix Letter(char c) => c switch
    {
        'X' => Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, 1 }, { 1, 0 } }),
        'Y' => Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 0, -Complex.ImaginaryOne }, { Complex.ImaginaryOne, 0 } }),
        'Z' => Matrix<Complex>.Build.DenseOfArray(new Complex[,] { { 1, 0 }, { 0, -1 } }),
        _ => Matrix<Complex>.Build.DenseIdentity(2),
    };

    /// <summary>The Pauli string with <paramref name="letters"/> on the consecutive sites
    /// starting at <paramref name="start"/> and the identity elsewhere.</summary>
    private static ComplexMatrix Term(int N, string letters, int start)
    {
        ComplexMatrix m = Letter(start == 0 ? letters[0] : 'I');
        for (int site = 1; site < N; site++)
        {
            int k = site - start;
            m = m.KroneckerProduct(Letter(k >= 0 && k < letters.Length ? letters[k] : 'I'));
        }
        return m;
    }

    private static ComplexMatrix SumOver(int N, string letters, double c)
    {
        var m = Matrix<Complex>.Build.Dense(1 << N, 1 << N);
        for (int start = 0; start + letters.Length <= N; start++)
            m += c * Term(N, letters, start);
        return m;
    }

    private static ComplexMatrix XYChain(int N) => SumOver(N, "XX", 1.0) + SumOver(N, "YY", 1.0);

    /// <summary>J_b·(XX + YY + Δ·ZZ) on each bond b of the chain.</summary>
    private static ComplexMatrix BondChain(int N, double[] J, double delta)
    {
        var m = Matrix<Complex>.Build.Dense(1 << N, 1 << N);
        for (int b = 0; b < N - 1; b++)
            m += J[b] * (Term(N, "XX", b) + Term(N, "YY", b) + delta * Term(N, "ZZ", b));
        return m;
    }

    /// <summary>Popcount-conserving chains, each named by what it adds to the XY chain.</summary>
    private static ComplexMatrix Chain(string name, int N) => name switch
    {
        "XY chain" => XYChain(N),
        "Heisenberg chain" => XYChain(N) + SumOver(N, "ZZ", 1.0),
        "XXZ chain, Δ = 0.6" => XYChain(N) + SumOver(N, "ZZ", 0.6),
        "XY chain + 0.3·(XYZ − YXZ)" => XYChain(N) + SumOver(N, "XYZ", 0.3) - SumOver(N, "YXZ", 0.3),
        "XY chain + 0.7·Z on site 0" => XYChain(N) + 0.7 * Term(N, "Z", 0),
        "XY chain + 0.4·Z on every site" => XYChain(N) + SumOver(N, "Z", 0.4),
        "XY chain + 0.5·(XY − YX) on bond 0" => XYChain(N) + 0.5 * (Term(N, "XY", 0) - Term(N, "YX", 0)),
        "XY chain + 0.4·(XZX + YZY)" => XYChain(N) + SumOver(N, "XZX", 0.4) + SumOver(N, "YZY", 0.4),
        "XY chain + 0.4·ZZZ" => XYChain(N) + SumOver(N, "ZZZ", 0.4),
        _ => throw new ArgumentException(name),
    };

    // (name, H commutes with X^⊗N, the pairing rule holds sector by sector). The XY − YX bond on
    // a chain commutes with X^⊗N only in a diagonal frame (the frame theorem of
    // experiments/THE_PALINDROME_AS_A_COLOURING.md): it pairs, yet fails the check, which is
    // sufficient and not necessary, so the engine solves its sectors one by one.
    public static TheoryData<string, bool, bool> Cases => new()
    {
        { "XY chain", true, true },
        { "Heisenberg chain", true, true },
        { "XXZ chain, Δ = 0.6", true, true },
        { "XY chain + 0.3·(XYZ − YXZ)", true, true },
        { "XY chain + 0.7·Z on site 0", false, false },
        { "XY chain + 0.4·Z on every site", false, false },
        { "XY chain + 0.5·(XY − YX) on bond 0", false, true },
        { "XY chain + 0.4·(XZX + YZY)", false, false },
        { "XY chain + 0.4·ZZZ", false, false },
    };

    // No exact route connects two eigensolves, so the spectral comparisons below are read
    // against an error model: a backward-stable eigensolver returns the eigenvalues of L + δL
    // with ‖δL‖ of order eps·‖L‖_F, which moves an eigenvalue of condition κ by about
    // κ·eps·‖L‖_F. Each deviation, the largest gap in a greedy one-to-one matching of the two
    // spectra (an upper bound on the optimal matching's), is divided by eps·‖L‖_F and read at
    // three scales of (H, γ) a hundredfold apart. Rounding keeps the ratio bounded and constant
    // across the scales; the wrong spectra below put it within a decade of 1/eps, constant as
    // well; an error of fixed size
    // would grow it as the scale shrinks, which the spread across the scales reads (a ratio
    // below 1, less than one rounding, counts as 1 there).
    private const double Eps = 2.220446049250313e-16; // 2^−52, the spacing of doubles at 1
    private static readonly double[] Scales = { 1e-2, 1.0, 1e2 };
    // Measured on this machine's MKL on every case, size and scale below: 0.92 to 8.5 where the spectra agree,
    // 4.9e14 to 7.7e14 where the sectors do not pair, a spread across the scales of at most 2.5.
    private const double RoundingRatioCeiling = 100.0;
    private const double StructuralRatioFloor = 1e10;
    private const double SpreadCeiling = 10.0;

    private static double[] Rates(int N) => Enumerable.Range(0, N).Select(l => 0.3 + 0.1 * l).ToArray();

    private void AssertLaw(IReadOnlyList<double> ratios, bool agree, string context)
    {
        _out.WriteLine($"{context}: deviation / (eps·‖L‖_F) = " +
                       string.Join(" | ", ratios.Select(r => r.ToString("G4"))) + " at the three scales");
        foreach (double r in ratios)
            Assert.True(agree ? r <= RoundingRatioCeiling : r >= StructuralRatioFloor, $"{context}: ratio {r:G4}");
        double spread = ratios.Max(r => Math.Max(r, 1.0)) / ratios.Min(r => Math.Max(r, 1.0));
        _out.WriteLine($"{context}: spread {spread:G4}");
        Assert.True(spread <= SpreadCeiling, $"{context}: spread {spread:G4} across the scales");
    }

    private static double GreedyMatchDistance(IReadOnlyList<Complex> a, IReadOnlyList<Complex> b)
    {
        Assert.Equal(a.Count, b.Count);
        var taken = new bool[b.Count];
        double worst = 0.0;
        foreach (var z in a)
        {
            int best = -1;
            double bestDistance = double.MaxValue;
            for (int j = 0; j < b.Count; j++)
            {
                if (taken[j]) continue;
                double distance = (z - b[j]).Magnitude;
                if (distance < bestDistance) { bestDistance = distance; best = j; }
            }
            taken[best] = true;
            worst = Math.Max(worst, bestDistance);
        }
        return worst;
    }

    [Theory]
    [MemberData(nameof(Cases))]
    public void CommutesWithXN_DecidesExactly(string name, bool commutes, bool _)
    {
        foreach (int N in new[] { 3, 4 })
            Assert.Equal(commutes, LiouvillianBlockSpectrum.CommutesWithXN(Chain(name, N), N));
    }

    // The production chain builders hand the engines Hamiltonians that pass the check
    // exactly, so the chain runs pair; a builder that
    // drifted by one rounding would leave the spectra right and solve every sector.
    [Theory]
    [InlineData(3)]
    [InlineData(4)]
    [InlineData(7)]
    public void CommutesWithXN_HoldsForTheProductionChains(int N)
    {
        Assert.True(LiouvillianBlockSpectrum.CommutesWithXN(PauliHamiltonian.XYChain(N, 0.7).ToMatrix(), N));
        Assert.True(LiouvillianBlockSpectrum.CommutesWithXN(PauliHamiltonian.HeisenbergChain(N, 0.7).ToMatrix(), N));
    }

    // The check compares entries exactly: the hopping entry |001⟩ ↔ |010⟩ of the XY chain
    // raised by one unit in the last place, with its Hermitian partner, keeps popcount
    // conserved and leaves its X^⊗N image |110⟩ ↔ |101⟩ one unit behind, and the check refuses
    // it, where a tolerance would accept it.
    [Fact]
    public void CommutesWithXN_RefusesAOneUlpAsymmetry()
    {
        const int N = 3;
        var H = XYChain(N);
        Assert.True(LiouvillianBlockSpectrum.CommutesWithXN(H, N));
        Assert.Equal(new Complex(2.0, 0.0), H[1, 2]);
        Assert.Equal(H[1, 2], H[6, 5]);
        var raised = new Complex(Math.BitIncrement(2.0), 0.0);
        H[1, 2] = raised;
        H[2, 1] = raised;
        Assert.False(LiouvillianBlockSpectrum.CommutesWithXN(H, N));
    }

    // The engine's spectrum against the dense eigensolver of the full L, at N = 3 and 4.
    [Theory]
    [MemberData(nameof(Cases))]
    public void ComputeSpectrumPerBlock_MatchesTheFullEigensolver(string name, bool _, bool __)
    {
        foreach (int N in new[] { 3, 4 })
        {
            var ratios = new List<double>();
            foreach (double s in Scales)
            {
                var H = s * Chain(name, N);
                var gammaPerSite = Rates(N).Select(g => s * g).ToArray();
                var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);
                var full = L.Evd().EigenValues.ToArray();
                var perBlock = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N);
                ratios.Add(GreedyMatchDistance(perBlock, full) / (Eps * L.FrobeniusNorm()));
            }
            AssertLaw(ratios, agree: true, $"{name}, N={N}");
        }
    }

    // The F71-refined engine returns the spectra of the F71-diagonal sub-blocks, which are L's
    // spectrum when H and the rates are symmetric under the site reflection: these cases are.
    [Theory]
    [InlineData("XY chain")]
    [InlineData("Heisenberg chain")]
    [InlineData("XY chain + 0.4·Z on every site")]
    [InlineData("XY chain + 0.4·(XZX + YZY)")]
    [InlineData("XY chain + 0.4·ZZZ")]
    public void F71Refined_MatchesTheFullEigensolver_ForReflectionSymmetricH(string name)
    {
        foreach (int N in new[] { 3, 4 })
        {
            var ratios = new List<double>();
            foreach (double s in Scales)
            {
                var H = s * Chain(name, N);
                var gammaPerSite = Enumerable.Range(0, N)
                    .Select(l => s * (0.3 + 0.1 * Math.Min(l, N - 1 - l))).ToArray();
                var L = PauliDephasingDissipator.BuildZ(H, gammaPerSite);
                var full = L.Evd().EigenValues.ToArray();
                var perBlock = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N);
                ratios.Add(GreedyMatchDistance(perBlock, full) / (Eps * L.FrobeniusNorm()));
            }
            AssertLaw(ratios, agree: true, $"{name}, N={N}");
        }
    }

    /// <summary>The largest distance of a follower sector's spectrum from the pairing rule
    /// applied to its orbit primary's, in an output whose sectors were all solved
    /// independently. That they were is checked first: an output in which every sector equals
    /// its X^⊗N partner bit for bit is a fill (at even N the Π-fixed sector keeps the reflected
    /// relation false, so the copied one carries this guard); a fill of single sectors would
    /// pass it.</summary>
    private static double WorstFollowerDistance(Complex[] spectrum, int N, double sumGamma)
    {
        Assert.Equal((false, false), BitRelations(spectrum, N, sumGamma));
        var decomp = JointPopcountSectorBuilder.Build(N);
        var bySector = new List<Complex[]>();
        int offset = 0;
        foreach (var sector in decomp.SectorRanges)
        {
            bySector.Add(spectrum.Skip(offset).Take(sector.Size).ToArray());
            offset += sector.Size;
        }
        var (_, followers) = F1PalindromeOrbitPairing.PartitionByPiOrbit(
            N, decomp.SectorRanges, s => (s.PCol, s.PRow));
        double worst = 0.0;
        foreach (var (followerIdx, follower) in followers)
        {
            var primary = bySector[follower.PrimaryIndex];
            var expected = follower.Kind == F1PalindromeOrbitPairing.F1FollowerKind.F1Reflect
                ? primary.Select(z => new Complex(-2.0 * sumGamma - z.Real, -z.Imaginary)).ToArray()
                : primary;
            worst = Math.Max(worst, GreedyMatchDistance(bySector[followerIdx], expected));
        }
        return worst;
    }

    // The pairing rule itself, read off independent eigensolves of every sector: it holds for
    // every Hamiltonian that commutes with X^⊗N, the complex one included, and for the XY − YX
    // bond, which commutes with it in a diagonal frame, and fails for the others. The failing
    // rows are the control that the reading can fail.
    [Theory]
    [MemberData(nameof(Cases))]
    public void PairingRule_HoldsWhereTheSectorsPair(string name, bool _, bool pairs)
    {
        const int N = 4;
        var ratios = new List<double>();
        foreach (double s in Scales)
        {
            var H = s * Chain(name, N);
            var gammaPerSite = Rates(N).Select(g => s * g).ToArray();
            double norm = PauliDephasingDissipator.BuildZ(H, gammaPerSite).FrobeniusNorm();
            var spectrum = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
                H, gammaPerSite, N, LiouvillianBlockSpectrum.SectorPairing.None);
            ratios.Add(WorstFollowerDistance(spectrum, N, gammaPerSite.Sum()) / (Eps * norm));
        }
        AssertLaw(ratios, pairs, name);
    }

    // The F71-refined engine returns the spectra of the F71-diagonal sub-blocks, which are not
    // L's spectrum when H or the rates break the site reflection; they pair all the same when H
    // commutes with X^⊗N, both maps commuting with the reflection. Uneven bonds and uneven
    // rates, with a Z field as the control.
    [Theory]
    [InlineData(0.0, 0.0, true)]
    [InlineData(0.6, 0.0, true)]
    [InlineData(0.0, 0.7, false)]
    public void F71Refined_PairingRule_HoldsWithoutReflectionSymmetry(double delta, double field, bool pairs)
    {
        const int N = 4;
        var ratios = new List<double>();
        foreach (double s in Scales)
        {
            var H = s * (BondChain(N, new[] { 1.0, 0.7, 1.3 }, delta) + field * Term(N, "Z", 0));
            var gammaPerSite = Rates(N).Select(g => s * g).ToArray();
            double norm = PauliDephasingDissipator.BuildZ(H, gammaPerSite).FrobeniusNorm();
            var spectrum = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(
                H, gammaPerSite, N, LiouvillianBlockSpectrum.SectorPairing.None);
            ratios.Add(WorstFollowerDistance(spectrum, N, gammaPerSite.Sum()) / (Eps * norm));
        }
        AssertLaw(ratios, pairs, $"bonds 1.0, 0.7, 1.3, Δ = {delta}, Z field {field} on site 0");
    }

    /// <summary>Which bit-for-bit relations an output carries, over every sector: whether each
    /// sector and its Π-image are related by the reflection λ ↦ −2σ − λ, entry for entry and
    /// in order, in one direction or the other, and whether each sector equals its X^⊗N
    /// partner. Phase 3 of both engines writes a filled follower exactly so, while two sectors
    /// solved by separate eigensolves agree only to rounding. The relations are read on the
    /// output alone, whichever sector of an orbit the engine solved, and it is the pair that
    /// tells the cases apart: a single reflected follower moved by one unit in the last place
    /// can leave the first relation standing, and the second catches it.</summary>
    private static (bool Reflected, bool Copied) BitRelations(Complex[] spectrum, int N, double sumGamma)
    {
        var decomp = JointPopcountSectorBuilder.Build(N);
        var slices = new Dictionary<(int PCol, int PRow), Complex[]>();
        int offset = 0;
        foreach (var sector in decomp.SectorRanges)
        {
            slices[(sector.PCol, sector.PRow)] = spectrum.Skip(offset).Take(sector.Size).ToArray();
            offset += sector.Size;
        }
        Complex[] Reflect(Complex[] a) =>
            a.Select(z => new Complex(-2.0 * sumGamma - z.Real, -z.Imaginary)).ToArray();
        bool reflected = true, copied = true;
        foreach (var ((pCol, pRow), slice) in slices)
        {
            var image = slices[(N - pRow, pCol)];
            reflected &= image.SequenceEqual(Reflect(slice)) || slice.SequenceEqual(Reflect(image));
            copied &= slices[(N - pCol, N - pRow)].SequenceEqual(slice);
        }
        return (reflected, copied);
    }

    // The requested pairing reaches both engines: the X^⊗N copy alone fills the copies and
    // solves the reflections, no pairing solves both, and a Z field, which fails the check,
    // is solved sector by sector whatever is requested. N = 5 has no Π-fixed sector.
    [Theory]
    [InlineData("XY chain", LiouvillianBlockSpectrum.SectorPairing.PiOrbit, true, true)]
    [InlineData("XY chain", LiouvillianBlockSpectrum.SectorPairing.XNCopy, false, true)]
    [InlineData("XY chain", LiouvillianBlockSpectrum.SectorPairing.None, false, false)]
    [InlineData("XY chain + 0.7·Z on site 0", LiouvillianBlockSpectrum.SectorPairing.PiOrbit, false, false)]
    [InlineData("XY chain + 0.7·Z on site 0", LiouvillianBlockSpectrum.SectorPairing.XNCopy, false, false)]
    public void RequestedPairing_ReachesBothEngines(
        string name, LiouvillianBlockSpectrum.SectorPairing pairing, bool reflected, bool copied)
    {
        const int N = 5;
        var H = Chain(name, N);
        var gammaPerSite = Rates(N);
        double sumGamma = gammaPerSite.Sum();
        var plain = LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N, pairing);
        Assert.Equal((reflected, copied), BitRelations(plain, N, sumGamma));
        var refined = F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N, pairing);
        Assert.Equal((reflected, copied), BitRelations(refined, N, sumGamma));
    }

    // The overloads without a SectorPairing argument, the default, the eigen-path one and the
    // explicit Z dephase letter, reach the Π-orbit pairing in both engines.
    [Fact]
    public void DefaultOverloads_UseThePiOrbitPairing()
    {
        const int N = 5;
        var H = XYChain(N);
        var gammaPerSite = Rates(N);
        double sumGamma = gammaPerSite.Sum();
        var filled = (true, true);
        Assert.Equal(filled, BitRelations(
            LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N), N, sumGamma));
        Assert.Equal(filled, BitRelations(LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.Auto), N, sumGamma));
        Assert.Equal(filled, BitRelations(LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, LiouvillianBlockSpectrum.EigenPath.Auto, PauliLetter.Z), N, sumGamma));
        Assert.Equal(filled, BitRelations(
            F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N), N, sumGamma));
        Assert.Equal(filled, BitRelations(F71MirrorBlockRefinement.ComputeSpectrumPerBlock(
            H, gammaPerSite, N, PauliLetter.Z), N, sumGamma));
    }

    // A Hamiltonian that does not conserve popcount is refused by both engines, rather than
    // returned as the spectrum of L's diagonal blocks.
    [Fact]
    public void BothEngines_RefuseAHamiltonianThatDoesNotConservePopcount()
    {
        const int N = 3;
        var H = XYChain(N) + SumOver(N, "YZ", 0.5);
        var gammaPerSite = Rates(N);
        var plain = Assert.Throws<ArgumentException>(
            () => LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N));
        Assert.StartsWith("H is not popcount-conserving", plain.Message);
        var refined = Assert.Throws<ArgumentException>(
            () => F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N));
        Assert.StartsWith("H is not popcount-conserving", refined.Message);
    }

    // The popcount check is exact as well: one entry of 2^−52 between |001⟩ and |011⟩
    // (popcounts 1 and 2), the size of the residue a dense basis rotation leaves, is refused.
    [Fact]
    public void BothEngines_RefuseAPopcountResidueOfOneRounding()
    {
        const int N = 3;
        var H = XYChain(N);
        Assert.Equal(Complex.Zero, H[1, 3]);
        H[1, 3] = Eps;
        H[3, 1] = Eps;
        var gammaPerSite = Rates(N);
        var plain = Assert.Throws<ArgumentException>(
            () => LiouvillianBlockSpectrum.ComputeSpectrumPerBlock(H, gammaPerSite, N));
        Assert.StartsWith("H is not popcount-conserving", plain.Message);
        var refined = Assert.Throws<ArgumentException>(
            () => F71MirrorBlockRefinement.ComputeSpectrumPerBlock(H, gammaPerSite, N));
        Assert.StartsWith("H is not popcount-conserving", refined.Message);
    }
}
