using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using RCPsiSquared.Core.Numerics;
using RCPsiSquared.Core.Pauli;
using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Diagnostics.F87;
using RCPsiSquared.Diagnostics.Foundation;
using Xunit;

namespace RCPsiSquared.Diagnostics.Tests.Foundation;

/// <summary>The live battery of F173: the claim's GF(2) readings agree with Theorem 2 of the complement connection,
/// walked on the bitstrings (Theorems D and B′ also against the balance), on every family the witness reads, and the
/// eighteen controls (the pieces of Theorem B swapped, the power dropped from Theorem C, Theorem D with the parity of N or
/// with i₀ at the greatest valuation, Theorem E's conditions each dropped, the tests of Corollaries E2 and E3 each
/// dropped, Theorem B′ read off the unreduced basis, the ℤ/4 carrier without its coset term, Corollary E3′ without N₁* or without its ordering term,
/// Theorem H with every rate read as one or with the rates in reverse site order) disagree, so the battery can fail.</summary>
public class SignedMaskBalanceWitnessTests
{
    [Fact]
    public void EveryFamilyAgreesWithTheoremTwo_AndTheControlsDisagree()
    {
        var w = new SignedMaskBalanceWitness();
        Assert.All(w.Families, f => Assert.True(f.Agrees, $"{f.Name}: {f.FirstMismatch}"));
        Assert.True(w.AllAgree);
        Assert.True(w.ControlsBreak, $"controls {w.SwappedWidthsMismatches}, {w.SwappedQuotientsMismatches}, {w.DroppedRingPowerMismatches}, " +
                                     $"{w.ParityOnlyMismatches}, {w.MostValuedMismatches}, {string.Join(",", w.DroppedConditionMismatches)}, {w.DroppedReadSetMismatches}, {w.DroppedOrderingMismatches}, {w.DroppedReadSumMismatches}, {w.DroppedE3OrderingMismatches}, {w.UnreducedBasisMismatches}, {w.DroppedCosetTermMismatches}, {w.DroppedE3PrimeReadSumMismatches}, {w.DroppedE3PrimeOrderingMismatches}, {w.RatesReadAwayMismatches}, {w.RatesReversedMismatches}, {w.PackingOverlapMismatches}, {w.CutWhereChainBreaksMismatches}, {w.CycleDoubleCountMismatches}");
        Assert.Equal(new[] { (864, 3024), (368, 1036), (19, 32), (49, 290), (1, 8), (1694, 9140), (1046, 28232), (416, 10418), (20, 100), (600, 9901), (600, 4800), (600, 600), (648, 648), (428, 3936), (16, 128), (428, 3936), (17, 108), (352, 2176), (352, 10560), (128, 768), (928, 1856), (32, 224), (1088, 2848), (500, 3000), (200, 1200), (80, 80), (864, 6912), (1292, 10408), (10, 30), (428, 192), (864, 15552), (432, 4528), (864, 3456), (24, 24) },
            w.Families.Select(f => (f.Cells, f.Readings)).ToArray());                                // measured
        Assert.Equal(new Dictionary<int, int> { [5] = 6, [6] = 8 }, w.F115FloorsAtK4);
        Assert.Equal(new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>                     // measured, quoted on the page
            {
                [SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail] = 928, [SignedMaskBalanceClaim.ChainShapeBranch.NoRelation] = 192,
                [SignedMaskBalanceClaim.ChainShapeBranch.ReadSumNonzero] = 928, [SignedMaskBalanceClaim.ChainShapeBranch.HolonomyPlus] = 400,
                [SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus] = 400,
            }, w.E3Branches);
        Assert.Equal((464, 576), (w.DroppedReadSumMismatches, w.DroppedE3OrderingMismatches));
        Assert.Equal((301, 468, 473), (w.UnreducedBasisMismatches, w.BPrimeBreaking, w.UnreducedBasisDeficient));   // measured
        Assert.Equal((1784, 192, 74), (w.Z4Paired, w.DroppedCosetTermMismatches, w.Z4FourStringPaired));
        Assert.Equal((4, 6), (w.DroppedE3PrimeReadSumMismatches, w.DroppedE3PrimeOrderingMismatches));            // measured
        Assert.Equal(new Dictionary<SignedMaskBalanceClaim.ChainShapeBranch, int>                     // measured, quoted on the page
            {
                [SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail] = 2566, [SignedMaskBalanceClaim.ChainShapeBranch.NoRelation] = 265,
                [SignedMaskBalanceClaim.ChainShapeBranch.ReadSumNonzero] = 9, [SignedMaskBalanceClaim.ChainShapeBranch.HolonomyPlus] = 398,
                [SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus] = 962,
            }, w.E3PrimeBranches);                       // measured
        Assert.Equal((2768, 2768, 1028), (w.StrongDephasingBroken, w.RatesReadAwayMismatches, w.RatesReversedMismatches));   // measured, quoted on the page
        Assert.True(w.StrongDephasingWorstRatio < 1, $"{w.StrongDephasingWorstRatio:F3} of dim·eps·‖𝓜‖");
        Assert.True(w.KernelZeroMax < 1e-13 && w.KernelGapMin > 1e-3, $"zeros to {w.KernelZeroMax:E2}, the least nonzero {w.KernelGapMin:E2}");   // the count's 10⁻⁹ between
        Assert.True(w.MostNegative > -1e-13, $"an eigenvalue of 𝓛 or 𝓜 at {w.MostNegative:E2}·‖·‖");                // positive semidefinite to rounding
        Assert.True(w.TrellisWorstRatio < 1, $"the trellis against the enumeration at {w.TrellisWorstRatio:F4} of 2N·eps");   // 0.138 measured
        Assert.Equal((588, 3752, 0, 2390), (w.PackingEqualities, w.PackingOverlapMismatches, w.CutEqualities, w.CutWhereChainBreaksMismatches));   // measured, quoted on the page
        // the cut family is not empty of breaks: rings broken while the chain pairs, λ up to 0.79 of twice the cut, measured
        Assert.Equal(128, w.CutBroken);
        Assert.Equal((216, 192, 0, 20), (w.CycleRelations, w.CyclePackingReadings, w.CyclePackingEqualities, w.CycleDoubleCountMismatches));   // measured, quoted on the page
        // the packing at most 0.82 of λ_min(𝓜) on Theorem E's cells (0.811 measured); the triangle with a potential: 27 readings
        // with moduli that differ, 18 of them with (3/2)·λ_* above 2·μ_min, which no cycle without a potential reaches
        Assert.True(w.CycleTightest > 0.8 && w.CycleTightest < 0.82, $"the packing at {w.CycleTightest:F4} of λ_min(𝓜)");
        Assert.Equal((27, 18, 0), (w.CyclePotentialReadings, w.CyclePotentialAboveTwoHops, w.CycleEqualModuliAboveTwoHops));
        Assert.True(w.CutTightest > 0.75 && w.CutTightest <= 1, $"λ at most {w.CutTightest:F4} of twice the cut");
    }

    // Corollary H2 by hand at rate 1 and coefficient 1: on XIY + XZX on the ring of six each triangle of flip sets {l, l + 2}
    // is a 3-cycle of hops of modulus 1 with holonomy −1, so λ_* = 2(1 − cos(π/3)) = 1, the bound (3/2)·1 per triangle and
    // 3 for the two, below λ_min(𝓜) = 4; on XZ + ZY on the chain of four the two end squares, walks a, b, a, b with s = 4,
    // give λ_min(𝓜) = 2(3 − √5) exactly; a walk that revisits a bitstring reads null.
    [Fact]
    public void CorollaryH2ReadsThePagesCycles()
    {
        int n = 6;
        var gamma = Enumerable.Repeat(1.0, n).ToArray();
        var placed = SignedMaskBalanceClaim.PlaceStrings(new[] { ("XIY", 1L), ("XZX", 1L) }, n, ring: true);
        var h = ComplementConnectionGraph.Columns(n, placed.Select(p => (new PauliMask(p.X, p.Z), (System.Numerics.BigInteger)p.C)));
        double lam = StrongDephasingEdges.Build(n, h, gamma).Far.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian).EigenValues.Min(z => z.Real);
        Assert.Equal(4.0, lam, 12);
        foreach (var tri in new[] { new[] { 0b000101, 0b010100, 0b010001 }, new[] { 0b001010, 0b101000, 0b100010 } })
        {
            var r = StrongDephasingEdges.CycleReading(n, h, gamma, tri)!.Value;
            Assert.Equal(1.0, r.LambdaStar, 12);
            Assert.Equal(1.0, r.HopMin, 12);
            Assert.Equal(System.Math.PI, r.ThetaMin, 12);
        }
        Assert.Null(StrongDephasingEdges.CycleReading(n, h, gamma, new[] { 0b000101, 0b000101, 0b000101, 0b000101 }));   // revisits

        int m = 4;
        var g4 = Enumerable.Repeat(1.0, m).ToArray();
        var p4 = SignedMaskBalanceClaim.PlaceStrings(new[] { ("XZ", 1L), ("ZY", 1L) }, m, ring: false);
        var h4 = ComplementConnectionGraph.Columns(m, p4.Select(p => (new PauliMask(p.X, p.Z), (System.Numerics.BigInteger)p.C)));
        double lam4 = StrongDephasingEdges.Build(m, h4, g4).Far.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian).EigenValues.Min(z => z.Real);
        double squares = StrongDephasingEdges.CycleReading(m, h4, g4, new[] { 1, 2, 1, 2 })!.Value.LambdaStar
                       + StrongDephasingEdges.CycleReading(m, h4, g4, new[] { 8, 4, 8, 4 })!.Value.LambdaStar;
        Assert.Equal(2 * (3 - System.Math.Sqrt(5)), lam4, 12);
        Assert.Equal(lam4, squares, 12);

        // with a potential: 10·XXI + 9·XXZ + 10·IXX + 9·ZXX + 10·XIY at rate 1, a = 0.9 scaled by ten: λ_min(𝓜) = 100(1 + a²) = 181,
        // λ_* = 90.5, the least hop 100(1 − a²)/2 = 9.5, holonomy −1; (3/2)·λ_* = 135.75 against 2·μ_min = 19
        var g3 = Enumerable.Repeat(1.0, 3).ToArray();
        var p3 = SignedMaskBalanceClaim.PlaceStrings(new[] { ("XXI", 10L), ("XXZ", 9L), ("IXX", 10L), ("ZXX", 9L), ("XIY", 10L) }, 3, ring: false);
        var h3 = ComplementConnectionGraph.Columns(3, p3.Select(p => (new PauliMask(p.X, p.Z), (System.Numerics.BigInteger)p.C)));
        Assert.Equal(181.0, StrongDephasingEdges.Build(3, h3, g3).Far.Evd(MathNet.Numerics.LinearAlgebra.Symmetricity.Hermitian).EigenValues.Min(z => z.Real), 9);
        var t3 = StrongDephasingEdges.CycleReading(3, h3, g3, new[] { 0b011, 0b110, 0b101 })!.Value;
        Assert.Equal(90.5, t3.LambdaStar, 9);
        Assert.Equal(9.5, t3.HopMin, 9);
        Assert.Equal(System.Math.PI, t3.ThetaMin, 12);
    }

    // Corollary H1 and the trellis by hand at rate 1 and coefficient 1 (the second outside reply's examples): the trellis
    // reads the enumeration's examples and rings far past its reach (the odd XY ring 1 at every length, XIY 2 exactly at
    // N ≡ 2 mod 4, both the weight over one cut, so the cut bound is met); +XIY + XXX on the chain grows in steps, at least
    // twice the least weight 1/3 of each translate of Theorem B's relation by 3 (pairwise disjoint), the steps measured;
    // and the trellis refuses a ring below its guard, two templates on one flip set and a diagonal template of odd weight.
    [Fact]
    public void CorollaryH1ReadsTheOutsideRepliesExamples()
    {
        double T(int n, bool ring, params (string T, long C)[] t) => SignedMaskBalanceClaim.StrongDephasingBreakTrellis(t, Enumerable.Repeat(1.0, n).ToArray(), ring);
        Assert.Equal(1.0, T(5, true, ("XY", 1)));
        Assert.Equal(1.0, T(33, true, ("XY", 1)));
        Assert.Equal(0.0, T(32, true, ("XY", 1)));
        Assert.Equal(2.0, T(6, true, ("XIY", 1)));
        Assert.Equal(2.0, T(30, true, ("XIY", 1)));
        Assert.Equal(0.0, T(28, true, ("XIY", 1)));
        Assert.Equal(0.0, T(4, false, ("XIY", 1), ("XXX", 1)));
        var steps = new[] { 1.0 / 3, 1.0 / 2, 1.0 / 2, 2.0 / 3, 1.0, 1.0, 7.0 / 6, 4.0 / 3 };          // measured, N = 5 to 12
        for (int n = 5; n <= 40; n++)
        {
            double lam = T(n, false, ("XIY", 1), ("XXX", 1));
            Assert.True(lam >= 2 * ((n - 5) / 3 + 1) / 3.0 - 1e-12, $"N = {n}: λ = {lam} below the packing bound");
            if (n <= 12) Assert.Equal(2 * steps[n - 5], lam, 12);
        }
        Assert.True(T(40, false, ("XIY", 1), ("XXX", 1)) > 8 - 1e-12);                                // the bound at 40, 2·12/3
        Assert.Throws<System.ArgumentException>(() => T(4, true, ("XIY", 1)));
        Assert.Throws<System.ArgumentException>(() => T(6, false, ("XY", 1), ("YX", 1)));
        Assert.Throws<System.ArgumentException>(() => T(6, false, ("XY", 1), ("Z", 1)));
        Assert.Throws<System.ArgumentException>(() => T(1, false, ("ZZ", 1)));                         // wider than N, though diagonal
        Assert.Throws<System.ArgumentException>(() => T(20, false, ("X" + new string('I', 15) + "Y", 1)));   // wider than 16
        Assert.Throws<System.ArgumentOutOfRangeException>(() => T(65, false, ("XY", 1)));
    }

    // Theorem H's closed form by hand at rate 1 on every site, every coefficient 1: the odd XY ring misses one bond of
    // weight 1/2, XIY on the ring of six one bond in each odd triangle, +XIY + XXX on the chain of five one XXX of weight
    // 1/3; a balanced cell reads 0; and the reader refuses two strings on one flip set and a diagonal string of odd weight.
    [Fact]
    public void TheoremHReadsThePagesExamples()
    {
        double At(int n, bool ring, params (string T, long C)[] t) =>
            SignedMaskBalanceClaim.StrongDephasingBreak(SignedMaskBalanceClaim.PlaceStrings(t, n, ring), Enumerable.Repeat(1.0, n).ToArray());
        Assert.Equal(1.0, At(5, true, ("XY", 1)));
        Assert.Equal(1.0, At(7, true, ("XY", 1)));
        Assert.Equal(0.0, At(6, true, ("XY", 1)));
        Assert.Equal(2.0, At(6, true, ("XIY", 1)));
        Assert.Equal(2.0 / 3, At(5, false, ("XIY", 1), ("XXX", 1)), 15);
        Assert.Equal(0.0, At(4, false, ("XIY", 1), ("XXX", 1)));
        Assert.Equal(2.0 / 4, SignedMaskBalanceClaim.StrongDephasingBreak(SignedMaskBalanceClaim.PlaceStrings(new[] { ("XY", 1L) }, 5, true), new[] { 1.0, 1.0, 1.0, 1.0, 3.0 }));   // a bond beside the rate 3 is missed, weight 1/4
        Assert.Throws<System.ArgumentException>(() => At(3, false, ("XZ", 1), ("XI", 1)));
        Assert.Throws<System.ArgumentException>(() => At(3, false, ("XY", 1), ("Z", 1)));
    }

    // The readings themselves, by hand: the page's worked examples.
    [Fact]
    public void TheClosedFormsReadThePagesExamples()
    {
        Assert.Equal((true, (int?)6), SignedMaskBalanceClaim.Chain("IXZX", "XIZX"));
        Assert.Equal((true, (int?)5), SignedMaskBalanceClaim.Chain("XIY", "XXX"));
        Assert.Equal((true, (int?)4), SignedMaskBalanceClaim.Chain("XIY", "XXZ"));
        Assert.Equal((false, (int?)null), SignedMaskBalanceClaim.Chain("XIX", "XXX"));
        Assert.Equal((true, (int?)3), SignedMaskBalanceClaim.Chain(new[] { "XX", "XZ", "XXX" }));                  // the triangle {1,2} ⊕ {0} = {0,1,2}
        Assert.Equal((true, (int?)5), SignedMaskBalanceClaim.Chain(new[] { "XIY", "XXX" }));                       // Theorem B′ at two templates is Theorem B
        Assert.Equal((false, (int?)null), SignedMaskBalanceClaim.Chain(new[] { "XY" }));                           // one template: no relation
        Assert.Equal(new ulong[] { 0b10, 0b11, 1 }, SignedMaskBalanceClaim.ChainRelationBasis(new[] { "XX", "XZ", "XIY" }).Single(b => b[2] != 0));   // (x, 1 + x, 1)
        Assert.True(SignedMaskBalanceClaim.RingBreaks("XY", 5));
        Assert.False(SignedMaskBalanceClaim.RingBreaks("XY", 6));
        Assert.True(SignedMaskBalanceClaim.RingBreaks("XIY", 6));
        Assert.False(SignedMaskBalanceClaim.RingBreaks("XIY", 8));
        Assert.True(SignedMaskBalanceClaim.RingBreaks("XIIIY", 12));
        Assert.False(SignedMaskBalanceClaim.RingBreaks("XIIIY", 16));
        Assert.False(SignedMaskBalanceClaim.RingBreaks("XYX", 7));
        Assert.False(SignedMaskBalanceClaim.RingBreaks("YY", 5));
        Assert.True(SignedMaskBalanceClaim.RingBreaks(new[] { "XXZ", "XZY" }, 7));     // v = 1, 2 at T = 1: the odd sign
        Assert.False(SignedMaskBalanceClaim.RingBreaks(new[] { "XXZ", "XZY" }, 6));    // T = 2: s_j + s_i₀·0 = 0
        Assert.True(SignedMaskBalanceClaim.RingBreaks(new[] { "XX", "XIY" }, 16));     // the larger valuation's sign is odd
        Assert.True(SignedMaskBalanceClaim.RingBreaks(new[] { "XX", "XIIY" }, 8));     // v = 1, 1: s₁ ≠ s₂
        Assert.False(SignedMaskBalanceClaim.RingBreaks(new[] { "XY", "XIIY" }, 8));    // v = 1, 1, both odd
        Assert.True(SignedMaskBalanceClaim.RingBreaks(new[] { "XY", "XIIY" }, 7));     // T = 1: every mask even, a sign odd
        Assert.True(SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnChain(new[] { "XXI", "XIX" }, 8)));
        Assert.False(SignedMaskBalanceClaim.Balanced(SignedMaskBalanceClaim.OnChain(new[] { "XXI", "XIY" }, 4)));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => SignedMaskBalanceClaim.OnRing(new[] { "XIY" }, 4).ToList());
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.OnChain(new[] { "XX", "XY" }, 3).ToList());
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.RingBreaks(System.Array.Empty<string>(), 8));
        Assert.Throws<System.ArgumentOutOfRangeException>(() => SignedMaskBalanceClaim.OnRing(new[] { "XY" }, 65).ToList());
    }

    // The character ratios are the solutions of the pairing equation. A flip set's entry is g(x) = Σ_z ĝ(z)·(−1)^⟨z, x⟩,
    // ĝ(z) = c·(−i)^#Y for its string with Z part z; g(x̄)/g(x) = ω·(−1)^⟨R, x⟩ at every x exactly when
    // (−1)^|z|·ĝ(z) = ω·ĝ(z ⊕ R) for every z, the characters being independent. On the same strings as the identity test,
    // wherever g vanishes nowhere: a solution (ω, R) exists exactly when ReadRatio reads a character, and it is ReadRatio's
    // (q, R), unique, a second route that never forms the ratio.
    [Fact]
    public void TheCharacterRatiosAreTheSolutionsOfThePairingEquation()
    {
        static GaussianInteger Turn(GaussianInteger a, int q) { for (int k = 0; k < q; k++) a = a * GaussianInteger.I; return a; }
        int nowhereZero = 0, characters = 0, nonCharacters = 0;
        void Check(ulong f, int w, IReadOnlyList<SignedMaskBalanceClaim.PlacedString> strings)
        {
            var hat = new Dictionary<ulong, GaussianInteger>();
            foreach (var s in strings) hat[s.Z] = Turn(new GaussianInteger(s.C, 0), (3 * BitOperations.PopCount(s.X & s.Z)) % 4);
            GaussianInteger Hat(ulong z) => hat.TryGetValue(z, out var v) ? v : GaussianInteger.Zero;
            for (ulong x = 0; x < 1UL << w; x++)
            {
                var g = GaussianInteger.Zero;
                foreach (var (z, v) in hat) g = g + ((BitOperations.PopCount(z & x) & 1) == 1 ? -v : v);
                if (g == GaussianInteger.Zero) return;
            }
            nowhereZero++;
            var solutions = new List<(int Q, ulong R)>();
            for (int q = 0; q < 4; q++)
                for (ulong r = 0; r < 1UL << w; r++)
                    if (hat.Keys.Concat(hat.Keys.Select(z => z ^ r)).Distinct().All(z =>
                            ((BitOperations.PopCount(z) & 1) == 1 ? -Hat(z) : Hat(z)) == Turn(Hat(z ^ r), q)))
                        solutions.Add((q, r));
            var read = SignedMaskBalanceClaim.ReadRatio(f, strings);
            if (read.Kind == SignedMaskBalanceClaim.RatioKind.Character)
            {
                characters++;
                Assert.Equal(new[] { (read.Quarter, read.Read) }, solutions.ToArray());
            }
            else
            {
                Assert.Empty(solutions);
                nonCharacters++;
            }
        }
        foreach (var (w, most) in new[] { (1, 4), (2, 4), (3, 4), (4, 3) })
            for (ulong f = 1; f < 1UL << w; f++)
            {
                var onMask = Enumerable.Range(0, 1 << w).Select(choice => (ulong)choice).ToArray();   // Y on the mask, Z off it
                void Grow(int start, List<SignedMaskBalanceClaim.PlacedString> chosen)
                {
                    if (chosen.Count > 0) Check(f, w, chosen);
                    if (chosen.Count == most) return;
                    for (int k = start; k < onMask.Length; k++)
                        foreach (long c in chosen.Count == 0 ? new[] { 1L } : new[] { 1L, -1L })
                        {
                            chosen.Add(new SignedMaskBalanceClaim.PlacedString(f, onMask[k], c));
                            Grow(k + 1, chosen);
                            chosen.RemoveAt(chosen.Count - 1);
                        }
                }
                Grow(0, new List<SignedMaskBalanceClaim.PlacedString>());
            }
        Assert.Equal((39892, 9828, 30064), (nowhereZero, characters, nonCharacters));                              // measured
    }

    // Open 1 at ω = ±i: ĝ(z) is real on the strings with an even number of Y and imaginary on the others, g = A + iB; at odd
    // q the pairing equation swaps the two kinds, so g(x) = A(x) ∓ i·(−1)^⟨R, x⟩·A(x̄) and |g(x)|² = A(x)² + A(x̄)²: g
    // vanishes nowhere exactly when A(x) and A(x̄) are never both zero. Read on every solution with odd q of
    // every set of one to four strings of one mask in a window of up to three sites, the solutions where g vanishes
    // somewhere included, the identity exactly in integers.
    [Fact]
    public void AtOddQuarterTheRealPartDecidesWhereTheEntryVanishes()
    {
        static GaussianInteger Turn(GaussianInteger a, int q) { for (int k = 0; k < q; k++) a = a * GaussianInteger.I; return a; }
        int solutions = 0, nowhereZero = 0;
        foreach (int w in new[] { 1, 2, 3 })
        {
            int full = (1 << w) - 1;
            for (ulong f = 1; f < 1UL << w; f++)
            {
                void Check(List<SignedMaskBalanceClaim.PlacedString> strings)
                {
                    var hat = strings.ToDictionary(s => s.Z, s => Turn(new GaussianInteger(s.C, 0), (3 * BitOperations.PopCount(s.X & s.Z)) % 4));
                    GaussianInteger Hat(ulong z) => hat.TryGetValue(z, out var v) ? v : GaussianInteger.Zero;
                    for (int q = 1; q < 4; q += 2)
                        for (ulong r = 0; r < 1UL << w; r++)
                        {
                            if (!hat.Keys.Concat(hat.Keys.Select(z => z ^ r)).Distinct().All(z =>
                                    ((BitOperations.PopCount(z) & 1) == 1 ? -Hat(z) : Hat(z)) == Turn(Hat(z ^ r), q))) continue;
                            solutions++;
                            var g = new GaussianInteger[full + 1];
                            var a = new BigInteger[full + 1];
                            for (int x = 0; x <= full; x++)
                            {
                                g[x] = GaussianInteger.Zero;
                                a[x] = BigInteger.Zero;
                                foreach (var (z, v) in hat)
                                {
                                    bool odd = (BitOperations.PopCount(z & (ulong)x) & 1) == 1;
                                    g[x] = g[x] + (odd ? -v : v);
                                    a[x] += odd ? -v.Re : v.Re;
                                }
                            }
                            Assert.True((q & 1) == (BitOperations.PopCount(r & f) & 1), $"q ≢ |R ∩ m| at mask {f}");   // the reality of the coefficients
                            for (int x = 0; x <= full; x++)
                            {
                                var turned = Turn(new GaussianInteger((BitOperations.PopCount(r & (ulong)x) & 1) == 1 ? -a[x ^ full] : a[x ^ full], BigInteger.Zero), 4 - q);
                                Assert.True(g[x] == new GaussianInteger(a[x], BigInteger.Zero) + turned, $"g ≠ A ∓ i·χ_R·A(x̄) at mask {f}, x = {x}");
                                Assert.True(g[x].Re * g[x].Re + g[x].Im * g[x].Im == a[x] * a[x] + a[x ^ full] * a[x ^ full], $"|g|² at mask {f}, x = {x}");
                            }
                            bool vanishesNowhere = g.All(v => v != GaussianInteger.Zero);
                            bool noZeroPair = Enumerable.Range(0, full + 1).All(x => !(a[x].IsZero && a[x ^ full].IsZero));
                            Assert.Equal(vanishesNowhere, noZeroPair);
                            if (vanishesNowhere) nowhereZero++;
                        }
                }
                void Grow(int start, List<SignedMaskBalanceClaim.PlacedString> chosen)
                {
                    if (chosen.Count > 0) Check(chosen);
                    if (chosen.Count == 4) return;
                    for (int z = start; z < 1 << w; z++)
                        foreach (long c in chosen.Count == 0 ? new[] { 1L } : new[] { 1L, -1L })
                        {
                            chosen.Add(new SignedMaskBalanceClaim.PlacedString(f, (ulong)z, c));
                            Grow(z + 1, chosen);
                            chosen.RemoveAt(chosen.Count - 1);
                        }
                }
                Grow(0, new List<SignedMaskBalanceClaim.PlacedString>());
            }
        }
        Assert.Equal((546, 346), (solutions, nowhereZero));                                                 // measured
    }

    // Theorem B's c is the girth ladder's syzygy polynomial B = q₂·Q₁ + q₁·Q₂ (Q the Y/Z mask) read at x = 1, since
    // Q(1) is the sign: on every pair of different shapes among the two-site and window templates.
    [Fact]
    public void TheoremBsCIsTheGirthLaddersSyzygyAtOne()
    {
        var templates = (from a in "XYZ" from b in "XYZ" select $"{a}{b}")
            .Concat(from a in "XYZ" from m in "IXYZ" from b in "XYZ" select $"{a}{m}{b}")
            .Where(t => SignedMaskBalanceClaim.Mask(t) != 0).ToList();
        static (ulong P, ulong Q) Masks(string t)
        {
            var (p, q, _) = WindowedObstructionScan.TemplateMasks(t.Select(PauliLetterExtensions.FromSymbol).ToArray());
            return (p, q);
        }
        int pairs = 0, odd = 0;
        for (int i = 0; i < templates.Count; i++)
            for (int j = i + 1; j < templates.Count; j++)
            {
                ulong ma = SignedMaskBalanceClaim.Mask(templates[i]), mb = SignedMaskBalanceClaim.Mask(templates[j]);
                if (ma >> BitOperations.TrailingZeroCount(ma) == mb >> BitOperations.TrailingZeroCount(mb)) continue;
                var (p1, q1) = Masks(templates[i]);
                var (p2, q2) = Masks(templates[j]);
                bool atOne = (BitOperations.PopCount(WindowedObstructionScan.SyzygyB(p1, q1, p2, q2)) & 1) == 1;
                Assert.True(SignedMaskBalanceClaim.Chain(templates[i], templates[j]).Breaks == atOne, $"{templates[i]} + {templates[j]}");
                pairs++;
                if (atOne) odd++;
            }
        Assert.Equal((648, 324), (pairs, odd));                                                              // measured
    }

    // Theorem E by hand: the read sets of two-string flip sets and the branches of Corollary E2.
    [Fact]
    public void TheCharacterConnectionReadsThePagesExamples()
    {
        var mid = SignedMaskBalanceClaim.ReadRatio(0b111, SignedMaskBalanceClaim.PlaceStrings(new[] { ("XXX", 1L), ("XYX", 1L) }, 3, ring: false));
        Assert.Equal(SignedMaskBalanceClaim.RatioKind.Character, mid.Kind);
        Assert.Equal(0b010UL, mid.Read);                       // X·(X + Y)·X reads its middle site
        Assert.Equal(1, mid.Quarter % 2);                      // ω = ±i, the middle site inside the mask
        var unequal = SignedMaskBalanceClaim.ReadRatio(0b111, SignedMaskBalanceClaim.PlaceStrings(new[] { ("XXX", 1L), ("XYX", 2L) }, 3, ring: false));
        Assert.Equal(SignedMaskBalanceClaim.RatioKind.NotCharacter, unequal.Kind);   // c₁² ≠ c₂²
        var flipFlop = SignedMaskBalanceClaim.ReadRatio(0b11, SignedMaskBalanceClaim.PlaceStrings(new[] { ("XX", 1L), ("YY", 1L) }, 2, ring: false));
        Assert.Equal(SignedMaskBalanceClaim.RatioKind.AbsentEdge, flipFlop.Kind);    // XX + YY vanishes where the bits agree
        Assert.Equal(SignedMaskBalanceClaim.RingShapeBranch.ReadSetNotDivisible, SignedMaskBalanceClaim.RingOneShape(new[] { ("XXX", 1L), ("XYX", 1L) }, 9).Branch);
        Assert.Equal(SignedMaskBalanceClaim.RingShapeBranch.NoRelation, SignedMaskBalanceClaim.RingOneShape(new[] { ("XXX", 1L), ("XYX", 1L) }, 8).Branch);
        Assert.Equal(SignedMaskBalanceClaim.RingShapeBranch.OnePlusXDoesNotDivide, SignedMaskBalanceClaim.RingOneShape(new[] { ("XXX", 1L), ("YYY", 1L) }, 6).Branch);
        Assert.Equal(SignedMaskBalanceClaim.RingShapeBranch.HolonomyMinus, SignedMaskBalanceClaim.RingOneShape(new[] { ("XIY", 1L), ("XZX", 1L) }, 6).Branch);
        Assert.Equal(SignedMaskBalanceClaim.RingShapeBranch.HolonomyPlus, SignedMaskBalanceClaim.RingOneShape(new[] { ("XIY", 1L), ("XZX", 1L) }, 8).Branch);
        var product = new[] { ("XX", 1L), ("XY", 1L), ("YX", -1L), ("YY", -1L) };          // (X − Y)⊗(X + Y): hol(g) = (−1)^(N − 2)
        Assert.True(SignedMaskBalanceClaim.RingOneShape(product, 7).Breaks);
        Assert.False(SignedMaskBalanceClaim.RingOneShape(product, 8).Breaks);
        Assert.False(SignedMaskBalanceClaim.RingOneShape(product, 7, orderingTerm: false).Breaks);
        Assert.Equal(SignedMaskBalanceClaim.CharacterVerdict.Breaks,
            SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(new[] { ("XZ", 1L), ("ZY", 1L) }, 6, ring: false)));   // the chain's end
        Assert.Equal(SignedMaskBalanceClaim.CharacterVerdict.Outside,
            SignedMaskBalanceClaim.DecideCharacter(SignedMaskBalanceClaim.PlaceStrings(new[] { ("XX", 1L), ("YY", 1L) }, 4, ring: false)));
    }

    // Corollary E3 by hand: two flip shapes on the chain, the read sum D and hol(u₀).
    [Fact]
    public void TheTwoShapeChainReadsThePagesExamples()
    {
        var field = new[] { ("X", 1L), ("Y", 1L) };
        var plus = new[] { ("XXX", 1L), ("YYY", 1L) };
        var minus = new[] { ("XXX", 1L), ("YYY", -1L) };
        Assert.Equal((3, 0UL), SignedMaskBalanceClaim.ChainTwoShapeData(field, plus));                  // q₂R₁ + q₁R₂ = f·1 + 1·f
        Assert.Equal((false, SignedMaskBalanceClaim.ChainShapeBranch.HolonomyPlus), SignedMaskBalanceClaim.ChainTwoShapes(field, plus, 5));
        Assert.Equal((true, SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus), SignedMaskBalanceClaim.ChainTwoShapes(field, minus, 5));
        Assert.True(SignedMaskBalanceClaim.ChainTwoShapes(field, plus, 5, orderingTerm: false).Breaks);   // the ordering term is −1 here
        var middle = new[] { ("XXX", 1L), ("XYX", 1L) };
        var wide = new[] { ("XIIX", 1L) };
        Assert.Equal((4, 0b110UL), SignedMaskBalanceClaim.ChainTwoShapeData(middle, wide));              // D = (1 + x)·x
        Assert.Equal(SignedMaskBalanceClaim.ChainShapeBranch.ReadSumNonzero, SignedMaskBalanceClaim.ChainTwoShapes(middle, wide, 4).Branch);
        Assert.Equal(SignedMaskBalanceClaim.ChainShapeBranch.SquaresFail, SignedMaskBalanceClaim.ChainTwoShapes(middle, wide, 5).Branch);
        Assert.Equal((false, SignedMaskBalanceClaim.ChainShapeBranch.NoRelation), SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XIY", 1L) }, new[] { ("XXX", 1L) }, 4));
        Assert.Equal((true, SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus), SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XIY", 1L) }, new[] { ("XXX", 1L) }, 5));   // Theorem B's c = 1 from N* = 5
        var fieldBesideProduct = SignedMaskBalanceClaim.ChainTwoShapes(field, new[] { ("XX", 1L), ("XY", 1L), ("YX", -1L), ("YY", -1L) }, 4);
        Assert.Equal((true, SignedMaskBalanceClaim.ChainShapeBranch.HolonomyMinus), fieldBesideProduct);  // ω_f²·ω_p = −1, the ordering term +1
        Assert.Equal(fieldBesideProduct, SignedMaskBalanceClaim.ChainTwoShapes(field, new[] { ("XX", 1L), ("XY", 1L), ("YX", -1L), ("YY", -1L) }, 4, orderingTerm: false));
        Assert.Equal((true, SignedMaskBalanceClaim.ChainShapeBranch.ModuliDiffer),
            SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("IX", 1L), ("ZX", 1L) }, new[] { ("XXX", 1L) }, 3));   // 1 + (−1)^x₀ vanishes on one side only
        Assert.Equal((false, SignedMaskBalanceClaim.ChainShapeBranch.Outside),
            SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XX", 1L), ("YY", 1L) }, new[] { ("X", 1L) }, 3));     // an absent edge: undecided
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.ChainTwoShapeData(new[] { ("XX", 1L), ("YY", 1L) }, new[] { ("X", 1L) }));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XX", 1L) }, new[] { ("XY", 1L) }, 4));
        Assert.Throws<System.ArgumentException>(() => SignedMaskBalanceClaim.ChainTwoShapes(new[] { ("XX", 1L), ("XIY", 1L) }, wide, 5));
    }

    // Two identities every character ratio obeys: r(x)·r(x̄) = 1 by definition gives ω² = (−1)^|R|, which ReadRatio
    // enforces itself by reading the character at every x (asserted here as a consistency check, it cannot fail on a
    // Character reading), and H Hermitian, r(x ⊕ m) = 1/r(x), gives ω² = (−1)^⟨R, m⟩, which ReadRatio does not use and
    // which tests its entry convention; together the read set has an even number of sites off its flip set. A reading
    // that misses the second is a fault of the reader, not of the physics. Read on every set of one to four strings of one mask in a window of up to three sites and one to
    // three strings in a window of four, signs ±1, and on the page's examples.
    [Fact]
    public void EveryCharacterRatioReadsAnEvenNumberOfSitesOffItsFlipSet()
    {
        static int Parity(ulong v) => BitOperations.PopCount(v) & 1;
        int characters = 0, oddQuarter = 0, offMask = 0, realBesideEvenRead = 0;
        void Check(ulong f, IReadOnlyList<SignedMaskBalanceClaim.PlacedString> strings)
        {
            var r = SignedMaskBalanceClaim.ReadRatio(f, strings);
            if (r.Kind != SignedMaskBalanceClaim.RatioKind.Character) return;
            characters++;
            Assert.True((r.Quarter & 1) == Parity(r.Read), $"ω² ≠ (−1)^|R| at mask {f}, read {r.Read}, quarter {r.Quarter}");
            Assert.True(Parity(r.Read & ~f) == 0, $"|R \\ m| odd at mask {f}, read {r.Read}");
            if ((r.Quarter & 1) == 1) oddQuarter++;
            if ((r.Read & ~f) != 0) offMask++;
            if (r.Quarter % 2 == 0 && r.Read != 0) realBesideEvenRead++;
        }
        foreach (var (w, most) in new[] { (1, 4), (2, 4), (3, 4), (4, 3) })
            for (ulong f = 1; f < 1UL << w; f++)
            {
                var onMask = Enumerable.Range(0, 1 << w).Select(choice =>
                {
                    ulong z = 0;
                    for (int j = 0; j < w; j++) if (((choice >> j) & 1) == 1) z |= 1UL << j;   // Y on the mask, Z off it
                    return z;
                }).Distinct().ToArray();
                void Grow(int start, List<SignedMaskBalanceClaim.PlacedString> chosen)
                {
                    if (chosen.Count > 0) Check(f, chosen);
                    if (chosen.Count == most) return;
                    for (int k = start; k < onMask.Length; k++)
                        foreach (long c in chosen.Count == 0 ? new[] { 1L } : new[] { 1L, -1L })
                        {
                            chosen.Add(new SignedMaskBalanceClaim.PlacedString(f, onMask[k], c));
                            Grow(k + 1, chosen);
                            chosen.RemoveAt(chosen.Count - 1);
                        }
                }
                Grow(0, new List<SignedMaskBalanceClaim.PlacedString>());
            }

        // the page's examples, with their read sets
        var xzzy = SignedMaskBalanceClaim.ReadRatio(0b010, new[] { new SignedMaskBalanceClaim.PlacedString(0b010, 0b100, 1), new SignedMaskBalanceClaim.PlacedString(0b010, 0b011, 1) });
        Assert.Equal((SignedMaskBalanceClaim.RatioKind.Character, 0b111UL), (xzzy.Kind, xzzy.Read));          // XZ ± ZY at {1}: R \ m = {0, 2}
        var cz = SignedMaskBalanceClaim.ReadRatio(0b001, new[] { new SignedMaskBalanceClaim.PlacedString(1, 0, 1), new SignedMaskBalanceClaim.PlacedString(1, 0b100, 1),
                                                                  new SignedMaskBalanceClaim.PlacedString(1, 0b010, 1), new SignedMaskBalanceClaim.PlacedString(1, 0b110, -1) });
        Assert.Equal((SignedMaskBalanceClaim.RatioKind.Character, 2, 0b110UL), (cz.Kind, cz.Quarter, cz.Read));   // X ⊗ 2·CZ: ω = −1, R \ m = {1, 2}
        var product = SignedMaskBalanceClaim.ReadRatio(0b11, SignedMaskBalanceClaim.PlaceStrings(new[] { ("XX", 15L), ("XY", 15L), ("YX", -15L), ("YY", -15L) }, 2, ring: false));
        Assert.Equal((SignedMaskBalanceClaim.RatioKind.Character, 0, 0b11UL), (product.Kind, product.Quarter, product.Read));   // (3X − 3Y)⊗(5X + 5Y): ω = 1, R = m
        var unequalProduct = SignedMaskBalanceClaim.ReadRatio(0b11, SignedMaskBalanceClaim.PlaceStrings(new[] { ("XX", 10L), ("XY", 14L), ("YX", 15L), ("YY", 21L) }, 2, ring: false));
        Assert.Equal(SignedMaskBalanceClaim.RatioKind.NotCharacter, unequalProduct.Kind);                    // (2X + 3Y)⊗(5X + 7Y): a² ≠ b², no character
        foreach (var e in new[] { xzzy, cz, product })
            Assert.True(Parity(e.Read & ~e.Mask) == 0 && (e.Quarter & 1) == Parity(e.Read));

        // anti-vacuity: both identities are read where they say something
        Assert.True(oddQuarter > 0 && offMask > 0 && realBesideEvenRead > 0, $"{characters}, {oddQuarter}, {offMask}, {realBesideEvenRead}");
        Assert.Equal((9828, 1370, 504, 148), (characters, oddQuarter, offMask, realBesideEvenRead));          // measured
    }
}
