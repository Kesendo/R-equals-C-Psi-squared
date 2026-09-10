using System.Globalization;
using System.Numerics;
using System.Text.Json;
using MathNet.Numerics.LinearAlgebra;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Diagnostics.Foundation;

/// <summary>F164, recomputed at inspect time: the N=5 Route-B A₂ locus whose coupling AND eigenvalue are
/// both real carries a SEMISIMPLE double eigenvalue, and the reason is a field argument rather than a
/// contour reading. This witness executes the parts of that argument that are self-contained structure,
/// and names the one input it consumes from elsewhere.
///
/// <para><b>Scope.</b> Of the N=5 inventory's 58 q-loci, all 24 imaginary-q loci are directly
/// Hermitian-certified. F164 certifies the complete 26-member R-odd orbit, overlapping 12 Hermitian loci
/// and adding the 2 real-q plus 12 nonreal-q loci. Thus 38 distinct loci are exact-certified semisimple;
/// only the 20 nonreal-q R-even loci remain numerical-only. The classifier still reads all 32 nonreal-q
/// loci. The artifact's empty <c>exactRankCertificates</c> array is classifier-local provenance, not the
/// current verdict. Not an all-N theorem, and no metrological claim follows from it.</para>
///
/// <para><b>What the proof does, and which half lives here.</b> Let w₀ be the positive real root of the
/// R-odd A₂ layer and t = √(−w). The proof's first step is arithmetic: A₂_O is irreducible over ℚ of
/// degree 13, so ℚ(w₀) is a degree-13 subfield of ℝ while t₀ = √(−w₀) is nonzero purely imaginary and lies
/// outside it, giving [ℚ(t₀) : ℚ] = 26 and a SINGLE Galois orbit containing all 26 lifts. That step is
/// exact integer arithmetic over the A₂_O coefficients and is certified by
/// <c>simulations/o2b_gcd_certificate.py</c>; this witness CONSUMES it rather than redoing it, and says so
/// on <see cref="OrbitBookkeeping"/>. The second step is the one recomputed here: at each of the 12 real t
/// the R-odd block L_O = D_O + t·K_O is real symmetric, hence semisimple, and rank is preserved by every
/// embedding because a minor is zero or nonzero before and after it. So the Hermitian-axis verdict travels
/// the orbit to the physical point.</para>
///
/// <para><b>The book is pinned, not assumed.</b> The repository carries two coupling conventions (the
/// glossary's q and Q, a factor of two apart), and a witness that picked one would be reading a different
/// operator half the time. <see cref="PinnedBook"/> scans the three candidate scalings and reports the one
/// that actually places the stored λ₀ in the R-odd spectrum. If none does, or more than one does, the
/// reading says so instead of choosing.</para>
///
/// <para><b>The direction law, and why it is here.</b> The gate beside the proof measures the coupling
/// response σ₂ = c·(δJ/J). The constant is direction-dependent and no parity fixes it: three PALINDROMIC
/// directions give three different c. What parity decides is whether a response exists at all, and that is
/// a theorem rather than a sample. L is affine in the bond couplings and the chain reflection carries
/// ∂L/∂J_b to ∂L/∂J_{N−2−b}, so an ANTI-palindromic pattern obeys R V R = −V at every base point.
/// For two R-odd vectors u, v that forces vᵀVu = −vᵀVu; hence P_- V P_- vanishes identically, at every
/// N and independently of the base point. Reflection symmetry of the base is needed only so the R-even
/// and R-odd subspaces are invariant under L.
/// <see cref="DirectionReach"/> recomputes it as a projection norm ratio, with no eigensolver.
/// It is a BLIND SPOT of the sector reading and not a null response: inside the sector PVP = 0 holds σ₂
/// on its floor, while the physical pair splits through the EVEN sector at second order. A reading that
/// stopped at the projection would call such a detune harmless.</para></summary>
public sealed record N6RealQExclusionReading(int LocusCount, int RealAxisTouchCount)
{
    public int ExcludedCount => LocusCount - RealAxisTouchCount;
    public bool AllExcluded => LocusCount == 266 && RealAxisTouchCount == 0;
    public string Summary => $"{ExcludedCount} of {LocusCount} certified exact rational tBox.real intervals " +
                             $"exclude zero; real-axis touches = {RealAxisTouchCount}";
}

public sealed class RouteBN5RealQSemisimpleWitness : IInspectable
{
    private const int Sites = 5;
    private const string N6AtlasCanonicalSha256 = "c4b19015a750e56a55b0c8fdba3d05534c0c133bc5355a74eb89b77020532d47";
    private const string N6Model = "open-uniform-XY-delta0-field0-gamma1-SEket-DEbra";
    private const string N6ParameterConvention = "t=i*qCSharp;qCSharp=-i*t";

    /// <summary>The stored physical eigenvalue at the locus (route_b_a2_n5.json, R-odd positive-real root).</summary>
    public const double Lambda0 = -4.7919603651796410;

    /// <summary>w₀, the positive real root of the R-odd A₂ layer (route_b_a2_n5.json).</summary>
    public const double W0 = 5.1008310208864164;

    /// <summary>|w| at one of the layer's six NEGATIVE real roots (route_b_a2_n5.json), whose t = √|w| is
    /// real: a sampled member of the twelve orbit points that sit on the Hermitian axis.</summary>
    public const double NegativeRootMagnitude = 2.8627663562428665;

    public string DisplayName => "F164: the N=5 real-q A₂ locus is semisimple by Galois conjugacy";

    public string Summary
    {
        get
        {
            var book = PinnedBook;
            var axis = HermitianAxisResidual;
            var n6 = N6RealQExclusion;
            return $"the only real-q points in the N=5/N=6 A2 inventories have settable parameter q = ±1.1292509708747671 " +
                   $"and resulting spectral eigenvalue λ = {Lambda0.ToString("0.0000000000", CultureInfo.InvariantCulture)}; " +
                   $"q is the settable parameter and lambda is the resulting spectral eigenvalue; " +
                   $"book pinned by spectrum/operator reconstruction to {book}; " +
                   $"on the Hermitian axis the R-odd block is real symmetric with residual " +
                   $"{axis.ToString("0.###e+00", CultureInfo.InvariantCulture)}, so the orbit carries " +
                   $"semisimplicity to every R-odd conjugate. N=5: 24 imaginary-q directly Hermitian-certified; " +
                   $"F164 certifies all 26 R-odd loci, so 38 distinct loci are exact-certified and only 20 " +
                   $"nonreal-q R-even loci are numerical-only; empty exactRankCertificates is " +
                   $"classifier-local provenance. N=6: {n6.Summary}; not an all-N theorem (Tier 1, derived)";
        }
    }

    /// <summary>Exact N=6 real-q exclusion from the committed atlas. Since t=i q, q is real iff
    /// Re(t)=0; every certified rational real-part interval must therefore exclude zero. No display
    /// midpoint participates in this reading.</summary>
    public N6RealQExclusionReading N6RealQExclusion
    {
        get
        {
            string path = Path.Combine(RepoRootLocator.Require(), "simulations", "results", "route_b_a2_n6.json");
            return ReadN6RealQExclusion(File.ReadAllText(path));
        }
    }

    /// <summary>Read the exact rational t-boxes from an N=6 atlas payload. Public so the gate can mutate
    /// one interval to touch zero and prove that the same decision path rejects it.</summary>
    public static N6RealQExclusionReading ReadN6RealQExclusion(string json)
    {
        using JsonDocument document = JsonDocument.Parse(json);
        JsonElement root = document.RootElement;
        if (root.ValueKind != JsonValueKind.Object)
            throw new InvalidDataException("the N=6 atlas root must be an object");
        if (ReadInt32(root, "schemaVersion") != 3 || ReadInt32(root, "n") != 6)
            throw new InvalidDataException("the real-q exclusion requires schemaVersion=3 and n=6");
        if (Required(root, "model", JsonValueKind.String).GetString() != N6Model)
            throw new InvalidDataException("the real-q exclusion requires the canonical N=6 model");
        JsonElement conventions = Required(root, "conventions", JsonValueKind.Object);
        if (Required(conventions, "parameter", JsonValueKind.String).GetString() != N6ParameterConvention
            || Required(conventions, "eigenvalue", JsonValueKind.String).GetString() != "Lambda=2*lambda"
            || Required(conventions, "coefficientOrder", JsonValueKind.String).GetString()
                != "lambda-lowest-first,t-lowest-first")
            throw new InvalidDataException("the real-q exclusion requires the canonical t/q/lambda conventions");

        JsonElement degrees = Required(root, "a2Degrees", JsonValueKind.Object);
        if (ReadInt32(degrees, "E") != 133 || ReadInt32(degrees, "O") != 133)
            throw new InvalidDataException("the N=6 atlas must declare A2 degrees E=133 and O=133");

        JsonElement loci = Required(root, "loci", JsonValueKind.Array);
        if (loci.GetArrayLength() != 266)
            throw new InvalidDataException("the N=6 atlas must contain exactly 266 certified loci");

        var ids = new HashSet<string>(StringComparer.Ordinal);
        int even = 0, odd = 0;
        int touching = 0;
        foreach (JsonElement locus in loci.EnumerateArray())
        {
            if (locus.ValueKind != JsonValueKind.Object)
                throw new InvalidDataException("every N=6 locus must be an object");
            string id = Required(locus, "id", JsonValueKind.String).GetString()!;
            if (string.IsNullOrWhiteSpace(id) || !ids.Add(id))
                throw new InvalidDataException("the N=6 atlas must contain 266 unique nonempty locus ids");

            string parity = Required(locus, "parity", JsonValueKind.String).GetString()!;
            if (parity == "E") even++;
            else if (parity == "O") odd++;
            else throw new InvalidDataException($"{id}: parity must be E or O");
            if (ReadInt32(locus, "algebraicMultiplicity") != 2)
                throw new InvalidDataException($"{id}: certified A2 locus must have algebraicMultiplicity=2");

            JsonElement tBox = Required(locus, "tBox", JsonValueKind.Object);
            JsonElement interval = Required(tBox, "real", JsonValueKind.Object);
            var lower = ReadRational(Required(interval, "lower", JsonValueKind.Object));
            var upper = ReadRational(Required(interval, "upper", JsonValueKind.Object));
            if (Compare(lower, upper) > 0)
                throw new InvalidDataException($"{id}: tBox.real interval has lower > upper");
            if (lower.Numerator.Sign <= 0 && upper.Numerator.Sign >= 0)
                touching++;
        }
        if (even != 133 || odd != 133)
            throw new InvalidDataException($"the N=6 atlas parity counts must be E=133 and O=133, got E={even}, O={odd}");
        if (RouteBN6UnfoldingCertificate.CanonicalUtf8LfSha256(json) != N6AtlasCanonicalSha256)
            throw new InvalidDataException("the N=6 boxes are detached from the semantically certified canonical carrier");
        return new N6RealQExclusionReading(loci.GetArrayLength(), touching);
    }

    private readonly record struct Rational(BigInteger Numerator, BigInteger Denominator);

    private static Rational ReadRational(JsonElement value)
    {
        BigInteger numerator = ReadBigInteger(value, "numerator");
        BigInteger denominator = ReadBigInteger(value, "denominator");
        if (denominator.Sign <= 0)
            throw new InvalidDataException("a rational box endpoint must have a positive denominator");
        return new Rational(numerator, denominator);
    }

    private static JsonElement Required(JsonElement parent, string name, JsonValueKind kind)
    {
        if (!parent.TryGetProperty(name, out JsonElement value) || value.ValueKind != kind)
            throw new InvalidDataException($"{name} must be a JSON {kind}");
        return value;
    }

    private static int ReadInt32(JsonElement parent, string name)
    {
        JsonElement value = Required(parent, name, JsonValueKind.Number);
        if (!value.TryGetInt32(out int result))
            throw new InvalidDataException($"{name} must be an exact 32-bit integer");
        return result;
    }

    private static BigInteger ReadBigInteger(JsonElement parent, string name)
    {
        string text = Required(parent, name, JsonValueKind.String).GetString()!;
        if (!BigInteger.TryParse(text, NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out BigInteger result))
            throw new InvalidDataException($"{name} must be an exact integer string");
        return result;
    }

    private static int Compare(Rational left, Rational right) =>
        (left.Numerator * right.Denominator).CompareTo(right.Numerator * left.Denominator);

    // ---------------------------------------------------------------- the block

    // Everything below is built by the repository's own (SE,DE) builder, XxzCoherenceBlock, rather than
    // by a second copy of it: the bond-weight parameter and the R = -1 sector were added there so this
    // witness and the rest of the repo cannot drift apart.

    private static double[] Uniform() => Enumerable.Repeat(1.0, Sites - 1).ToArray();

    /// <summary>The full 50-dimensional block at Δ=0, γ=1 per site, with per-bond hopping weights.</summary>
    public static Complex[,] BuildFull(Complex q, double[] bondWeights) =>
        XxzCoherenceBlock.BuildFull(Sites, q, 0.0, bondWeights);

    /// <summary>Orthonormal columns spanning the R = −1 sector.</summary>
    public static Matrix<Complex> OddSectorColumns() => XxzCoherenceBlock.BuildOddColumns(Sites);

    /// <summary>The ±1 orbit basis of the same sector. U<sup>T</sup>·M·U is an exact rearrangement there,
    /// so a projection that vanishes in exact arithmetic comes back as 0.0 rather than as 1e-17.</summary>
    public static Matrix<Complex> OddSectorColumnsExact() => XxzCoherenceBlock.BuildOddColumns(Sites, true);

    private static Matrix<Complex> OddBlock(Complex q, double[] bondWeights)
    {
        var u = OddSectorColumns();
        var full = Matrix<Complex>.Build.DenseOfArray(BuildFull(q, bondWeights));
        return u.ConjugateTranspose() * full * u;
    }

    // ---------------------------------------------------------------- readings

    /// <summary>(full, R-even, R-odd) dimensions of the block. The proof consumes the R-odd sector, so its
    /// size is the shape the whole argument sits on.</summary>
    public (int Full, int Even, int Odd) SectorDimensions
    {
        get
        {
            // Both sectors are BUILT, neither derived by subtraction, so the split is two independent
            // readings that happen to add up rather than one reading and an arithmetic consequence.
            int odd = OddSectorColumns().ColumnCount;
            int even = XxzCoherenceBlock.BuildSym(Sites, Complex.One, 0.0).ColumnCount;
            var full = BuildFull(Complex.One, Uniform());
            return (full.GetLength(0), even, odd);
        }
    }

    private static readonly (string Label, double Scale)[] Books =
    {
        ("q = qUnitHop", 1.0), ("q = 2·qUnitHop", 2.0), ("q = qUnitHop/2", 0.5),
    };

    /// <summary>Which coupling book actually places the stored λ₀ in the R-odd spectrum. Scanned, never
    /// assumed: the repository's q and Q differ by a factor of two and a wrong book reads a different
    /// operator. Reports the ambiguity rather than resolving it if the scan does not single one out.</summary>
    public string PinnedBook
    {
        get
        {
            var hits = MatchingBooks();
            return hits.Count == 1 ? hits[0]
                 : hits.Count == 0 ? "NONE: no candidate scaling places λ₀ in the R-odd spectrum"
                 : "AMBIGUOUS: " + string.Join(", ", hits);
        }
    }

    private static List<string> MatchingBooks()
    {
        var hits = new List<string>();
        foreach (var (label, scale) in Books)
        {
            var odd = OddBlock(new Complex(scale * Math.Sqrt(W0), 0), Uniform());
            double d = odd.Evd().EigenValues.Enumerate().Min(z => (z - new Complex(Lambda0, 0)).Magnitude);
            if (d < 1e-9) hits.Add(label);
        }
        return hits;
    }

    /// <summary>The coupling the pinned book names, or an exception if the scan did not single one out.
    /// It must not fall back on a first match: reporting an ambiguity in one property and quietly
    /// resolving it in another is how a witness comes to read a different operator than it advertises.</summary>
    private static Complex PinnedQ
    {
        get
        {
            var hits = MatchingBooks();
            if (hits.Count != 1)
                throw new InvalidOperationException(
                    $"the coupling book is not pinned ({hits.Count} candidates match); " +
                    "every reading below would be on an undetermined operator");
            return new Complex(Books.First(b => b.Label == hits[0]).Scale * Math.Sqrt(W0), 0);
        }
    }

    /// <summary>The physical coupling of one of the layer's six NEGATIVE real roots, in the pinned book:
    /// a sampled member of the twelve orbit points that sit on the Hermitian axis. The artifact stores
    /// the root's qUnitHop, and the book halves it, which is the conversion this witness exists to keep
    /// straight.</summary>
    private static double NegativeRootCoupling =>
        PinnedQ.Real / Math.Sqrt(W0) * Math.Sqrt(NegativeRootMagnitude);

    /// <summary>Step 2 of the proof, recomputed: on the Hermitian axis (a NEGATIVE real root of the layer,
    /// where t is real) the R-odd block is real symmetric, hence every one of its eigenvalues is
    /// semisimple. The residual is max(|Im L|, |L − Lᵀ|) and there is an exact route to it, so it is read
    /// against 0.0 and not against a tolerance.</summary>
    public double HermitianAxisResidual
    {
        get
        {
            // A negative real root of the layer gives −w > 0, so t = √(−w) is REAL; in this builder that
            // is the purely imaginary coupling. The block is affine in t with real coefficients,
            // L_O = D_O + t·K_O, so reading three real t (including t = 0, which isolates D_O, and one
            // taken from an actual negative root) certifies the structure at every real t and not only at
            // the sampled ones. Any residual here would be a defect in the construction, so it is compared
            // to 0.0 and not to a tolerance.
            double worst = 0;
            foreach (double t in new[] { 0.0, 1.0, NegativeRootCoupling })
            {
                var odd = OddBlock(new Complex(0, t), Uniform());
                for (int r = 0; r < odd.RowCount; r++)
                    for (int c = 0; c < odd.ColumnCount; c++)
                    {
                        worst = Math.Max(worst, Math.Abs(odd[r, c].Imaginary));
                        worst = Math.Max(worst, (odd[r, c] - odd[c, r]).Magnitude);
                    }
            }
            return worst;
        }
    }

    /// <summary>How far a bond detune direction reaches into the R-odd sector, as the largest entry of
    /// U_Oᵀ V U_O in the ±1 ORBIT basis, where the projection is an exact rearrangement. An
    /// anti-palindromic direction returns exactly 0.0 by R V R = −V; a palindromic one returns a number of
    /// the coupling's size. A direction that is neither returns something in between and is not covered by
    /// the theorem. Throws on a vanishing perturbation rather than returning a zero that would be
    /// indistinguishable from the theorem's.</summary>
    public static double DirectionReach(double[] pattern)
    {
        if (pattern is null) throw new ArgumentNullException(nameof(pattern));
        var q = PinnedQ;
        var build = Matrix<Complex>.Build;
        var baseFull = build.DenseOfArray(BuildFull(q, Uniform()));
        var moved = build.DenseOfArray(BuildFull(q, pattern.Select(p => 1.0 + p).ToArray()));
        var v = moved - baseFull;
        if (v.Enumerate().All(z => z == Complex.Zero))
            throw new ArgumentException("the perturbation vanishes, so there is no direction to read",
                                        nameof(pattern));
        var u = OddSectorColumnsExact();
        return (u.ConjugateTranspose() * v * u).Enumerate().Max(z => z.Magnitude);
    }

    /// <summary>‖R·L − L·R‖ at the pinned coupling on the given bond profile, with R the site reflection
    /// as a permutation. Both sides are the same float values rearranged, so this is an exact route and a
    /// nonzero value is a fact about the profile rather than about rounding. It is 0.0 on a palindromic
    /// profile, where the R-odd subspace is L-invariant and its restriction is a sub-spectrum, and nonzero
    /// off it, where the same restriction is only a compression. The theorem about anti-palindromic
    /// directions does NOT need this; the meaning of the restricted singular values does.</summary>
    public static double ReflectionCommutator(double[] bondWeights)
    {
        var build = Matrix<Complex>.Build;
        var l = build.DenseOfArray(BuildFull(PinnedQ, bondWeights ?? Uniform()));
        var r = XxzCoherenceBlock.ReflectionPermutation(Sites);
        return (r * l - l * r).Enumerate().Max(z => z.Magnitude);
    }

    /// <summary>Is a bond direction palindromic under the chain reflection, w[i] = w[n−1−i]?</summary>
    public static bool IsPalindromic(double[] pattern)
    {
        if (pattern is null) throw new ArgumentNullException(nameof(pattern));
        for (int i = 0; i < pattern.Length; i++)
            if (pattern[i] != pattern[pattern.Length - 1 - i]) return false;
        return true;
    }

    /// <summary>Is a bond direction ANTI-palindromic, w[i] = −w[n−1−i]? This is the predicate the theorem
    /// needs, and it is NOT the negation of <see cref="IsPalindromic"/>: a direction can be neither, and
    /// such a direction does reach the sector. [1,2,0,0] is neither and reaches at 0.67.</summary>
    public static bool IsAntiPalindromic(double[] pattern)
    {
        if (pattern is null) throw new ArgumentNullException(nameof(pattern));
        for (int i = 0; i < pattern.Length; i++)
            if (pattern[i] != -pattern[pattern.Length - 1 - i]) return false;
        return true;
    }

    /// <summary>The three smallest singular values of (L_O − λ₀·I), scaled by ‖L_O‖. Two below the floor
    /// and the third clear is the nullity-2 reading the character verdict rests on.</summary>
    public (double S1, double S2, double S3) SingularTriple
    {
        get
        {
            var odd = OddBlock(PinnedQ, Uniform());
            var shifted = odd - Matrix<Complex>.Build.DenseIdentity(odd.RowCount) * new Complex(Lambda0, 0);
            var s = shifted.Svd(false).S.Enumerate().Select(x => x.Real).OrderBy(x => x).Take(3).ToArray();
            double norm = odd.L2Norm();
            return (s[0] / norm, s[1] / norm, s[2] / norm);
        }
    }

    public InspectablePayload Payload => InspectablePayload.Empty;

    public IEnumerable<IInspectable> Children
    {
        get
        {
            var (full, even, odd) = SectorDimensions;
            yield return new InspectableNode("the sector the argument sits on",
                summary: $"the (1,2) coherence block of the open N=5 chain splits under the spatial " +
                         $"reflection into {full} = {even} even + {odd} odd; the A₂ layer the proof " +
                         $"consumes is the odd one, of degree 13.");

            yield return new InspectableNode("the coupling book, pinned by spectrum/operator reconstruction",
                summary: $"{PinnedBook}. Scanned over the three candidate scalings rather than assumed, " +
                         $"because the repository's q and Q differ by a factor of two.",
                    provenance: NodeProvenance.Live);

            var n6 = N6RealQExclusion;
            yield return new InspectableNode("N=6 real-q exclusion from exact rational boxes",
                summary: $"{n6.Summary}. Since t=i q, q is real iff Re(t)=0. " +
                         $"This exact interval decision does not inspect display midpoints.",
                provenance: NodeProvenance.Live);

            double axis = HermitianAxisResidual;
            yield return new InspectableNode("step 2: the orbit meets the Hermitian axis",
                summary: $"at a real t the R-odd block is real symmetric, so geometric multiplicity equals " +
                         $"algebraic multiplicity for every one of its eigenvalues. Residual " +
                         $"max(|Im L|, |L − Lᵀ|) = {axis.ToString("0.###e+00", CultureInfo.InvariantCulture)}, " +
                         $"read against 0.0 exactly because the route to it is exact. Twelve of the orbit's " +
                         $"26 members are real t, so the orbit meets this axis.",
                    provenance: NodeProvenance.Live);

            yield return new InspectableNode("step 1: the orbit bookkeeping (consumed, not recomputed here)",
                summary: "A₂_O is irreducible over ℚ of degree 13, so [ℚ(t₀) : ℚ] = 26 and the 26 lifts form " +
                         "one Galois orbit: 1 positive w gives 2 imaginary t, 6 negative w give 12 real t, " +
                         "6 nonreal w give 12 nonreal t. That is exact integer arithmetic on the A₂_O " +
                         "coefficients, certified by simulations/o2b_gcd_certificate.py; this witness takes it " +
                         "as input and recomputes the structure it is applied to.");

            var (s1, s2, s3) = SingularTriple;
            yield return new InspectableNode("the reading at the physical point",
                summary: $"the three smallest scaled singular values of (L_O − λ₀·I) are " +
                         $"{s1.ToString("0.###e+00", CultureInfo.InvariantCulture)}, " +
                         $"{s2.ToString("0.###e+00", CultureInfo.InvariantCulture)}, " +
                         $"{s3.ToString("0.###e+00", CultureInfo.InvariantCulture)}: two on the floor and the " +
                         $"third clear, a nullity of 2. This corroborates the field argument; it does not " +
                         $"carry it, and a defective point with a weak Jordan coupling would read the same way.",
                    provenance: NodeProvenance.Live);

            yield return new InspectableNode("the sector is invariant, not merely a compression",
                summary: $"max|RL − LR| = " +
                         $"{ReflectionCommutator(null).ToString("0.###e+00", CultureInfo.InvariantCulture)} on " +
                         $"the uniform profile, an exact route because both sides are the same floats " +
                         $"rearranged. That is what makes the restricted singular values a sub-spectrum. On " +
                         $"the profile [1, 1.3, 0.7, 1.9] it reads " +
                         $"{ReflectionCommutator(new[] { 1.0, 1.3, 0.7, 1.9 }).ToString("0.###", CultureInfo.InvariantCulture)}, " +
                         $"so there the same restriction would only be a compression.",
                provenance: NodeProvenance.Live);

            foreach (var pattern in new[]
                     {
                         new[] { 1.0, 1.0, 1.0, 1.0 }, new[] { 1.0, 0.0, 0.0, 1.0 },
                         new[] { 1.0, -1.0, -1.0, 1.0 }, new[] { 1.0, -1.0, 1.0, -1.0 },
                         new[] { 1.0, 0.0, 0.0, -1.0 }, new[] { 0.0, 1.0, -1.0, 0.0 },
                         new[] { 1.0, 2.0, 0.0, 0.0 },
                     })
            {
                string label = "[" + string.Join(",", pattern.Select(x => x.ToString("0", CultureInfo.InvariantCulture))) + "]";
                bool anti = IsAntiPalindromic(pattern);
                bool pal = IsPalindromic(pattern);
                string kind = anti ? "anti-palindromic" : pal ? "palindromic" : "neither";
                yield return new InspectableNode($"detune direction {label} ({kind})",
                    summary: $"largest entry of U_Oᵀ V U_O in the ±1 orbit basis = " +
                             $"{DirectionReach(pattern).ToString("0.###e+00", CultureInfo.InvariantCulture)}. " +
                             (anti ? "Exactly zero, by R V R = −V, independently of the base point. Reflection " +
                                     "symmetry of the base is needed only to make the parity subspaces invariant. " +
                                     "The physical pair still splits, through the even sector at second " +
                                     "order: a blind spot of the sector reading, not a null response."
                              : pal ? "Of the coupling's size, so the response here is linear, with a constant " +
                                      "of its own that no parity fixes."
                              : "Neither palindromic nor anti-palindromic, so the theorem says nothing about it " +
                                "and it does reach the sector. The two named classes are not a dichotomy."),
                    provenance: NodeProvenance.Live);
            }

            yield return new InspectableNode("what this witness cannot see",
                summary: "the commutator sign is a CONSTRUCTION choice here, ket −2iq and bra +2iq, and no " +
                         "reading below can check it: flipping it conjugates the R-odd spectrum, and λ₀ is " +
                         "real, so the pinned book, the singular triple and the axis residual all come back " +
                         "identical. The pin certifies the operator only up to that conjugation. Neither does " +
                         "anything here bound a Jordan coupling from below, so a defective point with a weak " +
                         "enough coupling would read this same nullity 2. The field argument rests on neither: " +
                         "its two inputs are the exact certificate, not this builder.",
                provenance: NodeProvenance.Live);
        }
    }
}
