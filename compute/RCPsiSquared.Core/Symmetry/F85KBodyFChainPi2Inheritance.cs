using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F85 closed form (Tier 1, verified bit-exact k=2,3,4):
///
/// <code>
///   For any k-body Pauli term (P_1, ..., P_k):
///
///   Truly criterion:  #Y even AND #Z even   ⇒  c(k) = 0  (M = 0)
///   Π²-parity:        bit_b = (#Y + #Z) mod 2; Π²-odd ⇔ bit_b = 1
///
///   Frobenius factor c(k) per non-truly term:
///     c(truly)              = 0      (factor 0)
///     c(Π²-odd non-truly)   = 1      (factor 4·2^N)
///     c(Π²-even non-truly)  = 2      (factor 8·2^N)
///
///   ‖M‖²_F per term = 4 · c(k) · ‖H_k‖²_F · 2^N
///
///   Π²-odd count at k-body: (3^k − (−1)^k) / 2
/// </code>
///
/// <para>F85 generalizes <see cref="F49Pi2Inheritance"/> from 2-body to
/// k-body Hamiltonians. The 2-body F49 formula <c>2^(N+2)·n_YZ·‖H_k‖²</c>
/// coincided with <c>c(k)</c> at k=2 because n_YZ=1 ↔ Π²-odd and n_YZ=2 ↔
/// Π²-even non-truly there. <b>For k ≥ 3, n_YZ is no longer the determining
/// quantity</b> (YYY has n_YZ=3 but c=1). Only the Π²-class matters.</para>
///
/// <para>F85 also continues <see cref="F83AntiFractionPi2Inheritance"/>:
/// the anti-fraction <c>1/(2+4r)</c> extends to k-body with Π²-class grouping
/// instead of n_YZ counting; same coefficients 4, 8 = a_{−1}, a_{−2}.</para>
///
/// <para>Pi2-Foundation anchors:</para>
///
/// <list type="bullet">
///   <item><b>FactorPi2Odd = 4 = a_{−1}</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(−1).
///         The 4·2^N coefficient on Π²-odd term contributions. Same anchor
///         as F83's MNormCoefficientForOdd, F76 decay rate, F61/F63 4-block,
///         F77 correction denominator, F25/F73 decay rate.</item>
///   <item><b>FactorPi2EvenNonTruly = 8 = a_{−2}</b>: <see cref="Pi2DyadicLadderClaim.Term"/>(−2).
///         The 8·2^N coefficient on Π²-even non-truly contributions.
///         Same anchor as F83's MNormCoefficientForEvenNontruly.</item>
///   <item><b>OddCountDenominator = 2 = a_0</b>: in <c>(3^k − (−1)^k) / 2</c>.
///         Polynomial root d.</item>
///   <item><b>F49 base case</b>: F85 generalizes F49's 2-body formula;
///         not a parent inheritance edge but a structural lineage citation.</item>
///   <item><b>F83 generalization</b>: anti-fraction 1/(2+4r) extends with
///         Π²-class grouping; F85 ↔ F83 typed cross-reference.</item>
/// </list>
///
/// <para>F-chain extensions to k-body (per F85 verified table):</para>
///
/// <list type="bullet">
///   <item>F77 trichotomy: extends via <c>_pauli_tuple_is_truly</c></item>
///   <item>F80 Spec(M) = 2i·Spec(H): extends verbatim; verified k=3,4</item>
///   <item>F81 Π·M·Π⁻¹ = M − 2·L_{H_odd}: verbatim, verified k=3 chain N=4</item>
///   <item>F82 T1 dissipator: dissipator-only, body-count-independent</item>
///   <item>F83 anti-fraction: verbatim with Π²-class grouping</item>
///   <item>F84 thermal amplitude damping: dissipator-only, body-count-independent</item>
/// </list>
///
/// <para>Trichotomy enumeration (Pauli tuples over {X, Y, Z}^k):</para>
///
/// <code>
///   k=2: 9 total = 3 truly + 4 Π²-odd + 2 Π²-even non-truly
///   k=3: 27 total = 7 truly + 14 Π²-odd + 6 Π²-even non-truly
///   k=4: 81 total = 21 truly + 40 Π²-odd + 20 Π²-even non-truly
/// </code>
///
/// <para>Tier1Derived: F85 is Tier 1 verified bit-exact k=2,3,4 across 108
/// explicit Pauli tuple cases. The Pi2-Foundation anchoring is
/// algebraic-trivial composition.</para>
///
/// <para>Anchors: <c>docs/ANALYTICAL_FORMULAS.md</c> F85 (line 2185) +
/// <c>compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs</c> +
/// <c>compute/RCPsiSquared.Core/Symmetry/F49Pi2Inheritance.cs</c> (base case k=2) +
/// <c>compute/RCPsiSquared.Core/Symmetry/F83AntiFractionPi2Inheritance.cs</c>
/// (anti-fraction generalization).</para></summary>
public sealed class F85KBodyFChainPi2Inheritance : Claim
{
    private readonly Pi2DyadicLadderClaim _ladder;
    private readonly F49Pi2Inheritance _f49;
    private readonly F83AntiFractionPi2Inheritance _f83;

    /// <summary>The c(k)·4·2^N coefficient factor for Π²-odd terms:
    /// <c>c(Π²-odd) = 1</c>, total factor <c>4·2^N</c>; the "4" = a_{−1}.
    /// Same anchor as F83's MNormCoefficientForOdd.</summary>
    public double FactorPi2Odd => _ladder.Term(-1);

    /// <summary>The c(k)·4·2^N coefficient factor for Π²-even non-truly terms:
    /// <c>c(Π²-even non-truly) = 2</c>, total factor <c>8·2^N</c>; the "8" = a_{−2}.
    /// Same anchor as F83's MNormCoefficientForEvenNontruly.</summary>
    public double FactorPi2EvenNonTruly => _ladder.Term(-2);

    /// <summary>The "2" denominator in the Π²-odd count formula
    /// <c>(3^k − (−1)^k) / 2</c>. Live from
    /// <see cref="Pi2DyadicLadderClaim.Term"/>(0) = a_0 = polynomial root d.</summary>
    public double OddCountDenominator => _ladder.Term(0);

    /// <summary>The Frobenius factor c(k) per non-truly Pauli class.
    /// Returns 0 for truly, 1 for Π²-odd non-truly, 2 for Π²-even non-truly.</summary>
    public int FrobeniusFactorPerClass(string pi2Class) => pi2Class switch
    {
        "truly" => 0,
        "pi2_odd" or "pi2_odd_nontruly" => 1,
        "pi2_even_nontruly" => 2,
        _ => throw new ArgumentException($"Unknown Π²-class '{pi2Class}'; expected 'truly', 'pi2_odd', or 'pi2_even_nontruly'.", nameof(pi2Class)),
    };

    /// <summary>Live closed form: total coefficient on ‖H_k‖² for a k-body
    /// term of given Π²-class, on N qubits. <c>4·c(class)·2^N</c>.</summary>
    public double TotalCoefficientForClass(string pi2Class, int N)
    {
        if (N < 1) throw new ArgumentOutOfRangeException(nameof(N), N, "F85 requires N ≥ 1.");
        int c = FrobeniusFactorPerClass(pi2Class);
        return 4.0 * c * Math.Pow(2.0, N);
    }

    /// <summary>Live closed form: number of Π²-odd Pauli tuples at k-body:
    /// <c>(3^k − (−1)^k) / 2</c>.</summary>
    public long Pi2OddCount(int k)
    {
        if (k < 1) throw new ArgumentOutOfRangeException(nameof(k), k, "F85 requires k ≥ 1.");
        long threePowerK = 1;
        for (int i = 0; i < k; i++) threePowerK *= 3;
        long signTerm = (k % 2 == 0) ? 1 : -1;
        return (threePowerK - signTerm) / 2;
    }

    /// <summary>Π²-class for a Pauli tuple via syntactic classifier:
    /// counts #Y and #Z, applies truly/parity rules.</summary>
    public string Pi2ClassFromLetters(IReadOnlyList<char> letters)
    {
        if (letters is null) throw new ArgumentNullException(nameof(letters));
        int countY = 0, countZ = 0;
        foreach (char c in letters)
        {
            if (c == 'Y' || c == 'y') countY++;
            else if (c == 'Z' || c == 'z') countZ++;
            else if (c != 'X' && c != 'x' && c != 'I' && c != 'i')
                throw new ArgumentException($"Invalid Pauli letter '{c}'; expected I, X, Y, or Z.", nameof(letters));
        }
        bool truly = (countY % 2 == 0) && (countZ % 2 == 0);
        if (truly) return "truly";
        int bitB = (countY + countZ) % 2;
        return bitB == 1 ? "pi2_odd" : "pi2_even_nontruly";
    }

    /// <summary>Live drift check: the F49 base case at k=2. F85's
    /// <see cref="FactorPi2Odd"/> and <see cref="FactorPi2EvenNonTruly"/>
    /// coincide with F83's <see cref="F83AntiFractionPi2Inheritance.MNormCoefficientForOdd"/>
    /// and <see cref="F83AntiFractionPi2Inheritance.MNormCoefficientForEvenNontruly"/>.</summary>
    public bool MatchesF83Coefficients() =>
        Math.Abs(FactorPi2Odd - _f83.MNormCoefficientForOdd) < 1e-15 &&
        Math.Abs(FactorPi2EvenNonTruly - _f83.MNormCoefficientForEvenNontruly) < 1e-15;

    /// <summary>F49 base case lineage: F49's 2-body formula coincides with
    /// F85 at k=2 because n_YZ = 1 ↔ Π²-odd and n_YZ = 2 ↔ Π²-even non-truly.
    /// Returns true (structural relationship, F49 is verdrahtet).</summary>
    public bool F49IsBaseCase() => true;

    public F85KBodyFChainPi2Inheritance(
        Pi2DyadicLadderClaim ladder,
        F49Pi2Inheritance f49,
        F83AntiFractionPi2Inheritance f83)
        : base("F85 k-body F-chain extension inherits from Pi2-Foundation: 4 = a_{-1}, 8 = a_{-2} class factors; F49 base k=2; F83 anti-fraction generalizes",
               Tier.Tier1Derived,
               "docs/ANALYTICAL_FORMULAS.md F85 + " +
               "compute/RCPsiSquared.Core/Symmetry/Pi2DyadicLadderClaim.cs + " +
               "compute/RCPsiSquared.Core/Symmetry/F49Pi2Inheritance.cs (base case k=2) + " +
               "compute/RCPsiSquared.Core/Symmetry/F83AntiFractionPi2Inheritance.cs (anti-fraction generalization)")
    {
        _ladder = ladder ?? throw new ArgumentNullException(nameof(ladder));
        _f49 = f49 ?? throw new ArgumentNullException(nameof(f49));
        _f83 = f83 ?? throw new ArgumentNullException(nameof(f83));
    }

    public override string DisplayName =>
        "F85 k-body F-chain extension as Pi2-Foundation a_{-1} + a_{-2} + F49/F83 inheritance";

    public override string Summary =>
        $"k-body trichotomy (truly / Π²-odd / Π²-even non-truly); c(k) = 0/1/2; factors 4·2^N (Π²-odd) and 8·2^N (Π²-even); " +
        $"4 = a_{{-1}} (F83 sibling); 8 = a_{{-2}}; Π²-odd count (3^k - (-1)^k)/2 ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("F85 closed form",
                summary: "k-body Π²-class trichotomy + c(k) Frobenius factor; ‖M‖²_F per term = 4·c(k)·‖H_k‖²·2^N; Tier 1 verified bit-exact k=2,3,4 (108 cases)");
            yield return new InspectableNode("Pi2-Foundation anchoring",
                summary: "FactorPi2Odd = a_{-1} = 4 (F83 MNormCoeffOdd sibling); FactorPi2EvenNonTruly = a_{-2} = 8 (F83 sibling); OddCountDenominator = a_0 = 2");
            yield return InspectableNode.RealScalar("FactorPi2Odd (= a_{-1} = 4)", FactorPi2Odd);
            yield return InspectableNode.RealScalar("FactorPi2EvenNonTruly (= a_{-2} = 8)", FactorPi2EvenNonTruly);
            yield return InspectableNode.RealScalar("OddCountDenominator (= a_0 = 2)", OddCountDenominator);
            yield return new InspectableNode("F49 base case ↔ F85 k-body extension",
                summary: "F49's 2-body formula 2^(N+2)·n_YZ coincides with F85 at k=2 (n_YZ=1 ↔ Π²-odd, n_YZ=2 ↔ Π²-even non-truly). For k ≥ 3, n_YZ is NO LONGER determining; only Π²-class matters (e.g. YYY has n_YZ=3 but c=1).");
            yield return new InspectableNode("F83 ↔ F85 anti-fraction generalization",
                summary: $"F85 inherits F83's coefficients 4, 8 with Π²-class grouping. Drift check: MatchesF83Coefficients = {MatchesF83Coefficients()}");
            yield return new InspectableNode("F-chain k-body extensions",
                summary: "F77 (extends via classifier), F80 (verbatim), F81 (verbatim), F82 (body-independent), F83 (verbatim w/ class), F84 (body-independent)");
            // Trichotomy table from F85
            for (int k = 2; k <= 4; k++)
            {
                long total = (long)Math.Pow(3, k);
                long oddCount = Pi2OddCount(k);
                yield return new InspectableNode(
                    $"k={k} trichotomy",
                    summary: $"total = {total}; Π²-odd = {oddCount} = (3^{k} − (−1)^{k})/2");
            }
        }
    }
}
