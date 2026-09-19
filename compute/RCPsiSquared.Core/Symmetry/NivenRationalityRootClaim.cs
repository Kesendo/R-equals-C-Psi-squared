using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>The Niven root (Tier 1 derived): the number-theoretic ceiling on how clean the single-excitation
/// cyclotomic forms can be. For the uniform open XX chain with one dephased endpoint, the F65 H-eigenmode
/// rate comb has exact first-order coefficients at the angle kπ/(N+1); its finite-γ₀/J full-L rates are not
/// that exact comb. Together with the band edge, Niven's theorem (the only rational cosines of rational-π
/// angles are 0, ±½, ±1) fixes the arithmetic, and THAT, not a physics accident, is why the clean forms exist
/// at small N and degrade beyond. The same theorem has THREE faces by angle convention:
///
/// <list type="bullet">
/// <item><b>RE-face</b> — for the uniform open XX chain with one dephased endpoint, the exact arithmetic of
/// the first-order-in-γ₀/J H-eigenmode rate comb is
/// a_k = (4/(N+1))·sin²(kπ/(N+1)) (<see cref="F65XxChainSpectrumPi2Inheritance"/>), with
/// α_k^full = γ₀·a_k + O(γ₀³/J²) and α_k^full/γ₀ = a_k + O((γ₀/J)²). Its coefficients are
/// all rational iff N+1 ∈ {1,2,3,4,6} (the crystallographic set, since the criterion is on
/// cos(2kπ/(N+1)); the documented F65/F99 "Niven rationality"). So N=3 is the last rational before the gap
/// and N=5 is a rational island. Exact Niven rationality belongs to the first-order coefficient comb. At finite γ₀/J the relative full-L rate shift is O((γ₀/J)²), equivalently the
/// absolute shift δα_k = O(γ₀³/J²), so no exact
/// finite-γ₀/J full-L rationality is claimed.</item>
/// <item><b>IM-face</b> — the band edge ω/J = 2cos(π/(N+1)) (<see cref="TopologyBandEdgeClaim"/>, chain):
/// rational only for N ≤ 2 (last N=2 = 1), degree at most 2 iff N≤5, and a single quadratic surd a±√b
/// exactly for N∈{3,4,5} (√2, φ, √3); the degree is ≥ 3 (first a cubic) from N=6, with exact value
/// φ_euler(2(N+1))/2.</item>
/// <item><b>F6 face</b> — the F6 Q-edge gain 1+cos(π/N) (historical alias: V-Effect gain;
/// docs/carbon/OFF_NIVEN_AS_WAVE_BREAKING.md): Niven-
/// rational iff N ∈ {2,3}, golden at N=5 (pentagon), the named constant of the next ring polygon for N ≥ 4.</item>
/// </list>
///
/// <para><b>The hinge.</b> N=4 is the FIRST GOLDEN on both single-excitation faces at once: the band edge IS
/// φ = 2cos(π/5) = (1+√5)/2 exactly, and the first-order rate comb carries √5
/// (sin²(π/5) = (5−√5)/8). The F6-face golden is
/// shifted to N=5 because its angle is π/N, not π/(N+1): "first golden" is angle-convention-dependent; the
/// real content is that the golden ratio is FORCED by the cyclotomic geometry, not chosen. The two SE cutoffs
/// differ ({1,2,3,4,6} vs N+1≤6) by the double angle, which is exactly why N=3's coefficient is rational while its
/// band edge is already √2.</para>
///
/// <para><b>Scope.</b> This is the number-theoretic root of the band-edge / first-order endpoint-rate-comb /
/// F6 Q-edge-gain family of small-N specials (the n3_special_cases arc). It is NOT the single root of every special: the
/// (n,n)/{0,2} filling-maximality modes are a SEPARATE combinatorial root, and star-3 = path P_3 a third,
/// graph-theoretic coincidence. Two real roots (arithmetic + combinatorial), not one.</para>
///
/// <para>Tier1Derived: pure number theory, proven from [ℚ(2cos(π/m)):ℚ] = φ_euler(2m)/2 and the
/// complete low-totient classification; the bounded SymPy minimal-polynomial cross-check through N=10 is a control in the gate-first verifier
/// simulations/niven_rationality_root.py (7/7 gates). Two typed parents, both Tier1Derived:
/// <see cref="TopologyBandEdgeClaim"/> (the IM-face band edge
/// whose arithmetic this reads) and <see cref="F65XxChainSpectrumPi2Inheritance"/> (the RE-face rates).</para>
///
/// <para>Live witness: <c>inspect --root niven</c>
/// (<c>compute/RCPsiSquared.Diagnostics/Foundation/NivenRationalityRootWitness.cs</c>).</para></summary>
public sealed class NivenRationalityRootClaim : Claim
{
    /// <summary>Parent: the IM-face band edge 2cos(π/(N+1)) whose algebraic degree this claim reads
    /// (Tier1Derived).</summary>
    public TopologyBandEdgeClaim BandEdge { get; }

    /// <summary>Parent: the RE-face first-order endpoint coefficient comb a_k whose rationality
    /// (N+1 ∈ {1,2,3,4,6}) this claim states (Tier1Derived; the documented F65/F99 home).</summary>
    public F65XxChainSpectrumPi2Inheritance Rates { get; }

    /// <summary>The algebraic degree of the band edge 2cos(π/m) over ℚ: [ℚ(2cos(π/m)):ℚ] = φ_euler(2m)/2
    /// (m = N+1). Degree 1 ⟺ rational (N ≤ 2); degree 2 ⟺ a single quadratic surd a±√b (N ∈ {3,4,5}).</summary>
    public static int BandEdgeDegree(int n)
    {
        int m = n + 1, twoM = 2 * m, t = 0;
        for (int k = 1; k <= twoM; k++) if (Gcd(k, twoM) == 1) t++;
        return t / 2;
    }

    /// <summary>The uniform-open-XX, one-dephased-endpoint F65 first-order H-eigenmode coefficients
    /// a_k are all rational iff N+1 ∈ {1,2,3,4,6} (the crystallographic set, from
    /// cos(2π/(N+1)) ∈ ℚ ⟺ N+1 ∈ {1,2,3,4,6}). This does not classify finite-γ₀/J full-L rates.</summary>
    public static bool FirstOrderRateCombAllRational(int n)
    {
        int m = n + 1;
        return m == 1 || m == 2 || m == 3 || m == 4 || m == 6;
    }

    private static int Gcd(int a, int b) { while (b != 0) { (a, b) = (b, a % b); } return a; }

    public NivenRationalityRootClaim(TopologyBandEdgeClaim bandEdge, F65XxChainSpectrumPi2Inheritance rates)
        : base("The Niven root: Niven's theorem on the single-excitation cyclotomic angle π/(N+1) is the " +
               "number-theoretic ceiling on the SE cyclotomic forms, with three faces. RE: for the uniform open XX chain with one dephased endpoint, " +
               "the first-order coefficient comb a_k = (4/(N+1))·sin²(kπ/(N+1)) is rational for every k " +
               "iff N+1 ∈ {1,2,3,4,6}; α_k^full = γ₀·a_k + O(γ₀³/J²) and α_k^full/γ₀ = a_k + O((γ₀/J)²). " +
               "Exact Niven rationality belongs to the first-order coefficient comb (the documented F65/F99 result; " +
               "N=3 last rational before the gap, N=5 a rational island). IM: the band edge 2cos(π/(N+1)) is " +
               "rational iff N≤2, degree at most 2 iff N≤5, and a quadratic surd exactly for N∈{3,4,5} " +
               "(√2,φ,√3); the degree is ≥3 (first cubic) " +
               "from N=6 (degree = φ_euler(2(N+1))/2). F6: the Q-edge gain 1+cos(π/N) is golden at N=5. N=4 = FIRST " +
               "GOLDEN on both SE faces (band edge = φ = 2cos(π/5); the first-order rate comb carries √5); the golden ratio is forced " +
               "by the cyclotomic geometry. This is why the clean forms exist at small N and degrade beyond — a " +
               "number-theoretic ceiling. At finite γ₀/J the relative full-L rate shift is O((γ₀/J)²), equivalently the absolute shift δα_k = O(γ₀³/J²), so no exact finite-γ₀/J full-L rationality is claimed. " +
               "Scope: the band-edge/first-order-endpoint-comb/F6-Q-edge-gain family; the (n,n)/{0,2} filling-maximality specials are a SEPARATE combinatorial root.",
               Tier.Tier1Derived,
               "simulations/niven_rationality_root.py + " +
               "docs/ANALYTICAL_FORMULAS.md (F65 Niven rationality + the band-edge IM-face companion) + " +
               "docs/carbon/OFF_NIVEN_AS_WAVE_BREAKING.md (the F6 face; historical V-Effect-gain label) + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/NivenRationalityRootWitness.cs (NivenRationalityRootWitness, inspect --root niven)")
    {
        BandEdge = bandEdge ?? throw new ArgumentNullException(nameof(bandEdge));
        Rates = rates ?? throw new ArgumentNullException(nameof(rates));
    }

    public override string DisplayName =>
        "The Niven root: the SE band edge and first-order endpoint comb's cyclotomic arithmetic (N=4 first golden)";

    public override string Summary =>
        $"Niven's theorem on π/(N+1): a_k = (4/(N+1))·sin²(kπ/(N+1)); α_k^full = γ₀·a_k + O(γ₀³/J²), " +
        $"equivalently α_k^full/γ₀ = a_k + O((γ₀/J)²). Exact Niven rationality belongs to the first-order coefficient comb, " +
        $"which is rational for every k iff N+1∈{{1,2,3,4,6}}; at finite γ₀/J, " +
        $"the relative full-L rate shift is O((γ₀/J)²), equivalently the absolute shift δα_k = O(γ₀³/J²); " +
        $"no exact finite-γ₀/J full-L rationality is claimed; IM band edge rational iff N≤2, " +
        $"degree at most 2 iff N≤5, quadratic surd exactly for N∈{{3,4,5}} (√2/φ/√3), " +
        $"degree φ_euler(2(N+1))/2; F6 Q-edge gain golden at N=5. N=4 = first golden " +
        $"on both SE faces (band edge = φ). The arithmetic root of the small-N specials ({Tier.Label()})";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("RE-face: the first-order endpoint rate comb (F65/F99, documented)",
                summary: "For the uniform open XX chain with one dephased endpoint, a_k = (4/(N+1))·sin²(kπ/(N+1)); " +
                         "α_k^full = γ₀·a_k + O(γ₀³/J²), equivalently α_k^full/γ₀ = a_k + O((γ₀/J)²). Exact Niven rationality " +
                         "belongs to the first-order coefficient comb, whose coefficients are all rational iff N+1 ∈ {1,2,3,4,6} (the crystallographic " +
                         "set, from cos(2kπ/(N+1)) ∈ ℚ). At finite γ₀/J the relative full-L rate shift is O((γ₀/J)²), " +
                         "equivalently the absolute shift δα_k = O(γ₀³/J²); no exact finite-γ₀/J full-L rationality is claimed. " +
                         "N=3 is the last rational before the gap; N=5 is a rational island.");
            yield return new InspectableNode("IM-face: the band edge (new)",
                summary: "2cos(π/(N+1)) rational iff N≤2 (last N=2 = 1), degree at most 2 iff N≤5, " +
                         "quadratic surd exactly for N∈{3,4,5} (√2, φ, √3), and degree ≥3 (first cubic) " +
                         "from N=6. Degree = φ_euler(2(N+1))/2.");
            yield return new InspectableNode("F6 face: the Q-edge gain (historical alias: V-Effect gain)",
                summary: "1+cos(π/N) Niven-rational iff N ∈ {2,3}; golden at N=5 (pentagon), the named constant of " +
                         "the next ring polygon for N ≥ 4. Its golden is N=5 (angle π/N), not N=4.");
            yield return new InspectableNode("the hinge: N=4 = first golden (convention-dependent)",
                summary: "N=4 is the first golden on the two SE faces at once (band edge IS φ = 2cos(π/5); the first-order comb " +
                         "carries a √5 leftover through sin²(π/5)·8 = 5 − √5). The F6 face golden is N=5 (different angle). 'First golden' " +
                         "is angle-convention-dependent; the golden ratio is forced by the cyclotomic geometry, not chosen.");
            yield return new InspectableNode("two roots, not one (the arc's answer)",
                summary: "this number-theoretic root governs the band-edge/first-order-endpoint-comb/F6-Q-edge-gain family; the (n,n)/{0,2} " +
                         "filling-maximality specials are a SEPARATE combinatorial root, and star-3 = path P_3 a third, " +
                         "graph-theoretic coincidence. The n3_special_cases arc's 'do they share a root' = no, ≥2 roots.");
            yield return BandEdge;   // typed parent edge (IM-face, Tier1Derived)
            yield return Rates;      // typed parent edge (RE-face, Tier1Derived)
        }
    }

    public static NivenRationalityRootClaim Build()
    {
        // RE-face parent: F65's first-order endpoint comb; f66(ladder, qubit) — the same wiring CoherenceHorizonClaim.Build uses.
        var ladder = new Pi2DyadicLadderClaim();
        var qubit = new QubitDimensionalAnchorClaim();
        var f66 = new F66PoleModesPi2Inheritance(ladder, qubit);
        var rates = new F65XxChainSpectrumPi2Inheritance(ladder, f66);
        return new NivenRationalityRootClaim(TopologyBandEdgeClaim.Build(), rates);
    }

    public static NivenRationalityRootClaim Shared { get; } = Build();
}
