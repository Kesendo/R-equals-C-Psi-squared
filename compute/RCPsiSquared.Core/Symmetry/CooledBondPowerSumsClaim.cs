using System.Numerics;
using RCPsiSquared.Core.Inspection;
using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>Where the pairing under cooling alone has to fail, and where it cannot be seen early
/// (<c>experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md</c>, the seventh power sum and the rings; F137's
/// index entry in <c>docs/ANALYTICAL_FORMULAS.md</c>).
///
/// <para>The setting: a chain of N sites, every bond the same element of the class
/// J_XX·XX + J_XY·XY + J_YX·YX + J_YY·YY + Δ·ZZ, and loss σ⁻ at a rate γ_l ≥ 0 on every site, nothing else;
/// c = Σγ_l. Write S = J_XX + J_YY, A = J_XY − J_YX and W = (J_XX − J_YY)² + (J_XY + J_YX)² (twice the turn
/// invariants s and a, four times |w|²). The power sums p_k = Tr((2L + c)^k) are homogeneous in the
/// couplings and the rates together, and a phase count on the turns about Z (every factor z, z̄, w, w̄ of a
/// bond gives each of its two sites a phase, which must cancel at every site) makes p₁, p₃ and p₅ vanish at
/// every N and every profile: a product that H ↦ H̄ changes has four such factors on two neighbouring bonds, beside an
/// odd power of Δ, and the rates enter at least once (a commutator's spectrum pairs about 0 on its own) and an even
/// number of times. The seventh is the first that can break, and it is</para>
///
/// <para>  p₇ = −7·2¹⁷·4^(N−3)·Δ·S·A·W·(γ₁γ₂ − γ_{N−1}γ_N),</para>
///
/// <para>a sum over three-site windows that telescopes to the two ends: the seventh reads the chain's ends
/// and nothing between them. So with Δ, S, A and W nonzero the spectrum does not pair about −c/2 wherever
/// γ₁γ₂ ≠ γ_{N−1}γ_N.</para>
///
/// <para>On a ring of N ≥ 4 sites with one bond of the class repeated, the parity of the factor counts is
/// the same on neighbouring bonds, so a product that H ↦ H̄ changes carries at least four X/Y factors
/// and the rates enter p_k at most k − 5 times; a polynomial of degree two or less in the rates that keeps
/// every rotation and changes sign under every reflection is zero, so p₃, p₅ and p₇ vanish at every rate
/// profile (on the triangle the two closing products Im z³ and Im z·|w|² carry numbers that two bonds with
/// independent pairs fix to zero). The ninth is the first that can read a ring's reflection.</para>
///
/// <para>Live witness: <c>compute/RCPsiSquared.Diagnostics/Foundation/CooledBondPowerSumsWitness.cs</c>,
/// <c>inspect --root cooledbonds</c>, which builds the generator in exact integer arithmetic on its own
/// road (column-major matrix units, no prime), beside the page's gate
/// <c>simulations/cooled_bonds_pairing.py</c>. F137 itself is not typed, so this claim has no typed parent;
/// its parent in the registry's prose is F137's criterion, the spectrum pairing exactly when exchanging
/// loss and gain leaves it unchanged. Adjacent, not parents: <see cref="PalindromeTwoEndCountClaim"/> (F158, the
/// palindrome criterion by the two ends' kernels, which notes F137's recentring under loss) and
/// <see cref="CorrelatedZDephasingPriceClaim"/> (moments of the spectrum about F137's centre).</para>
///
/// <para>The setting is F137's formal one: a bond that does not conserve the excitation number (J_XX ≠ J_YY or
/// J_XY ≠ −J_YX) admits no frame rotating with the qubits, and the generator is read as written, at zero
/// splitting; XXZ bonds with axial DM terms are read in the rotating frame.</para></summary>
public sealed class CooledBondPowerSumsClaim : Claim
{
    /// <summary>One bond of the class, integer couplings.</summary>
    public readonly record struct Bond(long XX, long XY, long YX, long YY, long ZZ)
    {
        public long S => XX + YY;
        public long A => XY - YX;
        public long W => (XX - YY) * (XX - YY) + (XY + YX) * (XY + YX);

        public override string ToString()
        {
            var parts = new List<string>();
            foreach (var (c, name) in new[] { (XX, "XX"), (XY, "XY"), (YX, "YX"), (YY, "YY"), (ZZ, "ZZ") })
                if (c != 0) parts.Add($"{(c < 0 ? "−" : parts.Count == 0 ? "" : "+")}{(Math.Abs(c) == 1 ? "" : Math.Abs(c).ToString())}{name}");
            return parts.Count == 0 ? "0" : string.Join(" ", parts).Replace(" −", " − ").Replace(" +", " + ");
        }
    }

    /// <summary>p₇ = Tr((2L + c)⁷) on the open chain, the closed form, N ≥ 3 (at N = 2 it is zero, the two
    /// ends being one).</summary>
    public static BigInteger SeventhPowerSum(Bond bond, IReadOnlyList<long> rates)
    {
        if (rates is null) throw new ArgumentNullException(nameof(rates));
        int n = rates.Count;
        if (n < 2) throw new ArgumentOutOfRangeException(nameof(rates), n, "a chain needs N ≥ 2");
        if (n == 2) return BigInteger.Zero;
        BigInteger ends = new BigInteger(rates[0]) * rates[1] - new BigInteger(rates[n - 1]) * rates[n - 2];
        return -7 * BigInteger.Pow(2, 17) * BigInteger.Pow(4, n - 3) * bond.ZZ * bond.S * bond.A * bond.W * ends;
    }

    public CooledBondPowerSumsClaim()
        : base("Under cooling alone a chain of one class bond has p₁ = p₃ = p₅ = 0 and Tr((2L + c)⁷) = " +
               "−7·2¹⁷·4^(N−3)·Δ·S·A·W·(γ₁γ₂ − γ_{N−1}γ_N), reading only the chain's ends; on every ring of one " +
               "repeated class bond p₃, p₅ and p₇ vanish at every rate profile",
               Tier.Tier1Derived,
               "experiments/WHICH_BONDS_PAIR_UNDER_COOLING.md + docs/ANALYTICAL_FORMULAS.md (F137, the index entry)")
    {
    }

    public override string DisplayName =>
        "Cooled bonds: the seventh power sum reads the chain's ends, rings are silent through the seventh";

    public override string Summary =>
        "chain of N sites, one bond J_XX·XX + J_XY·XY + J_YX·YX + J_YY·YY + Δ·ZZ on every edge, loss at γ_l on every site, " +
        "c = Σγ_l, S = J_XX + J_YY, A = J_XY − J_YX, W = (J_XX − J_YY)² + (J_XY + J_YX)²: p₁ = p₃ = p₅ = 0 by the phase " +
        "count, p₇ = −7·2¹⁷·4^(N−3)·Δ·S·A·W·(γ₁γ₂ − γ_{N−1}γ_N), so with Δ, S, A, W nonzero the spectrum breaks wherever " +
        "the end products differ; on every ring p₃ = p₅ = p₇ = 0 at every profile (Tier 1 derived; live witness: inspect --root cooledbonds).";

    protected override IEnumerable<IInspectable> ExtraChildren
    {
        get
        {
            yield return new InspectableNode("the phase count",
                summary: "a turn about Z by θ_l on site l multiplies z, z̄, w and w̄ of a bond by a phase from each of its two sites; " +
                         "a trace survives only if the phases cancel at every site, and the end sites of a chain force every bond to carry " +
                         "an even number of such factors; a product the reversal H ↦ H̄ changes needs at least four, beside an odd power of Δ, " +
                         "so the couplings enter at least five times; the rates enter at least once (a commutator's spectrum pairs about 0 " +
                         "on its own) and, the total degree being odd, an even number of times: none of p₁, p₃ and p₅ has room, so they vanish");
            yield return new InspectableNode("the seventh reads the ends",
                summary: "the seventh is the seventh cumulant, a sum over three-site windows, each contributing γ_{l+1}(γ_l − γ_{l+2}) " +
                         "times −7·2¹⁷·4^(N−3)·Δ·S·A·W; along the chain the windows telescope to γ₁γ₂ − γ_{N−1}γ_N");
            yield return new InspectableNode("the rings",
                summary: "round a ring the factor counts of neighbouring bonds have one parity, all even or all odd (all odd needs N factors), " +
                         "so from N = 4 a reversed product has four X/Y factors and p_k carries at most k − 5 rates; a rotation-invariant, " +
                         "reflection-odd polynomial of degree two or less is zero; on the triangle α·Im z³ + β·Im z·|w|² with α = β = 0 from " +
                         "two bonds with independent pairs");
            yield return new InspectableNode("what this does not say",
                summary: "where γ₁γ₂ = γ_{N−1}γ_N the seventh is silent and decides nothing; the ninth reads the whole chain there " +
                         "(its closed form is on the page, not typed here), and on a ring the ninth is the first sum that can break. " +
                         "For a bond that does not conserve the excitation number the generator is F137's formal one at zero splitting, " +
                         "read as written; the physical bath there has equal cooling and heating, where every H commuting with Z⊗N pairs");
        }
    }
}
