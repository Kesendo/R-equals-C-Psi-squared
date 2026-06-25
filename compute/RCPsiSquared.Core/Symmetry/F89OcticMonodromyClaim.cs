using RCPsiSquared.Core.Knowledge;

namespace RCPsiSquared.Core.Symmetry;

/// <summary>F89 path-3 octic monodromy generates S_8 (Tier 1 derived): the geometric/spectral route to
/// Gal(F_8) = S_8, the companion to the algebraic Frobenius route (<see cref="F89Path3OcticGaloisClaim"/>).
///
/// <para>As q = J/γ loops the complex plane the eight octic rates braid. A lasso around each genuine EP
/// (a simple zero of P_10, a defective √-branch point) reads a transposition of the 8 strands in one
/// common base labelling (base q0 = 2); the diabolic q_EP is silent (a double discriminant zero, loop =
/// identity). The transposition graph on the 8 strands is CONNECTED, so the transpositions generate the
/// full symmetric group: Gal(F_8) = S_8, reconstructed purely from eigenvalue braids, monodromy = Galois
/// from below, independent of the algebraic specialization + Dedekind + Jordan certificate.</para>
///
/// <para>Live witness <c>inspect --root galoismonodromy</c> (<c>GaloisMonodromyWitness</c>, gate G3:
/// components = 1, all 8 strands one orbit). Anchors: <c>experiments/F89_MONODROMY_MIRROR.md</c>.</para></summary>
public sealed class F89OcticMonodromyClaim : Claim
{
    // Parent-edge marker (the algebraic S_8 this geometric route reproduces).
    public F89Path3OcticGaloisClaim Galois { get; }
    // Parent-edge marker (the EP locations the lassos enclose; the silent diabolic).
    public F89Path3OcticEpClaim Ep { get; }

    public F89OcticMonodromyClaim(F89Path3OcticGaloisClaim galois, F89Path3OcticEpClaim ep)
        : base("F89 path-3 octic monodromy generates S_8: lassoing every genuine EP from a common base reads a transposition of the 8 strands; the transposition graph is connected ⟹ Gal(F_8) = S_8 reconstructed from eigenvalue braids (monodromy = Galois, from below), the geometric companion to the algebraic Frobenius certificate; the diabolic q_EP is silent (loop = identity)",
               Tier.Tier1Derived,
               "experiments/F89_MONODROMY_MIRROR.md + " +
               "compute/RCPsiSquared.Diagnostics/Foundation/GaloisMonodromyWitness.cs (inspect --root galoismonodromy, G3)")
    {
        Galois = galois ?? throw new ArgumentNullException(nameof(galois));
        Ep = ep ?? throw new ArgumentNullException(nameof(ep));
    }

    public override string DisplayName =>
        "F89 octic monodromy = S_8 from below: EP braids generate the Galois group (geometric route)";

    public override string Summary =>
        $"the octic's Galois group reconstructed from eigenvalue braids: every EP lassoed from a common base is a transposition, the graph on the 8 strands is connected ⟹ Gal(F_8) = S_8 (the diabolic is silent); monodromy = Galois, independent of the algebraic Frobenius route ({Tier.Label()})";
}
