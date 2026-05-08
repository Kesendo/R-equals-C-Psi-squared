using RCPsiSquared.Core.Symmetry;
using RCPsiSquared.Runtime.ObjectManager;

namespace RCPsiSquared.Runtime.PolarityArchitecture;

/// <summary>Schicht-1 wiring of <see cref="Pi2Z4KleinDistinctionClaim"/> — the typed
/// answer to Tom's 2026-05-08 question "are the Z₄ Π-sectors and the Klein four cells
/// the same partition?". Two parent edges:
///
/// <list type="bullet">
///   <item><see cref="Pi2I4MemoryLoopClaim"/>: provides the Z₄ side (Π itself, cyclic
///         order 4, eigenvalues {±1, ±i}).</item>
///   <item><see cref="KleinFourCellClaim"/>: provides the Z₂ × Z₂ side (Π²_Z, Π²_X
///         elementary abelian, two commuting involutions).</item>
/// </list>
///
/// <para>The distinction claim does not assert isomorphism — it asserts the opposite:
/// the two are the only two abelian groups of order 4 (up to isomorphism), and the
/// framework uses both. The squaring map Z₄ → Z₂ is the bridge that connects them
/// via Π² = Π²_Z (one of the Klein generators).</para>
///
/// <para>Tier consistency: all three claims (this and both parents) are Tier1Derived
/// (5 ≥ 5).</para>
///
/// <para>Requires: <see cref="Pi2FamilyRegistration.RegisterPi2Family"/> +
/// <see cref="Pi2DyadicLadderRegistration.RegisterPi2DyadicLadder"/> +
/// <see cref="Pi2I4MemoryLoopRegistration.RegisterPi2I4MemoryLoop"/>. The Klein cell
/// claim comes from the base Pi2 family.</para></summary>
public static class Pi2Z4KleinDistinctionRegistration
{
    public static ClaimRegistryBuilder RegisterPi2Z4KleinDistinction(
        this ClaimRegistryBuilder builder) =>
        builder.Register<Pi2Z4KleinDistinctionClaim>(b =>
        {
            _ = b.Get<Pi2I4MemoryLoopClaim>();    // Z₄ side
            _ = b.Get<KleinFourCellClaim>();       // Z₂ × Z₂ side
            return new Pi2Z4KleinDistinctionClaim();
        });
}
