namespace RCPsiSquared.Core.Decomposition.Views;

/// <summary>The four named modes of <see cref="FourModeBasis"/>: <c>|c_1⟩, |c_3⟩</c>
/// (channel-uniform; where the Dicke probe lives) and <c>|u_0⟩, |v_0⟩</c> (SVD-top inter-channel
/// coupling; where the EP partners live).
/// </summary>
public static class FourModeNames
{
    public const string C1 = "|c_1⟩";
    public const string C3 = "|c_3⟩";
    public const string U0 = "|u_0⟩";
    public const string V0 = "|v_0⟩";

    public static readonly IReadOnlyList<string> All = new[] { C1, C3, U0, V0 };
    public static readonly IReadOnlyList<string> ChannelUniform = new[] { C1, C3 };
    public static readonly IReadOnlyList<string> SvdTop = new[] { U0, V0 };
}
