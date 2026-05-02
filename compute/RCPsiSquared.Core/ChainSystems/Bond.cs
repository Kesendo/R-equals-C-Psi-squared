namespace RCPsiSquared.Core.ChainSystems;

/// <summary>A single bond between two qubits with a per-bond coupling.</summary>
public record Bond(int Site1, int Site2, double Coupling);
