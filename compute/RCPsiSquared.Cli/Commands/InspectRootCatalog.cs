using RCPsiSquared.Core.Inspection;

namespace RCPsiSquared.Cli.Commands;

/// <summary>The parsed inputs a root factory needs: the raw <see cref="ArgParser"/> plus the
/// common flags <see cref="InspectCommand"/> reads once up front. A factory pulls whatever
/// root-specific options it wants off <see cref="Parser"/>.</summary>
public sealed record InspectRootContext(
    ArgParser Parser,
    int N,
    bool WithQSweep,
    bool WithMeasured,
    int? QGridPoints);

/// <summary>One row of the root catalog: the <c>--root</c> name, a one-line description (shown
/// as the node summary under <c>--root world</c>), and the factory that turns a parsed
/// <see cref="InspectRootContext"/> into the root's <see cref="IInspectable"/>.</summary>
public sealed record RootCatalogEntry(
    string Name,
    string Description,
    Func<InspectRootContext, IInspectable> Factory,
    int DefaultDepth = 4,
    bool RequiresN = true,
    bool HonorsOptionalN = false);
