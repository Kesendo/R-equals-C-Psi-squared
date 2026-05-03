namespace RCPsiSquared.Cli;

/// <summary>Tiny option parser for <c>--name value</c> and <c>--name=value</c> flags.
///
/// Tracks which `--key` actually had a value so a missing value is reported clearly
/// rather than as a downstream parse failure.
/// </summary>
public sealed class ArgParser
{
    private readonly Dictionary<string, string> _options = new(StringComparer.Ordinal);
    private readonly List<string> _positional = new();

    public ArgParser(IReadOnlyList<string> args)
    {
        for (int i = 0; i < args.Count; i++)
        {
            string a = args[i];
            if (a.StartsWith("--"))
            {
                string body = a[2..];
                int eq = body.IndexOf('=');
                if (eq >= 0)
                {
                    _options[body[..eq]] = body[(eq + 1)..];
                    continue;
                }
                if (i + 1 < args.Count && !args[i + 1].StartsWith("--"))
                {
                    _options[body] = args[i + 1];
                    i++;
                }
                else
                {
                    _options[body] = ""; // flag with no value; downstream RequireX throws
                }
            }
            else
            {
                _positional.Add(a);
            }
        }
    }

    public IReadOnlyList<string> Positional => _positional;

    public void RequireNoPositional()
    {
        if (_positional.Count > 0)
            throw new ArgumentException($"unexpected positional arguments: {string.Join(", ", _positional)}");
    }

    public int RequireInt(string key)
    {
        string v = RequireValue(key);
        return int.Parse(v, System.Globalization.CultureInfo.InvariantCulture);
    }

    public double RequireDouble(string key)
    {
        string v = RequireValue(key);
        return double.Parse(v, System.Globalization.CultureInfo.InvariantCulture);
    }

    public string? OptionalString(string key) =>
        _options.TryGetValue(key, out var v) && v.Length > 0 ? v : null;

    /// <summary>True if <c>--key</c> was passed (with or without a value). Use for boolean
    /// flags where presence alone is the signal.</summary>
    public bool HasFlag(string key) => _options.ContainsKey(key);

    public double? OptionalDouble(string key) =>
        OptionalString(key) is { } v
            ? double.Parse(v, System.Globalization.CultureInfo.InvariantCulture)
            : null;

    private string RequireValue(string key)
    {
        if (!_options.TryGetValue(key, out var v))
            throw new ArgumentException($"missing required option --{key}");
        if (v.Length == 0)
            throw new ArgumentException($"option --{key} requires a value");
        return v;
    }
}
