"""Canonical UTF-8 provenance for Route-B Python and JSON artifacts."""
import ast
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def canonical_text_sha256(path):
    """SHA-256 of UTF-8 Python/JSON text with CRLF and lone CR mapped to LF."""
    path = Path(path)
    if path.suffix not in ('.py', '.json'):
        raise ValueError(f'artifact provenance requires .py or .json text: {path}')
    try:
        content = path.read_bytes().decode('utf-8')
    except UnicodeDecodeError as error:
        raise ValueError(f'artifact provenance requires UTF-8 text: {path}') from error
    if '\x00' in content:
        raise ValueError(f'artifact provenance requires text without NUL bytes: {path}')
    canonical = content.replace('\r\n', '\n').replace('\r', '\n')
    return hashlib.sha256(canonical.encode('utf-8')).hexdigest()


def artifact_provenance(script, source):
    """Hash a source and every transitive local import, including package startup.

    Static import traversal includes delayed imports without executing producers.
    All hashes normalize UTF-8 text to LF; binary inputs are outside this contract.
    """
    script = Path(script).resolve()
    pending, visited = [script], set()
    while pending:
        path = pending.pop().resolve()
        if path in visited:
            continue
        visited.add(path)
        for parent in path.parents:
            if parent == ROOT:
                break
            initializer = parent/'__init__.py'
            if initializer.is_file():
                pending.append(initializer)
        for node in ast.walk(ast.parse(path.read_bytes(), filename=str(path))):
            if isinstance(node, ast.Import):
                modules = [ROOT/'simulations'/name.name.replace('.', '/')
                           for name in node.names]
            elif isinstance(node, ast.ImportFrom):
                base = path.parent if node.level else ROOT/'simulations'
                for _ in range(max(0, node.level-1)):
                    base = base.parent
                module = base/(node.module or '').replace('.', '/')
                modules = [module]+[module/name.name for name in node.names]
            else:
                continue
            for module in modules:
                for candidate in (module/'__init__.py', module.with_suffix('.py')):
                    if candidate.is_relative_to(ROOT) and candidate.is_file():
                        pending.append(candidate)
                        break
    dependencies = sorted(path.relative_to(ROOT).as_posix()
                          for path in visited if path != script)
    return dict(script_sha256=canonical_text_sha256(script),
                source_sha256=canonical_text_sha256(source),
                dependencies={path: canonical_text_sha256(ROOT/path)
                              for path in dependencies})


def verify_artifact_provenance(artifact, script, source):
    """Reject absent/stale source, script or dependency hashes before writing."""
    for field, expected in artifact_provenance(script, source).items():
        if artifact.get(field) != expected:
            raise ValueError(f'artifact provenance mismatch: {field}')
