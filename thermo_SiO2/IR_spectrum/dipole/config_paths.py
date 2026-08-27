"""Path handling for YAML-configured IR-spectrum commands."""
from pathlib import Path


def resolve_config_path(in_file, value):
    """Resolve a configured path relative to its YAML file.

    Absolute paths are preserved and ``None`` remains ``None``.
    """
    if value is None:
        return None
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = Path(in_file).resolve().parent / path
    return str(path)

