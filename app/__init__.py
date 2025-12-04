__all__ = ["__version__"]

# Resolve version from the installed package metadata, with safe fallbacks.
try:
    from importlib.metadata import version, PackageNotFoundError  # Python 3.8+
except Exception:  # pragma: no cover
    from importlib_metadata import version, PackageNotFoundError  # backport

try:
    __version__ = version("hulk")
except PackageNotFoundError:  # editable installs during development, or no dist metadata
    try:
        # Optional: if you generate/set a version at build time
        from ._version import __version__  # type: ignore
    except Exception:
        __version__ = "0.0.0"
