from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("amrrules")
except PackageNotFoundError:
    __version__ = "unknown"
