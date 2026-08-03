"""``python3 -m mzml_utils.structure`` entry point."""
from .cli import _cli

if __name__ == "__main__":
    raise SystemExit(_cli())
