#!/usr/bin/env python3
"""NVD CLI entry point.

This module serves as a thin wrapper that delegates to the py_nvd.cli module.
The actual CLI implementation lives in lib/py_nvd/cli.py.
"""

from py_nvd.cli import main

if __name__ == "__main__":
    main()
