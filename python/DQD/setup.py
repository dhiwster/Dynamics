"""Compatibility wrapper for older imports.

Use :mod:`DQD.system` for new code. This module remains so notebooks and
scripts that import ``DQD.setup`` keep working.
"""

from .system import DQDsystem

__all__ = ["DQDsystem"]
