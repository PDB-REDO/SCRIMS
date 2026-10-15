"""
SCRIMS — statistical inference of functional interaction sites from structure-prediction screens.

Builds a confidence matrix (partners x bait residues) from AlphaFold2,
AlphaFold3 or AlphaBridge output, and calls the interaction sites that are
statistically supported across a screen.
"""

from scrims.screen import analyze, load_matrix, DEFAULTS, Result
from scrims.stats import call_patches
from scrims.matrix import ptm_func, bait_confidence_row, confidence_signal

__version__ = "0.1.0.dev0"
__all__ = ["analyze", "load_matrix", "call_patches", "DEFAULTS", "Result",
           "ptm_func", "bait_confidence_row", "confidence_signal"]
