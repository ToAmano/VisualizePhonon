"""
VASP Vibration Analyzer
=======================

A Python package for analyzing vibrational modes from VASP OUTCAR files.

Main features:
- Extract vibrational modes and frequencies from VASP OUTCAR files
- Handle both real and imaginary frequencies as complex numbers
- Export vibration modes to VESTA format
"""

__version__ = "0.1.0"

from VisualizePhonon.vibrational_analysis import VibrationalMode, VibrationAnalysis
__all__ = ["VibrationalMode", "VibrationAnalysis"]
