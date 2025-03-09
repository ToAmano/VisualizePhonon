import tempfile
import os
from typing import List, Tuple
import ase.atoms
import ase.io
import numpy as np


class VibrationalMode:
    """class to express a phonon mode"""

    def __init__(self, atoms: ase.Atoms, frequency: float, eigenvector: np.ndarray):
        assert eigenvector.shape == atoms.positions.shape
        self.atoms: ase.Atoms = atoms
        self.frequency: float = frequency  # frequency (THz)
        self.eigenvector: np.ndarray = eigenvector

    def __repr__(self):
        freq_type = "real" if self.frequency >= 0 else "imaginary"
        return f"VibrationalMode(frequency={self.frequency:.6f} THz, {freq_type})"


class VibrationAnalysis:
    """class to contain all the modes from VASP OUTCAR"""

    def __init__(self):
        self.modes = None

    def set_params(self, modes: list[VibrationalMode]):
        self.modes = modes

    def get_mode(self, index: int) -> VibrationalMode:
        """get the vibrational mode with the specified index"""
        if 0 <= index < len(self.modes):
            return self.modes[index]
        raise IndexError(
            f"Index {index} out of range for {len(self.modes)} modes")

    def get_real_modes(self) -> list[VibrationalMode]:
        """get modes with a positive frequency"""
        return [mode for mode in self.modes if mode.frequency >= 0]

    def get_imaginary_modes(self) -> list[VibrationalMode]:
        """get modes with a negative (=imaginary) frequency"""
        return [mode for mode in self.modes if mode.frequency < 0]

    def get_frequencies(self) -> np.ndarray:
        """get all the frequencies"""
        return np.array([mode.frequency for mode in self.modes])

    def get_eigenvectors(self) -> np.ndarray:
        """get all the eigenvectors"""
        return np.array([mode.eigenvector for mode in self.modes])

    def save_xsf(self, filename: str, index: int, scale: float = 1.0) -> str:
        """
        Write the position and eigenvector in XSF format.
        """
        vibmode = self.modes[index]
        vector = vibmode.eigenvector * scale
        atoms: ase.Atoms = vibmode.atoms
        chem_symbs = atoms.get_chemical_symbols()
        nions: int = len(chem_symbs)  # number of ions

        pos_vec = np.hstack((atoms.positions, vector))

        # XSF形式で文字列を構築
        lines = ["CRYSTAL", "PRIMVEC"]

        # 格子ベクトルを追加
        for vec in atoms.cell:
            lines.append(' '.join(['%21.16f' % a for a in vec]))

        # 原子座標セクション
        lines.extend(["PRIMCOORD", f"{nions:3d} {1:d}"])

        # 各原子の情報を追加
        for i in range(nions):
            atom_line = f"{chem_symbs[i]:3s}" + \
                ' '.join(['%21.16f' % a for a in pos_vec[i]])
            lines.append(atom_line)

        # 最終的な文字列を生成
        output_str = '\n'.join(lines)

        # ファイルに書き込み
        with open(filename, 'w') as output:
            output.write(output_str)

        #     output.write(line)
        return output_str
