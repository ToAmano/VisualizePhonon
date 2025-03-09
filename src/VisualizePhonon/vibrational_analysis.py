import ase.io
import numpy as np

class VibrationalMode:
    """振動モードを表すクラス"""
    def __init__(self, atoms:ase.Atoms,frequency:float, eigenvector:np.ndarray):
        self.atoms:ase.Atoms = atoms
        self.frequency:float = frequency  # frequency (THz)
        self.eigenvector:np.ndarray = eigenvector  # eigenvector (shape: (nions, 3))
    
    def __repr__(self):
        freq_type = "real" if self.frequency>=0 else "imaginary"
        return f"VibrationalMode(frequency={self.frequency:.6f} THz, {freq_type})"

class VibrationAnalysis:
    """VASPの振動解析結果を管理するクラス"""
    def __init__(self):
        self.modes = None

    def set_params(self,modes:list[VibrationalMode]):
        self.modes = modes
    
    def get_mode(self, index:int)->VibrationalMode:
        """get the vibrational mode with the specified index"""
        if 0 <= index < len(self.modes):
            return self.modes[index]
        raise IndexError(f"Index {index} out of range for {len(self.modes)} modes")
    
    def get_real_modes(self)->list[VibrationalMode]:
        """get modes with a positive frequency"""
        return [mode for mode in self.modes if mode.frequency >= 0]
    
    def get_imaginary_modes(self)->list[VibrationalMode]:
        """get modes with a negative (=imaginary) frequency"""
        return [mode for mode in self.modes if mode.frequency < 0]
    
    def get_frequencies(self)->np.ndarray:
        """get all the frequencies"""
        return np.array([mode.frequency for mode in self.modes])
    
    def get_eigenvectors(self)->np.ndarray:
        """get all the eigenvectors"""
        return np.array([mode.eigenvector for mode in self.modes])

    def save_xsf(self,filename:str, index:int, scale:float=1.0)->int:
        """
        Write the position and vector field in XSF format.
        """
        vibmode = self.modes[index]

        vector = np.asarray(vibmode.eigenvector, dtype=float) * scale
        atoms  = vibmode.atoms
        assert vector.shape == atoms.positions.shape
        pos_vec = np.hstack((atoms.positions, vector))
        nions = pos_vec.shape[0]
        chem_symbs = atoms.get_chemical_symbols()
        with open(filename, 'w') as out:
            line = "CRYSTAL\n"
            line += "PRIMVEC\n"
            line += '\n'.join([
                ' '.join(['%21.16f' % a for a in vec])
                for vec in atoms.cell
            ])
            line += "\nPRIMCOORD\n"
            line += "{:3d} {:d}\n".format(nions, 1)
            line += '\n'.join([
                '{:3s}'.format(chem_symbs[ii]) +
                ' '.join(['%21.16f' % a for a in pos_vec[ii]])
                for ii in range(nions)
            ])

            out.write(line)
        return 0
