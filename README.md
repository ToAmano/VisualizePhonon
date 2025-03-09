# VisualizePhonon

Phonon visualization tool from VASP OUTCAR

## Installation

We support both local and pypi installation.

For local installation:

"""
$git clone https://github.com/ToAmano/VisualizePhonon.git
$cd VisualizePhonon
$pip install .
"""

For pip installation:

"""bash
$pip install VisualizePhonon
"""

## Usage

After installation, we have "vaspvis.py" command to extract data from OUTCAR.

"""
%vaspvis generate --help
usage: vaspvis generate [-h] [--input INPUT] [--format {xsf,xyz,asy}]
                        [--mode MODE] [--scale SCALE]

options:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        input data. VASP OUTCAR is supported.
  --format {xsf,xyz,asy}, -f {xsf,xyz,asy}
                        Output data format. default to xsf
  --mode MODE, -m MODE  The vibration mode. -1 for all modes. default to -1.
  --scale SCALE, -s SCALE
                        Scale factor of the eigenvector. default to 1.
"""

The following example extracts all the modes to xsf format:

"""bash
vaspvis generate -i OUTCAR -f xsf -m -1 -s 1.0
"""

You can find an example in the "examples/" directory.

## Citation

If you find the tool useful and use it in the publication, please consider citing the following paper.

- [Transferability of the chemical bond-based machine learning model for dipole moment: the GHz to THz dielectric properties of liquid propylene glycol and polypropylene glycol](https://arxiv.org/abs/2410.22718)


## Other similar tools

There are several similar tools to extract phonon modes from VASP OUTCAR.

- [VaspVib2XSF](https://github.com/QijingZheng/VaspVib2XSF)
- [VASP-plot-modes](https://github.com/Stanford-MCTG/VASP-plot-modes)
- [Phonopy_VESTA](https://github.com/AdityaRoy-1996/Phonopy_VESTA)
- [p4vasp](https://github.com/orest-d/p4vasp)
