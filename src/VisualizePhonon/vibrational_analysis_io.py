import ase.io
import ase
import numpy as np
from VisualizePhonon import VibrationalMode, VibrationAnalysis


class CreateVibrationAnalysis:
    @classmethod
    def load_from_outcar(cls, outcar_file: str = 'OUTCAR', exclude_imag=False) -> list[VibrationalMode]:
        """load vib modes from OUTCAR"""
        atoms = ase.io.read(outcar_file)
        with open(outcar_file, 'r') as f:
            out = [line for line in f if line.strip()]

        # get the number of atoms in the system
        for line in out:
            if "NIONS =" in line:
                n_ions: int = int(line.split()[-1])
                break

        if n_ions == 0:
            raise ValueError("NIONS not found in OUTCAR")

        # search for vibrational analysis part
        ln = len(out)
        THz_index = []
        i_index = 0

        for ii in range(ln-1, 0, -1):
            if '2PiTHz' in out[ii]:
                THz_index.append(ii)
            if 'Eigenvectors and eigenvalues of the dynamical matrix' in out[ii]:
                i_index = ii + 2
                break

        if i_index == 0:
            raise ValueError(
                "Vibrational analysis section not found in OUTCAR")

        j_index = THz_index[0] + n_ions + 2

        # frequencies and eigenvectors
        freq_lines = [line for line in out[i_index:j_index]
                      if '2PiTHz' in line]
        mode_lines = [line for line in out[i_index:j_index]
                      if ('dx' not in line) and ('2PiTHz' not in line)]

        # construct eigenvectors from mode_lines
        modes = []
        eigenvectors = np.array([line.split()[3:6] for line in mode_lines],
                                dtype=float).reshape((-1, n_ions, 3))

        for i, line in enumerate(freq_lines):
            is_real = 'f/i' not in line
            frequency = float(line.split()[-4])
            # use negative frequency for imaginary modes
            if not is_real:
                frequency = -frequency

            mode = VibrationalMode(
                atoms=atoms,
                frequency=frequency,
                eigenvector=eigenvectors[i]
            )
            modes.append(mode)

        # remove imaginary mode (option)
        if exclude_imag:
            modes = [mode for mode in modes if mode.frequency >= 0]

        return modes


def read_file(filename: str = "OUTCAR", exclude_imag: bool = False) -> VibrationAnalysis:
    """read OUTCAR and return VibrationAnalysis"""

    vibrationalanalysis_instance = VibrationAnalysis()
    modes = CreateVibrationAnalysis.load_from_outcar(filename, exclude_imag)
    vibrationalanalysis_instance.set_params(modes)
    return vibrationalanalysis_instance


def save_xsf(filename: str, vibmode: VibrationalMode, scale: float = 1.0) -> str:
    """
    Write the position and eigenvector in XSF format.
    """
    vector = vibmode.eigenvector * scale
    atoms: ase.Atoms = vibmode.atoms
    chem_symbs: list[str] = atoms.get_chemical_symbols()
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

    return output_str

# Function to generate frames for vibration mode animation


def save_xyz(filename: str, vibmode: VibrationalMode, scale: float = 1.0, num_frames=20) -> list[ase.Atoms]:
    atoms_list = []

    positions = vibmode.atoms.get_positions()
    atom_types = vibmode.atoms.get_atomic_numbers()
    displacements = vibmode.eigenvector

    # Create a sinusoidal variation along the displacements
    t_vals = np.linspace(0, 2 * np.pi, num_frames)
    for t in t_vals:
        displaced_positions = positions + scale * np.sin(t) * displacements
        atoms = ase.Atoms(atom_types, positions=displaced_positions)
        atoms_list.append(atoms)

    ase.io.write(filename, atoms_list)
    return atoms_list


def generate_asymptote_phonon_code(vibmode: VibrationalMode,
                                   output_file="phonon_visualization.asy",
                                   sphere_radius=0.3,
                                   cyl_radius=0.1,
                                   figure_size=(4, 3),
                                   camera_position=(3, 1, 0.5)):
    """
    フォノンの振動モードをAsymptote+LaTeXで可視化するためのコードを生成します。

    Parameters:
    ----------
    atoms_coords : list of tuples
        原子の座標のリスト [(x1, y1, z1), (x2, y2, z2), ...]
    atoms_types : list of str
        原子の種類のリスト ['C', 'O', 'H', ...]
    vibration_vectors : list of tuples
        振動ベクトルのリスト [(dx1, dy1, dz1), (dx2, dy2, dz2), ...]
    output_file : str, optional
        出力ファイル名, デフォルトは "phonon_visualization.asy"
    sphere_radius : float, optional
        原子球の半径, デフォルトは 0.3
    cyl_radius : float, optional
        結合円柱の半径, デフォルトは 0.1
    figure_size : tuple, optional
        図のサイズ (幅, 高さ), デフォルトは (4, 3)
    camera_position : tuple, optional
        カメラ位置, デフォルトは (3, 1, 0.5)

    Returns:
    -------
    None. ファイルに書き込みます。
    """
    atoms_coords = vibmode.atoms.get_positions()
    atoms_types = vibmode.atoms.get_chemical_symbols()
    vibration_vectors = vibmode.eigenvector

    if len(atoms_coords) != len(atoms_types) or len(atoms_coords) != len(vibration_vectors):
        raise ValueError("原子の座標、種類、振動ベクトルの数は一致する必要があります。")

    # 元素ごとの色を定義
    element_colors = {
        'H': 'white',
        'C': 'gray',
        'N': 'blue',
        'O': 'red',
        'F': 'green',
        'P': 'orange',
        'S': 'yellow',
        'Cl': 'green',
        'Br': 'brown',
        'I': 'purple',
        'Ti': 'gray',
        # 他の元素も必要に応じて追加可能
    }

    # 元素ごとのデフォルト色
    default_color = 'gray'

    # コードの先頭部分を作成
    header = f"""\\documentclass[lualatex]{{standalone}}
\\usepackage{{asymptote}}

\\begin{{document}}
\\begin{{asy}}
import three;
settings.render=8;
settings.prc=false;
settings.outformat="pdf";
size({figure_size[0]}cm,{figure_size[1]}cm);

// the camera position
currentprojection = orthographic({camera_position});

// materials for different atoms
"""

    # マテリアル定義を追加
    materials = ""
    used_elements = set(atoms_types)
    for element in used_elements:
        color = element_colors.get(element, default_color)
        materials += f"material {element}_color = material(diffusepen={color}, specularpen=white);\n"

    materials += "material cylcolor = material(diffusepen=gray, emissivepen=gray);\n"

    # 関数定義
    functions = f"""
// cylinder radius
real cylRadius = {cyl_radius};
// point radius
real sphereRadius = {sphere_radius};

// draw rod(line)
void drawRod(triple a, triple b) {{
  surface rod = extrude(scale(cylRadius)*unitcircle, axis=length(b-a)*Z);
  triple orthovector = cross(Z, b-a);
  if (length(orthovector) > .01) {{
    real angle = aCos(dot(Z, b-a) / length(b-a));
    rod = rotate(angle, orthovector) * rod;
  }}
  draw(shift(a)*rod, surfacepen=cylcolor);
}}

// draw atoms
void drawAtom(triple center, material atomColor) {{
     draw(shift(center)*scale3(sphereRadius)*unitsphere, surfacepen=atomColor);
}}

// draw arrow
void drawArray(triple center, triple direction) {{
     triple end = center+direction;
     draw(center--end, green, Arrow3);
}}
"""

    # 原子定義部分を作成
    atoms_def = "\n// Atom positions\n"
    for i, (coord, atom_type) in enumerate(zip(atoms_coords, atoms_types)):
        atoms_def += f"triple {atom_type}_{i} = ({coord[0]},{coord[1]},{coord[2]});\n"

    # 振動ベクトル定義部分を作成
    vibrations_def = "\n// Vibration vectors\n"
    for i, (vib_vec, atom_type) in enumerate(zip(vibration_vectors, atoms_types)):
        vibrations_def += f"triple {atom_type}_arrow_{i} = ({vib_vec[0]},{vib_vec[1]},{vib_vec[2]});\n"

    # 描画命令部分を作成
    draw_commands = "\n// Draw atoms and vibrations\n"
    for i, atom_type in enumerate(atoms_types):
        draw_commands += f"drawAtom({atom_type}_{i}, {atom_type}_color);\n"
        draw_commands += f"drawArray({atom_type}_{i}, {atom_type}_arrow_{i});\n"

    # 結合描画のための部分（オプション）
    bonds = "\n// Bond connections - customize this part based on your molecule\n"
    bonds += "// Example: drawRod(C_0, O_1);\n"

    # ファイル終了部分
    footer = """
\\end{asy}
\\end{document}
"""

    # 全コードを結合
    full_code = header + materials + functions + atoms_def + \
        vibrations_def + draw_commands + bonds + footer

    # ファイルに書き込み
    with open(output_file, 'w') as f:
        f.write(full_code)

    print(f"{output_file}にコードを生成しました。")


def extract_data_from_asymptote(asymptote_code):
    """
    既存のAsymptoteコードから原子座標、種類、振動ベクトルのデータを抽出します。

    Parameters:
    ----------
    asymptote_code : str
        Asymptoteコード

    Returns:
    -------
    tuple:
        (atoms_coords, atoms_types, vibration_vectors)
    """
    lines = asymptote_code.split('\n')

    atoms_coords = []
    atoms_types = []
    atoms_indices = []
    vibration_vectors = []

    # 座標とタイプを抽出
    for line in lines:
        # 原子座標の行を検出
        if line.strip().startswith('triple') and '_arrow_' not in line:
            parts = line.split('=')
            if len(parts) < 2:
                continue

            atom_info = parts[0].strip().split()[-1]
            coord_str = parts[1].strip().rstrip(';')

            # 座標を抽出
            try:
                # (x,y,z)の形式から数値のタプルに変換
                coord_str = coord_str.strip('()')
                coords = tuple(float(x) for x in coord_str.split(','))

                # 原子タイプとインデックスを抽出
                if '_' in atom_info:
                    atom_type, idx = atom_info.split('_')
                    atoms_types.append(atom_type)
                    atoms_indices.append(int(idx))
                    atoms_coords.append(coords)
            except:
                continue

    # 振動ベクトルを抽出
    for line in lines:
        if line.strip().startswith('triple') and '_arrow_' in line:
            parts = line.split('=')
            if len(parts) < 2:
                continue

            arrow_info = parts[0].strip().split()[-1]
            vector_str = parts[1].strip().rstrip(';')

            # ベクトルを抽出
            try:
                # (x,y,z)の形式から数値のタプルに変換
                vector_str = vector_str.strip('()')
                vector = tuple(float(x) for x in vector_str.split(','))

                # インデックスを取得して正しい順序で追加
                if '_arrow_' in arrow_info:
                    idx = int(arrow_info.split('_')[-1])
                    if idx < len(vibration_vectors):
                        vibration_vectors[idx] = vector
                    else:
                        while len(vibration_vectors) < idx:
                            vibration_vectors.append((0, 0, 0))
                        vibration_vectors.append(vector)
            except:
                continue

    # インデックスでソートして元の順序を復元
    sorted_data = sorted(zip(atoms_indices, atoms_types, atoms_coords))
    atoms_types = [t for _, t, _ in sorted_data]
    atoms_coords = [c for _, _, c in sorted_data]

    return atoms_coords, atoms_types, vibration_vectors


def add_bonds_to_molecule(asymptote_file, bonds_list):
    """
    分子の結合情報を追加します。

    Parameters:
    ----------
    asymptote_file : str
        Asymptoteコードが含まれるファイル
    bonds_list : list of tuples
        結合のリスト。各結合は (atom1_type, atom1_idx, atom2_type, atom2_idx) の形式
        例: [('C', 1, 'O', 0), ('C', 1, 'C', 2)]

    Returns:
    -------
    None. ファイルを更新します。
    """
    with open(asymptote_file, 'r') as f:
        code = f.read()

    # 結合コードを生成
    bonds_code = "\n// Bond connections\n"
    for bond in bonds_list:
        atom1_type, atom1_idx, atom2_type, atom2_idx = bond
        bonds_code += f"drawRod({atom1_type}_{atom1_idx}, {atom2_type}_{atom2_idx});\n"

    # コード内の結合部分を探して置き換える
    if "// Bond connections" in code:
        # 既存の結合コードを置き換え
        parts = code.split("// Bond connections")
        if len(parts) > 1:
            # 次のメジャーセクションまでをスキップ
            next_part = parts[1].split("\\end{asy}", 1)
            if len(next_part) > 1:
                new_code = parts[0] + bonds_code + \
                    "\n\\end{asy}" + next_part[1]
                with open(asymptote_file, 'w') as f:
                    f.write(new_code)
                print(f"{asymptote_file}に結合情報を追加しました。")
                return

    # 結合セクションが見つからない場合、ファイル末尾に追加
    end_marker = "\\end{asy}"
    if end_marker in code:
        new_code = code.replace(end_marker, bonds_code + "\n" + end_marker)
        with open(asymptote_file, 'w') as f:
            f.write(new_code)
        print(f"{asymptote_file}に結合情報を追加しました。")
    else:
        print("ファイル形式が正しくないため、結合情報を追加できませんでした。")
