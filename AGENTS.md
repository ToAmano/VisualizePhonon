# Repository Guide

## 概要

VisualizePhononは、VASPの`OUTCAR`からphonon/vibrational modeを抽出し、
可視化用ファイルへ変換する小型Pythonパッケージです。主なユーザー向け入口は
`vaspvis`コンソールコマンドです。

```bash
vaspvis generate -i OUTCAR -f xsf -m -1 -s 1.0
```

CLIで指定できる出力形式は`xsf`、`xyz`、`asy`です。

## リポジトリ構造

- `src/VisualizePhonon/`
  - `vibrational_analysis.py`: 中核のデータ構造。`VibrationalMode`はASEの
    atomsオブジェクト、THz単位の周波数、固有ベクトルを保持します。
    `VibrationAnalysis`はmodeのリストを保持し、mode/frequency/eigenvector
    の取得とXSF出力メソッドを提供します。
  - `vibrational_analysis_io.py`: ファイルI/Oとエクスポータ。
    `CreateVibrationAnalysis.load_from_outcar()`と`read_file()`によるVASP
    `OUTCAR`の解析、`save_xsf()`、`save_xyz()`、
    `generate_asymptote_phonon_code()`、Asymptote補助関数があります。
  - `cmdline.py`: `argparse`による`vaspvis generate`のCLI実装。
  - `__init__.py`: `VibrationalMode`と`VibrationAnalysis`をパッケージから
    exportします。
- `test/`
  - `test_vibrationalanalysis.py`: core class、OUTCAR parser、mode filtering、
    accessor、eigenvector、XSF出力のpytest。
  - `OUTCAR`: テストで使う小さなH2Oの参照OUTCAR。
- `examples/`
  - `OUTCAR`と`POSCAR`: example用のVASP入出力。
  - `test.py`: `OUTCAR`を読み、modeを周波数順に並べ、XSF/XYZ/Asymptote
    LaTeX出力を作る例。
  - `mode_*`ファイルの多くは可視化の生成物です。明示的にexample出力の更新を
    求められていない限り、sourceではなく生成物として扱ってください。
- `pyproject.toml`: package metadata、console script、pytest設定、Black/isort設定。
- `setup.py`: 最小限のsetuptools shim。
- `README.md`: install方法、CLI usage、citation、関連ツール。

## 主要データフロー

1. `cmdline.py`の`vaspvis generate`が入力ファイルの存在を確認します。
2. `vibrational_analysis_io.py`の`read_file()`が`VibrationAnalysis`を作ります。
3. `CreateVibrationAnalysis.load_from_outcar()`がASEで構造を読み、`NIONS`を
   parseし、VASPのvibrational analysis sectionを探し、周波数と固有ベクトルを
   抽出して`VibrationalMode`を作ります。
4. Imaginary modeは負の周波数として表現されています。
5. 選択されたmodeをexportします。
   - `save_xsf()`: 静的構造とdisplacement vectorをXSFへ書き出します。
   - `save_xyz()`: sinusoidalに変位させたanimation frameをXYZへ書き出します。
   - `generate_asymptote_phonon_code()`: LaTeX/Asymptote sourceを書き出します。

## 開発コマンド

ローカルインストール:

```bash
pip install -e .
```

開発用extra込みのインストール:

```bash
pip install -e ".[dev]"
```

テスト実行:

```bash
python -m pytest
```

同梱fixtureに対するCLI実行例:

```bash
vaspvis generate -i test/OUTCAR -f xsf -m 0 -s 1
```

## 今後の作業時の注意

- codeは`ase`に依存します。parser/exporterやCIを触る場合は、`pyproject.toml`の
  dependencyとworkflowのinstall手順を同期させてください。
- `pyproject.toml`のversionは`0.0.2`、`src/VisualizePhonon/__init__.py`の
  `__version__`は`0.1.0`で不一致です。
- `cmdline.py`の`--scale`は`float`です。README例のように`1.0`を受け取れます。
- `vibrational_analysis.py`と`vibrational_analysis_io.py`の両方にXSF writer logic
  があります。出力形式を変える場合は挙動をそろえてください。
- CLIやexample実行後に`mode_*.xsf`、`mode_*.xyz`、`mode_*.tex`、`mode_*.pdf`、
  LaTeX auxiliary filesが増えることがあります。生成物を不用意にsource扱いしないで
  ください。
- 現在のtestは`test/OUTCAR`とVASP風のtext layoutに依存しています。parser変更時は
  このfixtureに加え、可能なら別の実OUTCARでも確認してください。
