from VisualizePhonon.vibrational_analysis_io import read_file,save_xsf,save_xyz,generate_asymptote_phonon_code

analysis = read_file('OUTCAR')

print(f"Total modes: {len(analysis.modes)}")
print(f"Real modes: {len(analysis.get_real_modes())}")
print(f"Imaginary modes: {len(analysis.get_imaginary_modes())}")

# 周波数の昇順にソート
sorted_modes = sorted(analysis.modes, key=lambda mode: mode.frequency)

# 最低の振動モードを表示
if sorted_modes:
    lowest_mode = sorted_modes[0]
    print(f"Lowest mode: {lowest_mode}")
    # export in xsf format
    save_xsf("mode_0000.xsf",sorted_modes[0])
    # export to xyz (for animation)
    save_xyz("mode_0000.xyz",sorted_modes[0])
    # export to asymptote    
    generate_asymptote_phonon_code(sorted_modes[26], "mode_0026.tex")
