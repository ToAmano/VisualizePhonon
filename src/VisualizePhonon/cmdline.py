import argparse
import os
from VisualizePhonon.vibrational_analysis_io import read_file, save_xsf, save_xyz, generate_asymptote_phonon_code


def main():
    parser = argparse.ArgumentParser(
        description='Visualize phonon modes calculated in VASP.')

    # subcommand. dest for help
    subparsers = parser.add_subparsers(dest='command')

    # generate :: generating files
    gen_parser = subparsers.add_parser(
        'generate', help='read data from VASP OUTCAR and extract phonon modes.')
    gen_parser.add_argument('--input', '-i', type=str,
                            help='input data. VASP OUTCAR is supported.')
    # gen_parser.add_argument('--output', '-o', type=str, default='phonon_visualization.asy', help='出力Asymptoteファイル')
    gen_parser.add_argument('--format', '-f', type=str, default='xsf', choices=[
                            'xsf', 'xyz', 'asy'], help='Output data format. default to xsf')
    gen_parser.add_argument('--mode', '-m', type=int, default=-1,
                            help='The vibration mode. -1 for all modes. default to -1.')
    gen_parser.add_argument('--scale', '-s', type=int, default=1,
                            help='Scale factor of the eigenvector. default to 1.')

    # sample :: generate a sample
    # sample_parser = subparsers.add_parser('sample', help='サンプルデータからコードを生成')
    # sample_parser.add_argument('--output', '-o', default='sample_phonon.asy', help='出力Asymptoteファイル')

    # get command line variables
    args = parser.parse_args()

    # print help
    if not args.command:
        parser.print_help()
        return

    # generate command
    if args.command == 'generate':
        if not os.path.isfile(args.input):
            raise FileNotFoundError(f"{args.input} does not exist!!")
        # read input file
        analysis = read_file(args.input)
        print(f"Total modes: {len(analysis.modes)}")
        print(f"Real modes: {len(analysis.get_real_modes())}")
        print(f"Imaginary modes: {len(analysis.get_imaginary_modes())}")

        if args.mode == -1:  # for all mode
            for index, mode in enumerate(analysis.modes):
                if args.format == "xsf":
                    save_xsf(f"mode_{index}.xsf", mode, args.scale)
                elif args.format == "xyz":  # export to xyz (for animation)
                    save_xyz(f"mode_{index}.xyz", mode, args.scale)
                elif args.format == "asy":
                    generate_asymptote_phonon_code(mode, f"mode_{index}.tex")
        else:  # for a single mode
            if args.format == "xsf":
                save_xsf(f"mode_{args.mode}.xsf",
                         analysis.modes[args.mode], args.scale)
            elif args.format == "xyz":
                save_xyz(f"mode_{args.mode}.xyz",
                         analysis.modes[args.mode], args.scale)
            elif args.format == "asy":
                generate_asymptote_phonon_code(
                    analysis.modes[args.mode], f"mode_{args.mode}.tex")
