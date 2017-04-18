from Vcf import Vcf
import argparse
import subprocess
import os, sys


def get_parser():
    parser = argparse.ArgumentParser(description='vcf filtering tool')
    parser.add_argument('--vcf', help='vcf file')
    parser.add_argument('-f', '--filter', help='filtering criteria')
    parser.add_argument('-p', '--operator', default='any', help='detail whether all or any samples need to pass the filtering criteria to pass')
    parser.add_argument('-n', '--no-print', help='do not print the filtered vcf')
    parser.add_argument('-o', '--out', help='name of filtered vcf file')
    parser.add_argument('-v', '--version', action='store_true', help='display the current version')
    parser.add_argument('--test', action='store_true', help='test installation')
    return parser

def cli():
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['test']:
        here = os.path.realpath(__file__)
        test = "/".join(here.split("/")[:-2]) + "/test/test_filters.py"
        subprocess.call(['python3', test])
        return

    if args['version']:
        print(__version__)
        return

    if args['operator'] not in ['any', 'all']:
        raise ValueError("--operator only accepts 'any' or 'all' as an argument") 
    
    # create and filter the Vcf object
    vcf = Vcf(args['vcf'], sys.argv)
    conditions = [x.strip(" ") for x in args['filter'].split(",")]
    vcf.filter_vcf(conditions, args['operator'])
    print("{} variants remain".format(vcf.count))

    if args['out']:
        vcf.write_filtered_vcf(args['out'])


if __name__ == '__main__':
    cli()
