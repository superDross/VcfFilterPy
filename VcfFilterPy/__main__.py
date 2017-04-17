from Vcf import Vcf
import argparse

def get_parser():
    parser = argparse.ArgumentParser(description='vcf filtering tool')
    parser.add_argument('--vcf', help='vcf file')
    parser.add_argument('-f', '--filter', help='filtering criteria')
    parser.add_argument('-p', '--operator', default='any', help='detail whether all or any samples need to pass the filtering criteria to pass')
    parser.add_argument('-n', '--no-print', help='do not print the filtered vcf')
    parser.add_argument('-o', '--out', help='name of filtered vcf file')
    parser.add_argument('-v', '--version', help='display the current version')
    return parser

def cli():
    parser = get_parser()
    args = vars(parser.parse_args())

    if args['version']:
        print(__version__)
        return

    if args['operator'] not in ['any', 'all']:
        raise ValueError("--operator only accepts 'any' or 'all' as an argument") 
    
    # create and filter the Vcf object
    vcf = Vcf(args['vcf'])
    conditions = [x.strip(" ") for x in args['filter'].split(",")]
    vcf.filter_vcf(conditions, args['operator'])
    print("{} variants remain".format(vcf.count))

    if args['out']:
        vcf.write_filtered_vcf(args['out'])


if __name__ == '__main__':
    cli()
