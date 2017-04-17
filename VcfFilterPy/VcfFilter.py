from vcf_dict import vcfline2dict
from custom_fields import custom_fields
from get_vcf_values import get_values
import operator
import re

ops = { 
    
    ">": operator.gt,
    "<": operator.lt,
    ">=": operator.ge,
    "<=": operator.le,
    "==": operator.eq,
    "=": operator.eq,
    "!=": operator.ne,
    "&": operator.and_,
    "|": operator.or_ 
}



def filter_vcf(vcf, conditions, how='any'):
    ''' Filter a VCF file based on a set of conditions.

    Args:
        vcf: a vcf file
        conditions: a list of filtering conditions
        how: specify whether all or any samples in the VCF need
             to meet the conditions for the VCF to pass filter.

    Returns:
        filtered vcf
    '''
    vcf = open(vcf)
    count = 0

    while True:

        line = vcf.readline()
        line = line.strip("\n")
        
        # when the file ends
        if not line:
            print(count)
            break
        
        # header
        if line.startswith("#"):
            handle_headers(line)
        
        # vcf line
        if not line.startswith("#"):
            fline = filter_line(line, conditions, how)
            count += 1 if fline else 0


def handle_headers(line):
    pass


def filter_line(line, conditions, how):
    ''' Filter a line of a vcf file with a set of conditions
    '''
    # store all vcf fields in a dict
    all_d = vcfline2dict(line)
    all_d = custom_fields(all_d, conditions)
   
    values = [get_values(conditions, x) for x in all_d]

    
    # test all genotype values in all samples for parsed conditions (returns Boolean list, each element is a sample in the VCF line)
    tested_samples =  [test_condition(conditions, x) for x in values]
    
    if how == 'any' and any(tested_samples):
        return line

    if how == 'all' and all(tested_samples):
        return line
 


def test_condition(conditions, sam_dict, combine="&"):
    ''' Test the filtering conditions with a given samples vcf field values.

    Args:
        conditions: list of filtering conditions
        sam_dict: dict containing fields for a sample in a vcf line 
                  (keys e.g. GT, CHROME) and their corresponding values 
        combine: all (&) or any (|) conditions must be met

    Returns:
        Boolean; True if all/any values in the sam_dict meet the given conditions
    '''
    if not sam_dict:
        return False
    
    tested_sample = []
    
    for cond in conditions:

        field, op_sign, val = cond.split(" ")
        op = ops[op_sign]
        if val.replace(".", "").isdigit():
            val =float(val)
        
        # get the samples value for the given field in the condition
        test = sam_dict.get(re.sub(r"\[.*\]", "", field))
        test = float(test) if str(test).replace(".", "").isdigit() else test
        
        # test the condition against the sample value
        if test:
            test_cond = op(test, val) 
            tested_sample.append(test_cond)
        else:
            tested_sample.append(False)

    if combine == "&":
        return all(tested_sample)

    elif combine == "|":
        return any(tested_sample)


# TODO: specify fields with commas and MANDATORY fields

f = '/home/david/projects/pdVCF/test/vcfs/testing.vcf'
filter_vcf(f, ['GT = 1/1', 'DP > 100', 'AB > 0']) # == 3
filter_vcf(f, ['GT = 1/1', 'DP > 100', 'AB > 0', 'AC < 50', 'DEPTH = 20']) # == 2

#f = '/home/david/projects/pdVCF/test/vcfs/testing3.vcf'
#filter_vcf(f, ['GT = 1/1', 'DP > 100']) # == 257
#filter_vcf(f, ['GT = 0/1', 'DP >= 50', 'GQ >= 30',  'DEPTH > 50']) # == 1896
#filter_vcf(f, ['GT = 0/1', 'DP >= 50', 'GQ >= 30', 'AC[0] > 20', 'CHROM = 1', 'POS > 2235893', 'ID != .']) # == 3
#
#f = '/home/david/projects/pdVCF/test/vcfs/testing4.vcf'
#filter_vcf(f, ['GT = 0/1', 'DP >= 50', 'AC[0] > 20', 'GQ >= 30']) # == 60
#filter_vcf(f, ['GT = 0/1', 'DP >= 50', 'AC[0] > 20', 'GQ >= 30', 'AD[1] >= 20']) # == 49


