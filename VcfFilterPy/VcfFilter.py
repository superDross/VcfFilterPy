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
    all_d = vcf_dict(line)
    all_d = custom_fields(all_d, conditions)
   
    values = [get_values(conditions, x) for x in all_d]

    
    # test all genotype values in all samples for parsed conditions (returns Boolean list, each element is a sample in the VCF line)
    tested_samples =  [check_condition(conditions, x) for x in values]
    
    if how == 'any' and any(tested_samples):
        return line

    if how == 'all' and all(tested_samples):
        return line
 

## VcfDict Object ##

def vcf_dict(line):
    ''' For every sample in a given vcf line, store the vcf fields
        in a dictionary.

    Returns:
        a list of dicts, one for each sample in the given vcf line.
    '''
    mand_dict = get_mand_dict(line)
    info_dict = get_info_dict(line)
    all_gt_dicts = get_genotype_dicts(line)
    # combine all 3 dicts, PYTHON 3.5 ONLY
    combined = [{**mand_dict, **info_dict, **x} for x in all_gt_dicts]
    return combined



def get_mand_dict(line):
    ''' Return a dict containing all the mandatory
        fields (keys) and their values.
    '''
    mandatory = line.split("\t")[:7]
    header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
    mand_dict = {x:y for x, y in zip(header, mandatory)}
    return mand_dict


def get_info_dict(line):
    ''' Return a dict containing all the info
        fields (keys) and their values.
    '''
    sline = line.split("\t")
    info = sline[7]
    sinfo = [x.split("=") for x in info.split(";")]
    # if an info field has no = sign then assign it an ''
    sinfo = [[x[0],""] if len(x) < 2 else x for x in sinfo ]
    info_dict = {x[0]:x[1] for x in sinfo}
    info_dict['DEPTH'] = info_dict.pop('DP')
    
    return info_dict


def get_genotype_dicts(line):
    ''' Get the genotype values for all samples in a vcf
        line and store as a dict
    '''
    sline = line.split("\t")
    # all sample genotype values
    gt_dict_store = []

    #check all samples genotype values in each line
    for n in range(9, len(sline)):

        sam_fields = sline[n].split(":")
        format_fields = sline[8].split(":")
        gt_dict = {x:y for x, y in zip(format_fields, sam_fields)}
        gt_dict_store.append(gt_dict)
     
    return gt_dict_store


def custom_fields(vcf_dict, conditions):
    ''' Determine whether custom fields need to be created
        in a given list of vcf dicts
    '''
    all_fields = [x.split(" ")[0] for x in conditions]

    if 'AB' in all_fields:
        vcf_dict = [calc_AB(x) for x in vcf_dict]
    
    return vcf_dict


def calc_AB(d):
    ''' Calculate the allele balance for each sample in 
        a list of vcf dicts.
    '''
    AD = d.get("AD")
    DP = d.get("DP")
    if not AD == '':
        AB = int(AD.split(",")[1]) / int(DP)
        d['AB'] = AB
    else:
        d['AB'] = ''
    return d


def get_values(conditions, vcf_dict):
    ''' Filter the info dict for fields present in the 
        given conditions.

    Args:
        conditions: list of filtering conditions
        vcf_dict: dictionary containing info fields and values for a vcf line

    Returns:
        filtered vcf_dict
    ''' 
    all_values = {}

    for cond in conditions:
        field = cond.split(" ")[0]
        alt_field = re.sub(r"\[.*\]", "", field)

        if alt_field in vcf_dict.keys():
            v = vcf_dict.get(alt_field) 
            all_values = assign_value(field, v, all_values)
    
    return all_values


def assign_value(field, v, d):
    ''' Assign a value to a key (field) and store in
        the given dict.

    Args:
        field: VCF field (AD, AC, CHROM etc.)
        v: value to assign to key
        d: dict used to store VCF fields and their values

    Returns:
        updated dict
    '''
    # vcftools fills some fields as '.' which breaks this func
    v = '' if v == '.' else v

    if '[' in str(field):
        index = int(field.split("[")[1].replace(']', ''))
        field = re.sub(r"\[.*\]", "", field)
        d[field] = v.split(",")[index] if "," in v else v
    elif ',' in str(v):
        d[field] = v.split(",")[0]
    else:
        d[field] = v
    
    return d




def check_condition(conditions, sam_dict, combine="&"):
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

f = '/home/david/projects/pdVCF/test/vcfs/testing3.vcf'
filter_vcf(f, ['GT = 1/1', 'DP > 100']) # == 257
filter_vcf(f, ['GT = 0/1', 'DP >= 50', 'GQ >= 30',  'DEPTH > 50']) # == 1896
filter_vcf(f, ['GT = 0/1', 'DP >= 50', 'GQ >= 30', 'AC[0] > 20', 'CHROM = 1', 'POS > 2235893', 'ID != .']) # == 3

f = '/home/david/projects/pdVCF/test/vcfs/testing4.vcf'
filter_vcf(f, ['GT = 0/1', 'DP >= 50', 'AC[0] > 20', 'GQ >= 30']) # == 60
filter_vcf(f, ['GT = 0/1', 'DP >= 50', 'AC[0] > 20', 'GQ >= 30', 'AD[1] >= 20']) # == 49


