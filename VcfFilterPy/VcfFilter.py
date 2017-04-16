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
        how: specify whether all o any of the conditions
             need to be met

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
        
        if line.startswith("#"):
            handle_headers(line)
        
        if not line.startswith("#"):
             mand_dict = get_mand_dict(line)
             mand_values = get_info_values(conditions, mand_dict)

             info_dict = get_info_dict(line)
             info_values = get_info_values(conditions, info_dict)
             
             gt_conditions = [cond for cond in conditions 
                              if not re.sub(r"\[.*\]", "", cond.split(" ")[0]) in info_dict.keys() and not cond.split(" ")[0] in mand_dict.keys()]
             
             # store the genotype values to be tested for each sample in the vcf within a nested list
             gt_values = get_genotype_value(line, gt_conditions)

             # combine gt and info values ONLY WORKS IN PYTHON 3.5
             values = [{**info_values, **mand_values, **x} for x in gt_values]
             
             # test all genotype values in all samples for parsed conditions
             tested_samples =  [check_condition(conditions, x) for x in values]
                
             if  how == 'any' and any(tested_samples):
                count += 1
                # print(line)


             if how == 'all' and all(tested_samples):
                count += 1
                # print(line)

 
def handle_headers(line):
    pass

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


def get_info_values(conditions, info_dict):
    ''' Filter the info dict for fields present in the 
        given conditions.

    Args:
        conditions: list of filtering conditions
        info_dict: dictionary containing info fields and values for a vcf line

    Returns:
        filtered info_dict
    ''' 
    all_values = {}

    for cond in conditions:
        field = cond.split(" ")[0]
        alt_field = re.sub(r"\[.*\]", "", field)

        if alt_field in info_dict.keys():
            v = info_dict.get(alt_field) 
            all_values = assign_value(field, v, all_values)
             
    #print(all_values)
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
    # for string indexing fields with ',' in cells
    if '[' in field:
        index = int(field.split("[")[1].replace(']', ''))
        field = field.split("[")[0]
        d[field] = v.split(",")[index]
    elif ',' in v:
        d[field] = v.split(",")[0]
    else:
        d[field] = v
    
    return d


def get_genotype_value(line, conditions):
    ''' Get the genotype values for the genotype fields 
        specified in a set of conditions from the parsed
        vcf line and do so for each sample in the vcf.

    Args:
        line: vcf line
        conditions: list of filtering conditions

    Returns:
        a nested list of genotype values where each list
        represents a sample in the vcf 
    '''

    sline = line.split("\t")

    all_value_store = []

    #check all samples genotype values in each line
    for n in range(9, len(sline)):

        all_values = {}
        all_value_store.append(all_values)

        # get all values in the given genotype fields
        for cond in conditions:
            field = cond.split(" ")[0]
            alt_field = re.sub("\[.*\]", "", field)

            all_fields = sline[n].split(":")
            format_fields = sline[8].split(":")

            if all_fields[0] != './.' and len(all_fields) > 1:
                
                if field == "AB" and all_fields[1] != "":
                    AD = all_fields[format_fields.index('AD')]
                    DP = all_fields[format_fields.index('DP')]
                    AB = int(AD.split(",")[1]) / int(DP)
                    #all_values.append(AB)
                    all_values['AB'] = AB

                if field != 'AB' and all_fields[format_fields.index(alt_field)] not in ['./.','']:

                    index = format_fields.index(alt_field)
                    value = all_fields[index]
                    
                    all_values = assign_value(field, value, all_values)

    return all_value_store


def check_condition(conditions, values, combine="&"):
    ''' Test the filtering conditions with the given values.

    Args:
        conditions: list of filtering conditions
        values: dict containing keys and values in which to test the conditions against
        combine: all (&) or any (|) conditions must be met

    Returns:
        Boolean; True if all values meet the given conditions

    Note:
        The conditions and values order must be the same e.g. 
        if the first element of conditions is for DP then the
        first element of values should refer to a DP value.
    '''
    if not values:
        return False
    
    tested_sample = []
    
    for cond in conditions:

        field, op_sign, val = cond.split(" ")
        op = ops[op_sign]
        if val.replace(".", "").isdigit():
            val =float(val)
        
        # test the condition against the value
        test = values.get(re.sub(r"\[.*\]", "", field))
        test = float(test) if str(test).replace(".", "").isdigit() else test
        
        if test:
            #print(field, test, val)
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


