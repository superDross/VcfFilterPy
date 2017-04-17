''' Allows one to filter a vcf dict line for fields present in a list of filtering conditions.
'''
import re

def filter_dict_keys(conditions, vcf_dict):
    ''' Filter the vcf dict for fields present in the 
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


