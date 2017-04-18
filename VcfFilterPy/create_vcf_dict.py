''' Allows one to create a list of dictionarys for every sample in a vcf line.
'''
def vcfline2dict(line):
    ''' For every sample in a given vcf line, store the vcf fields
        in a dictionary.
    
    Args:
        line: a vcf line

    Returns:
        a list of dicts, one for each sample in the given vcf line.
    '''
    mand_dict = get_mand_dict(line)
    info_dict = get_info_dict(line)
    all_gt_dicts = get_genotype_dicts(line)
    # combine all 3 dicts, PYTHON 3.5 ONLY
    #combined = [{**mand_dict, **info_dict, **x} for x in all_gt_dicts]
    combined = [merge_dicts(mand_dict, info_dict, x) for x in all_gt_dicts]
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


def merge_dicts(*dict_args):
    '''Given any number of dicts, shallow copy and merge into a new dict,
       precedence goes to key value pairs in latter dicts.

    Note:
        in python 3.5; merged = {**dict1, **dict2, **dict3}
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


