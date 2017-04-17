''' Allows one to define custom fields to filter.
'''

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


