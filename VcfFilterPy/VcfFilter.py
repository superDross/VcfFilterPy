import operator

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

             # get_info_values(line, conditions)

             # store the genotype values to be tested for each sample in the vcf within a nested list
             values = get_genotype_value(line, conditions)
             
             # test all genotype values in all samples for parsed conditions
             tested_samples =  [check_condition(conditions, x) for x in values]
                
             if  how == 'any' and any(tested_samples):
                count += 1
                # print(line)

             if how == 'all' and all(tested_samples):
                count += 1
                # print(line)

    
def handle_headers(line):
    if line.startswith("##INFO"):
        pass


#def get_info_values(line, conditions):
#
#    sline = line.split("\t")
#    info = sline[7]
#    sinfo = [x.split("=") for x in info.split(";")]
#    info_dict = {x[0]:x[1] for x in sinfo}
#    info_dict['DEPTH'] = info_dict.pop('DP')
#   
#    for cond in conditions:
#       field = cond.split(" ")[0]
#
#       if field in info_dict.keys():
#          v = info_dict.get(field) 
#          
#    import sys
#    sys.exit()

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

    all_field_store = []

    #check all samples genotype values in each line
    for n in range(9, len(sline)):

        field_store = []
        all_field_store.append(field_store)

        # get all values in the given genotype fields
        for cond in conditions:

            field = cond.split(" ")[0]
            sam = sline[n]
            all_fields = sam.split(":")
            format_fields = sline[8].split(":")

            if all_fields[0] != './.' and len(all_fields) > 1:
                
                if field == "AB" and all_fields[1] != "":
                    AD = all_fields[format_fields.index('AD')]
                    DP = all_fields[format_fields.index('DP')]
                    AB = int(AD.split(",")[1]) / int(DP)
                    field_store.append(AB)

                if field != 'AB' and all_fields[format_fields.index(field)] not in ['./.','']:

                    index = format_fields.index(field)
                    value = all_fields[index]
                    field_store.append(value)
    
    return all_field_store

            

def check_condition(conditions, values, combine="&"):
    ''' Test the filtering conditions with the given values.

    Args:
        conditions: list of filtering conditions
        values: values in which to test the conditions against
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
    
    values = [float(x) if str(x).replace(".","").isdigit() else x for x in values ]
    
    tested_sample = []

    for cond, value in zip(conditions, values):

        field, op_sign, val = cond.split(" ")
        op = ops[op_sign]
        if val.replace(".", "").isdigit():
            val =float(val)
        # test the condition against the value
        test_cond = op(value, val) 
        tested_sample.append(test_cond)


    if combine == "&":
        return all(tested_sample)

    elif combine == "|":
        return any(tested_sample)



f = '/home/david/projects/pdVCF/test/vcfs/testing.vcf'
filter_vcf(f, ['GT = 1/1', 'DP > 100']) # == 3
#filter_vcf(f, ['GT = 1/1', 'DP > 100', 'AC = 2'])


f = '/home/david/projects/pdVCF/test/vcfs/testing3.vcf'
filter_vcf(f, ['GT = 1/1', 'DP > 100']) # == 257


