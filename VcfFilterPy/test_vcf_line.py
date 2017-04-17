''' Test a VcfLine dict against a set of conditions.
'''
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



def test_vcf_line(conditions, sam_dict, combine="&"):
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


