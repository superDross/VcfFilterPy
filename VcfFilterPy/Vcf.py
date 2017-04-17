from VcfLine import VcfLine

class Vcf(object):

    def __init__(self, vcf):
        self.vcf = open(vcf)
        self.filtered_vcf = []
        self.count = None

    def filter_vcf(self, conditions, how='any'):
        ''' Filter a Vcf based on a set of conditions.

        Args:
            vcf: a vcf file
            conditions: a list of filtering conditions
            how: specify whether all or any samples in the VCF need
                 to meet the conditions for the VCF to pass filter.
        '''
        count = 0

        while True:
            line = self.vcf.readline().strip("\n")
            
            if not line:
                break
            
            if line.startswith("#"):
                self.filtered_vcf.append(line)
           
            if line.startswith("#CHROM"):
                # Add the command used to filer just above #CHROM line
                pass

            if not line.startswith("#"):
                vline = VcfLine(line)
                fline = vline.filter_line(conditions, how)
                if fline:
                    self.filtered_vcf.append(fline)
                    count += 1 if fline else 0
        
        self.count = count

    def write_filtered_vcf(self, out):
        with open(out, 'w') as f:
            for l in self.filtered_vcf:
                f.write(l+"\n")


def lazy_test(f, conditions):
    vcf = Vcf(f)
    vcf.filter_vcf(conditions)
    print(vcf.count)

# TODO: specify fields with commas and MANDATORY fields

#f = '/home/david/projects/pdVCF/test/vcfs/testing.vcf'
#lazy_test(f, ['GT = 1/1', 'DP > 100', 'AB > 0']) # == 3
#lazy_test(f, ['GT = 1/1', 'DP > 100', 'AB > 0', 'AC < 50', 'DEPTH = 20']) # == 2

#f = '/home/david/projects/pdVCF/test/vcfs/testing3.vcf'
#lazy_test(f, ['GT = 1/1', 'DP > 100']) # == 257
#lazy_test(f, ['GT = 0/1', 'DP >= 50', 'GQ >= 30',  'DEPTH > 50']) # == 1896
#lazy_test(f, ['GT = 0/1', 'DP >= 50', 'GQ >= 30', 'AC[0] > 20', 'CHROM = 1', 'POS > 2235893', 'ID != .']) # == 3
#
#f = '/home/david/projects/pdVCF/test/vcfs/testing4.vcf'
#lazy_test(f, ['GT = 0/1', 'DP >= 50', 'AC[0] > 20', 'GQ >= 30']) # == 60
#lazy_test(f, ['GT = 0/1', 'DP >= 50', 'AC[0] > 20', 'GQ >= 30', 'AD[1] >= 20']) # == 49


