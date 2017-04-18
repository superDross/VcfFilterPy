from VcfLine import VcfLine

class Vcf(object):

    def __init__(self, vcf, command=None):
        self.vcf = open(vcf)
        self.filtered_vcf = []
        self.command = command # command run at the command line
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
            
            if line.startswith("#CHROM"):
                if self.command:
                    self.filtered_vcf.append('##source=' + ' '.join(self.command))

            if line.startswith("#"):
                self.filtered_vcf.append(line)
           
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



