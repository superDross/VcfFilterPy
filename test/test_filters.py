import sys, os
# perform relative import (fuck python3)
path = os.path.dirname(os.path.abspath(__file__))[:-4]+"VcfFilterPy"
sys.path.append(path)
from VcfFilterPy.Vcf import Vcf
import unittest


class TestMandatory(unittest.TestCase):

    def test_pos_ranges(self):
        ''' Filter for variants between range

        Tested:
            SnpSift filter "( CHROM = '1') & ( POS > 2000 ) & ( POS < 2235601 )" testing2.vcf 
        '''
        v = Vcf('vcfs/testing2.vcf')
        v.filter_vcf(['CHROM = 1', 'POS > 2000', 'POS < 2235601'])

        self.assertEqual(v.count, 3)

    def test_loads(self):
        ''' Test numerous mandatory fields

        Tested:
            SnpSift filter "( FILTER != 'PASS' ) & ( ALT = 'G' ) & ( QUAL > 200 )" testing2.vcf
        '''
        v = Vcf('vcfs/testing2.vcf')
        v.filter_vcf(['FILTER != PASS', 'ALT = G', 'QUAL > 200'])

        self.assertEqual(v.count, 1)


class TestFilterGenotype(unittest.TestCase):
    
    def test_gt_filter(self):
        ''' Simple genotype filter

        Tested:
            SnpSift filter "( GEN[*].GT == '1/1' ) & ( GEN[*].DP > 100 )" testing.vcf
        '''
        v = Vcf('vcfs/testing.vcf')
        v.filter_vcf(['GT = 1/1', 'DP > 100'], how='any')

        self.assertEqual(v.count, 3)

    def test_gt_big_file(self):
        ''' Same test as above but with a larger vcf file

        Tested:
            SnpSift filter "( GEN[*].GT == '1/1' ) & ( GEN[*].DP > 100 )" testing3.vcf | grep -v '^#' | wc -l
        '''
        v = Vcf('vcfs/testing3.vcf')
        v.filter_vcf(['GT = 1/1', 'DP > 100'], how='any')

        self.assertEqual(v.count, 257)


    def test_AB(self):
        ''' Filter AB < 0.09

        Tested:
            I cannot find a tool to filter Vcfs by allele balance
        '''
        v = Vcf('vcfs/testing2.vcf')
        v.filter_vcf(['AB > 0', 'AB < 0.09'])

        self.assertEqual(1, v.count)



class TestFilterInfo(unittest.TestCase):

    def test_info_depth(self):
        ''' Testing filtering the INFO DP field

        Tested:
            SnpSift filter "( DP > 100 )" testing.vcf
        '''
        v = Vcf('vcfs/testing.vcf')
        v.filter_vcf(['DEPTH > 100'])

        self.assertEqual(v.count, 2)

    def test_info_float_int(self):
        ''' Filter int and float fields

        Tested:
            SnpSift filter "( AF[0] < 0.1 ) & ( AC[0] > 2 )" testing2.vcf 
        '''
        v = Vcf('vcfs/testing2.vcf')
        v.filter_vcf(['AF[0] < 0.1', 'AC[0] > 2'])

        self.assertEqual(v.count, 1)

    def test_string(self):
        ''' Filter str fields

        Tested:
            SnpSift filter "( MEOW != 'ON' )" testing2.vcf 
        '''
        v = Vcf('vcfs/testing2.vcf')
        v.filter_vcf(['MEOW != ON'])

        self.assertEqual(v.count, 3)


class TestFilterArgs(unittest.TestCase):

    def test_all_genotype(self):
        ''' Test filtering using the all flag in genotype fields

        Tested:
            SnpSift filter "( GEN[?].DP >= 50 ) & ( GEN[?].GQ >= 30 ) " testing.vcf
        '''
        v = Vcf('vcfs/testing.vcf')
        v.filter_vcf(['DP >= 50', 'GQ >= 30'], how='all')

        self.assertEqual(v.count, 1)

    def test_any_filtering(self):
        ''' Test filtering using the any flag

        Tested:
            SnpSift filter "( GEN[*].DP >= 50 ) & ( GEN[*].GQ >= 30 ) " testing.vcf
        '''
        v = Vcf('vcfs/testing.vcf')
        v.filter_vcf(['DP >= 50', 'GQ >= 30'], how='any')

        self.assertEqual(v.count, 6)

    #def test_or_operator(self):
    #    ''' Testing or operator

    #    Tested:
    #        SnpSift filter "( MEOW != 'ON' ) | ( LOLZ > 200 )" testing2.vcf
    #    '''
    #    v = Vcf('vcfs/testing2.vcf')
    #    v.filter_vcf(['MEOW != ON', 'LOLZ > 200'], op="|")

    #    self.assertEqual(v.count, 4)


class TestUltimateFilter(unittest.TestCase):
    
    def test_ultimate(self):
        ''' Test limits of filtering

        Tested:
         SnpSift filter "( GEN[*].GT == '0/1' ) & ( GEN[*].DP >= 50 ) & ( GEN[*].GQ >= 30 ) & ( AC[0] > 20 ) & ( CHROM = '1' ) & ( POS > 2235893 ) & ( ID =~ 'rs' )" testing3.vcf
        '''
        m = Vcf('vcfs/testing3.vcf')
        m.filter_vcf(['GT = 0/1', 'DP >= 50', 'GQ >= 30', 'AC[0] > 20', 'CHROM = 1', 'POS > 2235893', 'ID != .'])

        answer = ['1:218519928-A/AAAAC', '1:218578726-ACTCT/A,ACTCTCT,ACT', 
                  '1:218607557-G/T'] 

        self.assertEqual(m.count, 3)



