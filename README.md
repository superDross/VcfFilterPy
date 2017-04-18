# VcfFilterPy
Perform complex filtering upon a VCF file.

## Install
To install and test VcfFilterPy
```bash
git clone https://github.com/superDross/VcfFilterPy
cd VcfFilterPy/
python3 VcfFilterPy --test
```
## Why another VCF filtering tool?
This project was born out of the inability to perform complex filtering upon certain fields such as allele balance with existing vcf filtering tools and it has been created to allow one to easily calculate custom fields and filter them. This can be done by altering the create_custom_fields.py file.

## Example Usage
Filtering can be done at the python interpretor or at the command line. The below example shows one how to do this at both by filtering a vcf for variants in which at least one sample has an allele balance above 0.3 and an ALT allele depth above 30.

At the python interpreter:
```python
from VcfFilterPy.Vcf import Vcf

# initilise a Vcf object and filter
vcf = VCF("test.vcf")
vcf.filter_vcf(['AB > 0.3', 'AD[1] > 30'], how='any')

# get the number of variants remaining after filtering
vcf.write_filtered_vcf("out.vcf")

# print the number variants remianing after filtering
vcf.count
```

At the command line:
```bash
python3 VcfFilterPy --vcf test.vcf \
                    --filter 'AB > 0.3, AD[1] > 30' \
                    --operator 'any' \
                    --out out.vcf
```

Extremely complex filtering can be performed such as; filter for variants in chromosome 1, above positions 2235893, has an reference allele count above 20 and in which at least one sample has a genotype 0/1, genotype depth of at least 50 and genotype quality of at least 30:
```bash
python3 VcfFilterPy \
--vcf test.py \
--filter 'CHROM = 1, POS > 2235893, AC[0] > 20, ID != ., GT = 0/1, DP >= 50, GQ >= 30' \
--operator 'any'
--out out.vcf
```
