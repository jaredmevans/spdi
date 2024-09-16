# SPDI
Generate the normalized SPDI Canonical allele representation for genomic variants in VCF style notation.

## Example

Command:

`python spdi.py --chr chr2 --pos 47641510 --ref T --alt TAT --genome GRCh37 --fasta hs37d5.fa`

Output:

`NC_000002.11:47641509:TATAT:TATATAT`

## Requirements:
- Python 3+
- pysam

## References
This code implements the algorithms described in the following publication.

Holmes, J. B., Moyer, E., Phan, L., Maglott, D., & Kattman, B. (2020). SPDI: data model for variants and applications at NCBI. *Bioinformatics (Oxford, England)*, 36(6), 1902–1907. https://doi.org/10.1093/bioinformatics/btz856