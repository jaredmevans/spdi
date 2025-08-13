# SPDI
Convert genomic variant coordinates to their normalized Canonical SPDI format.

## Install

Install latest SPDI version from github:

`pip install git+https://github.com/jaredmevans/spdi.git`

## Usage

Call SPDI function from your python scripts (recommended):

```
import spdi
import pysam

# load reference genome fasta
fasta = pysam.FastaFile("hs37d5.fa")

# get SPDI variant format
spdi_variant = spdi.convert("chr2", 47641510, "T", "TAT", fasta)
print(spdi_variant)
```

Output:

`NC_000002.11:47641509:TATAT:TATATAT`

Alternatively, the SPDI script can be run standalone for one-off conversions:

`python spdi.py --chr chr2 --pos 47641510 --ref T --alt TAT --genome GRCh37 --fasta hs37d5.fa`


## References
This code implements the algorithms described in the following publication.

Holmes, J. B., Moyer, E., Phan, L., Maglott, D., & Kattman, B. (2020). SPDI: data model for variants and applications at NCBI. *Bioinformatics (Oxford, England)*, 36(6), 1902–1907. https://doi.org/10.1093/bioinformatics/btz856