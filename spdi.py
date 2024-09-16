"""
Convert VCF style genomic coordinates and alleles to Canonical SPDI notation
"""

import argparse
import logging
import pysam


_CONTIG_ALIASES = {'GRCh38': {"chr1": "NC_000001.11", "chr2": "NC_000002.12", "chr3": "NC_000003.12",
                              "chr4": "NC_000004.12", "chr5": "NC_000005.10", "chr6": "NC_000006.12",
                              "chr7": "NC_000007.14", "chr8": "NC_000008.11", "chr9": "NC_000009.12",
                              "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
                              "chr13": "NC_000013.11", "chr14": "NC_000014.9", "chr15": "NC_000015.10",
                              "chr16": "NC_000016.10", "chr17": "NC_000017.11", "chr18": "NC_000018.10",
                              "chr19": "NC_000019.10", "chr20": "NC_000020.11", "chr21": "NC_000021.9",
                              "chr22": "NC_000022.11", "chrX": "NC_000023.11", "chrY": "NC_000024.10",
                              "chrM": "NC_012920.1"},
                   'GRCh37': {"chr1": "NC_000001.10", "chr2": "NC_000002.11", "chr3": "NC_000003.11",
                              "chr4": "NC_000004.11", "chr5": "NC_000005.9", "chr6": "NC_000006.11",
                              "chr7": "NC_000007.13", "chr8": "NC_000008.10", "chr9": "NC_000009.11",
                              "chr10": "NC_000010.10", "chr11": "NC_000011.9", "chr12": "NC_000012.11",
                              "chr13": "NC_000013.10", "chr14": "NC_000014.8", "chr15": "NC_000015.9",
                              "chr16": "NC_000016.9", "chr17": "NC_000017.10", "chr18": "NC_000018.9",
                              "chr19": "NC_000019.9", "chr20": "NC_000020.10", "chr21": "NC_000021.8",
                              "chr22": "NC_000022.10", "chrX": "NC_000023.10", "chrY": "NC_000024.9",
                              "chrM": "NC_001807.1"}}


def get_args():
    """
    commandline arguments
    :return: parsed args
    """
    parser = argparse.ArgumentParser(description='Convert VCF style genomic coordinates and alleles to Canonical SPDI')
    parser.add_argument('-c', '--chr', help='Chromosome', required=True)
    parser.add_argument('-p', '--pos', help='Position (1-based)', required=True, type=int)
    parser.add_argument('-r', '--ref', help='Ref Allele', required=True)
    parser.add_argument('-a', '--alt', help='Alt Allele', required=True)
    parser.add_argument('-f', '--fasta', help='Fasta for reference assembly', required=True)
    parser.add_argument('-g', '--genome', help='Genome build', default='GRCh38', choices=_CONTIG_ALIASES.keys())
    return parser.parse_args()


def right_trim_alleles(ref, alt):
    """
    Trim unnecessary bases off the right side of ref and alt alleles
    :param ref: reference allele
    :param alt: alternate allele
    :return: trimmed alleles tuple
    """
    new_ref = None
    new_alt = None
    for i in range(min(len(ref), len(alt)) - 1):
        if ref[len(ref) - 1 - i] == alt[len(alt) - 1 - i]:
            new_ref = ref[0:(len(ref) - 1 - i)]
            new_alt = alt[0:(len(alt) - 1 - i)]
        else:
            break
    if new_ref is not None:
        ref = new_ref
        alt = new_alt

    return ref, alt


def left_trim_alleles(pos, ref, alt, over_trim=False):
    """
    Trim unnecessary bases off the left side of the ref and alt alleles. Adjust position as well.
    :param pos: 1-based coordinate of variant
    :param ref: reference allele
    :param alt: alternate allele
    :param over_trim: Turn on over-trimming, which trims one base further than the VCF standard of keeping the left base matching for ref/alt indels
    :return: trimmed pos and alleles tuple
    """
    new_pos = None
    new_ref = None
    new_alt = None
    offset = 1
    if over_trim:
        offset = 0
    for i in range(min(len(ref), len(alt)) - offset):
        if ref[i] == alt[i]:
            new_ref = ref[(i + 1):]
            new_alt = alt[(i + 1):]
            new_pos = pos + i + 1
        else:
            break
    if new_ref is not None:
        ref = new_ref
        alt = new_alt
        pos = new_pos

    return pos, ref, alt


def grow_alleles(contig, pos, ref, alt, fasta):
    """
    As part of NCBI's VOCA algorithm, grow the alleles in both directions to adjust overly precise variant notation
    :param contig: contig
    :param pos: 1-based start position of variant
    :param ref: reference allele
    :param alt: alternate allele
    :param fasta: pysam fasta file object
    :return: grown position and alleles tuple
    """
    bud_seq = None
    if len(ref) == 0:
        bud_seq = alt
    elif len(alt) == 0:
        bud_seq = ref
    if bud_seq is not None:
        pos_offset = 0
        if len(ref) != 0:
            # use bud seq len for deletions only
            pos_offset = len(bud_seq)
        # try growing right side
        right_counter = 1
        right_adjacent = fasta.fetch(contig, pos+pos_offset-1, pos+pos_offset+len(bud_seq)-1)
        right_full_flank = right_adjacent
        while right_adjacent == bud_seq:
            right_counter += 1
            right_adjacent = fasta.fetch(contig, pos + pos_offset + (len(bud_seq) * (right_counter-1)) - 1,
                                         pos + pos_offset + (len(bud_seq) * right_counter) - 1)
            right_full_flank = right_full_flank + right_adjacent
        # trim any bases that don't match from the last full bud seq interval
        for i in range(len(bud_seq)):
            if bud_seq[i] != right_full_flank[(len(bud_seq) * (right_counter - 1)) + i]:
                right_full_flank = right_full_flank[0:(len(bud_seq) * (right_counter - 1)) + i]
                break
        # try growing left side
        left_counter = 1
        left_adjacent = fasta.fetch(contig, pos-len(bud_seq)-1, pos-1)
        left_full_flank = left_adjacent
        while left_adjacent == bud_seq:
            left_counter += 1
            left_adjacent = fasta.fetch(contig, (pos - (len(bud_seq) * left_counter) - 1),
                                        (pos - (len(bud_seq) * (left_counter - 1)) - 1))
            left_full_flank = left_adjacent + left_full_flank
        pos = pos - (len(bud_seq) * left_counter)
        # trim any bases that don't match from last full bud seq interval
        for i in range(len(bud_seq)-1, -1, -1):
            if bud_seq[i] != left_full_flank[i]:
                left_full_flank = left_full_flank[i+1:]
                pos += i + 1
                break
        ref = left_full_flank + ref + right_full_flank
        alt = left_full_flank + alt + right_full_flank
    else:
        logging.error("Can't find necessary bud sequence to start allele growing")
    return pos, ref, alt


def spdi_normalize(contig, pos, ref, alt, fasta):
    """
    Apply SPDI's Variant Overprecision Correction Algorithm (VOCA) to normalize variant
    :param contig: variant contig
    :param pos: coordinate of variant (1-based)
    :param ref: Reference Allele
    :param alt: Alternate Allele
    :param fasta: reference fasta as pysam FastaFile object
    :return: normalized variant tuple. Note: pos is now 0-based
    """
    new_pos = None
    new_ref = None
    new_alt = None
    if ref != alt:
        # Variant Overprecision Correction Algorithm (VOCA)
        # 1. Trim right and left
        ref, alt = right_trim_alleles(ref, alt)
        pos, ref, alt = left_trim_alleles(pos, ref, alt, over_trim=True)
        bud_seq = None
        if len(ref) == 0:
            bud_seq = alt
        elif len(alt) == 0:
            bud_seq = ref
        if bud_seq is not None:
            # 2. Grow sequence in both directions to see if current representation is overly precise
            pos, ref, alt = grow_alleles(contig, pos, ref, alt, fasta)
        new_pos = pos - 1
        new_ref = ref
        new_alt = alt
    return new_pos, new_ref, new_alt


def main():
    # set logging format
    logging.basicConfig(format='%(levelname)s\t%(asctime)s - %(message)s', level=logging.INFO)
    # parse arguments
    args = get_args()

    fasta = pysam.FastaFile(args.fasta)
    if args.chr not in fasta.references:
        logging.error("Chromosome not found in reference FASTA: " + args.chr)
    elif args.chr not in _CONTIG_ALIASES[args.genome]:
        logging.error("Chromosome not in alias list: " + args.chr)
    else:
        # convert to Canonical SPDI notation
        ncbi_alias = _CONTIG_ALIASES[args.genome][args.chr]
        new_pos, new_ref, new_alt = spdi_normalize(args.chr, args.pos, args.ref, args.alt, fasta)
        print(":".join([ncbi_alias, str(new_pos), new_ref, new_alt]))


if __name__ == '__main__':
    main()
