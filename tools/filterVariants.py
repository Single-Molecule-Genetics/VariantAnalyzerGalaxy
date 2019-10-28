"""filterVariants.py

Author -- Monika Heinzl
Contact -- monika.heinzl@edumail.at

Takes a tabular file with mutations and a vcf file from a Variant Caller and prints
all positions where an alternative allele is present to a user specified output file.
Creates tabular file of mutations.

=======  ==========  =================  ================================
Version  Date        Author             Description
0.1.1    2019-03-26  Monika Heinzl      -
=======  ==========  =================  ================================

USAGE: python filterVariants.py -variantAnnotatorFile VariantAnnotator_Mutations.tabular -vcfFile NVC_Variants.vcf -outputTabular DCS_Mutations.tabular
"""

import os
import vcf
import numpy as np
import argparse
import sys


def make_argparser():
    parser = argparse.ArgumentParser(description='Takes a tabular file with positions of mutations and a VCF file as input and prints all positions with an alternative allele. Creates a tabular file of mutations.')
    parser.add_argument('-variantAnnotatorFile',
                        help='TABULAR file with all positions of variants from the Variant Annotator.')
    parser.add_argument('-vcfFile',
                        help='VCF file from a Variant Caller.')
    parser.add_argument('-no_homozygousAlt', action="store_true",
                        help='If set to False (by default), the homozygous alternative alleles will not be filtered out.')
    parser.add_argument('-outputTabular',
                        help='Output TABULAR file of positions with an alternative allele.')

    return parser


def filterVariants(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])

    infile = args.variantAnnotatorFile
    in_vcf = args.vcfFile
    homozygousAlt = args.no_homozygousAlt
    outfile = args.outputTabular

    if os.path.isfile(infile) is False:
        print("Error: Could not find '{}'".format(infile))
        exit(0)

    if os.path.isfile(in_vcf) is False:
        print("Error: Could not find '{}'".format(in_vcf))
        exit(1)

    # read vcf file
    vcf_file = vcf.VCFReader(open(in_vcf, 'rb'))
    list_vcf = []
    for rec in vcf_file:
        list_vcf.append(np.array([rec.POS, rec.REF]))
    list_vcf = np.reshape(np.array(list_vcf), (len(list_vcf), 2))

    # filter tabular file with mutations
    with open(outfile, 'w') as out:
        header = "\t".join(
            ["#SAMPLE", "CHR", "POS", "A", "C", "G", "T", "CVRG", "ALLELES", "MAJOR", "MINOR", "MAF", "BIAS"])
        out.write(header)
        out.write("\n")
        with open(infile, 'r') as mutations_file:
            next(mutations_file)

            for line in mutations_file:
                line = line.rstrip('\n')
                splits = line.split('\t')
                index_vcf = np.where(list_vcf[:, 0] == splits[2])[0]
                ref = list_vcf[index_vcf, 1][0]

                if splits[10] in ["A", "T", "G", "C"]:
                    out.write("\t".join(splits))
                    out.write("\n")

                if homozygousAlt is False:
                    if str(splits[9]) != str(ref) and len(ref) == 1:
                        out.write("\t".join(splits[0:9]))
                        out.write("\t{}\t{}\t".format(ref, splits[9]))
                        out.write("\t".join(splits[11:12]))
                        out.write("\n")


if __name__ == '__main__':
    sys.exit(filterVariants(sys.argv))

