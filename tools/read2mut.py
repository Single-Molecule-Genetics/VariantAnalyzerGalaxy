#!/usr/bin/env python

"""read2mut.py

Author -- Gundula Povysil
Contact -- povysil@bioinf.jku.at

Looks for reads with mutation at known
positions and calculates frequencies and stats.

=======  ==========  =================  ================================
Version  Date        Author             Description
0.2.2    2019-10-27  Gundula Povysil    -
=======  ==========  =================  ================================


USAGE: python read2mut.py --mutFile DCS_Mutations.tabular --bamFile Interesting_Reads.trim.bam
                          --inputJson tag_count_dict.json --sscsJson SSCS_counts.json
                          --outputFile mutant_reads_summary_short_trim.xlsx --thresh 10 --phred 20 --trim 10 --chimera_correction

"""

from __future__ import division

import argparse
import csv
import itertools
import json
import operator
import os
import re
import sys

import numpy as np
import pysam
import xlsxwriter
from cyvcf2 import VCF


def make_argparser():
    parser = argparse.ArgumentParser(description='Takes a VCF file with mutations, a BAM file and JSON files as input and prints stats about variants to a user specified output file.')
    parser.add_argument('--mutFile',
                        help='VCF file with DCS mutations.')
    parser.add_argument('--bamFile',
                        help='BAM file with aligned raw reads of selected tags (FASTQ created by mut2read.py - trimming with Trimmomatic - alignment with bwa).')
    parser.add_argument('--inputJson',
                        help='JSON file with data collected by mut2read.py.')
    parser.add_argument('--sscsJson',
                        help='JSON file with SSCS counts collected by mut2sscs.py.')
    parser.add_argument('--outputFile',
                        help='Output xlsx file with summary of mutations.')
    parser.add_argument('--outputFile_csv',
                        help='Output csv file with summary of mutations.')
    parser.add_argument('--outputFile2',
                        help='Output xlsx file with allele frequencies of mutations.')
    parser.add_argument('--outputFile3',
                        help='Output xlsx file with examples of the tier classification.')
    parser.add_argument('--thresh', type=int, default=0,
                        help='Integer threshold for displaying mutations. Only mutations occuring less than thresh times are displayed. Default of 0 displays all.')
    parser.add_argument('--phred', type=int, default=20,
                        help='Integer threshold for Phred score. Only reads higher than this threshold are considered. Default 20.')
    parser.add_argument('--trim', type=int, default=10,
                        help='Integer threshold for assigning mutations at start and end of reads to lower tier. Default 10.')
    parser.add_argument('--chimera_correction', action="store_true",
                        help='Count chimeric variants and correct the variant frequencies')
    parser.add_argument('--softclipping_dist',  type=int, default=15,
                        help='Count mutation as an artifact if mutation lies within this parameter away from the softclipping part of the read.')
    parser.add_argument('--reads_threshold',  type=float, default=1.0,
                        help='Float number which specifies the minimum percentage of softclipped reads in a family to be considered in the softclipping tiers. Default: 1.0, means all reads of a family have to be softclipped.')
    return parser


def safe_div(x, y):
    if y == 0:
        return None
    return x / y


def read2mut(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])
    file1 = args.mutFile
    file2 = args.bamFile
    json_file = args.inputJson
    sscs_json = args.sscsJson
    outfile = args.outputFile
    outfile2 = args.outputFile2
    outfile3 = args.outputFile3
    outputFile_csv = args.outputFile_csv
    thresh = args.thresh
    phred_score = args.phred
    trim = args.trim
    chimera_correction = args.chimera_correction
    thr = args.softclipping_dist
    threshold_reads = args.reads_threshold

    if os.path.isfile(file1) is False:
        sys.exit("Error: Could not find '{}'".format(file1))
    if os.path.isfile(file2) is False:
        sys.exit("Error: Could not find '{}'".format(file2))
    if os.path.isfile(json_file) is False:
        sys.exit("Error: Could not find '{}'".format(json_file))
    if thresh < 0:
        sys.exit("Error: thresh is '{}', but only non-negative integers allowed".format(thresh))
    if phred_score < 0:
        sys.exit("Error: phred is '{}', but only non-negative integers allowed".format(phred_score))
    if trim < 0:
        sys.exit("Error: trim is '{}', but only non-negative integers allowed".format(thresh))
    if thr <= 0:
        sys.exit("Error: trim is '{}', but only non-negative integers allowed".format(thr))

    # load dicts
    with open(json_file, "r") as f:
        (tag_dict, cvrg_dict) = json.load(f)

    with open(sscs_json, "r") as f:
        (mut_pos_dict, ref_pos_dict) = json.load(f)

    # read bam file
    # pysam.index(file2)
    bam = pysam.AlignmentFile(file2, "rb")

    # create mut_dict
    mut_dict = {}
    mut_read_pos_dict = {}
    mut_read_dict = {}
    reads_dict = {}
    mut_read_cigar_dict = {}
    real_start_end = {}
    i = 0
    mut_array = []

    for count, variant in enumerate(VCF(file1)):
        chrom = variant.CHROM
        stop_pos = variant.start
        ref = variant.REF
        if len(variant.ALT) == 0:
            continue
        else:
            alt = variant.ALT[0]
        chrom_stop_pos = str(chrom) + "#" + str(stop_pos) + "#" + ref + "#" + alt

        if len(ref) == len(alt):
            mut_array.append([chrom, stop_pos, ref, alt])
            i += 1
            mut_dict[chrom_stop_pos] = {}
            mut_read_pos_dict[chrom_stop_pos] = {}
            reads_dict[chrom_stop_pos] = {}
            mut_read_cigar_dict[chrom_stop_pos] = {}
            real_start_end[chrom_stop_pos] = {}

            for pileupcolumn in bam.pileup(chrom, stop_pos - 1, stop_pos + 1, max_depth=100000000):
                if pileupcolumn.reference_pos == stop_pos:
                    count_alt = 0
                    count_ref = 0
                    count_indel = 0
                    count_n = 0
                    count_other = 0
                    count_lowq = 0
                    n = 0
                    for pileupread in pileupcolumn.pileups:
                        n += 1
                        if not pileupread.is_del and not pileupread.is_refskip:
                            tag = pileupread.alignment.query_name
                            nuc = pileupread.alignment.query_sequence[pileupread.query_position]
                            phred = ord(pileupread.alignment.qual[pileupread.query_position]) - 33
                            # if read is softclipped, store real position in reference
                            if "S" in pileupread.alignment.cigarstring:
                                # spftclipped at start
                                if re.search(r"^[0-9]+S", pileupread.alignment.cigarstring):
                                    start = pileupread.alignment.reference_start - int(pileupread.alignment.cigarstring.split("S")[0])
                                    end = pileupread.alignment.reference_end
                                # softclipped at end
                                elif re.search(r"S$", pileupread.alignment.cigarstring):
                                    end = pileupread.alignment.reference_end + int(re.split("[A-Z]", str(pileupread.alignment.cigarstring))[-2])
                                    start = pileupread.alignment.reference_start
                            else:
                                end = pileupread.alignment.reference_end
                                start = pileupread.alignment.reference_start

                            if phred < phred_score:
                                nuc = "lowQ"
                            if tag not in mut_dict[chrom_stop_pos]:
                                mut_dict[chrom_stop_pos][tag] = {}
                            if nuc in mut_dict[chrom_stop_pos][tag]:
                                mut_dict[chrom_stop_pos][tag][nuc] += 1
                            else:
                                mut_dict[chrom_stop_pos][tag][nuc] = 1
                            if tag not in mut_read_pos_dict[chrom_stop_pos]:
                                mut_read_pos_dict[chrom_stop_pos][tag] = [pileupread.query_position + 1]
                                reads_dict[chrom_stop_pos][tag] = [len(pileupread.alignment.query_sequence)]
                                mut_read_cigar_dict[chrom_stop_pos][tag] = [pileupread.alignment.cigarstring]
                                real_start_end[chrom_stop_pos][tag] = [(start, end)]
                            else:
                                mut_read_pos_dict[chrom_stop_pos][tag].append(pileupread.query_position + 1)
                                reads_dict[chrom_stop_pos][tag].append(len(pileupread.alignment.query_sequence))
                                mut_read_cigar_dict[chrom_stop_pos][tag].append(pileupread.alignment.cigarstring)
                                real_start_end[chrom_stop_pos][tag].append((start, end))
                            if nuc == alt:
                                count_alt += 1
                                if tag not in mut_read_dict:
                                    mut_read_dict[tag] = {}
                                    mut_read_dict[tag][chrom_stop_pos] = (alt, ref)
                                else:
                                    mut_read_dict[tag][chrom_stop_pos] = (alt, ref)
                            elif nuc == ref:
                                count_ref += 1
                            elif nuc == "N":
                                count_n += 1
                            elif nuc == "lowQ":
                                count_lowq += 1
                            else:
                                count_other += 1
                        else:
                            count_indel += 1

    mut_array = np.array(mut_array)
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            pure_tag = read.query_name[:-5]
            nuc = "na"
            for key in tag_dict[pure_tag].keys():
                if key not in mut_dict:
                    mut_dict[key] = {}
                if read.query_name not in mut_dict[key]:
                    mut_dict[key][read.query_name] = {}
                if nuc in mut_dict[key][read.query_name]:
                    mut_dict[key][read.query_name][nuc] += 1
                else:
                    mut_dict[key][read.query_name][nuc] = 1
    bam.close()

    # create pure_tags_dict
    pure_tags_dict = {}
    for key1, value1 in sorted(mut_dict.items()):
        i = np.where(np.array(['#'.join(str(i) for i in z)
                               for z in zip(mut_array[:, 0], mut_array[:, 1], mut_array[:, 2], mut_array[:, 3])]) == key1)[0][0]
        ref = mut_array[i, 2]
        alt = mut_array[i, 3]
        pure_tags_dict[key1] = {}
        for key2, value2 in sorted(value1.items()):
            for key3, value3 in value2.items():
                pure_tag = key2[:-5]
                if key3 == alt:
                    if pure_tag in pure_tags_dict[key1]:
                        pure_tags_dict[key1][pure_tag] += 1
                    else:
                        pure_tags_dict[key1][pure_tag] = 1

    # create pure_tags_dict_short with thresh
    if thresh > 0:
        pure_tags_dict_short = {}
        for key, value in sorted(pure_tags_dict.items()):
            if len(value) < thresh:
                pure_tags_dict_short[key] = value
    else:
        pure_tags_dict_short = pure_tags_dict

    csv_data = open(outputFile_csv, "w")
    csv_writer = csv.writer(csv_data, delimiter=",")

    # output summary with threshold
    workbook = xlsxwriter.Workbook(outfile)
    workbook2 = xlsxwriter.Workbook(outfile2)
    workbook3 = xlsxwriter.Workbook(outfile3)
    ws1 = workbook.add_worksheet("Results")
    ws2 = workbook2.add_worksheet("Allele frequencies")
    ws3 = workbook3.add_worksheet("Tiers")

    format1 = workbook.add_format({'bg_color': '#BCF5A9'})  # green
    format2 = workbook.add_format({'bg_color': '#FFC7CE'})  # red
    format3 = workbook.add_format({'bg_color': '#FACC2E'})  # yellow

    format12 = workbook2.add_format({'bg_color': '#BCF5A9'})  # green
    format22 = workbook2.add_format({'bg_color': '#FFC7CE'})  # red
    format32 = workbook2.add_format({'bg_color': '#FACC2E'})  # yellow

    format13 = workbook3.add_format({'bg_color': '#BCF5A9'})  # green
    format23 = workbook3.add_format({'bg_color': '#FFC7CE'})  # red
    format33 = workbook3.add_format({'bg_color': '#FACC2E'})  # yellow

    header_line = ('variant ID', 'tier', 'tag', 'mate', 'read pos.ab', 'read pos.ba', 'read median length.ab',
                   'read median length.ba', 'DCS median length',
                   'FS.ab', 'FS.ba', 'FSqc.ab', 'FSqc.ba', 'ref.ab', 'ref.ba', 'alt.ab', 'alt.ba',
                   'rel. ref.ab', 'rel. ref.ba', 'rel. alt.ab', 'rel. alt.ba',
                   'na.ab', 'na.ba', 'lowq.ab', 'lowq.ba', 'trim.ab', 'trim.ba',
                   'SSCS alt.ab', 'SSCS alt.ba', 'SSCS ref.ab', 'SSCS ref.ba',
                   'in phase', 'chimeric tag')
    ws1.write_row(0, 0, header_line)
    csv_writer.writerow(header_line)
    counter_tier11 = 0
    counter_tier12 = 0
    counter_tier21 = 0
    counter_tier22 = 0
    counter_tier23 = 0
    counter_tier24 = 0
    counter_tier31 = 0
    counter_tier32 = 0
    counter_tier25 = 0
    counter_tier4 = 0
    counter_tier51 = 0
    counter_tier52 = 0
    counter_tier53 = 0
    counter_tier54 = 0
    counter_tier55 = 0
    counter_tier6 = 0
    counter_tier7 = 0

    row = 1
    tier_dict = {}
    chimera_dict = {}
    for key1, value1 in sorted(mut_dict.items()):
        counts_mut = 0
        chimeric_tag_list = []
        chimeric_tag = {}
        if key1 in pure_tags_dict_short.keys():
            change_tier_after_print = []
            i = np.where(np.array(['#'.join(str(i) for i in z)
                                   for z in zip(mut_array[:, 0], mut_array[:, 1], mut_array[:, 2], mut_array[:, 3])]) == key1)[0][0]
            ref = mut_array[i, 2]
            alt = mut_array[i, 3]
            dcs_median = cvrg_dict[key1][2]
            whole_array = list(pure_tags_dict_short[key1].keys())

            tier_dict[key1] = {}
            values_tier_dict = [("tier 1.1", 0), ("tier 1.2", 0), ("tier 2.1", 0), ("tier 2.2", 0), ("tier 2.3", 0), ("tier 2.4", 0), ("tier 2.5", 0),
                                ("tier 3.1", 0), ("tier 3.2", 0), ("tier 4", 0), ("tier 5.1", 0), ("tier 5.2", 0), ("tier 5.3", 0), ("tier 5.4", 0), ("tier 5.5", 0),
                                ("tier 6", 0), ("tier 7", 0)]
            for k, v in values_tier_dict:
                tier_dict[key1][k] = v

            used_keys = []
            if 'ab' in mut_pos_dict[key1].keys():
                sscs_mut_ab = mut_pos_dict[key1]['ab']
            else:
                sscs_mut_ab = 0
            if 'ba' in mut_pos_dict[key1].keys():
                sscs_mut_ba = mut_pos_dict[key1]['ba']
            else:
                sscs_mut_ba = 0
            if 'ab' in ref_pos_dict[key1].keys():
                sscs_ref_ab = ref_pos_dict[key1]['ab']
            else:
                sscs_ref_ab = 0
            if 'ba' in ref_pos_dict[key1].keys():
                sscs_ref_ba = ref_pos_dict[key1]['ba']
            else:
                sscs_ref_ba = 0
            for key2, value2 in sorted(value1.items()):
                add_mut14 = ""
                add_mut23 = ""
                if (key2[:-5] in pure_tags_dict_short[key1].keys()) and (key2[:-5] not in used_keys) and (key1 in tag_dict[key2[:-5]].keys()):
                    if key2[:-5] + '.ab.1' in mut_dict[key1].keys():
                        total1 = sum(mut_dict[key1][key2[:-5] + '.ab.1'].values())
                        if 'na' in mut_dict[key1][key2[:-5] + '.ab.1'].keys():
                            na1 = mut_dict[key1][key2[:-5] + '.ab.1']['na']
                            # na1f = na1/total1
                        else:
                            # na1 = na1f = 0
                            na1 = 0
                        if 'lowQ' in mut_dict[key1][key2[:-5] + '.ab.1'].keys():
                            lowq1 = mut_dict[key1][key2[:-5] + '.ab.1']['lowQ']
                            # lowq1f = lowq1 / total1
                        else:
                            # lowq1 = lowq1f = 0
                            lowq1 = 0
                        if ref in mut_dict[key1][key2[:-5] + '.ab.1'].keys():
                            ref1 = mut_dict[key1][key2[:-5] + '.ab.1'][ref]
                            ref1f = ref1 / (total1 - na1 - lowq1)
                        else:
                            ref1 = ref1f = 0
                        if alt in mut_dict[key1][key2[:-5] + '.ab.1'].keys():
                            alt1 = mut_dict[key1][key2[:-5] + '.ab.1'][alt]
                            alt1f = alt1 / (total1 - na1 - lowq1)
                        else:
                            alt1 = alt1f = 0
                        total1new = total1 - na1 - lowq1
                        if (key2[:-5] + '.ab.1') in mut_read_dict.keys():
                            k1 = mut_read_dict[(key2[:-5] + '.ab.1')].keys()
                            add_mut1 = len(k1)
                            if add_mut1 > 1:
                                for k, v in mut_read_dict[(key2[:-5] + '.ab.1')].items():
                                    if k != key1:
                                        new_mut = str(k).split("#")[0] + "-" + str(int(str(k).split("#")[1]) + 1) + "-" + v[1] + "-" + v[0]
                                        if len(add_mut14) == 0:
                                            add_mut14 = new_mut
                                        else:
                                            add_mut14 = add_mut14 + ", " + new_mut
                        else:
                            k1 = []
                    else:
                        total1 = total1new = na1 = lowq1 = 0
                        ref1 = alt1 = ref1f = alt1f = 0
                        k1 = []

                    if key2[:-5] + '.ab.2' in mut_dict[key1].keys():
                        total2 = sum(mut_dict[key1][key2[:-5] + '.ab.2'].values())
                        if 'na' in mut_dict[key1][key2[:-5] + '.ab.2'].keys():
                            na2 = mut_dict[key1][key2[:-5] + '.ab.2']['na']
                            # na2f = na2 / total2
                        else:
                            # na2 = na2f = 0
                            na2 = 0
                        if 'lowQ' in mut_dict[key1][key2[:-5] + '.ab.2'].keys():
                            lowq2 = mut_dict[key1][key2[:-5] + '.ab.2']['lowQ']
                            # lowq2f = lowq2 / total2
                        else:
                            # lowq2 = lowq2f = 0
                            lowq2 = 0
                        if ref in mut_dict[key1][key2[:-5] + '.ab.2'].keys():
                            ref2 = mut_dict[key1][key2[:-5] + '.ab.2'][ref]
                            ref2f = ref2 / (total2 - na2 - lowq2)
                        else:
                            ref2 = ref2f = 0
                        if alt in mut_dict[key1][key2[:-5] + '.ab.2'].keys():
                            alt2 = mut_dict[key1][key2[:-5] + '.ab.2'][alt]
                            alt2f = alt2 / (total2 - na2 - lowq2)
                        else:
                            alt2 = alt2f = 0
                        total2new = total2 - na2 - lowq2
                        if (key2[:-5] + '.ab.2') in mut_read_dict.keys():
                            k2 = mut_read_dict[(key2[:-5] + '.ab.2')].keys()
                            add_mut2 = len(k2)
                            if add_mut2 > 1:
                                for k, v in mut_read_dict[(key2[:-5] + '.ab.2')].items():
                                    if k != key1:
                                        new_mut = str(k).split("#")[0] + "-" + str(int(str(k).split("#")[1]) + 1) + "-" + v[1] + "-" + v[0]
                                        if len(add_mut23) == 0:
                                            add_mut23 = new_mut
                                        else:
                                            add_mut23 = add_mut23 + ", " + new_mut
                        else:
                            k2 = []
                    else:
                        total2 = total2new = na2 = lowq2 = 0
                        ref2 = alt2 = ref2f = alt2f = 0
                        k2 = []

                    if key2[:-5] + '.ba.1' in mut_dict[key1].keys():
                        total3 = sum(mut_dict[key1][key2[:-5] + '.ba.1'].values())
                        if 'na' in mut_dict[key1][key2[:-5] + '.ba.1'].keys():
                            na3 = mut_dict[key1][key2[:-5] + '.ba.1']['na']
                            # na3f = na3 / total3
                        else:
                            # na3 = na3f = 0
                            na3 = 0
                        if 'lowQ' in mut_dict[key1][key2[:-5] + '.ba.1'].keys():
                            lowq3 = mut_dict[key1][key2[:-5] + '.ba.1']['lowQ']
                            # lowq3f = lowq3 / total3
                        else:
                            # lowq3 = lowq3f = 0
                            lowq3 = 0
                        if ref in mut_dict[key1][key2[:-5] + '.ba.1'].keys():
                            ref3 = mut_dict[key1][key2[:-5] + '.ba.1'][ref]
                            ref3f = ref3 / (total3 - na3 - lowq3)
                        else:
                            ref3 = ref3f = 0
                        if alt in mut_dict[key1][key2[:-5] + '.ba.1'].keys():
                            alt3 = mut_dict[key1][key2[:-5] + '.ba.1'][alt]
                            alt3f = alt3 / (total3 - na3 - lowq3)
                        else:
                            alt3 = alt3f = 0
                        total3new = total3 - na3 - lowq3
                        if (key2[:-5] + '.ba.1') in mut_read_dict.keys():
                            add_mut3 = len(mut_read_dict[(key2[:-5] + '.ba.1')].keys())
                            if add_mut3 > 1:
                                for k, v in mut_read_dict[(key2[:-5] + '.ba.1')].items():
                                    if k != key1 and k not in k2:
                                        new_mut = str(k).split("#")[0] + "-" + str(int(str(k).split("#")[1]) + 1) + "-" + v[1] + "-" + v[0]
                                        if len(add_mut23) == 0:
                                            add_mut23 = new_mut
                                        else:
                                            add_mut23 = add_mut23 + ", " + new_mut
                    else:
                        total3 = total3new = na3 = lowq3 = 0
                        ref3 = alt3 = ref3f = alt3f = 0

                    if key2[:-5] + '.ba.2' in mut_dict[key1].keys():
                        total4 = sum(mut_dict[key1][key2[:-5] + '.ba.2'].values())
                        if 'na' in mut_dict[key1][key2[:-5] + '.ba.2'].keys():
                            na4 = mut_dict[key1][key2[:-5] + '.ba.2']['na']
                            # na4f = na4 / total4
                        else:
                            # na4 = na4f = 0
                            na4 = 0
                        if 'lowQ' in mut_dict[key1][key2[:-5] + '.ba.2'].keys():
                            lowq4 = mut_dict[key1][key2[:-5] + '.ba.2']['lowQ']
                            # lowq4f = lowq4 / total4
                        else:
                            # lowq4 = lowq4f = 0
                            lowq4 = 0
                        if ref in mut_dict[key1][key2[:-5] + '.ba.2'].keys():
                            ref4 = mut_dict[key1][key2[:-5] + '.ba.2'][ref]
                            ref4f = ref4 / (total4 - na4 - lowq4)
                        else:
                            ref4 = ref4f = 0
                        if alt in mut_dict[key1][key2[:-5] + '.ba.2'].keys():
                            alt4 = mut_dict[key1][key2[:-5] + '.ba.2'][alt]
                            alt4f = alt4 / (total4 - na4 - lowq4)
                        else:
                            alt4 = alt4f = 0
                        total4new = total4 - na4 - lowq4
                        if (key2[:-5] + '.ba.2') in mut_read_dict.keys():
                            add_mut4 = len(mut_read_dict[(key2[:-5] + '.ba.2')].keys())
                            if add_mut4 > 1:
                                for k, v in mut_read_dict[(key2[:-5] + '.ba.2')].items():
                                    if k != key1 and k not in k1:
                                        new_mut = str(k).split("#")[0] + "-" + str(int(str(k).split("#")[1]) + 1) + "-" + v[1] + "-" + v[0]
                                        if len(add_mut14) == 0:
                                            add_mut14 = new_mut
                                        else:
                                            add_mut14 = add_mut14 + ", " + new_mut
                    else:
                        total4 = total4new = na4 = lowq4 = 0
                        ref4 = alt4 = ref4f = alt4f = 0

                    read_pos1 = read_pos2 = read_pos3 = read_pos4 = -1
                    read_len_median1 = read_len_median2 = read_len_median3 = read_len_median4 = 0
                    cigars_dcs1 = cigars_dcs2 = cigars_dcs3 = cigars_dcs4 = []
                    pos_read1 = pos_read2 = pos_read3 = pos_read4 = []
                    end_read1 = end_read2 = end_read3 = end_read4 = []
                    if key2[:-5] + '.ab.1' in mut_read_pos_dict[key1].keys():
                        read_pos1 = np.median(np.array(mut_read_pos_dict[key1][key2[:-5] + '.ab.1']))
                        read_len_median1 = np.median(np.array(reads_dict[key1][key2[:-5] + '.ab.1']))
                        cigars_dcs1 = mut_read_cigar_dict[key1][key2[:-5] + '.ab.1']
                        pos_read1 = mut_read_pos_dict[key1][key2[:-5] + '.ab.1']
                        end_read1 = reads_dict[key1][key2[:-5] + '.ab.1']
                        ref_positions1 = real_start_end[key1][key2[:-5] + '.ab.1']
                    if key2[:-5] + '.ab.2' in mut_read_pos_dict[key1].keys():
                        read_pos2 = np.median(np.array(mut_read_pos_dict[key1][key2[:-5] + '.ab.2']))
                        read_len_median2 = np.median(np.array(reads_dict[key1][key2[:-5] + '.ab.2']))
                        cigars_dcs2 = mut_read_cigar_dict[key1][key2[:-5] + '.ab.2']
                        pos_read2 = mut_read_pos_dict[key1][key2[:-5] + '.ab.2']
                        end_read2 = reads_dict[key1][key2[:-5] + '.ab.2']
                        ref_positions2 = real_start_end[key1][key2[:-5] + '.ab.2']
                    if key2[:-5] + '.ba.1' in mut_read_pos_dict[key1].keys():
                        read_pos3 = np.median(np.array(mut_read_pos_dict[key1][key2[:-5] + '.ba.1']))
                        read_len_median3 = np.median(np.array(reads_dict[key1][key2[:-5] + '.ba.1']))
                        cigars_dcs3 = mut_read_cigar_dict[key1][key2[:-5] + '.ba.1']
                        pos_read3 = mut_read_pos_dict[key1][key2[:-5] + '.ba.1']
                        end_read3 = reads_dict[key1][key2[:-5] + '.ba.1']
                        ref_positions3 = real_start_end[key1][key2[:-5] + '.ba.1']
                    if key2[:-5] + '.ba.2' in mut_read_pos_dict[key1].keys():
                        read_pos4 = np.median(np.array(mut_read_pos_dict[key1][key2[:-5] + '.ba.2']))
                        read_len_median4 = np.median(np.array(reads_dict[key1][key2[:-5] + '.ba.2']))
                        cigars_dcs4 = mut_read_cigar_dict[key1][key2[:-5] + '.ba.2']
                        pos_read4 = mut_read_pos_dict[key1][key2[:-5] + '.ba.2']
                        end_read4 = reads_dict[key1][key2[:-5] + '.ba.2']
                        ref_positions4 = real_start_end[key1][key2[:-5] + '.ba.2']

                    used_keys.append(key2[:-5])
                    counts_mut += 1
                    if (alt1f + alt2f + alt3f + alt4f) > 0.5:
                        total1new_trim, total2new_trim, total3new_trim, total4new_trim = total1new, total2new, total3new, total4new
                        if total1new == 0:
                            ref1f = alt1f = None
                            alt1ff = -1
                            alt1ff_trim = -1
                        else:
                            alt1ff = alt1f
                            alt1ff_trim = alt1f
                        if total2new == 0:
                            ref2f = alt2f = None
                            alt2ff = -1
                            alt2ff_trim = -1
                        else:
                            alt2ff = alt2f
                            alt2ff_trim = alt2f
                        if total3new == 0:
                            ref3f = alt3f = None
                            alt3ff = -1
                            alt3ff_trim = -1
                        else:
                            alt3ff = alt3f
                            alt3ff_trim = alt3f
                        if total4new == 0:
                            ref4f = alt4f = None
                            alt4ff = -1
                            alt4ff_trim = -1
                        else:
                            alt4ff = alt4f
                            alt4ff_trim = alt4f

                        beg1 = beg4 = beg2 = beg3 = 0

                        details1 = (total1, total4, total1new, total4new, ref1, ref4, alt1, alt4, ref1f, ref4f, alt1f, alt4f, na1, na4, lowq1, lowq4, beg1, beg4)
                        details2 = (total2, total3, total2new, total3new, ref2, ref3, alt2, alt3, ref2f, ref3f, alt2f, alt3f, na2, na3, lowq2, lowq3, beg2, beg3)

                        trimmed = False
                        contradictory = False
                        softclipped_mutation_allMates = False
                        softclipped_mutation_oneOfTwoMates = False
                        softclipped_mutation_oneOfTwoSSCS = False
                        softclipped_mutation_oneOfTwoSSCS_diffMates = False
                        softclipped_mutation_oneMate = False
                        softclipped_mutation_oneMateOneSSCS = False

                        trimmed_actual_high_tier = False

                        dist_start_read1 = dist_start_read2 = dist_start_read3 = dist_start_read4 = []
                        dist_end_read1 = dist_end_read2 = dist_end_read3 = dist_end_read4 = []
                        ratio_dist_start1 = ratio_dist_start2 = ratio_dist_start3 = ratio_dist_start4 = False
                        ratio_dist_end1 = ratio_dist_end2 = ratio_dist_end3 = ratio_dist_end4 = False

                        # mate 1 - SSCS ab
                        softclipped_idx1 = [True if re.search(r"^[0-9]+S", string) or re.search(r"S$", string) else False for string in cigars_dcs1]
                        safe_div_result = safe_div(sum(softclipped_idx1), float(len(softclipped_idx1)))
                        if (safe_div_result == None):
                            ratio1 = False
                        else:
                            ratio1 = safe_div_result >= threshold_reads
                        if any(ij is True for ij in softclipped_idx1):
                            softclipped_both_ends_idx1 = [True if (re.search(r"^[0-9]+S", string) and re.search(r"S$", string)) else False for string in cigars_dcs1]
                            softclipped_start1 = [int(string.split("S")[0]) if re.search(r"^[0-9]+S", string) else -1 for string in cigars_dcs1]
                            softclipped_end1 = [int(re.split("[A-Z]", str(string))[-2]) if re.search(r"S$", string) else -1 for string in cigars_dcs1]
                            dist_start_read1 = [(pos - soft) if soft != -1 else thr + 1000 for soft, pos in zip(softclipped_start1, pos_read1)]
                            dist_end_read1 = [(length_read - pos - soft) if soft != -1 else thr + 1000 for soft, pos, length_read in zip(softclipped_end1, pos_read1, end_read1)]
                            # if read at both ends softclipped --> select end with smallest distance between mut position and softclipping
                            if any(ij is True for ij in softclipped_both_ends_idx1):
                                for nr, indx in enumerate(softclipped_both_ends_idx1):
                                    if indx:
                                        if dist_start_read1[nr] <= dist_end_read1[nr]:
                                            dist_end_read1[nr] = thr + 1000  # use dist of start and set start to very large number
                                        else:
                                            dist_start_read1[nr] = thr + 1000  # use dist of end and set start to very large number
                            ratio_dist_start1 = safe_div(sum([True if x <= thr else False for x in dist_start_read1]), float(sum(softclipped_idx1))) >= threshold_reads
                            ratio_dist_end1 = safe_div(sum([True if x <= thr else False for x in dist_end_read1]), float(sum(softclipped_idx1))) >= threshold_reads

                        # mate 1 - SSCS ba
                        softclipped_idx4 = [True if re.search(r"^[0-9]+S", string) or re.search(r"S$", string) else False for string in cigars_dcs4]
                        safe_div_result = safe_div(sum(softclipped_idx4), float(len(softclipped_idx4)))
                        if (safe_div_result == None):
                            ratio4 = False
                        else:
                            ratio4 = safe_div_result >= threshold_reads
                        if any(ij is True for ij in softclipped_idx4):
                            softclipped_both_ends_idx4 = [True if (re.search(r"^[0-9]+S", string) and re.search(r"S$", string)) else False for string in cigars_dcs4]
                            softclipped_start4 = [int(string.split("S")[0]) if re.search(r"^[0-9]+S", string) else -1 for string in cigars_dcs4]
                            softclipped_end4 = [int(re.split("[A-Z]", str(string))[-2]) if re.search(r"S$", string) else -1 for string in cigars_dcs4]
                            dist_start_read4 = [(pos - soft) if soft != -1 else thr + 1000 for soft, pos in zip(softclipped_start4, pos_read4)]
                            dist_end_read4 = [(length_read - pos - soft) if soft != -1 else thr + 1000 for soft, pos, length_read in zip(softclipped_end4, pos_read4, end_read4)]
                            # if read at both ends softclipped --> select end with smallest distance between mut position and softclipping
                            if any(ij is True for ij in softclipped_both_ends_idx4):
                                for nr, indx in enumerate(softclipped_both_ends_idx4):
                                    if indx:
                                        if dist_start_read4[nr] <= dist_end_read4[nr]:
                                            dist_end_read4[nr] = thr + 1000  # use dist of start and set start to very large number
                                        else:
                                            dist_start_read4[nr] = thr + 1000  # use dist of end and set start to very large number
                            ratio_dist_start4 = safe_div(sum([True if x <= thr else False for x in dist_start_read4]), float(sum(softclipped_idx4))) >= threshold_reads
                            ratio_dist_end4 = safe_div(sum([True if x <= thr else False for x in dist_end_read4]), float(sum(softclipped_idx4))) >= threshold_reads

                        # mate 2 - SSCS ab
                        softclipped_idx2 = [True if re.search(r"^[0-9]+S", string) or re.search(r"S$", string) else False for string in cigars_dcs2]
                        safe_div_result = safe_div(sum(softclipped_idx2), float(len(softclipped_idx2)))
                        if (safe_div_result == None):
                            ratio2 = False
                        else:
                            ratio2 = safe_div_result >= threshold_reads
                        if any(ij is True for ij in softclipped_idx2):
                            softclipped_both_ends_idx2 = [True if (re.search(r"^[0-9]+S", string) and re.search(r"S$", string)) else False for string in cigars_dcs2]
                            softclipped_start2 = [int(string.split("S")[0]) if re.search(r"^[0-9]+S", string) else -1 for string in cigars_dcs2]
                            softclipped_end2 = [int(re.split("[A-Z]", str(string))[-2]) if re.search(r"S$", string) else -1 for string in cigars_dcs2]
                            dist_start_read2 = [(pos - soft) if soft != -1 else thr + 1000 for soft, pos in zip(softclipped_start2, pos_read2)]
                            dist_end_read2 = [(length_read - pos - soft) if soft != -1 else thr + 1000 for soft, pos, length_read in zip(softclipped_end2, pos_read2, end_read2)]
                            # if read at both ends softclipped --> select end with smallest distance between mut position and softclipping
                            if any(ij is True for ij in softclipped_both_ends_idx2):
                                for nr, indx in enumerate(softclipped_both_ends_idx2):
                                    if indx:
                                        if dist_start_read2[nr] <= dist_end_read2[nr]:
                                            dist_end_read2[nr] = thr + 1000  # use dist of start and set start to very large number
                                        else:
                                            dist_start_read2[nr] = thr + 1000  # use dist of end and set start to very large number
                            ratio_dist_start2 = safe_div(sum([True if x <= thr else False for x in dist_start_read2]), float(sum(softclipped_idx2))) >= threshold_reads
                            ratio_dist_end2 = safe_div(sum([True if x <= thr else False for x in dist_end_read2]), float(sum(softclipped_idx2))) >= threshold_reads

                        # mate 2 - SSCS ba
                        softclipped_idx3 = [True if re.search(r"^[0-9]+S", string) or re.search(r"S$", string) else False for string in cigars_dcs3]
                        safe_div_result = safe_div(sum(softclipped_idx3), float(len(softclipped_idx3)))
                        if (safe_div_result == None):
                            ratio3 = False
                        else:
                            ratio3 = safe_div_result >= threshold_reads
                        if any(ij is True for ij in softclipped_idx3):
                            softclipped_both_ends_idx3 = [True if (re.search(r"^[0-9]+S", string) and re.search(r"S$", string)) else False for string in cigars_dcs3]
                            softclipped_start3 = [int(string.split("S")[0]) if re.search(r"^[0-9]+S", string) else -1 for string in cigars_dcs3]
                            softclipped_end3 = [int(re.split("[A-Z]", str(string))[-2]) if re.search(r"S$", string) else -1 for string in cigars_dcs3]
                            dist_start_read3 = [(pos - soft) if soft != -1 else thr + 1000 for soft, pos in zip(softclipped_start3, pos_read3)]
                            dist_end_read3 = [(length_read - pos - soft) if soft != -1 else thr + 1000 for soft, pos, length_read in zip(softclipped_end3, pos_read3, end_read3)]
                            # if read at both ends softclipped --> select end with smallest distance between mut position and softclipping
                            if any(ij is True for ij in softclipped_both_ends_idx3):
                                for nr, indx in enumerate(softclipped_both_ends_idx3):
                                    if indx:
                                        if dist_start_read3[nr] <= dist_end_read3[nr]:
                                            dist_end_read3[nr] = thr + 1000  # use dist of start and set start to a larger number than thresh
                                        else:
                                            dist_start_read3[nr] = thr + 1000  # use dist of end and set start to very large number
                            ratio_dist_start3 = safe_div(sum([True if x <= thr else False for x in dist_start_read3]), float(sum(softclipped_idx3))) >= threshold_reads
                            ratio_dist_end3 = safe_div(sum([True if x <= thr else False for x in dist_end_read3]), float(sum(softclipped_idx3))) >= threshold_reads

                        if ((all(float(ij) >= 0.5 for ij in [alt1ff, alt4ff]) &  # contradictory variant
                            all(float(ij) == 0. for ij in [alt2ff, alt3ff])) |
                            (all(float(ij) >= 0.5 for ij in [alt2ff, alt3ff]) &
                             all(float(ij) == 0. for ij in [alt1ff, alt4ff]))):
                            alt1ff = 0
                            alt4ff = 0
                            alt2ff = 0
                            alt3ff = 0
                            trimmed = False
                            contradictory = True
                        # softclipping tiers
                        # information of both mates available --> all reads for both mates and SSCS are softclipped
                        elif (ratio1 & ratio4 & ratio2 & ratio3 &
                              (ratio_dist_start1 | ratio_dist_end1) & (ratio_dist_start4 | ratio_dist_end4) & (ratio_dist_start2 | ratio_dist_end2) & (ratio_dist_start3 | ratio_dist_end3) &
                              all(float(ij) > 0. for ij in [alt1ff, alt2ff, alt3ff, alt4ff])):  # all mates available
                            # if distance between softclipping and mutation is at start or end of the read smaller than threshold
                            softclipped_mutation_allMates = True
                            softclipped_mutation_oneOfTwoMates = False
                            softclipped_mutation_oneOfTwoSSCS = False
                            softclipped_mutation_oneOfTwoSSCS_diffMates = False
                            softclipped_mutation_oneMate = False
                            softclipped_mutation_oneMateOneSSCS = False
                            alt1ff = 0
                            alt4ff = 0
                            alt2ff = 0
                            alt3ff = 0
                            trimmed = False
                            contradictory = False
                        # information of both mates available --> only one mate softclipped
                        elif (((ratio1 & ratio4 & (ratio_dist_start1 | ratio_dist_end1) & (ratio_dist_start4 | ratio_dist_end4)) |
                               (ratio2 & ratio3 & (ratio_dist_start2 | ratio_dist_end2) & (ratio_dist_start3 | ratio_dist_end3))) &
                              all(float(ij) > 0. for ij in [alt1ff, alt2ff, alt3ff, alt4ff])):  # all mates available
                            # if distance between softclipping and mutation is at start or end of the read smaller than threshold
                            min_start1 = min(min([ij[0] for ij in ref_positions1]), min([ij[0] for ij in ref_positions4]))  # red
                            min_start2 = min(min([ij[0] for ij in ref_positions2]), min([ij[0] for ij in ref_positions3]))  # blue
                            max_end1 = max(max([ij[1] for ij in ref_positions1]), max([ij[1] for ij in ref_positions4]))  # red
                            max_end2 = max(max([ij[1] for ij in ref_positions2]), max([ij[1] for ij in ref_positions3]))  # blue
                            if (min_start1 > min_start2) or (max_end1 > max_end2):  # if mate1 is red and mate2 is blue
                                softclipped_mutation_oneOfTwoMates = False
                                # blue mate at beginning softclipped
                                if min_start1 > min_start2:
                                    n_spacer_barcode = min_start1 - min_start2
                                    read_pos2 = read_pos2 - n_spacer_barcode
                                    read_pos3 = read_pos3 - n_spacer_barcode
                                    read_len_median2 = read_len_median2 - n_spacer_barcode
                                    read_len_median3 = read_len_median3 - n_spacer_barcode
                                # red mate at end softclipped
                                if max_end1 > max_end2:
                                    n_spacer_barcode = max_end1 - max_end2
                                    read_len_median1 = read_len_median1 - n_spacer_barcode
                                    read_len_median4 = read_len_median4 - n_spacer_barcode
                            elif (min_start1 < min_start2) or (max_end1 < max_end2):  # if mate1 is blue and mate2 is red
                                softclipped_mutation_oneOfTwoMates = False
                                if min_start1 < min_start2:
                                    n_spacer_barcode = min_start2 - min_start1
                                    read_pos1 = read_pos1 - n_spacer_barcode
                                    read_pos4 = read_pos4 - n_spacer_barcode
                                    read_len_median1 = read_len_median1 - n_spacer_barcode
                                    read_len_median4 = read_len_median4 - n_spacer_barcode
                                if max_end1 < max_end2:  # if mate1 ends after mate 2 starts
                                    n_spacer_barcode = max_end2 - max_end1
                                    read_len_median2 = read_len_median2 - n_spacer_barcode
                                    read_len_median3 = read_len_median3 - n_spacer_barcode
                            else:
                                softclipped_mutation_oneOfTwoMates = True
                                alt1ff = 0
                                alt4ff = 0
                                alt2ff = 0
                                alt3ff = 0
                                trimmed = False
                                contradictory = False
                            softclipped_mutation_allMates = False
                            softclipped_mutation_oneOfTwoSSCS = False
                            softclipped_mutation_oneMate = False
                            softclipped_mutation_oneMateOneSSCS = False

                            if softclipped_mutation_oneOfTwoMates is False:  # check trimming tier
                                if ((read_pos1 >= 0) and ((read_pos1 <= trim) | (abs(read_len_median1 - read_pos1) <= trim))):
                                    beg1 = total1new
                                    total1new = 0
                                    alt1ff = 0
                                    alt1f = 0
                                    trimmed = True

                                if ((read_pos4 >= 0) and ((read_pos4 <= trim) | (abs(read_len_median4 - read_pos4) <= trim))):
                                    beg4 = total4new
                                    total4new = 0
                                    alt4ff = 0
                                    alt4f = 0
                                    trimmed = True

                                if ((read_pos2 >= 0) and ((read_pos2 <= trim) | (abs(read_len_median2 - read_pos2) <= trim))):
                                    beg2 = total2new
                                    total2new = 0
                                    alt2ff = 0
                                    alt2f = 0
                                    trimmed = True

                                if ((read_pos3 >= 0) and ((read_pos3 <= trim) | (abs(read_len_median3 - read_pos3) <= trim))):
                                    beg3 = total3new
                                    total3new = 0
                                    alt3ff = 0
                                    alt3f = 0
                                    trimmed = True
                                details1 = (total1, total4, total1new, total4new, ref1, ref4, alt1, alt4, ref1f, ref4f, alt1f, alt4f, na1, na4, lowq1, lowq4, beg1, beg4)
                                details2 = (total2, total3, total2new, total3new, ref2, ref3, alt2, alt3, ref2f, ref3f, alt2f, alt3f, na2, na3, lowq2, lowq3, beg2, beg3)

                        # information of both mates available --> only one mate softclipped
                        elif (((ratio1 & (ratio_dist_start1 | ratio_dist_end1)) | (ratio4 & (ratio_dist_start4 | ratio_dist_end4))) &
                              ((ratio2 & (ratio_dist_start2 | ratio_dist_end2)) | (ratio3 & (ratio_dist_start3 | ratio_dist_end3))) &
                              all(float(ij) > 0. for ij in [alt1ff, alt2ff, alt3ff, alt4ff])):  # all mates available
                            # if distance between softclipping and mutation is at start or end of the read smaller than threshold
                            softclipped_mutation_allMates = False
                            softclipped_mutation_oneOfTwoMates = False
                            softclipped_mutation_oneOfTwoSSCS = True
                            softclipped_mutation_oneOfTwoSSCS_diffMates = False
                            softclipped_mutation_oneMate = False
                            softclipped_mutation_oneMateOneSSCS = False
                            alt1ff = 0
                            alt4ff = 0
                            alt2ff = 0
                            alt3ff = 0
                            trimmed = False
                            contradictory = False

                        # information of one mate available --> all reads of one mate are softclipped
                        elif ((ratio1 & ratio4 & (ratio_dist_start1 | ratio_dist_end1) & (ratio_dist_start4 | ratio_dist_end4) &
                              all(float(ij) < 0. for ij in [alt2ff, alt3ff]) & all(float(ij) > 0. for ij in [alt1ff, alt4ff])) |
                              (ratio2 & ratio3 & (ratio_dist_start2 | ratio_dist_end2) & (ratio_dist_start3 | ratio_dist_end3) &
                              all(float(ij) < 0. for ij in [alt1ff, alt4ff]) & all(float(ij) > 0. for ij in [alt2ff, alt3ff]))):  # all mates available
                            # if distance between softclipping and mutation is at start or end of the read smaller than threshold
                            softclipped_mutation_allMates = False
                            softclipped_mutation_oneOfTwoMates = False
                            softclipped_mutation_oneOfTwoSSCS = False
                            softclipped_mutation_oneOfTwoSSCS_diffMates = False
                            softclipped_mutation_oneMate = True
                            softclipped_mutation_oneMateOneSSCS = False
                            alt1ff = 0
                            alt4ff = 0
                            alt2ff = 0
                            alt3ff = 0
                            trimmed = False
                            contradictory = False
                        # information of one mate available --> only one SSCS is softclipped
                        elif ((((ratio1 & (ratio_dist_start1 | ratio_dist_end1)) | (ratio4 & (ratio_dist_start4 | ratio_dist_end4))) &
                              (all(float(ij) < 0. for ij in [alt2ff, alt3ff]) & all(float(ij) > 0. for ij in [alt1ff, alt4ff]))) |
                              (((ratio2 & (ratio_dist_start2 | ratio_dist_end2)) | (ratio3 & (ratio_dist_start3 | ratio_dist_end3))) &
                              (all(float(ij) < 0. for ij in [alt1ff, alt4ff]) & all(float(ij) < 0. for ij in [alt2ff, alt3ff])))):  # all mates available
                            # if distance between softclipping and mutation is at start or end of the read smaller than threshold
                            softclipped_mutation_allMates = False
                            softclipped_mutation_oneOfTwoMates = False
                            softclipped_mutation_oneOfTwoSSCS = False
                            softclipped_mutation_oneOfTwoSSCS_diffMates = False
                            softclipped_mutation_oneMate = False
                            softclipped_mutation_oneMateOneSSCS = True
                            alt1ff = 0
                            alt4ff = 0
                            alt2ff = 0
                            alt3ff = 0
                            trimmed = False
                            contradictory = False

                        else:
                            if ((read_pos1 >= 0) and ((read_pos1 <= trim) | (abs(read_len_median1 - read_pos1) <= trim))):
                                beg1 = total1new
                                total1new = 0
                                alt1ff = 0
                                alt1f = 0
                                trimmed = True

                            if ((read_pos4 >= 0) and ((read_pos4 <= trim) | (abs(read_len_median4 - read_pos4) <= trim))):
                                beg4 = total4new
                                total4new = 0
                                alt4ff = 0
                                alt4f = 0
                                trimmed = True

                            if ((read_pos2 >= 0) and ((read_pos2 <= trim) | (abs(read_len_median2 - read_pos2) <= trim))):
                                beg2 = total2new
                                total2new = 0
                                alt2ff = 0
                                alt2f = 0
                                trimmed = True

                            if ((read_pos3 >= 0) and ((read_pos3 <= trim) | (abs(read_len_median3 - read_pos3) <= trim))):
                                beg3 = total3new
                                total3new = 0
                                alt3ff = 0
                                alt3f = 0
                                trimmed = True
                            details1 = (total1, total4, total1new, total4new, ref1, ref4, alt1, alt4, ref1f, ref4f, alt1f, alt4f, na1, na4, lowq1, lowq4, beg1, beg4)
                            details2 = (total2, total3, total2new, total3new, ref2, ref3, alt2, alt3, ref2f, ref3f, alt2f, alt3f, na2, na3, lowq2, lowq3, beg2, beg3)

                        # assign tiers
                        if ((all(int(ij) >= 3 for ij in [total1new, total4new]) &
                             all(float(ij) >= 0.75 for ij in [alt1ff, alt4ff])) |
                            (all(int(ij) >= 3 for ij in [total2new, total3new]) &
                             all(float(ij) >= 0.75 for ij in [alt2ff, alt3ff]))):
                            tier = "1.1"
                            counter_tier11 += 1
                            tier_dict[key1]["tier 1.1"] += 1

                        elif (all(int(ij) >= 1 for ij in [total1new, total2new, total3new, total4new]) &
                              any(int(ij) >= 3 for ij in [total1new, total4new]) &
                              any(int(ij) >= 3 for ij in [total2new, total3new]) &
                              all(float(ij) >= 0.75 for ij in [alt1ff, alt2ff, alt3ff, alt4ff])):
                            tier = "1.2"
                            counter_tier12 += 1
                            tier_dict[key1]["tier 1.2"] += 1

                        elif ((all(int(ij) >= 1 for ij in [total1new, total4new]) &
                               any(int(ij) >= 3 for ij in [total1new, total4new]) &
                               all(float(ij) >= 0.75 for ij in [alt1ff, alt4ff])) |
                              (all(int(ij) >= 1 for ij in [total2new, total3new]) &
                               any(int(ij) >= 3 for ij in [total2new, total3new]) &
                               all(float(ij) >= 0.75 for ij in [alt2ff, alt3ff]))):
                            tier = "2.1"
                            counter_tier21 += 1
                            tier_dict[key1]["tier 2.1"] += 1

                        elif (all(int(ij) >= 1 for ij in [total1new, total2new, total3new, total4new]) &
                              all(float(ij) >= 0.75 for ij in [alt1ff, alt2ff, alt3ff, alt4ff])):
                            tier = "2.2"
                            counter_tier22 += 1
                            tier_dict[key1]["tier 2.2"] += 1

                        elif ((all(int(ij) >= 1 for ij in [total1new, total4new]) &
                               any(int(ij) >= 3 for ij in [total2new, total3new]) &
                               all(float(ij) >= 0.75 for ij in [alt1ff, alt4ff]) &
                               any(float(ij) >= 0.75 for ij in [alt2ff, alt3ff])) |
                              (all(int(ij) >= 1 for ij in [total2new, total3new]) &
                               any(int(ij) >= 3 for ij in [total1new, total4new]) &
                               all(float(ij) >= 0.75 for ij in [alt2ff, alt3ff]) &
                               any(float(ij) >= 0.75 for ij in [alt1ff, alt4ff]))):
                            tier = "2.3"
                            counter_tier23 += 1
                            tier_dict[key1]["tier 2.3"] += 1

                        elif ((all(int(ij) >= 1 for ij in [total1new, total4new]) &
                               all(float(ij) >= 0.75 for ij in [alt1ff, alt4ff])) |
                              (all(int(ij) >= 1 for ij in [total2new, total3new]) &
                               all(float(ij) >= 0.75 for ij in [alt2ff, alt3ff]))):
                            tier = "2.4"
                            counter_tier24 += 1
                            tier_dict[key1]["tier 2.4"] += 1

                        elif ((len(pure_tags_dict_short[key1]) > 1) &
                              (all(float(ij) >= 0.5 for ij in [alt1ff, alt4ff]) |
                               all(float(ij) >= 0.5 for ij in [alt2ff, alt3ff]))):
                            tier = "3.1"
                            counter_tier31 += 1
                            tier_dict[key1]["tier 3.1"] += 1

                        elif ((all(int(ij) >= 1 for ij in [total1new, total4new]) &
                               all(float(ij) >= 0.5 for ij in [alt1ff, alt4ff])) |
                              (all(int(ij) >= 1 for ij in [total2new, total3new]) &
                               all(float(ij) >= 0.5 for ij in [alt2ff, alt3ff]))):
                            tier = "3.2"
                            counter_tier32 += 1
                            tier_dict[key1]["tier 3.2"] += 1

                        elif (trimmed):
                            tier = "4"
                            counter_tier4 += 1
                            tier_dict[key1]["tier 4"] += 1

                            # assign tiers
                            if ((all(int(ij) >= 3 for ij in [total1new_trim, total4new_trim]) &
                                 all(float(ij) >= 0.75 for ij in [alt1ff_trim, alt4ff_trim])) |
                                (all(int(ij) >= 3 for ij in [total2new_trim, total3new_trim]) &
                                 all(float(ij) >= 0.75 for ij in [alt2ff_trim, alt3ff_trim]))):
                                trimmed_actual_high_tier = True
                            elif (all(int(ij) >= 1 for ij in [total1new_trim, total2new_trim, total3new_trim, total4new_trim]) &
                                  any(int(ij) >= 3 for ij in [total1new_trim, total4new_trim]) &
                                  any(int(ij) >= 3 for ij in [total2new_trim, total3new_trim]) &
                                  all(float(ij) >= 0.75 for ij in [alt1ff_trim, alt2ff_trim, alt3ff_trim, alt4ff_trim])):
                                trimmed_actual_high_tier = True
                            elif ((all(int(ij) >= 1 for ij in [total1new_trim, total4new_trim]) &
                                   any(int(ij) >= 3 for ij in [total1new_trim, total4new_trim]) &
                                   all(float(ij) >= 0.75 for ij in [alt1ff_trim, alt4ff_trim])) |
                                  (all(int(ij) >= 1 for ij in [total2new_trim, total3new_trim]) &
                                   any(int(ij) >= 3 for ij in [total2new_trim, total3new_trim]) &
                                   all(float(ij) >= 0.75 for ij in [alt2ff_trim, alt3ff_trim]))):
                                trimmed_actual_high_tier = True
                            elif (all(int(ij) >= 1 for ij in [total1new_trim, total2new_trim, total3new_trim, total4new_trim]) &
                                  all(float(ij) >= 0.75 for ij in [alt1ff_trim, alt2ff_trim, alt3ff_trim, alt4ff_trim])):
                                trimmed_actual_high_tier = True
                            elif ((all(int(ij) >= 1 for ij in [total1new_trim, total4new_trim]) &
                                   any(int(ij) >= 3 for ij in [total2new_trim, total3new_trim]) &
                                   all(float(ij) >= 0.75 for ij in [alt1ff_trim, alt4ff_trim]) &
                                   any(float(ij) >= 0.75 for ij in [alt2ff_trim, alt3ff_trim])) |
                                  (all(int(ij) >= 1 for ij in [total2new_trim, total3new_trim]) &
                                   any(int(ij) >= 3 for ij in [total1new_trim, total4new_trim]) &
                                   all(float(ij) >= 0.75 for ij in [alt2ff_trim, alt3ff_trim]) &
                                   any(float(ij) >= 0.75 for ij in [alt1ff_trim, alt4ff_trim]))):
                                trimmed_actual_high_tier = True
                            elif ((all(int(ij) >= 1 for ij in [total1new_trim, total4new_trim]) &
                                   all(float(ij) >= 0.75 for ij in [alt1ff_trim, alt4ff_trim])) |
                                  (all(int(ij) >= 1 for ij in [total2new_trim, total3new_trim]) &
                                   all(float(ij) >= 0.75 for ij in [alt2ff_trim, alt3ff_trim]))):
                                trimmed_actual_high_tier = True
                            else:
                                trimmed_actual_high_tier = False

                        elif softclipped_mutation_allMates:
                            tier = "5.1"
                            counter_tier51 += 1
                            tier_dict[key1]["tier 5.1"] += 1

                        elif softclipped_mutation_oneOfTwoMates:
                            tier = "5.2"
                            counter_tier52 += 1
                            tier_dict[key1]["tier 5.2"] += 1

                        elif softclipped_mutation_oneOfTwoSSCS:
                            tier = "5.3"
                            counter_tier53 += 1
                            tier_dict[key1]["tier 5.3"] += 1

                        elif softclipped_mutation_oneMate:
                            tier = "5.4"
                            counter_tier54 += 1
                            tier_dict[key1]["tier 5.4"] += 1

                        elif softclipped_mutation_oneMateOneSSCS:
                            tier = "5.5"
                            counter_tier55 += 1
                            tier_dict[key1]["tier 5.5"] += 1

                        elif (contradictory):
                            tier = "6"
                            counter_tier6 += 1
                            tier_dict[key1]["tier 6"] += 1

                        else:
                            tier = "7"
                            counter_tier7 += 1
                            tier_dict[key1]["tier 7"] += 1

                        chrom, pos, ref_a, alt_a = re.split(r'\#', key1)
                        var_id = '-'.join([chrom, str(int(pos)+1), ref, alt])
                        sample_tag = key2[:-5]
                        array2 = np.unique(whole_array)  # remove duplicate sequences to decrease running time
                        # exclude identical tag from array2, to prevent comparison to itself
                        same_tag = np.where(array2 == sample_tag)
                        index_array2 = np.arange(0, len(array2), 1)
                        index_withoutSame = np.delete(index_array2, same_tag)  # delete identical tag from the data
                        array2 = array2[index_withoutSame]
                        if len(array2) != 0:  # only perform chimera analysis if there is more than 1 variant
                            array1_half = sample_tag[0:int(len(sample_tag) / 2)]  # mate1 part1
                            array1_half2 = sample_tag[int(len(sample_tag) / 2):int(len(sample_tag))]  # mate1 part 2
                            array2_half = np.array([ii[0:int(len(ii) / 2)] for ii in array2])  # mate2 part1
                            array2_half2 = np.array([ii[int(len(ii) / 2):int(len(ii))] for ii in array2])  # mate2 part2

                            min_tags_list_zeros = []
                            chimera_tags = []
                            for mate_b in [False, True]:
                                i = 0  # counter, only used to see how many HDs of tags were already calculated
                                if mate_b is False:  # HD calculation for all a's
                                    half1_mate1 = array1_half
                                    half2_mate1 = array1_half2
                                    half1_mate2 = array2_half
                                    half2_mate2 = array2_half2
                                elif mate_b is True:  # HD calculation for all b's
                                    half1_mate1 = array1_half2
                                    half2_mate1 = array1_half
                                    half1_mate2 = array2_half2
                                    half2_mate2 = array2_half
                                # calculate HD of "a" in the tag to all "a's" or "b" in the tag to all "b's"
                                dist = np.array([sum(itertools.imap(operator.ne, half1_mate1, c)) for c in half1_mate2])
                                min_index = np.where(dist == dist.min())  # get index of min HD
                                # get all "b's" of the tag or all "a's" of the tag with minimum HD
                                min_tag_half2 = half2_mate2[min_index]
                                min_tag_array2 = array2[min_index]  # get whole tag with min HD
                                min_value = dist.min()
                                # calculate HD of "b" to all "b's" or "a" to all "a's"
                                dist_second_half = np.array([sum(itertools.imap(operator.ne, half2_mate1, e))
                                                             for e in min_tag_half2])
                                dist2 = dist_second_half.max()
                                max_index = np.where(dist_second_half == dist_second_half.max())[0]  # get index of max HD
                                max_tag = min_tag_array2[max_index]

                                # tags which have identical parts:
                                if min_value == 0 or dist2 == 0:
                                    min_tags_list_zeros.append(tag)
                                    chimera_tags.append(max_tag)
                                i += 1
                            chimera_tags = [x for x in chimera_tags if x != []]
                            chimera_tags_new = []
                            for i in chimera_tags:
                                if len(i) > 1:
                                    for t in i:
                                        chimera_tags_new.append(t)
                                else:
                                    chimera_tags_new.extend(i)
                            chimera = ", ".join(chimera_tags_new)
                        else:
                            chimera_tags_new = []
                            chimera = ""

                        if len(chimera_tags_new) > 0:
                            chimera_tags_new.append(sample_tag)
                            key_chimera = ",".join(sorted(chimera_tags_new))
                            if key_chimera in chimeric_tag.keys():
                                chimeric_tag[key_chimera].append(float(tier))
                            else:
                                chimeric_tag[key_chimera] = [float(tier)]

                        if (read_pos1 == -1):
                            read_pos1 = read_len_median1 = None
                        if (read_pos4 == -1):
                            read_pos4 = read_len_median4 = None
                        if (read_pos2 == -1):
                            read_pos2 = read_len_median2 = None
                        if (read_pos3 == -1):
                            read_pos3 = read_len_median3 = None
                        line = (var_id, tier, key2[:-5], 'ab1.ba2', read_pos1, read_pos4, read_len_median1, read_len_median4, dcs_median) + details1 + (sscs_mut_ab, sscs_mut_ba, sscs_ref_ab, sscs_ref_ba, add_mut14, chimera)
                        # ws1.write_row(row, 0, line)
                        # csv_writer.writerow(line)
                        line2 = ("", "", key2[:-5], 'ab2.ba1', read_pos2, read_pos3, read_len_median2, read_len_median3, dcs_median) + details2 + (sscs_mut_ab, sscs_mut_ba, sscs_ref_ab, sscs_ref_ba, add_mut23, chimera)
                        # ws1.write_row(row + 1, 0, line2)
                        # csv_writer.writerow(line2)

                        # ws1.conditional_format('L{}:M{}'.format(row + 1, row + 2),
                        #                       {'type': 'formula',
                        #                        'criteria': '=OR($B${}="1.1", $B${}="1.2")'.format(row + 1, row + 1),
                        #                        'format': format1,
                        #                        'multi_range': 'L{}:M{} T{}:U{} B{}'.format(row + 1, row + 2, row + 1, row + 2, row + 1, row + 2)})
                        # ws1.conditional_format('L{}:M{}'.format(row + 1, row + 2),
                        #                       {'type': 'formula',
                        #                        'criteria': '=OR($B${}="2.1", $B${}="2.2", $B${}="2.3", $B${}="2.4", $B${}="2.5")'.format(row + 1, row + 1, row + 1, row + 1, row + 1),
                        #                        'format': format3,
                        #                        'multi_range': 'L{}:M{} T{}:U{} B{}'.format(row + 1, row + 2, row + 1, row + 2, row + 1, row + 2)})
                        # ws1.conditional_format('L{}:M{}'.format(row + 1, row + 2),
                        #                       {'type': 'formula',
                        #                        'criteria': '=$B${}>="3"'.format(row + 1),
                        #                        'format': format2,
                        #                        'multi_range': 'L{}:M{} T{}:U{} B{}'.format(row + 1, row + 2, row + 1, row + 2, row + 1, row + 2)})
                        change_tier_after_print.append((row, line, line2, trimmed_actual_high_tier))
                        row += 3

            if chimera_correction:
                chimeric_dcs_high_tiers = 0
                chimeric_dcs = 0
                for keys_chimera in chimeric_tag.keys():
                    tiers = chimeric_tag[keys_chimera]
                    chimeric_dcs += len(tiers) - 1
                    high_tiers = sum(1 for t in tiers if t < 3.)
                    if high_tiers == len(tiers):
                        chimeric_dcs_high_tiers += high_tiers - 1
                    else:
                        chimeric_dcs_high_tiers += high_tiers
                chimera_dict[key1] = (chimeric_dcs, chimeric_dcs_high_tiers)

            # write to file
            # move tier 4 counts to tier 2.5 if there other mutations with tier <= 2.4
            sum_highTiers = sum([tier_dict[key1][ij] for ij in list(sorted(tier_dict[key1].keys()))[:6]])
            correct_tier = False
            if tier_dict[key1]["tier 4"] > 0 and sum_highTiers > 0:
                tier_dict[key1]["tier 2.5"] = tier_dict[key1]["tier 4"]
                tier_dict[key1]["tier 4"] = 0
                correct_tier = True

            for sample in change_tier_after_print:
                row_number = sample[0]
                line1 = sample[1]
                line2 = sample[2]
                actual_high_tier = sample[3]
                current_tier = list(line1)[1]

                if correct_tier and (current_tier == "4") and actual_high_tier:
                    line1 = list(line1)
                    line1[1] = "2.5"
                    line1 = tuple(line1)
                    counter_tier25 += 1
                    counter_tier4 -= 1
                ws1.write_row(row_number, 0, line1)
                csv_writer.writerow(line1)
                ws1.write_row(row_number + 1, 0, line2)
                csv_writer.writerow(line2)

                ws1.conditional_format('L{}:M{}'.format(row_number + 1, row_number + 2),
                                       {'type': 'formula',
                                        'criteria': '=OR($B${}="1.1", $B${}="1.2")'.format(row_number + 1, row_number + 1),
                                        'format': format1,
                                        'multi_range': 'L{}:M{} T{}:U{} B{}'.format(row_number + 1, row_number + 2, row_number + 1, row_number + 2, row_number + 1, row_number + 2)})
                ws1.conditional_format('L{}:M{}'.format(row_number + 1, row_number + 2),
                                       {'type': 'formula',
                                        'criteria': '=OR($B${}="2.1", $B${}="2.2", $B${}="2.3", $B${}="2.4", $B${}="2.5")'.format(row_number + 1, row_number + 1, row_number + 1, row_number + 1, row_number + 1),
                                        'format': format3,
                                        'multi_range': 'L{}:M{} T{}:U{} B{}'.format(row_number + 1, row_number + 2, row_number + 1, row_number + 2, row_number + 1, row_number + 2)})
                ws1.conditional_format('L{}:M{}'.format(row_number + 1, row_number + 2),
                                       {'type': 'formula',
                                        'criteria': '=$B${}>="3"'.format(row_number + 1),
                                        'format': format2,
                                        'multi_range': 'L{}:M{} T{}:U{} B{}'.format(row_number + 1, row_number + 2, row_number + 1, row_number + 2, row_number + 1, row_number + 2)})
    # sheet 2
    if chimera_correction:
        header_line2 = ('variant ID', 'cvrg', 'AC alt (all tiers)', 'AF (all tiers)', 'cvrg (tiers 1.1-2.5)', 'AC alt (tiers 1.1-2.5)', 'AF (tiers 1.1-2.5)', 'chimera-corrected cvrg (tiers 1.1-2.5)', 'chimeras in AC alt (tiers 1.1-2.5)', 'chimera-corrected AF (tiers 1.1-2.5)', 'AC alt (orginal DCS)', 'AF (original DCS)',
                        'tier 1.1', 'tier 1.2', 'tier 2.1', 'tier 2.2', 'tier 2.3', 'tier 2.4', 'tier 2.5',
                        'tier 3.1', 'tier 3.2', 'tier 4', 'tier 5.1', 'tier 5.2', 'tier 5.3', 'tier 5.4', 'tier 5.5', 'tier 6', 'tier 7', 'AF 1.1-1.2', 'AF 1.1-2.1', 'AF 1.1-2.2',
                        'AF 1.1-2.3', 'AF 1.1-2.4', 'AF 1.1-2.5', 'AF 1.1-3.1', 'AF 1.1-3.2', 'AF 1.1-4', 'AF 1.1-5.1', 'AF 1.1-5.2', 'AF 1.1-5.3', 'AF 1.1-5.4', 'AF 1.1-5.5', 'AF 1.1-6', 'AF 1.1-7')
    else:
        header_line2 = ('variant ID', 'cvrg', 'AC alt (all tiers)', 'AF (all tiers)', 'cvrg (tiers 1.1-2.5)', 'AC alt (tiers 1.1-2.5)', 'AF (tiers 1.1-2.5)', 'AC alt (orginal DCS)', 'AF (original DCS)',
                        'tier 1.1', 'tier 1.2', 'tier 2.1', 'tier 2.2', 'tier 2.3', 'tier 2.4', 'tier 2.5',
                        'tier 3.1', 'tier 3.2', 'tier 4', 'tier 5.1', 'tier 5.2', 'tier 5.3', 'tier 5.4', 'tier 5.5', 'tier 6', 'tier 7', 'AF 1.1-1.2', 'AF 1.1-2.1', 'AF 1.1-2.2',
                        'AF 1.1-2.3', 'AF 1.1-2.4', 'AF 1.1-2.5', 'AF 1.1-3.1', 'AF 1.1-3.2', 'AF 1.1-4', 'AF 1.1-5.1', 'AF 1.1-5.2', 'AF 1.1-5.3', 'AF 1.1-5.4', 'AF 1.1-5.5', 'AF 1.1-6', 'AF 1.1-7')

    ws2.write_row(0, 0, header_line2)
    row = 0

    for key1, value1 in sorted(tier_dict.items()):
        if key1 in pure_tags_dict_short.keys():
            i = np.where(np.array(['#'.join(str(i) for i in z)
                                   for z in zip(mut_array[:, 0], mut_array[:, 1], mut_array[:, 2], mut_array[:, 3])]) == key1)[0][0]
            ref = mut_array[i, 2]
            alt = mut_array[i, 3]
            chrom, pos, ref_a, alt_a = re.split(r'\#', key1)
            ref_count = cvrg_dict[key1][0]
            alt_count = cvrg_dict[key1][1]
            cvrg = ref_count + alt_count

            var_id = '-'.join([chrom, str(int(pos)+1), ref, alt])
            lst = [var_id, cvrg]
            used_tiers = []
            cum_af = []
            for key2, value2 in sorted(value1.items()):
                # calculate cummulative AF
                used_tiers.append(value2)
                if len(used_tiers) > 1:
                    cum = safe_div(sum(used_tiers), cvrg)
                    cum_af.append(cum)
            if sum(used_tiers) == 0:  # skip mutations that are filtered by the VA in the first place
                continue
            lst.extend([sum(used_tiers), safe_div(sum(used_tiers), cvrg)])
            lst.extend([(cvrg - sum(used_tiers[-10:])), sum(used_tiers[0:7]), safe_div(sum(used_tiers[0:7]), (cvrg - sum(used_tiers[-10:])))])
            if chimera_correction:
                chimeras_all = chimera_dict[key1][1]
                new_alt = sum(used_tiers[0:7]) - chimeras_all
                fraction_chimeras = safe_div(chimeras_all, float(sum(used_tiers[0:7])))
                if fraction_chimeras is None:
                    fraction_chimeras = 0.
                new_cvrg = (cvrg - sum(used_tiers[-10:])) * (1. - fraction_chimeras)
                lst.extend([new_cvrg, chimeras_all, safe_div(new_alt, new_cvrg)])
            lst.extend([alt_count, safe_div(alt_count, cvrg)])
            lst.extend(used_tiers)
            lst.extend(cum_af)
            lst = tuple(lst)
            ws2.write_row(row + 1, 0, lst)
            if chimera_correction:
                ws2.conditional_format('M{}:N{}'.format(row + 2, row + 2), {'type': 'formula', 'criteria': '=$M$1="tier 1.1"', 'format': format12, 'multi_range': 'M{}:N{} M1:N1'.format(row + 2, row + 2)})
                ws2.conditional_format('O{}:S{}'.format(row + 2, row + 2), {'type': 'formula', 'criteria': '=$O$1="tier 2.1"', 'format': format32, 'multi_range': 'O{}:S{} O1:S1'.format(row + 2, row + 2)})
                ws2.conditional_format('T{}:AC{}'.format(row + 2, row + 2), {'type': 'formula', 'criteria': '=$T$1="tier 3.1"', 'format': format22, 'multi_range': 'T{}:AC{} T1:AC1'.format(row + 2, row + 2)})
            else:
                ws2.conditional_format('J{}:K{}'.format(row + 2, row + 2), {'type': 'formula', 'criteria': '=$J$1="tier 1.1"', 'format': format12, 'multi_range': 'J{}:K{} J1:K1'.format(row + 2, row + 2)})
                ws2.conditional_format('L{}:P{}'.format(row + 2, row + 2), {'type': 'formula', 'criteria': '=$L$1="tier 2.1"', 'format': format32, 'multi_range': 'L{}:P{} L1:P1'.format(row + 2, row + 2)})
                ws2.conditional_format('Q{}:Z{}'.format(row + 2, row + 2), {'type': 'formula', 'criteria': '=$Q$1="tier 3.1"', 'format': format22, 'multi_range': 'Q{}:Z{} Q1:Z1'.format(row + 2, row + 2)})
            row += 1

    # sheet 3
    sheet3 = [("tier 1.1", counter_tier11), ("tier 1.2", counter_tier12), ("tier 2.1", counter_tier21),
              ("tier 2.2", counter_tier22), ("tier 2.3", counter_tier23), ("tier 2.4", counter_tier24), ("tier 2.5", counter_tier25),
              ("tier 3.1", counter_tier31), ("tier 3.2", counter_tier32), ("tier 4", counter_tier4),
              ("tier 5.1", counter_tier51), ("tier 5.2", counter_tier52),
              ("tier 5.3", counter_tier53), ("tier 5.4", counter_tier54), ("tier 5.5", counter_tier55), ("tier 6", counter_tier6), ("tier 7", counter_tier7)]

    header = ("tier", "count")
    ws3.write_row(0, 0, header)

    for i in range(len(sheet3)):
        ws3.write_row(i + 1, 0, sheet3[i])
        ws3.conditional_format('A{}:B{}'.format(i + 2, i + 2),
                               {'type': 'formula',
                                'criteria': '=OR($A${}="tier 1.1", $A${}="tier 1.2")'.format(i + 2, i + 2),
                                'format': format1})
        ws3.conditional_format('A{}:B{}'.format(i + 2, i + 2),
                               {'type': 'formula',
                                'criteria': '=OR($A${}="tier 2.1", $A${}="tier 2.2", $A${}="tier 2.3", $A${}="tier 2.4", $A${}="tier 2.5")'.format(i + 2, i + 2, i + 2, i + 2, i + 2),
                                'format': format3})
        ws3.conditional_format('A{}:B{}'.format(i + 2, i + 2),
                               {'type': 'formula',
                                'criteria': '=$A${}>="3"'.format(i + 2),
                                'format': format2})

    description_tiers = [("Tier 1.1", "both ab and ba SSCS present (>75% of the sites with alternative base) and minimal FS>=3 for both SSCS in at least one mate"), ("", ""),
                         ("Tier 1.2", "both ab and ba SSCS present (>75% of the sites with alt. base) and mate pair validation (min. FS=1) and minimal FS>=3 for at least one of the SSCS"),
                         ("Tier 2.1", "both ab and ba SSCS present (>75% of the sites with alt. base) and minimal FS>=3 for at least one of the SSCS in at least one mate"),
                         ("Tier 2.2", "both ab and ba SSCS present (>75% of the sites with alt. base) and mate pair validation (min. FS=1)"),
                         ("Tier 2.3", "both ab and ba SSCS present (>75% of the sites with alt. base) and minimal FS=1 for both SSCS in one mate and minimal FS>=3 for at least one of the SSCS in the other mate"),
                         ("Tier 2.4", "both ab and ba SSCS present (>75% of the sites with alt. base) and minimal FS=1 for both SSCS in at least one mate"),
                         ("Tier 2.5", "variants at the start or end of the read (ignoring variant position tier 1.1-2.4) and recurring mutation on this position in tier 1.1-2.4"),
                         ("Tier 3.1", "both ab and ba SSCS present (>50% of the sites with alt. base) and recurring mutation on this position"),
                         ("Tier 3.2", "both ab and ba SSCS present (>50% of the sites with alt. base) and minimal FS>=1 for both SSCS in at least one mate"),
                         ("Tier 4", "variants at the start or end of the reads"),
                         ("Tier 5.1", "variant is close to softclipping in both mates and SSCS"),
                         ("Tier 5.2", "variant is close to softclipping in one of the mates but both SSCS"),
                         ("Tier 5.3", "variant is close to softclipping in one of the SSCS of both mates"),
                         ("Tier 5.4", "variant is close to softclipping in one mate and both SSCS (no information of second mate)"),
                         ("Tier 5.5", "variant is close to softclipping in one of the SSCS (no information of the second mate)"),
                         ("Tier 6", "mates with contradictory information"),
                         ("Tier 7", "remaining variants")]
    examples_tiers = [[("chr5-11068-C-G", "1.1", "AAAAAGATGCCGACTACCTT", "ab1.ba2", "254", "228", "287", "288", "289",
                        "3", "6", "3", "6", "0", "0", "3", "6", "0", "0", "1", "1", "0", "0", "0", "0", "0", "0",
                        "4081", "4098", "5", "10", "", ""),
                       ("", "", "AAAAAGATGCCGACTACCTT", "ab2.ba1", None, None, None, None,
                        "289", "0", "0", "0", "0", "0", "0", "0", "0", None, None, None, None,
                        "0", "0", "0", "0", "0", "0", "4081", "4098", "5", "10", "", "")],
                      [("chr5-11068-C-G", "1.1", "AAAAATGCGTAGAAATATGC", "ab1.ba2", "254", "228", "287", "288", "289",
                        "33", "43", "33", "43", "0", "0", "33", "43", "0", "0", "1", "1", "0", "0", "0", "0", "0",
                        "0", "4081", "4098", "5", "10", "", ""),
                       ("", "", "AAAAATGCGTAGAAATATGC", "ab2.ba1", "268", "268", "270", "288", "289",
                        "11", "34", "10", "27", "0", "0", "10", "27", "0", "0", "1", "1", "0", "0", "1",
                        "7", "0", "0", "4081", "4098", "5", "10", "", "")],
                      [("chr5-10776-G-T", "1.2", "CTATGACCCGTGAGCCCATG", "ab1.ba2", "132", "132", "287", "288", "290",
                        "4", "1", "4", "1", "0", "0", "4", "1", "0", "0", "1", "1", "0", "0", "0", "0",
                        "0", "0", "1", "6", "47170", "41149", "", ""),
                       ("", "", "CTATGACCCGTGAGCCCATG", "ab2.ba1", "77", "132", "233", "200", "290",
                        "4", "1", "4", "1", "0", "0", "4", "1", "0", "0", "1", "1", "0", "0", "0", "0",
                        "0", "0", "1", "6", "47170", "41149", "", "")],
                      [("chr5-11068-C-G", "2.1", "AAAAAAACATCATACACCCA", "ab1.ba2", "246", "244", "287", "288", "289",
                        "2", "8", "2", "8", "0", "0", "2", "8", "0", "0", "1", "1", "0", "0", "0", "0", "0", "0",
                        "4081", "4098", "5", "10", "", ""),
                       ("", "", "AAAAAAACATCATACACCCA", "ab2.ba1", None, None, None, None,
                        "289", "0", "0", "0", "0", "0", "0", "0", "0", None, None, None, None, "0", "0",
                        "0", "0", "0", "0", "4081", "4098", "5", "10", "", "")],
                      [("chr5-11068-C-G", "2.2", "ATCAGCCATGGCTATTATTG", "ab1.ba2", "72", "72", "217", "288", "289",
                        "1", "1", "1", "1", "0", "0", "1", "1", "0", "0", "1", "1", "0", "0", "0", "0", "0", "0",
                        "4081", "4098", "5", "10", "", ""),
                       ("", "", "ATCAGCCATGGCTATTATTG", "ab2.ba1", "153", "164", "217", "260", "289",
                        "1", "1", "1", "1", "0", "0", "1", "1", "0", "0", "1", "1", "0", "0", "0", "0", "0", "0",
                        "4081", "4098", "5", "10", "", "")],
                      [("chr5-11068-C-G", "2.3", "ATCAATATGGCCTCGCCACG", "ab1.ba2", None, None, None, None,
                        "289", "0", "5", "0", "5", "0", "0", "0", "5", None, None, None, "1", "0",
                        "0", "0", "0", "0", "0", "4081", "4098", "5", "10", "", ""),
                       ("", "", "ATCAATATGGCCTCGCCACG", "ab2.ba1", "202", "255", "277", "290", "289",
                        "1", "3", "1", "3", "0", "0", "1", "3", "0", "0", "1", "1", "0", "0", "0", "0",
                        "0", "0", "4081", "4098", "5", "10", "", "")],
                      [("chr5-11068-C-G", "2.4", "ATCAGCCATGGCTATTTTTT", "ab1.ba2", "72", "72", "217", "288", "289",
                        "1", "1", "1", "1", "0", "0", "1", "1", "0", "0", "1", "1", "0", "0", "0", "0", "0", "0", "4081",
                        "4098", "5", "10", "", ""),
                       ("", "", "ATCAGCCATGGCTATTTTTT", "ab2.ba1", "153", "164", "217", "260", "289",
                        "1", "1", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "1", "1", "0", "0", "0", "0", "4081",
                        "4098", "5", "10", "", "")],
                      [("chr5-11068-C-G", "2.5", "ATTGAAAGAATAACCCACAC", "ab1.ba2", "1", "100", "255", "276", "269",
                        "5", "6", "0", "6", "0", "0", "5", "6", "0", "0", "0", "1", "0", "0", "0", "0", "5", "0", "1", "1", "5348", "5350", "", ""),
                       ("", "", "AAAAAAAGAATAACCCACAC", "ab2.ba1", None, None, None, None,
                        "269", "0", "0", "0", "0", "0", "0", "0", "0", None, None, None, None, "0",
                        "0", "0", "0", "0", "0", "1", "1", "5348", "5350", "", "")],
                      [("chr5-10776-G-T", "3.1", "ATGCCTACCTCATTTGTCGT", "ab1.ba2", "46", "15", "287", "288", "290",
                        "3", "3", "3", "2", "3", "1", "0", "1", "1", "0.5", "0", "0.5", "0", "0", "0", "1",
                        "0", "0", "3", "3", "47170", "41149", "", ""),
                       ("", "", "ATGCCTACCTCATTTGTCGT", "ab2.ba1", None, "274", None,
                        "288", "290", "0", "3", "0", "2", "0", "1", "0", "1", None, "0.5", None, "0.5",
                        "0", "0", "0", "1", "0", "0", "3", "3", "47170", "41149", "", "")],
                      [("chr5-11315-C-T", "3.2", "ACAACATCACGTATTCAGGT", "ab1.ba2", "197", "197", "240", "255", "271",
                        "2", "3", "2", "3", "0", "1", "2", "2", "0", "0.333333333333333", "1",
                        "0.666666666666667", "0", "0", "0", "0", "0", "0", "1", "1", "6584", "6482", "", ""),
                       ("", "", "ACAACATCACGTATTCAGGT", "ab2.ba1", "35", "35", "240", "258", "271",
                        "2", "3", "2", "3", "0", "1", "2", "2", "0", "0.333333333333333", "1",
                        "0.666666666666667", "0", "0", "0", "0", "0", "0", "1", "1", "6584", "6482", "", "")],
                      [("chr5-13983-G-C", "4", "AAAAAAAGAATAACCCACAC", "ab1.ba2", "1", "100", "255", "276", "269",
                        "5", "6", "0", "6", "0", "0", "5", "6", "0", "0", "0", "1", "0", "0", "0", "0", "5", "0", "1", "1", "5348", "5350", "", ""),
                       ("", "", "AAAAAAAGAATAACCCACAC", "ab2.ba1", None, None, None, None,
                        "269", "0", "0", "0", "0", "0", "0", "0", "0", None, None, None, None, "0",
                        "0", "0", "0", "0", "0", "1", "1", "5348", "5350", "", "")],
                      [("" * 34), ("" * 34)], [("" * 34), ("" * 34)], [("" * 34), ("" * 34)], [("" * 34), ("" * 34)], [("" * 34), ("" * 34)],
                      [("chr5-13963-T-C", "6", "TTTTTAAGAATAACCCACAC", "ab1.ba2", "38", "38", "240", "283", "263",
                        "110", "54", "110", "54", "0", "0", "110", "54", "0", "0", "1", "1", "0", "0", "0",
                        "0", "0", "0", "1", "1", "5348", "5350", "", ""),
                       ("", "", "TTTTTAAGAATAACCCACAC", "ab2.ba1", "100", "112", "140", "145", "263",
                        "7", "12", "7", "12", "7", "12", "0", "0", "1", "1", "0",
                        "0", "0", "0", "0", "0", "0", "0", "1", "1", "5348", "5350", "", "")],
                      [("chr5-13983-G-C", "7", "ATGTTGTGAATAACCCACAC", "ab1.ba2", None, "186", None, "276", "269",
                        "0", "6", "0", "6", "0", "0", "0", "6", "0", "0", "0", "1", "0", "0", "0", "0", "0",
                        "0", "1", "1", "5348", "5350", "", ""),
                       ("", "", "ATGTTGTGAATAACCCACAC", "ab2.ba1", None, None, None, None,
                        "269", "0", "0", "0", "0", "0", "0", "0", "0", None, None, None, None, "0",
                        "0", "0", "0", "0", "0", "1", "1", "5348", "5350", "", "")]]

    start_row = 20
    ws3.write(start_row, 0, "Description of tiers with examples")
    ws3.write_row(start_row + 1, 0, header_line)
    row = 0
    for i in range(len(description_tiers)):
        ws3.write_row(start_row + 2 + row + i + 1, 0, description_tiers[i])
        ex = examples_tiers[i]
        for k in range(len(ex)):
            ws3.write_row(start_row + 2 + row + i + k + 2, 0, ex[k])
        ws3.conditional_format('L{}:M{}'.format(start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3), {'type': 'formula', 'criteria': '=OR($B${}="1.1", $B${}="1.2")'.format(start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 2), 'format': format13, 'multi_range': 'L{}:M{} T{}:U{} B{}'.format(start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3, start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3, start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3)})
        ws3.conditional_format('L{}:M{}'.format(start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3),
                               {'type': 'formula', 'criteria': '=OR($B${}="2.1",$B${}="2.2", $B${}="2.3", $B${}="2.4", $B${}="2.5")'.format(start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 2),
                                'format': format33,
                                'multi_range': 'L{}:M{} T{}:U{} B{}'.format(start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3, start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3, start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3)})
        ws3.conditional_format('L{}:M{}'.format(start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3),
                               {'type': 'formula',
                                'criteria': '=$B${}>="3"'.format(start_row + 2 + row + i + k + 2),
                                'format': format23,
                                'multi_range': 'L{}:M{} T{}:U{} B{}'.format(start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3, start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3, start_row + 2 + row + i + k + 2, start_row + 2 + row + i + k + 3)})
        row += 3
    workbook.close()
    workbook2.close()
    workbook3.close()
    csv_data.close()


if __name__ == '__main__':
    sys.exit(read2mut(sys.argv))
