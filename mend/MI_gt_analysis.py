#!/usr/bin/env python3
"""
Mendelian inheritance analysis for a single chrom

Usage: ./MI-anaylsis.py <pedigree> <chrom> <vcf>

Example:
./MI-analysis.py /gymreklab-tscc/helia/TR_1000G/1000G.ped chr21 /gymreklab-tscc/helia/ensembl/ensemble_out/merged_chr21_sorted_ver2.vcf.gz
"""

import pandas as pd
from cyvcf2 import VCF
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from collections import defaultdict
import sys
import csv


# Function used for filtering VNTRs
def get_normalized_sequence_similarity(a, b):
    # Compute Hamming distance between two short sequences.
    sequences_len = len(a)
    if len(a) != len(b):
        print("Warning: Computing sequence similarity between two sequences " +\
              "of different lengths: {} and {}".format(len(a), len(b)))
        sequences_len = min(len(a), len(b))
    similarity_score = 0
    for idx in range(sequences_len):
        # "M" character represent the masked character that matches
        # any other character for the purpose of this filter.
        if a[idx] == b[idx] or a[idx] == "M" or b[idx] == "M":
            similarity_score += 1
    return similarity_score / sequences_len

# Function used for filtering VNTRs
def get_motif_complexity_score(motif):
    self_match_score = 0
    for idx_in_mask in range(len(motif)):
        # Masking a single character at a time and computing the
        # max similarity score among all masked motifs.
        # Masked character matches any other character.
        masked_motif = motif[:idx_in_mask] + "M" + motif[idx_in_mask + 1:]
        # Creating a rolling window to compare the masked motif with itself.
        motif_window = masked_motif + masked_motif
        for idx in range(1, len(masked_motif)):
            end_idx_in_window = idx + len(masked_motif)
            # Compute the max score among all possible positions
            # when sliding the motif within the motif window to compare.
            self_match_score = max(self_match_score,
                                   get_normalized_sequence_similarity(masked_motif,
                                                         motif_window[idx:end_idx_in_window]))

    return self_match_score

# Function used for filtering VNTRs
def isValidVNTR(motif, motif_complexity_threshold):
    if get_motif_complexity_score(motif) > motif_complexity_threshold:
        # VNTR is very much STR like. Skip this VNTR.
        #print(motif)
        return False
    return True


try:
    pedfile = sys.argv[1]
    chrom = sys.argv[2]
    vcfpath = sys.argv[3]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

pedigree = pd.read_csv(pedfile, delim_whitespace=True)
trios = pedigree[(pedigree['FatherID'] != "0") & (pedigree['MotherID'] != "0")]
child_in_trios = set(trios['SampleID'].to_list())
mother_in_trios = set(trios['MotherID'].to_list())
father_in_trios = set(trios['FatherID'].to_list())
all_ids = (child_in_trios.union(mother_in_trios)).union(father_in_trios)

IDs = trios[['SampleID', 'MotherID', 'FatherID']]
trio_IDs = []
for index, row in IDs.iterrows():
    trio_IDs.append((row['SampleID'], row['MotherID'], row['FatherID']))

sys.stderr.write("Found %s trios...\n"%(len(trio_IDs)))

def CheckMI(sample_gt, mother_gt, father_gt):
    if sample_gt[0] in mother_gt[0:2] and sample_gt[1] in father_gt[0:2]:
        return True
    if sample_gt[1] in mother_gt[0:2] and sample_gt[0] in father_gt[0:2]:
        return True
    return False

def IsRef(gt):
    return gt[0]==0 and gt[1]==0

def printMIStats(stats_map):
    for counter in stats_map:
        print("{}: {}".format(counter, stats_map[counter]))

vcf = VCF(vcfpath, samples = list(all_ids))
samples = vcf.samples
# To store the stats on filtered or distarded variants.
stats_map = defaultdict(int)

for variant in vcf:
    n_families = 0  # Number of families with full call
    motif = variant.INFO.get('RU')
    stats_map["num all vntrs"] += 1
    if len(variant.ALT) == 0: # If there is no alt allele here
        stats_map["num vntrs with no alt allele"] += 1
        continue
    if not isValidVNTR(motif, motif_complexity_threshold=0.7):
        # VNTR is filtered out due to being str-like.
        stats_map["num str like vntrs"] += 1
        continue
    for family in trio_IDs:
        if "HG02567" in family: continue
        sample_index = samples.index(family[0])
        mother_index = samples.index(family[1])
        father_index = samples.index(family[2])
        fam_indices = [sample_index, mother_index, father_index]
        sample_GT = variant.genotypes[sample_index]
        mother_GT = variant.genotypes[mother_index]
        father_GT = variant.genotypes[father_index]
        if sample_GT[0] == -1 or mother_GT[0] == -1 or father_GT[0] == -1: # No call
            stats_map["num calls that are no call"] += 1
            continue
        if IsRef(sample_GT) and IsRef(mother_GT) and IsRef(father_GT): # all homozygous ref
            stats_map["num calls that are all homozygous ref"] += 1
            continue
        gbs=("%s,%s,%s"%(variant.format("GB")[sample_index],
                          variant.format("GB")[mother_index],
                          variant.format("GB")[father_index]))
        MI_val = CheckMI(sample_GT, mother_GT, father_GT)
        min_score_gt = np.min([variant.format('SCORE')[ind] for ind in fam_indices])
        items = [variant.CHROM, variant.POS, variant.INFO["RU"], family[0], variant.INFO["METHODS"], gbs, \
                 MI_val, min_score_gt]
        items = [variant.CHROM, variant.POS, variant.INFO["RU"], family[0], \
                 MI_val]
        sys.stdout.write("\t".join([str(item) for item in items])+"\n")

printMIStats(stats_map)
