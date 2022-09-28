
# Modified from a script written by J Guhlin

from cyvcf2 import VCF
from cyvcf2 import Writer

import numpy as np
import pandas as pd
import sys


# From: https://stackoverflow.com/questions/25699439/how-to-iterate-over-consecutive-chunks-of-pandas-dataframe-efficiently
def chunks(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


vcf_template = sys.argv[1]
chr = sys.argv[2]
cutoff = float(sys.argv[3])

haps = pd.read_csv(chr + "_output_ncycles20.haps", sep="\s+", header=None)


samples_chunked = list(chunks(haps, 4))


for sample_id_num in range(len(samples_chunked)):
    phasing = list()
    z = samples_chunked[sample_id_num]
    sample_name = z[0].iloc[0]
    if sample_name in ('DUMAL001','John-girl', 'Pegasus'):
        print(sample_name + " is being skipped because sample is not sequenced.")
    else:
        allele_probs = z.drop(columns=0)
        nsnps = allele_probs.shape[1]

        for i in range(nsnps):
            probs = allele_probs.iloc[:, i]
            if np.any(probs >= cutoff):
                if probs.iloc[0] >= cutoff:
                    phasing.append(0) # 0 is aa 
                if probs.iloc[1] >= cutoff:
                    phasing.append(1) # 1 is aA
                if probs.iloc[2] >= cutoff:
                    phasing.append(2) # 2 is Aa
                if probs.iloc[3] >= cutoff:
                    phasing.append(3) # 3 is AA
            else:
                phasing.append(False)
        print(sample_name, "Nsnps:", nsnps, "Phasing:", len(phasing))
        writer = Writer("phased_" + sample_name + "_" + chr +".vcf", 
                        tmpl=VCF(vcf_template, samples = sample_name))
        writer.write_header()

        for i, variant in enumerate(VCF(vcf_template, samples = sample_name)(chr)):
            if phasing[i]:
                phase = phasing[i]
                ## aa, aA, Aa, AA
                if phase == 0:
                    variant.genotypes[0] = [0, 0, True]
                elif phase == 1:
                    variant.genotypes[0] = [0, 1, True]
                elif phase == 2:
                    variant.genotypes[0] = [1, 0, True]
                elif phase == 3:
                    variant.genotypes[0] = [1, 1, True]
            variant.genotypes = variant.genotypes
            writer.write_record(variant)
        del(writer)
        print("Phased vcf written for ", sample_name, " ID number ", sample_id_num)