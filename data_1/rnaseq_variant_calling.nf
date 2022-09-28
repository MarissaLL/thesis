#!/usr/bin/env nextflow

/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */

//Also requires tabix and bgzip which I should possibly put in a container


genome = file("data/small_assembly.fasta")
genome_fai = file("data/small_assembly.fasta.fai")
deepvariant_vcf = file("data/deepvariant.r1.vcf")


Channel
.fromFilePairs('../rna-seq/files/HWLNGBCX2-4836-0*-49-01_S*_L00*_R{1,2}_001.fastq.gz', size: -1)              
    .set { raw_rnaseq_all }


//raw_rnaseq_all.subscribe {println it}


//raw_rnaseq_all.subscribe {println it}

//Decompress and merge all the forward reads together and all the reverse reads together.
/*Normalization based on lane shouldn't be quite such a big deal as for differential expression data but
might still need to be done */


process convert_rnaseq_data {
    //publishDir './output/05_rnaseq_parentage_assignment/', mode: 'copy'
    
    input:
        set lane_etc_id, file(raw_reads) from raw_rnaseq_all
    output:
        file("merged_R1.fastq") into (merged_R1_ch, merged_R1_ch2, merged_R1_ch3)
        file("merged_R2.fastq") into (merged_R2_ch, merged_R2_ch2)

    """
    gunzip ${raw_reads[0]} -c > merged_R1.fastq
    gunzip ${raw_reads[1]} -c > merged_R2.fastq
    """
    }


// STAR index the genome
//remove genomesaindexnbases later - this is just for the subset file
//genome_index_1 dir needs to be made beforehand if im using  ${baseDir}/output/05_rnaseq_parentage_assignment/genome_index_1
process index_genome {
    input:
        file(genome_file) from genome
    output:
        set file("chrStart.txt"),
            file("SA"),
            file("chrName.txt"),
            file("chrLength.txt"),
            file("chrNameLength.txt"),
            file("SAindex"),
            file("Genome"),
            file("genomeParameters.txt") into genome_index_1_ch2 

    """
    STAR --runMode genomeGenerate \
    --genomeDir . \
    --genomeFastaFiles ${genome_file} \
    --runThreadN 5 \
    --genomeSAindexNbases 8
    """
}


//genome_index_1_ch.collect().subscribe{println it}

process alignment_pass_one {
    publishDir './output/05_rnaseq_parentage_assignment/', mode: 'copy'

    input:
        set file(chrStart),
            file(sa),
            file(chrName),
            file(chrLength),
            file(chrNameLength),
            file(saindex),
            file(genome),
            file(genomeParameters) from genome_index_1_ch2
        file(fwd_reads) from merged_R1_ch2.collectFile() 
        file(rev_reads) from merged_R2_ch.collectFile()
    output:
        file("Aligned.out.sam") into aligned_one_sam_ch
        file("SJ.out.tab") into sjout_ch

    """
    STAR --genomeDir . \
    --readFilesIn ${fwd_reads} ${rev_reads} \
    --runThreadN 5
    """
}



// --sjdbOverhang is the length of the donor/acceptor sequence on each side of the junctions, 
// ideally = (mate_length - 1). Set to 199 as the average mapped length was 200


process create_new_sj_index {
    input:
        file(genome_file) from genome
        file(sam) from aligned_one_sam_ch
        file(sjout) from sjout_ch
    output:
        set file("sjdbList.out.tab"),
            file("sjdbInfo.txt"),
            file("chrStart.txt"),
            file("SA"),
            file("chrName.txt"),
            file("chrLength.txt"),
            file("chrNameLength.txt"),
            file("SAindex"),
            file("Genome"),
            file("genomeParameters.txt") into index_2_ch

    """
    STAR  --runMode genomeGenerate \
    --genomeDir .   \
    --genomeFastaFiles ${genome_file} \
    --sjdbFileChrStartEnd ${sjout} \
    --sjdbOverhang 199 \
    --runThreadN 5
    """
}

process create_final_alignments {
    publishDir './output/05_rnaseq_parentage_assignment/', mode: 'copy'
    input:
        set file(sjdbList),
            file(sjdbInfo),
            file(chrStart),
            file(sa),
            file(chrName),
            file(chrLength),
            file(chrNameLength),
            file(saindex),
            file(genome),
            file(genomeParameters) from index_2_ch
        file(fwd_reads) from merged_R1_ch3.collectFile() 
        file(rev_reads) from merged_R2_ch2.collectFile()
    output:
        file("Aligned.out.sam") into aligned_two_sam_ch

    """
    STAR --genomeDir . \
    --readFilesIn ${fwd_reads} ${rev_reads} \
    --runThreadN 5
    """
}


// the flags for these are sort order, read group id, read group library, platform, machine and sample.
// These are required. The important bit is that it is sorting by coordinate here though 

process coord_sort_sam {
    conda 'bioconda::picard=2.20.3-0'

    input:
        file(sam) from aligned_two_sam_ch
    output:
        file("rg_added_sorted.bam") into sorted_bam_ch

   """
    picard AddOrReplaceReadGroups \
    I=${sam} \
    O=rg_added_sorted.bam \
    SO=coordinate \
    RGID=id \
    RGLB=library \
    RGPL=platform \
    RGPU=machine \
    RGSM=nora-1b
    """
    }




process  mark_duplicates {
    conda 'bioconda::picard=2.20.3-0'

    input:
        file(sorted_bam) from sorted_bam_ch
    output:
        file("dedupped.bam") into dedupped_bam_ch
    """
    picard MarkDuplicates \
    I=${sorted_bam} \
    O=dedupped.bam  \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    M=output.metrics 
    """
}



// at this point, should consider doing indel realignment and base recalibration



process make_seq_dict {
    conda 'bioconda::picard=2.20.3-0'
    input:
        file genome
    output:
        file("small_assembly.dict") into (dict_ch, dict_ch2)
    """
    picard CreateSequenceDictionary \
    R= ${genome} \
    O=small_assembly.dict
    """
}


// Have to run this section on gatk 3.8 because gatk 4 has either removed or renamed these tools and I can't find t
// the fai and dict don't need to be called directly but they do need to be in the same directory

process reassign_map_qual {
    input:
        file genome
        file genome_fai
        file(dict) from dict_ch
        file(dedupped_bam) from dedupped_bam_ch
    output:
        set file("new_map_quals.bam"),
            file("new_map_quals.bai") into new_map_quals_bam_ch

    """
    java -jar /usr/GenomeAnalysisTK.jar \
    -T SplitNCigarReads \
    -R ${genome} \
    -I ${dedupped_bam} \
    -o new_map_quals.bam \
    -rf ReassignOneMappingQuality \
    -RMQF 255 \
    -RMQT 60 \
    -U ALLOW_N_CIGAR_READS
    """
}



process variant_calling {
    publishDir './output/05_rnaseq_parentage_assignment/', mode: 'copy'
    input:
        file genome
        file genome_fai
        file genome_dict from dict_ch2
        set file(final_bam),
            file(final_bai) from new_map_quals_bam_ch
    output:
        file("nora-1-b_rnaseq.vcf") into (full_vcf_ch, full_vcf_ch2)

    """
     java -jar /usr/GenomeAnalysisTK.jar \
     -T HaplotypeCaller \
     -R ${genome} \
     -I ${final_bam} \
     -dontUseSoftClippedBases \
     -stand_call_conf 20.0 \
     -o nora-1-b_rnaseq.vcf
     """
}


process find_shared_variant_positions {
    publishDir './output/05_rnaseq_parentage_assignment/', mode: 'copy'
    input:
        file(full_vcf_nora1b) from full_vcf_ch
        file deepvariant_vcf
    output:
        file("shared_snps.txt") into shared_snps_ch
    """
    grep -v "#" ${full_vcf_nora1b} | cut -f 1,2 > nora1b_list_of_snps.txt
    grep -v "#" ${deepvariant_vcf} | cut -f 1,2 > pop_list_of_snps.txt
    grep -Fwf pop_list_of_snps.txt nora1b_list_of_snps.txt > shared_snps.txt
    """
}


process subset_convert_vcfs {
    conda 'bioconda::vcftools'

    input:
        file(nora1b_vcf) from full_vcf_ch2
        file deepvariant_vcf
        file(shared_snps) from shared_snps_ch
    output:
        set file("nora1b_subset.recode.vcf.gz"),
            file("nora1b_subset.recode.vcf.gz.tbi") into nora1b_converted_vcf_ch
        set file("deepvar_subset.recode.vcf.gz"),
            file("deepvar_subset.recode.vcf.gz.tbi") into deepvariant_converted_vcf_ch
    
    """
    vcftools --vcf ${nora1b_vcf} \
    --positions ${shared_snps} \
    --recode \
    --out nora1b_subset
    bgzip nora1b_subset.recode.vcf
    tabix -p vcf nora1b_subset.recode.vcf.gz
    vcftools --vcf ${deepvariant_vcf} \
    --positions ${shared_snps} \
    --recode \
    --out deepvar_subset
    bgzip deepvar_subset.recode.vcf
    tabix -p vcf deepvar_subset.recode.vcf.gz
    """
}

process merge_vcfs {
    publishDir './output/05_rnaseq_parentage_assignment/', mode: 'copy'
    conda 'bioconda::vcftools'

    input:
        set file(nora1b_vcfgz),
            file(nora1b_tbi) from nora1b_converted_vcf_ch
        set file(deepvar_vcfgz),
            file(deepvar_tbi) from deepvariant_converted_vcf_ch
    output:
        file("merged.vcf") into merged_vcf_ch
    """
    vcf-merge ${nora1b_vcfgz} ${deepvar_vcfgz} > merged.vcf
    """
}


process filter_merged_vcf {
    publishDir './output/05_rnaseq_parentage_assignment/', mode: 'copy'
    conda 'bioconda::vcftools'

    input:
        file(merged_vcf) from merged_vcf_ch
    output:
        file("merged_filtered.recode.vcf")
    """
    vcftools \
    --vcf ${merged_vcf} \
    --minQ 3 \
    --max-missing 0.8 \
    --recode \
    --out merged_filtered
    """

}

// Now run check_uncertain_parentage.nf with the merged_filtered.vcf as input
/**/