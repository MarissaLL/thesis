#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */

pyscript = file("src/haploblocks_to_vcf.py")
chromosomes = Channel.from(1,2,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26)

// Just wrote this in using nano, from the header of the gwas set vcf
header_info = file("output/31_haploblocker/header_tmpl.txt")

chromosomes.map{ chr -> tuple("${chr}", 
				file("output/31_haploblocker/beagle_coverage_95/haplotype_startend_S${chr}.txt"),
				file("output/31_haploblocker/beagle_coverage_95/haplotype_blocks_S${chr}.txt"))}. set {test_ch}


process convert_genotypes {
	input:
		set val(chr),
			file(startend),
			file(blocks) from test_ch
	output:
		file("haploblocks_S${chr}.vcf") into headerless_vcf_ch
	"""
	python3 ${pyscript} ${chr}
	"""
}

process add_header_info {
	input:
		file(header_info)
		file(headerless_vcf) from headerless_vcf_ch
	output:
		file("${headerless_vcf.baseName}.full.vcf") into full_vcf_ch
	"""
	cat ${header_info} ${headerless_vcf} > ${headerless_vcf.baseName}.full.vcf
	"""

}

process sort_vcf {
	input:
		file(vcf) from full_vcf_ch
	output:
		file("${vcf.simpleName}.sorted.vcf") into sorted_vcf_ch
	"""
	bcftools sort ${vcf} > ${vcf.simpleName}.sorted.vcf
	"""
}


process cat_vcfs {
	
	input:
		file(vcf) from sorted_vcf_ch.collect()
		// file(vcf) from gz_ch.collect()
		// file(tbi) from tbi_ch.collect()
	output:
		file("all_haploblocks.vcf") into cat_vcf_ch
	"""
	bcftools concat ${vcf} > all_haploblocks.vcf
	"""
}

process sort_vcf_again {
	publishDir './output/31_haploblocker/beagle_coverage_95/', mode: 'copy'
	input:
		file(vcf) from cat_vcf_ch
	output:
		file("coverage_95_haploblocks.vcf")
	"""
	bcftools sort ${vcf} > coverage_95_haploblocks.vcf
	"""
 }
