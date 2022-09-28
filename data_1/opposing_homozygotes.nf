#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */
oh_script = file("src/opposing_homozygotes.py")

 variants_ch = Channel.from(tuple('original', file("data_stored/deepvariant.r1.vcf")),
 						tuple('new', file("data_stored/gwas_set_numeric_maf_0.05_f_missing_0.2.bcf")))


process make_vcf_all_chr {
 	publishDir './output/27_opposing_homozygotes/', mode: 'copy'
		input:
			set version,
				file(variants) from variants_ch
		output:
			set version,
				file("opposing_sites_${version}.txt") into (oh_loci_ch1, oh_loci_ch2)
		

	"""
	bcftools view ${variants} | grep "0/0" | grep "1/1" > opposing_sites.vcf
	bcftools view ${variants} | grep "#CHROM" > header.txt
	cat header.txt opposing_sites.vcf > opposing_sites_${version}.txt
	"""
}
oh_loci_ch1.map{version, variants -> tuple("$version", "all", file("$variants"))}
		.set {all_chr_ch}


process remove_sex_chr {
	publishDir './output/27_opposing_homozygotes/', mode: 'copy'
		input:
			set version,
				file(opposing_sites) from oh_loci_ch2
		output:
			set version,
				file("opposing_sites_${version}_autosomes.txt") into autosome_ch1
	"""
	grep -v "S3" ${opposing_sites} | grep -v "S13" | \
	grep -v "CM013773.1" | grep -v "CM013763.1" > "opposing_sites_${version}_autosomes.txt"
	"""
}

autosome_ch1.map{version, variants -> tuple("$version", "autosome", file("$variants"))}
		.set {autosome_ch}

sites_files_ch = all_chr_ch.concat(autosome_ch)

process calculate_opposing_homozygotes {
 	publishDir './output/27_opposing_homozygotes/', mode: 'copy'
 		input:
 			file oh_script
 			set version,
 				chrs, 
 				file(opposing_sites) from sites_files_ch
 		output:
 			file("opposing_homozygotes_${version}_${chrs}.txt")
 	"""
 	python3 ${oh_script} ${opposing_sites} ${version} ${chrs}
 	"""

 }






