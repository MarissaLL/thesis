#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */


vcf = file("data_stored/Trained.bcf")
bird_data = file("data/bird_info_rships_nov20.csv")
// chrs_ch = Channel.from("CM013761.1", "CM013762.1", "CM013764.1", "CM013765.1", 
// 	"CM013766.1", "CM013767.1", "CM013768.1", "CM013769.1", "CM013770.1", 
// 	"CM013771.1", "CM013772.1", "CM013774.1", "CM013775.1", "CM013776.1", 
// 	"CM013777.1", "CM013778.1", "CM013779.1", "CM013780.1", "CM013781.1", 
// 	"CM013782.1", "CM013783.1", "CM013784.1", "CM013785.1", "CM013786.1", 
// 	"CM013773.1", "CM013763.1")



chrs_ch = Channel.from(1..26).map{num -> "S${num}"}
seq_birds = file("data/all_seq_birdnames.txt")


fimpute_R_script = file("src/make_fimpute3_input_snps.R")
decoding_R_script = file("src/make_decoded_haplo_files.R")
plot_unique_haplos_script = file("src/count_haplos_founders_vs_new.R")
// fimpute_to_vcf_R_script = file("src/fimpute_to_vcf.R") // no longer used in favour of fimpute_to_vcf.py
fimpute_to_vcf_py_script = file("src/fimpute_to_vcf.py")

approach = "make_vcf"


process missing_filter_vcf {
	//publishDir './output/10_fimpute_phasing', mode: 'copy'
	
	input:
		file vcf
	output:
		file("miss08_dv.vcf") into vcf_filtered_ch
		
	"""
	bcftools view ${vcf} -m2 -M2 -v snps -e "QUAL<20 & N_MISSING > 0.1" > miss08_dv.vcf
	"""
}

process bgzip_vcf {
	input:
		file(vcf) from vcf_filtered_ch
	output:
		file("${vcf}.gz") into bgzip_vcf_ch
	"""
	bgzip ${vcf}
	"""
}

process tabix_vcf {
	input:
		file(vcf_gz) from bgzip_vcf_ch
	output:
		set file(vcf_gz),
			file("${vcf_gz}.tbi") into tbi_ch
	"""
	tabix ${vcf_gz}
	"""
}


process split_by_chr {
	// publishDir './output/10_fimpute_phasing', mode: 'copy'

	input:
		each chr from chrs_ch
		file(vcf_filtered) from tbi_ch
	output:
		set val(chr), 
			file("filtered_chr_${chr}.vcf") into (split_by_chr_ch, split_by_chr_ch2)

	"""
	bcftools view ${vcf_filtered} -r ${chr} > filtered_chr_${chr}.vcf
	"""
}

process bgzip_chr_vcfs {
	// publishDir './output/10_fimpute_phasing', mode: 'copy'
	input:
		set val(chr),
			file(vcf) from split_by_chr_ch2
	output:
		set val(chr),
			file("${vcf}.gz"),
			file("${vcf}.gz.tbi") into (tabixed_indiv_vcf_ch, tabixed_indiv_vcf_ch2)
	"""
	bgzip ${vcf}
	tabix ${vcf}.gz
	"""
}



process make_fimpute3_inputs {	
	publishDir './output/10_fimpute_phasing/input_files/', mode: 'copy'
	tag { "${chr}" }

	input:
		file fimpute_R_script
		file bird_data
		set val(chr),
			file(split_vcf) from split_by_chr_ch
	output:
		set val(chr), 
			file("chr_${chr}.txt") into genotypes_ch
		set val(chr), 
			file("snps_info_${chr}.txt") into snps_info_ch
		file("fimpute_pedigree.txt") into fimpute_pedigree_ch
		set val(chr),
			file("fimpute_config_${chr}.txt") into fimpute_config_ch

	"""
	Rscript ${fimpute_R_script} ${chr} ${bird_data}
	"""
}




process tidy_genotypes_files {
	publishDir './output/10_fimpute_phasing/input_files/', mode: 'copy'
	tag { "${chr}" }

	input:
		set val(chr), 
			file(genotypes) from genotypes_ch
	output:
		set val(chr),
			file("fimpute_${genotypes.baseName}.txt") into fimpute_genotypes_ch

	"""
	tr ':' '\t' < ${genotypes} > fimpute_${genotypes.baseName}.txt
	"""
}

fimpute_config_ch.join(fimpute_genotypes_ch).join(snps_info_ch).set{input_files_ch}
//input_files_ch.subscribe{println it}


process run_fimpute {
	publishDir './output/10_fimpute_phasing/', mode: 'copy'
	tag { "${chr}" }

 	input:
 		set val(chr),
 			file(config),
 			file(genotypes),
 			file(snps_info) from input_files_ch
 		file(pedigree) from fimpute_pedigree_ch
 	output:
 		set val(chr),
 			file("genotypes_imp_${chr}.txt") into genotypes_imp_ch
 		file("report_${chr}.txt")

 	"""
 	~/Downloads/FImpute3 ${config} -o
 	mv output_${chr}/genotypes_imp.txt genotypes_imp_${chr}.txt
 	mv output_${chr}/report.txt report_${chr}.txt
 	"""
}

tabixed_indiv_vcf_ch.join(genotypes_imp_ch).set{chrs_haplos_ch}


process make_decoded_haplo_files {
	publishDir'./output/10_fimpute_phasing/', mode: 'copy'
	queuesize = 1
	memory = 24.GB
	tag {"${chr}"}

 	input:
 		file decoding_R_script
 	 	file seq_birds
 		set val(chr),
 			file(vcf),
 			file(tbi),
 			file(genotypes_phased) from chrs_haplos_ch
 		
 	output:
 		set val(chr),
 			file("decoded_haplos_vcflike_${chr}.txt") into decoded_haplos_ch
 		file("decoded_haplos_${chr}.txt")
 	"""
 	Rscript ${decoding_R_script} ${chr} ${seq_birds} ${genotypes_phased} ${vcf}
 	"""
}

tabixed_indiv_vcf_ch2.join(decoded_haplos_ch).set{for_py_vcf_ch}

if(approach == "make_vcf"){

	process make_vcf {
		input:
			file fimpute_to_vcf_py_script
			set val(chr),
				file(vcf_gz),
				file(vcf_tbi),
				file(decoded_haplo_vcflike) from for_py_vcf_ch
		output:
			set val(chr),
				file("fimpute_*_${chr}.vcf") into indiv_vcfs_ch mode flatten
	"""
	python3 ${fimpute_to_vcf_py_script} ${chr}
	"""
	}

	process bgzip_indiv_vcfs {
		input:
			set val(chr),
				file(vcf) from indiv_vcfs_ch
		output:
			set val(chr),
				file("${vcf}.gz") into (bgzip_indiv_ch, bgzip_indiv_ch2)
		"""
		bgzip ${vcf}
		"""
	}

	process tabix_indiv_vcfs {
		input:
			set val(chr),
				file(vcf_gz) from bgzip_indiv_ch
		output:
			set val(chr),
				file("${vcf_gz}.tbi") into tbx_indiv_ch
		"""
		tabix ${vcf_gz}
		"""
	}

	bgzip_indiv_ch2.groupTuple(by:0)
	.set{grouped_vcfs_ch}

	tbx_indiv_ch.groupTuple(by:0)
	.set{grouped_tbis_ch}

	// Put the vcf.gz and vcf.gz.tbi files in the same channel
	grouped_vcfs_ch.join(grouped_tbis_ch)
	.set{all_by_chr_ch}


	process combine_vcfs {
		publishDir './output/10_fimpute_phasing/', mode: 'copy'
		queuesize = 1
		input:
			set val(chr),
				file(vcf_gz),
				file(tbi) from all_by_chr_ch
		output:
			file("fimpute_phased_${chr}.vcf") into per_chr_vcf_ch
		"""
		bcftools merge ${vcf_gz} --output fimpute_phased_${chr}.vcf
		"""
	}

	// Make sure indivs are in the same order in every vcf
	process order_samples {
		input:
			file(vcf) from per_chr_vcf_ch
		output:
			file("sorted_${vcf}") into sorted_per_chr_vcf_ch
		"""
		bcftools query -l ${vcf} | sort > ordered_samples.txt
		bcftools view -S ordered_samples.txt ${vcf} > sorted_${vcf}
		"""
	}

	// Bgzip again
	process bgzip_per_chr {
		input:
			file(vcf) from sorted_per_chr_vcf_ch
		output:
			file("${vcf}.gz") into (bgzip_per_chr_ch1, bgzip_per_chr_ch2)
		"""
		bgzip ${vcf}
		"""
	}

	// Tabix again
	process tbi_per_chr {
		input:
			file(vcf_gz) from bgzip_per_chr_ch1
		output:
			file("${vcf_gz}.tbi") into tbi_per_chr_ch
		"""
		tabix ${vcf_gz}
		"""
	}

	process combine_chromosomes {
		publishDir './output/10_fimpute_phasing/', mode: 'copy'
		input:
			file(vcf_gz) from bgzip_per_chr_ch2.collect()
			file(tbi) from tbi_per_chr_ch.collect()
		output:
			file("fimpute_phased_allchr.vcf")
		"""
		bcftools concat ${vcf_gz} --output fimpute_phased_allchr.vcf
		"""
	}
}


if(approach == "make_plots"){
	process plot_unique_haplos {
		publishDir './output/10_fimpute_phasing', mode: 'copy'

		input:
			file plot_unique_haplos_script
			set val(chr),
				file(decoded_haplos) from decoded_haplos_ch
			file bird_data
		output:
			file("chr_*_unique_haplos.png")
			file("summary_${chr}_haplos.txt")
			file("ath_all_${chr}.txt")
			file("ath_founders_${chr}.txt")

		"""
		Rscript ${plot_unique_haplos_script} ${decoded_haplos} "chr_${chr}_unique_haplos.png" ${bird_data}
		"""
	}
}




