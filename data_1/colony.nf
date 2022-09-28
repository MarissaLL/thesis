#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */

publish_path_base = "/media/drive_6tb/projects/kakapo-genomics/output/39_colony_retry/"

generate_input_src = file("src/make_colony_input_file.R")
ped_file = file("data/bird_info_rships_nov20.csv")
n_snp_ch = Channel.from(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
rep_ch = Channel.from(1,2,3,4,5,6,7,8,9,10)
dataset_ch = Channel.from("deepvariant","gwas_set")

n_snp_ch.combine(rep_ch).combine(dataset_ch)
	.map{ n_snp,rep,dataset -> tuple("${n_snp}",
				"${rep}",
				"${dataset}", 
					file("/media/drive_6tb/projects/kakapo-genomics/output/37_sequoia_pt2/sample_${n_snp}_rep${rep}_${dataset}.raw"))}
					.set{snps_ch}

process make_colony_input {
	publishDir "${publish_path_base}${dataset}/", mode: 'copy'
	input:
		file(generate_input_src)
		file(ped_file)
		set val(n_snp),
			val(rep),
			val(dataset),
			file(gts_file) from snps_ch
	output:
		set val(n_snp),
			val(rep),
			val(dataset),
			file("colony_sample_${n_snp}_rep${rep}.dat") into colony_dat_ch
	"""
	Rscript ${generate_input_src} kakapo_colony_pedigree_test ${n_snp} ${rep} ${dataset}
	"""
} 


process run_colony {
	publishDir "${publish_path_base}${dataset}", mode: 'copy'
	input:
		set val(n_snp),
			val(rep),
			val(dataset),
			file(colony_dat) from colony_dat_ch
	output:
		file("colony_out_sample${n_snp}_rep${rep}.BestConfig") optional true
		file("colony_out_sample${n_snp}_rep${rep}.log")
	"""
	~/Downloads/COLONY2/colony2s.ifort.out IFN:${colony_dat} > colony_out_sample${n_snp}_rep${rep}.log
	"""
}
