#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */

// Scripts
py_script = file("t_add_alphapeel_haplos.py")
filter_script = file("src/alphapeel_assess_phase_qual.R")
changepoint_script = file("src/alphapeel_recombination_points.R")
convert_index_script = file("src/alphapeel_convert_snp_index.R")

// Files
seq_birds = file("data/all_seq_birdnames.txt")
vcf_for_alphapeel = file("output/21_alphapeel/snps_more_filtering.vcf")

// Values
ncycle_values = Channel.from(20, 100,200,1000)
Channel.from(1,2,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26)
	.into{chromosomes; chromosomes2}
parent = Channel.from("mat", "pat")

//Options
add_haplos = "yes"



if(add_haplos == "yes"){
	ncycle_values.combine(chromosomes)
			.map {ncycles, chr -> tuple("$ncycles", "$chr", file("output/21_alphapeel/raw_alphapeel_files/S${chr}_output_ncycles${ncycles}.seg"))}
			.set{seg_file_ch}




 	process reformat_sum_alphapeel_output {
 		publishDir './output/25_tidy_alphapeel/added_haplos/', mode: 'copy'
		input:
			file py_script
			set ncycles, chr, file(seg_file) from seg_file_ch
		output:
			set ncycles, chr, file("added_haplos_S${chr}_ncycles${ncycles}.txt") into (added_haplo_ch_1, added_haplo_ch_2)
		

		"""
		python3 ${py_script} S${chr} ${ncycles}
		"""
	}
}

if(add_haplos == "no"){
// // Excluding adding step (because it has already run separately and it is slow)
	ncycle_values.combine(chromosomes)
		.map{ncycles, chr -> tuple("$ncycles", "$chr", file("output/25_tidy_alphapeel/added_haplos/added_haplos_S${chr}_ncycles${ncycles}.txt"))}
		.set {added_haplo_ch}

	added_haplo_ch.into{added_haplo_ch_1; added_haplo_ch_2}
}

// Exclude chromosomes from individuals if they are phased poorly
process filter_on_phase_qual {
	publishDir './output/25_tidy_alphapeel/phasing_filter/', mode: 'copy'
	tag { "S${chr}_${ncycles}_${parent}" }
	input:
		file filter_script
		file seq_birds
		each parent
		set ncycles, chr, file(added_file) from added_haplo_ch_1
	output:
		set ncycles, chr, parent, file("indivs_to_keep_S${chr}_ncycles_${ncycles}_${parent}.txt") into indivs_to_keep_ch
		set ncycles, chr, parent, file("filter_res_S${chr}_ncycles_${ncycles}_${parent}.txt") into filter_res_ch

	"""
	Rscript ${filter_script} S${chr} ${ncycles} ${parent}_only ${seq_birds}
	"""

}

indivs_to_keep_ch.combine(added_haplo_ch_2, by: [0,1])
	.set{for_recombination_ch}


// Call recombinations based on changepoints in the mean of inheritance values
//The script here can also be run as a standalone to produce plots
process call_recombinations {
	// publishDir './output/25_tidy_alphapeel/changepoint_files/', mode: 'copy', pattern: 'changepoints_S*'
	tag { "S${chr}_${ncycles}_${parent}" }
	input:
		file changepoint_script
		set ncycles, 
			chr, 
			parent, 
			file(indivs_to_keep),
			file(added_haplos) from for_recombination_ch
	output:
		set chr, 
			ncycles, 
			parent, 
			file("cpt_pos_*_S${chr}_ncycles${ncycles}_${parent}.txt") into cpt_pos_ch
		set chr,
			ncycles,
			parent,
			file("changepoints_S${chr}_ncycles${ncycles}_${parent}.txt") into changepoints_ch

	"""
	Rscript ${changepoint_script} S${chr} ${ncycles} ${parent}_only
	"""

}

changepoints_ch.groupTuple(by:[0,1])
	.set{changepoints_grouped_ch}

process concat_changepoint_files {
	publishDir './output/25_tidy_alphapeel/changepoint_files/', mode: 'copy'
	input:
		set chr,
			ncycles,
			parents,
			file(changepoints_files) from changepoints_grouped_ch
	output:
		file("changepoints_S${chr}_ncycles${ncycles}.txt")
	"""
	cat ${changepoints_files} > "changepoints_S${chr}_ncycles${ncycles}.txt"
	"""
}

// Concatenate recombination point output for all individuals into files per chromosome
process concatenate_cpt_pos_files {
	publishDir './output/25_tidy_alphapeel/recombination_positions_for_100/', mode: 'copy'
	errorStrategy 'ignore' // Don't stop if files are empty
	tag { "S${chr}_${ncycles}_${parent}" }
	input:
		set chr,
			ncycles,
			parent,
			file(cpt_pos) from cpt_pos_ch
	output:
		set chr, 
		ncycles,
		parent,
		file("cpt_pos_S${chr}_ncycles_${ncycles}_${parent}.txt") into full_cpt_pos_ch
	"""
	cat ${cpt_pos} | grep -v "name" > cpt_pos_S${chr}_ncycles_${ncycles}_${parent}.txt
	"""

}

// Extract positions of snps per chromosome in the vcf used as input for alphapeel
process extract_vcf_positions {
	tag { "S${chr}" }
	input:
		val chr from chromosomes2
		file vcf_for_alphapeel
	output:
		file("S${chr}_pos.txt") into vcf_pos_ch
	"""
	grep -w "S${chr}" ${vcf_for_alphapeel} | cut -f1,2 | grep -v "##" > S${chr}_pos.txt
	"""
}

// Convert snp indicies in the alphapeel output back to chromosomal position
// This doesn't work properly for ncycles100 because it was run with different alphapeel settings 
// so results for ncycles100 are currently in recombination_positions_for_100
process convert_snp_index_to_pos {
	publishDir './output/25_tidy_alphapeel/recombination_positions/', mode: 'copy'
	tag { "S${chr}_${ncycles}_${parent}" }
	input:
		file convert_index_script
		file vcf_pos from vcf_pos_ch.collect()
		set chr, 
			ncycles,
			parent,
			file(cpt_pos) from full_cpt_pos_ch
	output:
		set ncycles,
			file("recombination_S${chr}_ncycles_${ncycles}_${parent}.txt") into recombo_pos_ch
	"""
	Rscript ${convert_index_script} ${chr} ${ncycles} ${parent}
	"""
}

// Group all the chromosome files for each ncycles value
recombo_pos_ch.groupTuple()
	.set {grouped_recombo_ch}

// Concatenate recombination by chromosome into files per ncycles. Ncycles 20 is the one I plan to use.
process concatenate_recombo_by_ncycles {
	publishDir './output/25_tidy_alphapeel/recombination_positions/', mode: 'copy'
	input:
		set ncycles,
			file(files) from grouped_recombo_ch
	output:
		file("all_recombination_ncycles_${ncycles}.txt")

	"""
	cat ${files} >  recombo.txt
	head -n 1 recombo.txt > all_recombination_ncycles_${ncycles}.txt
	grep -v "snp_index" recombo.txt >> all_recombination_ncycles_${ncycles}.txt
	"""
}

