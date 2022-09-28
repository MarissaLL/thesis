#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */

/*NOTES
This is intended to be run after step.nf - before that it wasn't known that the original four dubious chicks needed checking.
Turns out there might actually be five - once the spelling is corrected so that it is included, 
 Elliott also seems to have errors
 Also used this script to check the no assigned paternity chicks. need to comment and uncomment lines in the R
 script to run it on different individuals

 */

// CONFIG FILE IS cupn_nf.config 
// either use dubious_ped_alts.R, check_sibling_parent_swap_me.R 
tidy_wiki_birdlist = file("data/bird_info_rships.csv")

use = "swaps" //"dubious_ped_alts" or "swaps"
make_plink_instruction_swaps_src = file("src/check_sibling_parent_swap_me.R")
make_plink_instruction_src = file("src/dubious_ped_alts.R")
n_variants = Channel.value(200000)
use_bcf = true // Using the gwas set vcf file doesn't currently work. Nor does this vcftools command
bcf = file("data_stored/Trained.bcf")
vcf = file("data_stored/gwas_set.vcf")
// outdir = "output/03_mendelian_error_distribution/"
// outdir = "output/04_test_dubious_trios/"
outdir = "output/041_test_dubious_trios_2020/"


founders = ["Alice", "Bella", "Cyndy", "Flossie", "Jean", "Lisa", "Margaret-Maree", "Nora", "Ruth", "Solstice",
"Sue", "Suzanne", "Arab", "Basil", "Blades", "Bonus", "Boss", "Felix", "Gumboots", "Joe", "Bill", "Fuchsia",
"Gunner", "Lee", "Luke", "Maggie", "Merv", "Nog", "Ox", "Sandra", "Barnard", "Merty", "Piripi", "Ralph", "Rangi",
"Sarah", "Wendy", "Sass", "Richard_Henry", "Whiskas", "Waynebo", "Ben", "Lionel", "Smoko", "Jimmy"]


if (use_bcf == true) {
	chr_handling = "--allow-extra-chr"

// Ideally I'd replace this with a bcftools command instead
process quality_missing_filter_vcf {
	conda 'bioconda::vcftools'
	input:
		file bcf
	output:
		file("qual3_missing_filter.recode.vcf") into (vcf_filtered_ch, vcf_filtered_ch2)
		
	"""
	vcftools --bcf ${bcf} --minQ 3 --max-missing 0.8 --not-chr S3 --not-chr S13 --recode \
	--out qual3_missing_filter 
	"""
}
}

if(use_bcf == false){
	Channel.from(vcf).into{vcf_filtered_ch; vcf_filtered_ch2}
	chr_handling = "--chr-set 89"
}

//CAUTION - conversion vcf to ped seems to remove underscores from the names! 
// Hence this step has got very convoluted with the sed commands
// Make the vcf into a ped file for plink to use. 
//Take only n_variants number of SNPs (e.g. 200,000 SNPs), drawn randomly
process convert_vcf_to_ped {
	input:
		file(q3_vcf) from vcf_filtered_ch
		val n_variants
		chr_handling
	output:
		file("part_two.ped") into ped_file_ch
		file("part_two.map") into map_file_ch
		file("subset_snps.txt")
		file("list_of_snps.txt")
		

	"""
	sed -i 's/-1/./g' ${q3_vcf}
	sed -i 's/_/~~/g' ${q3_vcf}
	plink --vcf ${q3_vcf} ${chr_handling} --vcf-half-call missing --recode --out part_one
	sed -i 's/~~/_/g' part_one.ped
	sed -i 's/~~/_/g' part_one.map
	cat part_one.map | cut -f 2 > list_of_snps.txt
	awk -F"\t" '!uniq[\$0]++' list_of_snps.txt > list_of_uniq_snps
	shuf list_of_uniq_snps | head -n $n_variants > subset_snps.txt
	plink --file part_one ${chr_handling} --extract subset_snps.txt \\
	--recode --out part_two
	"""
}


//Makes files to update id, parents and sex information in the ped file so plink can run mendelian error tests
// Currently the names of these individuals that are part of the dubious trios are hard coded into the R script 
process make_plink_instruction_files {
	publishDir "./${outdir}", mode: 'copy', pattern: 'candidate_trios_*'
	input:
		file make_plink_instruction_swaps_src
		file make_plink_instruction_src
		file tidy_wiki_birdlist
		file(vcf) from vcf_filtered_ch2
		founders.flatten()
	output:
		file("family_*_keep_info.txt") into keep_info_ch
		file("family_*_id_info.txt") into id_info_ch
		file("family_*_sex_info.txt") into sex_info_ch
		file("family_*_parent_info.txt") into parent_info_ch
		file("candidate_trios_${vcf.baseName}.txt") into all_run_trios_ch
		file("make_plink_instruction_files.log") into make_plink_instruction_files_ch
	script:
		if(use == "dubious_ped_alts")
			"""
			grep "CHROM" -m 1 ${vcf} | cut -f 10- > seq_birds.txt
			Rscript ${make_plink_instruction_src}  ${tidy_wiki_birdlist} seq_birds.txt candidate_trios_${vcf.baseName}.txt \
			chick_in_pedigree ${founders} \
			> make_plink_instruction_files.log
			"""
		else if(use == "swaps")
			"""
			grep "CHROM" -m 1 ${vcf} | cut -f 10- > seq_birds.txt
			Rscript ${make_plink_instruction_swaps_src}  ${tidy_wiki_birdlist} seq_birds.txt candidate_trios_${vcf.baseName}.txt \
			> make_plink_instruction_files.log
			"""	

    	else
        	error "Unknown mode: ${use}. Use one of dubious_ped_alts or swaps."

}

make_plink_instruction_files_ch.collectFile(storeDir: "${outdir}/logs/")
all_run_trios_ch.collectFile(storeDir: "${outdir}/")

// Update individual and family IDs and output a ped file for each trio, containing only those individuals
process update_ped_ids {
	//publishDir "${outdir}", mode: 'copy'
	executor 'local'
	

	input:
		file(keep_info) from keep_info_ch.flatMap().toSortedList().flatMap() 
		file(id_info) from id_info_ch.flatMap().toSortedList().flatMap()
		file(ped) from ped_file_ch
		file(map) from map_file_ch
	output:
		file("${keep_info.baseName}.ped") into ped_ids
		file("${keep_info.baseName}.map") into map_ids
		file("update_ped_ids.log") into update_ped_ids_ch

	"""
	plink --file ${ped.baseName} \
	--update-ids ${id_info} \
	--keep ${keep_info} \
	--allow-extra-chr \
	--recode --out ${keep_info.baseName} > update_ped_ids.log
	"""
}

update_ped_ids_ch.collectFile(storeDir: "${outdir}/logs/")



// Add information about relationship and sex to the PED file
process update_other_details {
	executor 'local'
	input:
		file(sex_info) from sex_info_ch.flatMap().toSortedList().flatMap()
		file(parent_info) from parent_info_ch.flatMap().toSortedList().flatMap()
		file(ped) from ped_ids.toSortedList {it.baseName}.flatMap()
		file(map) from map_ids.toSortedList {it.baseName}.flatMap()
		
			// file(map) from ped_ids.toSortedList({ a,b -> a[0] <=> b[0]  })
			
	output:
		file("${sex_info.baseName}_all.ped") into ped_all_info
		file("${sex_info.baseName}_all.map") into map_all_info
		file("update_other_details.log") into update_other_details_log_ch
		

	"""
	plink --file ${ped.baseName} \
	--update-parents ${parent_info} \
	--update-sex ${sex_info} \
	--allow-extra-chr \
	--recode \
	--out ${sex_info.baseName}_all > update_other_details.log
	"""
}

update_other_details_log_ch.collectFile(storeDir: "${outdir}/logs/")

// Calculate mendelian errors for each trio from the ped file
process calculate_mendelian_errors {
	executor 'local'

	input:
		file(ped_all) from ped_all_info
		file(map_all) from map_all_info
	output:
		file("${ped_all.baseName}.fmendel") into fmendel_ch optional true
		file("${ped_all.baseName}.imendel") optional true
		file("${ped_all.baseName}.lmendel") optional true
		file("${ped_all.baseName}.mendel") optional true
		file("mendel.log") into mendel_log_ch optional true

	"""
	plink --file ${ped_all.baseName} \
	--allow-extra-chr \
	--mendel \
	--out ${ped_all.baseName} > mendel.log
	"""
}

mendel_log_ch.collectFile(storeDir: "${outdir}/logs/")

// grep "SNPs with no founder genotypes" mendel.log | head -n 5000 > set_to_zero.txt


// Merge the separate fmendel outputs for each trio into one file
process merge_fmendel_outputs {
	executor 'local'

	input:
		each file(fmendel) from fmendel_ch
	output:
		file("all_candidates_fmendel.txt") into all_candidates_fmendel_ch

	"""
	sed '2q;d' ${fmendel} > all_candidates_fmendel.txt
	"""
}


all_candidates_fmendel_ch.collectFile(storeDir: "${outdir}")

