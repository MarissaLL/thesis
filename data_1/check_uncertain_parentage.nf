#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */

/*NOTES
This is intended to be run after step.nf - before that it wasn't known that these four needed checking.
Turns out there might actually be five - once the spelling is corrected so that it is included, 
 Elliott also seems to have errors
 Also used this script to check the no assigned paternity chicks. need to comment and uncomment lines in the R
 script to run it on different individuals

 */


//tidy_wiki_birdlist = file("data/konini.csv")
// tidy_wiki_birdlist = file("data/bird_info_rships.csv")
tidy_wiki_birdlist = file("data/bird_info_rships_nov20.csv")
n_variants = Channel.value(200000)
vcf = file("data_stored/gwas_set_numeric_maf_0.05_f_missing_0.2.bcf")
// vcf = file("data/deepvariant.r1.vcf")
//vcf = file("data/simple_36.vcf")
//contig_map = file("output/02_assess_vcfs/mendelian_errors/contig_map.txt")
contig_map = file("contig_map.txt")


process quality_missing_filter_vcf {
	//publishDir './output/03_mendelian_error_distribution', mode: 'copy'
	conda 'bioconda::vcftools'

	input:
		file vcf
	output:
		file("qual3_missing_filter.recode.vcf") into vcf_filtered_ch
		
	"""
	vcftools --vcf ${vcf} --minQ 3 --max-missing 0.8 --recode \
	--out qual3_missing_filter 
	"""
}

// Make the vcf into a ped file for plink to use. 
//Take only n_variants number of SNPs (e.g. 200,000 SNPs), drawn randomly
process convert_vcf_to_ped {
	//publishDir './output/04_test_dubious_trios', mode: 'copy'
	conda 'bioconda::vcftools'

	input:
		file(q3_vcf) from vcf_filtered_ch
		file contig_map
		val n_variants
	output:
		file("me_test.ped") into ped_file_ch
		file("me_test.map") into map_file_ch
		file("subset_snps.txt")
		file("list_of_snps.txt")
		

	"""
	grep -v "#" ${q3_vcf} | cut -f 1,2 > list_of_snps.txt
	awk -F"\t" '!uniq[\$0]++' list_of_snps.txt > list_of_uniq_snps
	shuf list_of_uniq_snps | head -n $n_variants > subset_snps.txt
	vcftools --vcf ${q3_vcf} --positions subset_snps.txt --plink \
	--chrom-map ${contig_map} --out me_test
	"""
}



/*Makes files to update id, parents and sex information in the ped file so plink can run mendelian error tests
Currently the names of these individuals that are part of the dubious trios are hard coded into the R script */
process make_plink_instruction_files_konini_16 {
	//publishDir './output/04_test_dubious_trios', mode: 'copy'
	input:
		file tidy_wiki_birdlist
		file vcf
	output:
		file("family_*_keep_info_2020.txt") into keep_info_ch
		file("family_*_id_info_2020.txt") into id_info_ch
		file("family_*_sex_info_2020.txt") into sex_info_ch
		file("family_*_parent_info_2020.txt") into parent_info_ch
		file("candidate_trios_${vcf.baseName}_2020.txt") into all_run_trios_ch
		file("make_plink_instruction_files.log") into make_plink_instruction_files_ch

	"""
	grep "CHROM" -m 1 ${vcf} | cut -f 10- > seq_birds.txt
	Rscript ${baseDir}/src/dubious_ped_alts.R ${tidy_wiki_birdlist} seq_birds.txt candidate_trios_${vcf.baseName}_2020.txt \
	> make_plink_instruction_files.log
	"""

}

make_plink_instruction_files_ch.collectFile(storeDir: 'output/logs/041_test_dubious_trios_2020')
all_run_trios_ch.collectFile(storeDir: 'output/041_test_dubious_trios_2020')

// Update individual and family IDs and output a ped file for each trio, containing only those individuals
process update_ped_ids {
	//publishDir './output/04_test_dubious_trios', mode: 'copy'
	executor 'local'
	
	conda 'bioconda::plink'

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
	--recode --out ${keep_info.baseName} > update_ped_ids.log
	"""
}

update_ped_ids_ch.collectFile(storeDir: 'output/logs/041_test_dubious_trios_2020')



// Add information about relationship and sex to the PED file
process update_other_details {
	//publishDir './output/04_test_dubious_trios', mode: 'copy'
	executor 'local'
	
	conda 'bioconda::plink'

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
	--recode \
	--out ${sex_info.baseName}_all > update_other_details.log
	"""
}

update_other_details_log_ch.collectFile(storeDir: 'output/logs/041_test_dubious_trios_2020')

// Calculate mendelian errors for each trio from the ped file
process calculate_mendelian_errors {
	//publishDir './output/04_test_dubious_trios', mode: 'copy'
	executor 'local'

	conda 'bioconda::plink'

	input:
		file(ped_all) from ped_all_info
		file(map_all) from map_all_info
	output:
		file("${ped_all.baseName}.fmendel") into fmendel_ch
		file("${ped_all.baseName}.imendel")
		file("${ped_all.baseName}.lmendel")
		file("${ped_all.baseName}.mendel")
		file("mendel.log") into mendel_log_ch

	"""
	plink --file ${ped_all.baseName} \
	--mendel \
	--out ${ped_all.baseName} > mendel.log
	"""
}

mendel_log_ch.collectFile(storeDir: 'output/logs/041_test_dubious_trios_2020')

// grep "SNPs with no founder genotypes" mendel.log | head -n 5000 > set_to_zero.txt


// Merge the separate fmendel outputs for each trio into one file
process merge_fmendel_outputs {
	//publishDir './output/04_test_dubious_trios', mode: 'copy'
	executor 'local'

	conda 'bioconda::plink'

	input:
		each file(fmendel) from fmendel_ch
	output:
		file("all_candidates_fmendel.txt") into all_candidates_fmendel_ch

	"""
	sed '2q;d' ${fmendel} > all_candidates_fmendel.txt
	"""
}


all_candidates_fmendel_ch.collectFile(storeDir: 'output/041_test_dubious_trios_2020')
      
// Use plot_me_dists.R for plotting the fmendel results 