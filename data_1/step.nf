#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */

/*NOTES
If using partially this and partially the m_err.sh bash script, remember to change the number of runs of the 
for loop to reflect the number of trios being used 

Don't currently have logs for the vcftools stages - they don't seem to print to stdout */


tidy_wiki_birdlist = file("old_files/formatted_bird_list_manual.csv")
n_variants = Channel.value(10000)
vcf = file("data/deepvariant.r1.vcf")
//vcf = file("data/simple_36.vcf")
contig_map = file("output/02_assess_vcfs/mendelian_errors/contig_map.txt")
outdir = "output/03_mendelian_error_distribution/"


// Filter the vcf for quality > 3 (can't go higher because it is using deepvariant quality scores)
// Filter for loci genotyped in at least 80% of individuals
process quality_missing_filter_vcf {
	publishDir "./${outdir}", mode: 'copy'
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
//Take only n_variants number of SNPs (e.g. 10,000 SNPs), drawn randomly
process convert_vcf_to_ped {
	publishDir "./${outdir}", mode: 'copy'
	conda 'bioconda::vcftools'

	input:
		file vcf_filtered from vcf_filtered_ch
		file contig_map
		val n_variants
	output:
		file("me_dist.ped") into ped_file_ch
		file("me_dist.map") into map_file_ch
		file("subset_snps.txt")
		file("list_of_snps.txt")
		

	"""
	grep -v "#" ${vcf_filtered} | cut -f 1,2 > list_of_snps.txt
	awk -F"\t" '!uniq[\$0]++' list_of_snps.txt > list_of_uniq_snps
	shuf list_of_uniq_snps | head -n $n_variants > subset_snps.txt
	vcftools --vcf ${vcf_filtered} --positions subset_snps.txt --plink \
	--chrom-map ${contig_map} --out me_dist
	"""
}



/*Makes files to update id, parents and sex information in the ped file so plink can run mendelian error tests*/
process make_plink_instruction_files {
	//publishDir "./${outdir}", mode: 'copy'
	input:
		file tidy_wiki_birdlist
		file vcf
	output:
		file("family_*_keep_info.txt") into keep_info_ch
		file("family_*_id_info.txt") into id_info_ch
		file("family_*_sex_info.txt") into sex_info_ch
		file("family_*_parent_info.txt") into parent_info_ch
		file("${vcf.baseName}_all_step_trios.txt") into all_run_trios_ch
		file("make_plink_instruction_files.log") into make_plink_instruction_files_ch

	"""
	grep "CHROM" -m 1 ${vcf} | cut -f 10- > seq_birds.txt
	Rscript ${baseDir}/src/stepwise_false_trios.R ${tidy_wiki_birdlist} seq_birds.txt ${vcf.baseName}_all_step_trios.txt \
	> make_plink_instruction_files.log
	"""

}

make_plink_instruction_files_ch.collectFile(storeDir: "${outdir}/logs/")
all_run_trios_ch.collectFile(storeDir:  "${outdir}")

// Update individual and family IDs and output a ped file for each trio, containing only those individuals
process update_ped_ids {
	//publishDir './output/03_mendelian_error_distribution', mode: 'copy'
	executor 'local'
	
	conda 'bioconda::plink'

	input:
		file(keep_info) from keep_info_ch.flatMap().toSortedList().flatMap() 
		file(id_info) from id_info_ch.flatMap().toSortedList().flatMap()
		file(ped) from ped_file_ch
		file(map) from map_file_ch
	output:
		file("me_dist_${keep_info.baseName}.ped") into ped_ids
		file("me_dist_${keep_info.baseName}.map") into map_ids
		file("update_ped_ids.log") into update_ped_ids_ch

	"""
	plink --file ${ped.baseName} \
	--update-ids ${id_info} \
	--keep ${keep_info} \
	--recode --out me_dist_${keep_info.baseName} > update_ped_ids.log
	"""
}

update_ped_ids_ch.collectFile(storeDir: 'output/logs/03_mendelian_error_distribution')



// Add information about relationship and sex to the PED file
process update_other_details {
	//publishDir './output/03_mendelian_error_distribution', mode: 'copy'
	executor 'local'
	
	conda 'bioconda::plink'

	input:
		file(sex_info) from sex_info_ch.flatMap().toSortedList().flatMap()
		file(parent_info) from parent_info_ch.flatMap().toSortedList().flatMap()
		file(ped) from ped_ids.toSortedList {it.baseName}.flatMap()
		file(map) from map_ids.toSortedList {it.baseName}.flatMap()
		
			// file(map) from ped_ids.toSortedList({ a,b -> a[0] <=> b[0]  })
			
	output:
		file("me_dist_${sex_info.baseName}_all.ped") into ped_all_info
		file("me_dist_${sex_info.baseName}_all.map") into map_all_info
		file("update_other_details.log") into update_other_details_log_ch
		

	"""
	plink --file ${ped.baseName} \
	--update-parents ${parent_info} \
	--update-sex ${sex_info} \
	--recode \
	--out me_dist_${sex_info.baseName}_all > update_other_details.log
	"""
}

update_other_details_log_ch.collectFile(storeDir: 'output/logs/03_mendelian_error_distribution')

// Calculate mendelian errors for each trio from the ped file
process calculate_mendelian_errors {
	publishDir './output/03_mendelian_error_distribution', mode: 'copy'
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

mendel_log_ch.collectFile(storeDir: 'output/logs/03_mendelian_error_distribution')

// grep "SNPs with no founder genotypes" mendel.log | head -n 5000 > set_to_zero.txt


// Merge the separate fmendel outputs for each trio into one file
process merge_fmendel_outputs {
	//publishDir './output/03_mendelian_error_distribution', mode: 'copy'
	executor 'local'

	conda 'bioconda::plink'

	input:
		each file(fmendel) from fmendel_ch
	output:
		file("all_stepwise_fmendel.txt") into all_stepwise_fmendel_ch

	"""
	sed '2q;d' ${fmendel} > all_stepwise_fmendel.txt
	"""
}


all_stepwise_fmendel_ch.collectFile(storeDir: 'output/03_mendelian_error_distribution')
      
// Use plot_me_dists.R for plotting the fmendel results (for now, this needs replacing with a tidier one)