#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */


vcf = file("data/deepvariant.r1.vcf")
r_script = file("src/make_true_trio_plink_files.R")
bird_rships = file("data/bird_info_rships.csv")
contig_map = file("scaffold_chrom_map.txt")
all_true_trios_ped = file("output/07_offspring_m_err/all_true_trios.ped")
all_true_trios_map = file("output/07_offspring_m_err/all_true_trios.map")
all_true_trio_parents = file("output/07_offspring_m_err/all_true_trio_parents.txt")

// Depending on system memory, change the upper value here to split the data into smaller pieces for plink
//At the moment running with one piece per trio
num_subsets = Channel.value(121)
data_subsets = Channel.from(1..121) 
family_numbers = Channel.from(1..16)


	// bcftools query --format '%CHROM\t%POS\t[%GQ\t]' ${trio_vcf} > ${vcf.baseName}.GQ.FORMAT 
	//Could also use but sample names are not preserved



// 01 Makes files to update id, parents and sex information in the ped file so plink can run mendelian error tests
// Make the key value pairs thingy here, then carry it through
// Would also be cool to make the script flexible to the number of groups you want to have

process make_plink_instruction_files {
	publishDir './output/07_offspring_m_err/', mode: 'copy'
	input:
		file bird_rships
		file vcf
		file r_script
		val num_subsets
		each num from data_subsets
	output:
		set val(num), file("part_${num}_keep_indivs.txt") into keep_info_ch
		set val(num), file("part_${num}_names.txt") into id_info_ch
		set val(num), file("part_${num}_update_sex.txt") into sex_info_ch
		set val(num), file("part_${num}_parents.txt") into parent_info_ch
		file("make_plink_instruction_files.log")
		file("all_run_trios_part_num.txt") into trios_run_file_ch

	"""
	grep "CHROM" -m 1 ${vcf} | cut -f 10- > seq_birds.txt
	Rscript ${r_script} ${bird_rships} seq_birds.txt ${num_subsets} > make_plink_instruction_files.log
	"""
}



/*file("pt2_plink_keep_indivs.txt"), 
			file("pt3_plink_keep_indivs.txt") into keep_info_ch
		set file("pt1_plink_names.txt"),
			file("pt2_plink_names.txt"),
			file("pt3_plink_names.txt") into id_info_ch
		set file("pt1_plink_update_sex.txt"),
			file("pt2_plink_update_sex.txt"),
			file("pt3_plink_update_sex.txt") into sex_info_ch
		set file("pt1_plink_parents.txt"),
			file("pt2_plink_parents.txt"),
			file("pt3_plink_parents.txt") into parent_info_ch*/



// 02 Make a vcf of only the trio individuals. First started with Hoki,  Zephyr, and Felix, now all 
process select_indivs {
	publishDir './output/07_offspring_m_err', mode: 'copy'
	
	input:
		file vcf
		set val(num), file(keep_indivs) from keep_info_ch
	output:
		set val(num), file("${keep_indivs.baseName}.vcf") into (subset_indivs_vcf_ch, subset_indivs_vcf_ch2, subset_indivs_vcf_ch3)

	"""
	bcftools view --samples-file ${keep_indivs} ${vcf} > ${keep_indivs.baseName}.vcf
	"""
}




//subset_indivs_vcf_ch3.subscribe{println it}



// 03 Extract genotype qualities for all the relevant individuals
process extract_GQ {
	publishDir './output/07_offspring_m_err', mode: 'copy'
	
	input:
		set val(num), file(subset_vcf) from subset_indivs_vcf_ch
	output:
		set val(num), file("${subset_vcf.baseName}.GQ.FORMAT") into gq_file_ch
	"""
	vcftools --vcf ${subset_vcf} --extract-FORMAT-info GQ --out ${subset_vcf.baseName}
	"""
}



// 04 Convert VCF to ped. This requires vcftools 0.1.17
process convert_vcf_to_ped {
	//publishDir './output/07_offspring_m_err', mode: 'copy'
	
	input:
		set val(num), file(vcf_filtered) from subset_indivs_vcf_ch2
		file contig_map
	output:
		set val(num), file("${vcf_filtered.baseName}.ped") into ped_file_ch
		set val(num), file("${vcf_filtered.baseName}.map") into map_file_ch
		
	"""
	vcftools --vcf ${vcf_filtered} --plink \
	--chrom-map ${contig_map} --out ${vcf_filtered.baseName}
	"""
}

id_info_ch.join(ped_file_ch).join(map_file_ch).into{update_ped_ids_ch; update_ped_ids_ch2}

//update_ped_ids_ch2.subscribe{println it}




// 05 Update individual and family IDs and output a ped file for each trio, containing only those individuals
process update_ped_ids {
	//publishDir './output/07_offspring_m_err', mode: 'copy'
	executor 'local'
	memory = '16 GB'
	
	input:
		set val(num), 
			file(id_info),
			file(ped),
			file(map) from update_ped_ids_ch
	output:
		set val(num), file("${ped.baseName}_ids.ped") into ped_ids
		set val(num), file("${ped.baseName}_ids.map") into map_ids
		file("update_ped_ids_${num}.log") into update_ped_ids_log_ch

	"""
	plink --file ${ped.baseName} \
	--update-ids ${id_info} \
	--recode --out ${ped.baseName}_ids > update_ped_ids_${num}.log
	"""
}

update_ped_ids_log_ch.collectFile(storeDir: 'output/logs/07_offspring_m_err')

sex_info_ch.join(parent_info_ch).join(ped_ids).join(map_ids).into{update_other_details_ch; update_other_details_ch2}

update_other_details_ch2.subscribe{println it}


// 06 Add information about relationship and sex to the PED file
process update_other_details {
	//publishDir './output/07_offspring_m_err', mode: 'copy'
	executor 'local'
	memory = '16 GB'
	
	input:
		set val(num),
			file(sex_info),
			file(parent_info),
			file(ped),
			file(map) from update_other_details_ch
	output:
		set val(num), file("all_details_${num}.ped") into ped_all_info
		set val(num), file("all_details_${num}.map") into map_all_info
		file("update_other_details.log") into update_other_details_log_ch
		

	"""
	plink --file ${ped.baseName} \
	--update-parents ${parent_info} \
	--update-sex ${sex_info} \
	--recode \
	--out all_details_${num} > update_other_details.log
	"""
}

update_other_details_log_ch.collectFile(storeDir: 'output/logs/07_offspring_m_err')

ped_all_info.join(map_all_info).into{calculate_mendelian_errors_ch; calculate_mendelian_errors_ch2}

//calculate_mendelian_errors_ch2.subscribe{println it}


// 07 Calculate mendelian errors for each trio from the ped file
process calculate_mendelian_errors {
	publishDir './output/07_offspring_m_err', mode: 'copy'
	executor 'local'
	memory = '16 GB'

	conda 'bioconda::plink'

	input:
		set val(num),
			file(ped_all),
			file(map_all) from calculate_mendelian_errors_ch
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

mendel_log_ch.collectFile(storeDir: 'output/logs/07_offspring_m_err')

/*

// 08 Run the mega-script
// Need to pass it the mendel files and vcf to get it working
process mega_script {
	publishDir './output/mega_script/', mode: 'copy'

	input:
		file all_true_trio_parents
		file(all_trios_file) from trios_run_file_ch
		each family_num from family_numbers
	output: 
		file("chick_error_attributions_genotypes_family_${family_num}.txt")
	"""
	Rscript ~/projects/kakapo-genomics/src/find_multisib_me.R ${all_true_trio_parents} ${all_trios_file} ${family_num}
	"""
}






/*
// 05 Update individual and family IDs and output a ped file for each trio, containing only those individuals
process update_ped_ids {
	//publishDir './output/07_offspring_m_err', mode: 'copy'
	executor 'local'
	
	input:
		set val(num), file(id_info) from id_info_ch
		set val(num), file(ped) from ped_file_ch
		set val(num), file(map) from map_file_ch
	output:
		file("${ped.baseName}_ids.ped") into ped_ids
		file("${ped.baseName}_ids.map") into map_ids
		file("update_ped_ids_${num}.log") into update_ped_ids_ch

	"""
	plink --file ${ped.baseName} \
	--update-ids ${id_info} \
	--recode --out ${ped.baseName}_ids > update_ped_ids_${num}.log
	"""
}

//--allow-extra-chr \


update_ped_ids_log_ch.collectFile(storeDir: 'output/logs/07_offspring_m_err')

// 06 Add information about relationship and sex to the PED file
process update_other_details {
	//publishDir './output/07_offspring_m_err', mode: 'copy'
	executor 'local'
	input:
		file(sex_info) from sex_info_ch
		file(parent_info) from parent_info_ch
		file(ped) from ped_ids
		file(map) from map_ids
	output:
		file("all_true_trios_info.ped") into ped_all_info
		file("all_true_trios_info.map") into map_all_info
		file("update_other_details.log") into update_other_details_log_ch
		

	"""
	 plink --file ${ped.baseName} \
	--update-parents ${parent_info} \
	--update-sex ${sex_info} \
	--allow-extra-chr \
	--recode \
	--out all_true_trios_info > update_other_details.log
	"""
}

update_other_details_log_ch.collectFile(storeDir: 'output/logs/07_offspring_m_err')

// 07 Calculate mendelian errors for each trio from the ped file
process calculate_mendelian_errors {
	publishDir './output/07_offspring_m_err', mode: 'copy'
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
	--allow-extra-chr \
	--mendel \
	--out ${ped_all.baseName} > mendel.log
	"""
}

mendel_log_ch.collectFile(storeDir: 'output/logs/07_offspring_m_err')
*/