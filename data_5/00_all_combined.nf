/*
marissa.lelec@gmail.com
github.com/MarissaLL
*/

//requires --sample, --dir and --parents_dir  

samples = Channel.from("${params.sample}" )
samples.map{samples -> tuple("${params.sample}", file("${params.dir}/${params.sample}.fastq"))}. set{fastq}
fastq.into{fastq_ch; fastq_ch2}
ref = file("data/NCBI_New_Reference.fasta")
ref_fai = file("data/NCBI_New_Reference.fasta.fai")
run_name = "${params.dir}/processed"
main_ds_1_indiv = file('data/full_ds_single_indiv.vcf.gz')
main_ds_1_indiv_tbi = file('data/full_ds_single_indiv.vcf.gz.tbi')
main_ds_pos = file('data/snp_pos_full_ds.txt')


parent_combos = Channel.fromPath("/nesi/nobackup/uoo03425/${params.parents_dir}/*.ped")
parent_combos.map{ped -> tuple(ped.baseName.split('__')[0], 
				ped.baseName.split('__')[1],
				ped)}.set{labelled_parent_combos}
wrangle_fmendels_src = file("src/wrangle_fmendel.R")
sex_test_src = file("src/assign_sex3.R")

train_data_file = file('data/training_data.tsv')
classifier_script = file("src/calc_parentage_probs_gaussian_bayes.py")



process minimap {
	publishDir "./${run_name}", mode: 'copy'
	cpus 4
  	queue 'large'
  	memory '24 GB' 
  	time '4h'
	module 'minimap2'
        input:
		set val(file_name),
			file(chopped_fastq) from fastq_ch
		file(ref)
	output:
		file("${file_name}.sam") into sam_ch
	"""
	minimap2 -ax map-ont ${ref} ${chopped_fastq} -t 10 > ${file_name}.sam
	"""
}

process sam_to_bam {
	cpus 4
  	queue 'large'
  	memory '8 GB' 
  	time '2h'
   	module 'SAMtools'
        input:
		file(sam) from sam_ch
	output:
		file("${sam.baseName}.bam") into unsorted_bam_ch
	"""
	samtools view -S -b ${sam} > ${sam.baseName}.bam
	"""
}

process sort_and_index_bam {
	cpus 4
  	queue 'large'
  	memory '8 GB' 
  	time '2h'
	module 'SAMtools'
	publishDir "./${run_name}", mode: 'copy'
        input:
		file(bam) from unsorted_bam_ch
	output:
		set val("${bam.baseName}"),
			file("${bam.baseName}_sorted.bam"),
			file("${bam.baseName}_sorted.bam.bai") into (idx_bam_ch, idx_bam_ch2) 
	"""
	samtools sort ${bam} -o ${bam.baseName}_sorted.bam
	samtools index ${bam.baseName}_sorted.bam
	"""
}


process idxstats {
        module 'SAMtools'
        input:
                set val(bam_name),
			file(bam),
			file(bai) from idx_bam_ch2
        output:
                file("${bam.baseName}.idxstats") into idxstats_ch
        
        """
        samtools idxstats ${bam} > ${bam.baseName}.idxstats
        """
}

process sex_test {
        module 'R/4.1.0-gimkl-2020a'
        input:
                file(idxstats) from idxstats_ch
                file(sex_test_src)
        output:
                file("coverage_stats.txt")
                file("assigned_sex.txt")
		stdout sex_res_ch
        """
        Rscript ${sex_test_src} 
	cat assigned_sex.txt | cut -f2 | tail -n 1 | tr -d '\n'
        """
}


process longshot {
        cpus 10
  	queue 'large'
  	memory '16 GB' 
  	time '8h'
	module 'Singularity'
	publishDir "./${run_name}", mode: 'copy'
        input:
		set val(file_name),
			file(bam),
			file(bai) from idx_bam_ch
		file(ref)
		file(ref_fai)
	output:
		set val("${file_name}"),
			file("${file_name}_sample.vcf") into vcf_ch
	"""
	longshot --bam ${bam} --ref ${ref} --strand_bias_pvalue_cutoff 0.01 --out ${file_name}_sample.vcf
	"""
}
process rename_samples_in_vcf {
        cpus 2
        queue 'large'
        memory '4 GB' 
        time '0.1h'
        input:
                set val(file_name),
                        file(vcf) from vcf_ch
        output:
                set val("${file_name}"),
			file("${file_name}.vcf") into (np_vcf_ch, np_vcf_ch2)
        """
        samplename=\$(echo ${vcf} | cut -d'_' -f 1)
        sed "s/SAMPLE/\$samplename/g" ${vcf} > ${file_name}.vcf
        """
} 

process bgzip {
	module 'tabix'
	input:
		set val(name),
			file(vcf) from np_vcf_ch.collect()
	output:
		set val("${name}"),
			file("${vcf}.gz"),
			file("${vcf}.gz.tbi") into np_gz_ch
	"""
	bgzip ${vcf}
	tabix ${vcf}.gz
	"""
}


process add_missing_positions {
	module 'BCFtools'
	maxForks 1
	input:
		file(main_ds_1_indiv)
		file(main_ds_1_indiv_tbi)
		file(main_ds_pos)
		set val(sample),
			file(nnp_vcfgz),
			file(nnp_tbi) from np_gz_ch
	output:
		set val("${sample}"),
			file("main_shared_${sample}.vcf") into plinku
		
	"""
	bcftools merge ${main_ds_1_indiv} ${nnp_vcfgz} --output ${sample}_allpos.vcf
	bgzip ${sample}_allpos.vcf
	tabix ${sample}_allpos.vcf.gz
	bcftools view -R ${main_ds_pos} --samples ^Adelaide ${sample}_allpos.vcf.gz > main_shared_${sample}.vcf
	"""
}
// Converting from and back to underscores because otherwise plink incorrectly
// splits names into FID and IID at the underscores
process convert_to_ped {
        module 'PLINK/1.09b6.16'
	maxForks 1
	publishDir './output/plink_files', mode: 'copy'
	input:
                set val(key),
                        file(vcf) from plinku
        output:
                set val("${key}"),
                        file("shared_${key}.ped"),
                        file("shared_${key}.map") into ped_original_ch
        """
        sed 's/_/!/g' ${vcf} > no_underscores.vcf 
        plink --vcf no_underscores.vcf --allow-extra-chr --vcf-half-call 'missing' --recode --out shared_${key}
        sed 's/!/_/g' -i shared_${key}.ped 
        sed 's/!/_/g' -i shared_${key}.map 
        """
}



labelled_parent_combos.combine(ped_original_ch).set{to_concat_chick_ch}


process cat_chick {
	maxForks 500
	queue 'large'
        cpus 1
        time '0.1h'
        memory '0.8GB'
	input:
		set val(father),
			val(mother),
			file(parent_combos), 
			val(name),
			file(ped),
			file(map) from to_concat_chick_ch
	output:
		set val("${name}"),
			file("${name}_${parent_combos.baseName}.ped"),
			file("${map}") into plink_check_ch
	"""
	cut -f7- -d" " ${ped} > genos.txt
	echo "kakapo_parentage ${name} ${father} ${mother} 0 -9" > info.txt
	paste -d" " info.txt genos.txt > edited_ped.txt
	cat edited_ped.txt ${parent_combos} > ${name}_${parent_combos.baseName}.ped
	"""
}	

process mendel {
	module 'PLINK/1.09b6.16'
	maxForks 100
        queue 'large'
        cpus 1
	time '0.1h'
        memory '16GB'
	errorStrategy 'ignore'
	//publishDir './fmendels/', mode: 'move'
	input:
              	set val(key),
                        file(ped),
                        file(map) from plink_check_ch
        output:
               	file("${ped.baseName}.fmendel")
		file("${ped.baseName}.log") into mendel_log_ch
        """
	mv ${map} ${ped.baseName}.map
	plink --file ${ped.baseName} --mendel --aec --out ${ped.baseName}
        """
}

process combine_results {
	queue 'large'
	cpus 4
	time '0.1h'
	memory '2GB'
	publishDir './fmendels/', mode: 'copy'
	input:
		params.sample
		file(logs) from mendel_log_ch.collect()
	output:
		file("results_${params.sample}.fmendel") into fmendels_ch
	"""
	grep "rate" -A 1 *.log > results_${params.sample}.fmendel
	"""
}

process wrangle_fmendels {
	queue 'large'
        cpus 4
	time '0.1h'
        memory '2GB'
        module 'R'
	publishDir './fmendels/', mode: 'copy'
        input:
              	file(fmendel) from fmendels_ch
		file(wrangle_fmendels_src)
        output:
		file("${fmendel.baseName}.txt") into putative_parents_ch
	"""
	Rscript ${wrangle_fmendels_src} ${fmendel}
	"""
}

process get_putative_parents {
	input:
		file(putative_parents) from putative_parents_ch
	output:
		stdout parents_ch 
	"""
	head -n 2 ${putative_parents} | tail -n 1 | cut -f 1 | tr -d '\n'
	"""
}

process get_passed_snps {
	input:
		set val(name),
			file(vcf) from np_vcf_ch2
	output:
		stdout passed_snps_ch
	"""
	grep "PASS" ${vcf} | wc -l | tr -d '\n'
	"""
}

process get_yield {
	input:
		set val(name),
			file(fastq) from fastq_ch2
	output:
		stdout yield_ch
	"""
	cat ${fastq} | paste - - - - | cut -f 2 | tr -d '\n' | wc -c | tr -d '\n'
	"""
}

yield_ch.combine(passed_snps_ch).set{predictors}
predictors.into{predictors_ch; predictors_ch2}

process combine {
	input:
		set val(yield),
			val(passed_snps) from predictors_ch
	output:
		file("file.test") into for_classification
	"""
	printf "%s \n" "${yield},${passed_snps}\n0,0" > file.test
	"""
}

process classify_naive_bayes {
	input:
		file(res) from for_classification
		file(classifier_script)
		file(train_data_file)
	output:
		file("probs.txt") into probability_ch
	"""
	python3 ${classifier_script} --test_file ${res} --training_file ${train_data_file}
	head -n 1 prob_list.txt > probs.txt
	"""
}		

parents_ch.combine(probability_ch).set{final_res_ch}

process produce_output {
	publishDir './results/', mode: 'copy'
	input:
		set val(parents),
			file(probs) from final_res_ch
		set val(yield),
			val(passed_snps) from predictors_ch2
		val(sex) from sex_res_ch
	output:
		file("results_${params.sample}.txt")
	"""
	echo "Parents Sex Sequencing_yield SNPs_filtered Assignment_probability" > results_${params.sample}.txt
	printf "%s " "${parents} ${sex} ${yield} ${passed_snps}" >> results_${params.sample}.txt
	cat ${probs} >> results_${params.sample}.txt
	"""
}
