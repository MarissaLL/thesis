#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */

// Inputs and parameters
new_var_vcf = file("/media/drive_6tb/projects/kakapo-genomics/data_stored/gwas_set.vcf")
old_var_vcf = file("/media/drive_6tb/projects/kakapo-genomics/data_stored/deepvariant.r1.vcf")
bird_details = file("/media/drive_6tb/projects/kakapo-genomics/output/37_sequoia_pt2/birds.txt")

outdir = "output/37_sequoia_pt2/"
sampling_size = Channel.from(100,200,300,400,500,600,700,800,900,1000)
replicate = Channel.from(1,2,3,4,5,6,7,8,9,10)


//Scripts
shuf_src = file("src/shuffle.awk")
sequoia_src = file("src/sequoia_parentage.R")


//Main
vcfs = Channel.from(new_var_vcf, old_var_vcf)
sampling_size.combine(replicate)
	.map{size, replicate -> tuple(size, replicate, size+6)}
	.set{size_x_replicate} 




process filter_vcfs {
	input:
		file(vcf) from vcfs
	output:
		file("${vcf.simpleName}.filtered.vcf") into filtered_vcfs
	"""
	bcftools view -i 'F_MISSING<0.1' --min-af 0.3:minor ${vcf} > ${vcf.simpleName}.filtered.vcf
	"""
}


process convert_to_ped {
	publishDir "${outdir}", mode: "copy"
	conda 'bioconda::plink'

	input:
		file filtered_vcf from filtered_vcfs 
	output:
		file("${filtered_vcf.simpleName}.raw") into full_raw
	"""
	sed 's/_/!/g' ${filtered_vcf} > no_underscores.vcf 
	plink --vcf no_underscores.vcf --indep 50 5 2 --vcf-half-call 'missing' \\
	--chr-set 89 --allow-extra-chr --recode A --out ${filtered_vcf.simpleName}
	sed 's/!/_/g' -i ${filtered_vcf.simpleName}.raw
	"""
}

process subsample {
	publishDir "${outdir}", mode: "copy"
	input:
		each file(raw) from full_raw
		file(shuf_src)
		set val(sampling_size),
			val(replicate),
			val(tot_filesize) from size_x_replicate
	output:
		file("sample_${sampling_size}_rep${replicate}_${raw.baseName}.raw") into ds_raw
	"""
	cat ${raw} | awk -f ${shuf_src} -v ncols=${tot_filesize} > sample_${sampling_size}_rep${replicate}_${raw.simpleName}.raw
	"""

}

process sequoia {
	publishDir "${outdir}", mode: "copy"
	input:
		file(raw) from ds_raw
		file(bird_details)
		file(sequoia_src)
	output:
		file("sequoia_${raw.simpleName}.txt") optional true
	"""
	Rscript ${sequoia_src} ${raw} ${bird_details}
	"""
}

