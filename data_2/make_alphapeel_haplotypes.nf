#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */

make_vcf_script = file("src/alphapeel_haps_to_vcf.py")

// haps_file = file("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/raw_alphapeel_files/S8_output_ncycles20.haps")
alphapeel_input_vcf = file("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/snps_used_for_alphapeel.vcf.gz")
chrs = Channel.from(1,2,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26)
cutoff = 0.95

chrs.map{chr -> tuple("$chr", file("/media/drive_6tb/projects/kakapo-genomics/output/21_alphapeel/raw_alphapeel_files/S${chr}_output_ncycles20.haps"))}
	.set{haps_ch}

// Tabix the main vcf
process tabix_vcf {
	input:
		file alphapeel_input_vcf
	output:
		set file(alphapeel_input_vcf),
			file("${alphapeel_input_vcf}.tbi") into vcf_ch
	"""
	tabix ${alphapeel_input_vcf}
	"""
}

// Convert the haplotypes from alphapeel into phased vcf format, using the alphapeel_input_vcf as a template to
// convert the genotypes. Produces one vcf per chromosome per individual
process make_indiv_vcfs{
	queuesize = 1
	input:
		file make_vcf_script
		set chr,
			file(haps_file) from haps_ch
		set file(vcf),
			file(tbi) from vcf_ch
		cutoff
	output:
		set chr,
			file("phased_*_S${chr}.vcf") into indiv_vcf_ch mode flatten
	"""
	python3 ${make_vcf_script} ${vcf} S${chr} ${cutoff}
	"""

}


// Bgzip all the vcfs so they can be tabix indexed
process bgzip_indiv_vcfs {
	input:
		set chr,
			file(vcf) from indiv_vcf_ch
	output:
		set chr,
			file("${vcf}.gz") into (indiv_vcf_gz_ch1, indiv_vcf_gz_ch2)
	"""
	bgzip ${vcf}
	"""
}

// Tabix index all the vcfs so that they can be merged back together
process tabix_indiv_vcfs {
	queuesize = 1
	input:
		set chr,
			file(vcf_gz) from indiv_vcf_gz_ch1
	output:
		set chr,
			file("${vcf_gz}.tbi") into indiv_tbi_ch
	"""
	tabix ${vcf_gz}
	"""
}

// Group all the vcf.gz and vcf.gz.tbi files by chromosome so that they can be merged
indiv_vcf_gz_ch2.groupTuple(by:0)
	.set{grouped_vcfs_ch}

indiv_tbi_ch.groupTuple(by:0)
	.set{grouped_tbis_ch}

// Put the vcf.gz and vcf.gz.tbi files in the same channel
grouped_vcfs_ch.join(grouped_tbis_ch)
	.set{all_by_chr_ch}


// Merge all the individuals back together into one vcf per chromosome
process merge_indiv_vcfs {
	publishDir './output/25_tidy_alphapeel/haplotypes/', mode: 'copy'
	queuesize = 1
		input:
			set chr,
				file(vcf_gz), 
				file(tbi) from all_by_chr_ch
		output:
			file("alphapeel_phased_S${chr}.vcf") into per_chr_vcf_ch
	"""
	bcftools merge ${vcf_gz} --output alphapeel_phased_S${chr}.vcf
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

// Different sample names between vcfs. Look into this.
process combine_chromosomes {
	publishDir './output/25_tidy_alphapeel/haplotypes/', mode: 'copy'
	input:
		file(vcf_gz) from bgzip_per_chr_ch2.collect()
		file(tbi) from tbi_per_chr_ch.collect()
	output:
		file("alphapeel_phased_allchr.vcf")
	"""
	bcftools concat ${vcf_gz} --output alphapeel_phased_allchr.vcf
	"""
}
