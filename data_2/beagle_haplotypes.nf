#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */


vcf = file("data_stored/gwas_set.vcf")
plot_unique_haplos_script = file("src/reformat_beagle_haplos.R")
bird_rships = file("data/bird_info_rships_nov20.csv")



process make_beagle_haplos {
	publishDir './output/11_beagle_phasing', mode: 'copy'
	
	input:
		//set val(filtered_chr),
		//	file(vcf) from chr_files_ch
		file(vcf)
	output:
		file("phasing_${vcf.baseName}.vcf.gz") into phased_chr_ch
		file("beagle.log")
		
	"""
	java -jar ~/Downloads/beagle.25Nov19.28d.jar gt=${vcf} out=phasing_${vcf.baseName} > beagle.log
	"""
}


/*
process plot_unique_haplos {
	publishDir './output/11_beagle_phasing', mode: 'copy'

	input:
		file plot_unique_haplos_script
		file(phased_vcf) from phased_chr_ch
		file seq_bird_names
		file bird_rships
	output:
		file("chr_*_beagle_unique_haplos.png")
		file("summary_beagle_*_haplos.txt")
		file("ath_all_*_beagle.txt")
		file("ath_founders_*_beagle.txt")

	"""
	Rscript ${plot_unique_haplos_script} ${phased_vcf} ${seq_bird_names} ${bird_rships}
	"""
}
*/
