#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */
//Channel.from(1,2,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26).into{ chr; chr2; chr3}

n_vcfs = 4

Channel.from (file('output/25_tidy_alphapeel/haplotypes/alphapeel_phased_allchr_ordered.vcf'), file('output/11_beagle_phasing/phased_gwas_set_autosomes.vcf'),
	file('output/10_fimpute_phasing/fimpute_phased_allchr_ordered.vcf'), file('output/13_whatshap_phasing/whatshap_phased_allchr.vcf'))
	.set{vcfs}


// chr.map{chr -> tuple(file("whatshap_nesi/S${chr}.phased.vcf"), 
// 					 file("output/25_tidy_alphapeel/haplotypes/alphapeel_phased_S${chr}.vcf"),
// 					 file("output/10_fimpute_phasing/fimpute_phased_S${chr}.vcf"))}
// 	.flatten()
// 	.set{vcfs}

parse_wc_script = file("src/parse_whatshap_compare_files.R")

indiv = Channel.from( "Adelaide")

// , "Al", "Alice", "Aparima", "Arab", "Aranga", "Ariki", 
// "Atareta", "Attenborough", "Aumaria", "Awarua", "Awhero", "Barnard", 
// "Basil", "Bella", "Ben", "Bill", "Blades", "Blake", "Bluster_Murphy", 
// "Bonus", "Boomer", "Boss", "Clout", "Cyndy", "Dobbie", "Doc", 
// "Dusky", "Egilsay", "Ellie", "Elliott", "Elwin", "Esperance", 
// "Evohe", "Faulkner", "Felix", "Flossie", "Fuchsia", "George", 
// "Gertrude", "Guapo", "Gulliver", "Gumboots", "Gunner", "Hakatere", 
// "Hananui", "Hauturu", "Heather", "Henry", "Hera", "Hillary", 
// "Hine_Taumai", "Hinemoa", "Hoki", "Hokonui", "Horton", "Hugh", 
// "Huhana", "Hurihuri", "Ian", "Ihi", "JEM", "Jack", "Jamieson", 
// "Jean", "Jemma", "Jester", "Jimmy", "Joe", "Juanma", "Kanawera", 
// "Kokoto", "Komaru", "Konini", "Konini_3-4-16", "Kotiu", "Kuia", 
// "Kuihi", "Kumi", "Kura", "Lee", "Lionel", "Lisa", "Luke", "Maggie", 
// "Mahli", "Makorea", "Manu", "Marama", "Margaret-Maree", "Marian", 
// "Matangi", "Merty", "Merv", "Mila", "Millie", "Moorhouse", "Morehu", 
// "Moss", "Mukeke", "Ngatapa", "Ninihi", "Nog", "Nora", "Oraka", 
// "Ox", "Paddy", "Palmersan", "Pearl", "Percy", "Piripi", "Pounamu", 
// "Pura", "Queenie", "Ra", "Rakiura", "Ralph", "Rangi", "Richard_Henry", 
// "Rimu", "Robbie", "Roha", "Ruapuke", "Ruggedy", "Ruth", "Sandra", 
// "Sarah", "Sass", "Scratch", "Sinbad", "Sirocco", "Smoko", "Solstice", 
// "Stella", "Stumpy", "Sue", "Suzanne", "Taeatanga", "Takitimu", 
// "Tamahou", "Tau_Kuhurangi", "Te_Atapo", "Te_Awa", "Te_Here", 
// "Te_Kingi", "Tia", "Tiaka", "Titapu", "Tiwai", "Tiwhiri", "Tohu", 
// "Toitiiti", "Trevor", "Tukaha", "Tumeke", "Tuterangi", "Tutoko", 
// "Waa", "Waihopai", "Waikawa", "Waynebo", "Weheruatanga-o-te-po", 
// "Wendy", "Wharetutu", "Whiskas", "Wiremu", "Wolf", "Yasmine", 
// "Zephyr")

// indiv = Channel.from( "Adelaide", "Al")


process bgzip_vcfs {
	input:
		each file(vcf) from vcfs
	output:
		file("${vcf}.gz") into bgzipped_ch
	"""
	bgzip ${vcf}
	"""
}

process tabix_vcfs {
	input:
		file(vcf_gz) from bgzipped_ch
	output:
		file(vcf_gz) into gz_ch
		file("${vcf_gz}.tbi") into tabix_ch
	"""
	tabix ${vcf_gz}
	"""
}


process compare_chr {
	publishDir "./output/29_whatshap_compare/by_indiv_${chr_num}/", mode: 'copy'
	maxForks 1
	// tag { "${indiv} ${chr_num}" }
	input:
		// val(chr_num) from chr2.first()
		each indiv
		file(vcf_gz) from gz_ch.collect()
		file(tbi) from tabix_ch.collect()
	output:
		file("whatshap_compare_${indiv}.txt") into whatshap_compare_ch
	"""
	whatshap compare ${vcf_gz} --sample ${indiv} --tsv-multiway whatshap_compare_${indiv}.txt
	"""
}

/*
process combine_files {
	publishDir './output/29_whatshap_compare/', mode: 'copy'
	input:
		val n_vcfs
		// val(chr_num) from chr3.first()
		file parse_wc_script
		file(whatshap) from whatshap_compare_ch.collect()
	output:
		file("all_whatshap_compare_S${chr_num}.txt")
	"""
	echo ${whatshap} > filelist.txt
	Rscript ${parse_wc_script} filelist.txt ${n_vcfs}
	mv all_whatshap_compare.txt all_whatshap_compare_S${chr_num}.txt
	"""
}
*/