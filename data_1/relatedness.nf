#!/usr/bin/env nextflow


/*
 * marissa.lelec@gmail.com
 * github.com/MarissaLL
 */

vcf_ch = Channel.from(tuple('original', 'vcf', file("data_stored/deepvariant.r1.vcf")),
                                          tuple('new', 'vcf', file("data_stored/gwas_set.vcf")))


process filter_dataset {
	input:
		set val(dataset),
			val(type),
			file(unfiltered) from vcf_ch
	output:
		set val(dataset),
			file("${unfiltered.baseName}_intermediate.bcf") into filtered_bcf
	"""
	bcftools view -v snps -m2 -M2 -i 'QUAL>40 & FMT/GQ>10 & FMT/DP>3 & FMT/DP<100 & F_MISSING<0.1' -q 0.05 ${unfiltered} -Ob -o ${unfiltered.baseName}_intermediate.bcf
	"""
}

process ld_prune {
        input:
		set val(dataset),
			file(bcf) from filtered_bcf
        output:
		set val(dataset),
			file("${bcf.baseName}_ld.bcf") into (ld_filtered_bcf, ld_filtered_bcf2, ld_filtered_bcf3)
        """
	bcftools +prune --max 0.8 -w 1000 ${bcf} -Ob -o ${bcf.baseName}_ld.bcf
	"""
}

process index {
	input:
		set val(dataset),
			file(bcf) from ld_filtered_bcf
	output:
		set file("${bcf}"),
			file("${bcf}.csi") into (bcf_for_rel, bcf_for_rel2)
	"""
	bcftools index ${bcf}
	"""
}

process get_sample_names {
	publishDir 'output/28_relatedness', mode: 'copy'
	input:
		set file(bcf),
			file(csi) from bcf_for_rel2
	output:
		 file("${bcf.baseName}.samplenames")
	"""
	bcftools query -l ${bcf} > ${bcf.baseName}.samplenames
	"""
}	

process calc_rxy_rel_matrix {
	publishDir 'output/28_relatedness', mode: 'copy'

	input:
		set file(bcf),
			file(csi) from bcf_for_rel
	output:
		file("${bcf.baseName}.rxy")
	"""
	ngsRelate -h ${bcf} -O ${bcf.baseName}.rxy
	"""
}


process make_chr_set {
        input:
                set version,
                    file(variants) from ld_filtered_bcf2
        output:
                set version,
                        stdout,
                        file(variants) into chrset_ch
        when:
                version == 'new'
        """ 
        bcftools view ${variants} | head -n 400 | grep "##contig" | wc -l 
        """
}


process bcf_to_plink2_aec {
        input:
                set version,
			file(bcf) from ld_filtered_bcf3
        output:
                set file("${version}.psam"),
                        file("${version}.pvar"),
                        file("${version}.pgen") into plink2_format_ch1
	when:
		version == "original"
        """
        plink2 --bcf ${bcf} -aec --vcf-half-call 'missing' --out ${version}
        """
}

process bcf_to_plink2_chrset {
	input:
		set version,
			chrset,
			file(variants) from chrset_ch
	output:
		set file("${version}.psam"),
			file("${version}.pvar"),
			file("${version}.pgen") into plink2_format_ch2

	"""
	plink2 --bcf ${variants} --out ${version} -aec --vcf-half-call 'missing' --chr-set ${chrset} 
	"""
}

plink2_format_ch1.concat(plink2_format_ch2).set{plink_files_ch}


process make_king_table {
	 publishDir './output/28_relatedness/', mode: 'copy'
        input:
                set file(psam),
                        file(pvar),
                        file(pgen) from plink_files_ch
        output:
                file("${pvar.baseName}.kin0")
                file("${pvar.baseName}.log")

        """
        plink2 --pfile ${pvar.baseName} -aec --make-king-table --out ${pvar.baseName}
        """
}
