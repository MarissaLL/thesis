input_convert_src = file("src/make_franz_input.R")
sex_data_file = file("data/bird_info_rships_nov20.csv")
sampling_size = Channel.from(100,200,300,400,500,600,700,800,900,1000)
replicate = Channel.from(1,2,3,4,5,6,7,8,9,10)

sampling_size.combine(replicate)
	.map{size, replicate -> tuple(size, replicate, file("output/37_sequoia_pt2/sample_${size}_rep${replicate}_deepvariant.raw"))}
	.set{sequoia_format_ch}
//sequoia_format_ch.view()



process make_franz_format {
	input:
		set val(size),
			val(replicate),
			file(sequoia_format) from sequoia_format_ch
		file(input_convert_src)
		file(sex_data_file)
	output:
		file("${sequoia_format.baseName}.franz") into franz_input_ch
	"""
	Rscript ${input_convert_src} ${sequoia_format}
	echo "1 ${size} / ${sequoia_format.baseName}" > ${sequoia_format.baseName}.franz_final
	echo "169 location" >> ${sequoia_format.baseName}.franz_final
	cat ${sequoia_format.baseName}.franz >> ${sequoia_format.baseName}.franz_final
	mv ${sequoia_format.baseName}.franz_final ${sequoia_format.baseName}.franz
	sed -i 's/"//g' ${sequoia_format.baseName}.franz
	"""
}

process run_franz {
	publishDir './output/42_franz/', mode: 'copy'
	input:
		file(franz_input) from franz_input_ch
	output:
		file("${franz_input.baseName}.parentage")
	"""
	FRANz ${franz_input}
	mv parentage.csv ${franz_input.baseName}.parentage
	"""
}
