process formatFuscia {
	label 'process_tiny'
	publishDir "${params.out_dir}", mode: 'copy'

	input:
	path(fuscia_files)
	val(output_file)

	output:
	path("${output_file}")

	script:
	def input_files = fuscia_files.collect { f -> f.name }.join(' ')
	"""
	format_barcodes.py --type fuscia --output '${output_file}' ${input_files}
	"""
}

process formatFlexiplex {
	label 'process_tiny'
	publishDir "${params.out_dir}", mode: 'copy'

	input:
	path(barcode_files)
	val(output_file)

	output:
	path("${output_file}")

	script:
	def input_files = barcode_files.collect { f -> f.name }.join(' ')
	"""
	format_barcodes.py --type flexiplex --output '${output_file}' ${input_files}
	"""
}

process formatArriba {
	label 'process_tiny'
	publishDir "${params.out_dir}", mode: 'copy'

	input:
	path(barcode_files)
	val(output_file)

	output:
	path("${output_file}")

	script:
	def input_files = barcode_files.collect { f -> f.name }.join(' ')
	"""
	format_barcodes.py --type arriba --output '${output_file}' ${input_files}
	"""
}
