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
	val(protocol)

	output:
	path("${output_file}")

	script:
	def input_files = barcode_files.collect { f -> f.name }.join(' ')

	if (protocol == "10x-3prime-visiumHD") {
	"""
		format_barcodes.py --type flexiplex_hd --output '${output_file}' ${input_files}
	"""
	} else {
	"""
		format_barcodes.py --type flexiplex --output '${output_file}' ${input_files}
	"""
	}
}

process formatArriba {
	label 'process_tiny'
	publishDir "${params.out_dir}", mode: 'copy'

	input:
	path(barcode_files)
	val(output_file)
	val(protocol)

	output:
	path("${output_file}")

	script:
	def input_files = barcode_files.collect { f -> f.name }.join(' ')
	if (protocol == "10x-3prime-visiumHD") {
	"""
	format_barcodes.py --type arriba_hd --output '${output_file}' ${input_files}
	"""
	} else {
	"""
	format_barcodes.py --type arriba --output '${output_file}' ${input_files}
	"""
	}
}

process combineFusionCalls {
    tag "combine fusion calls"
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path arriba_csv
    path flexiplex_csv
    path fuscia_csv

    output:
    path "combined_fusions.csv"

    script:
    """
    combine_fusions.py \
    	   --arriba $arriba_csv \
  	   --flexiplex $flexiplex_csv \
  	   --fuscia $fuscia_csv \
  	   -o combined_fusions.csv
  """
}
process convertToSpatialBarcodes {
    label 'process_tiny'
    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path combined_csv
    path whitelist
    val  bin_size

    output:
    path "combined_fusions_spatial.csv"

    script:
    """
    convert_barcodes_spatial.py \
        --input $combined_csv \
        --whitelist $whitelist \
        --bin-size $bin_size \
        --output combined_fusions_spatial.csv
    """
}
