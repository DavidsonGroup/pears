nextflow.enable.dsl=2
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'
include { downloadReferences } from './modules/download_references.nf'
include { genFusionTargets ; mergeFusionTargetLists } from './modules/gen_fusion_targets.nf'
include { prepareIncludeList } from './modules/prepare_include_list.nf'
include { calculateReadLength } from './modules/calculate_read_length.nf'
include { buildSTARIndex; runSTARSolo } from './modules/star_solo.nf'
include { runFuscia } from './modules/fuscia.nf'
include { getFusionReadsFlexiplex; getBarcodesFlexiplex } from './modules/flexiplex.nf'
include { formatFuscia; formatFlexiplex; formatArriba ; combineFusionCalls } from './modules/formatting.nf'
include { runArriba ; getFusionReadsArriba; getBarcodesArriba ; get_novel_fusions } from './modules/arriba.nf'

// Calculate barcode length from first line of barcode file (handles gzipped files)
def getBarcodeLength(barcode_path) {
	def f = file(barcode_path)
	def isGz = barcode_path.toString().endsWith('.gz')
	def first_line = null
	if (isGz) {
		f.withInputStream { stream ->
			new java.util.zip.GZIPInputStream(stream).withReader { reader ->
				first_line = reader.readLine()
			}
		}
	} else {
		f.withReader { reader ->
			first_line = reader.readLine()
		}
	}

	if (!first_line) {
		error "No valid barcode line found in file: ${barcode_path}"
	}
	return first_line.length()
}

workflow {
	// Validate parameters against schema
	validateParameters()
	log.info paramsSummaryLog(workflow)

	// 10x Chromium protocol definitions
	// Maps protocol name to [barcode_file, umi_length]
	def protocol_config = [
		// 3' Gene Expression chemistries
		'10x-3prime-v2': ['737K-august-2016.txt.gz', 10],
		'10x-3prime-v3': ['3M-february-2018.txt.gz', 12],
		'10x-3prime-v4': ['3M-3pgex-may-2023.txt.gz', 12],
		// 5' Gene Expression chemistries
		'10x-5prime-v2': ['737K-august-2016.txt.gz', 10],
		'10x-5prime-v3': ['3M-5pgex-jan-2023.txt.gz', 12],
		//spatial
		'10x-3prime-visiumHD': ['visium-hd.txt.gz', 9]
	]

	def umi_length = null
	def barcode_file = null
	// Resolve protocol to barcode file and UMI length
	if (params.protocol) {
		def config = protocol_config[params.protocol]
		barcode_file = params.barcode_include_list ?: "${projectDir}/assets/${config[0]}"
		umi_length = params.umi_len ?: config[1]
		log.info "Protocol ${params.protocol}: barcode_file=${barcode_file}, umi_len=${umi_length}"
	} else if (params.barcode_include_list && params.umi_len) {
		// Manual configuration
		barcode_file = params.barcode_include_list
		umi_length = params.umi_len
	} else {
		error "Either --protocol or both --barcode_include_list and --umi_len must be specified"
	}

	// Build default flexiplex demultiplex options if not provided
	def barcode_length = getBarcodeLength(barcode_file)
	def barcode_pattern = "?" * barcode_length.toInteger()
	def umi_pattern = "?" * umi_length.toInteger()
	def default_flexiplex_opts = "-b \"${barcode_pattern}\" -u \"${umi_pattern}\" -e 1 -f 0"
	flexiplex_demultiplex_options = params.flexiplex_demultiplex_options ?: default_flexiplex_opts
	log.info "Flexiplex demultiplex options: ${flexiplex_demultiplex_options} (barcode_len=${barcode_length}, umi_len=${umi_length})"

	// Use pre-built references if all three are provided, otherwise download
	if (params.ref_fasta && params.ref_gtf && params.star_genome_index) {
		log.info "Using pre-built references: skipping download and index building"
		ref_fasta = channel.value(file(params.ref_fasta))
		ref_gtf = channel.value(file(params.ref_gtf))
		star_index = channel.value(file(params.star_genome_index))
	} else {
		// Download reference genome and annotation
		references = downloadReferences(params.genome_version)
		ref_fasta = references.fasta
		ref_gtf = references.gtf
		star_index = null  // Will be built below
	}

	// Auto-enable discovery when no known fusion list is supplied
	discover_fusions = params.discover_fusions || !params.known_fusions_list

	log.info "discover_fusions = ${discover_fusions}"
	log.info "known_fusions_list = ${params.known_fusions_list}"

	// Prepare barcode include list (decompress if gzipped)
	include_list  = prepareIncludeList(file(barcode_file))


	// Build STAR index if not already provided via pre-built references
	if (!params.star_genome_index) {
		// Calculate R2 read length for STAR index generation
		r2_files = channel.fromPath(params.fastq_r2).collect()
		read_length = calculateReadLength(r2_files)

		star_index = buildSTARIndex(
			ref_fasta,
			ref_gtf,
			read_length
		)
	}

	star_solo_result = runSTARSolo(
		channel.fromPath(params.fastq_r1).collect(),
		channel.fromPath(params.fastq_r2).collect(),
		star_index,
		include_list,
		umi_length,
		params.protocol
	)

	arriba_output = runArriba(star_solo_result.bam, ref_fasta, ref_gtf)

	final_target_list = params.known_fusions_list ? file(params.known_fusions_list) : null

	if( discover_fusions ) {
    	    extra_targets = get_novel_fusions(arriba_output)
    	    final_target_list = final_target_list ? mergeFusionTargetLists(final_target_list, extra_targets) : extra_targets
	}
	
	fusion_targets = genFusionTargets(
		final_target_list,
		ref_gtf,
		ref_fasta,
		params.flexiplex_searchlen,
		params.fuscia_up,
		params.fuscia_down
	)

	fusion_target_rows = fusion_targets \
		| splitCsv(header:true) \
		| map { row ->
			tuple(
				row.fusion_genes,
				row.chrom1,
				row.gene1,
				row.base1,
				row.sequence1,
				row.chrom2,
				row.gene2,
				row.base2,
				row.sequence2
			)
		}

	fuscia_result = runFuscia(fusion_target_rows, star_solo_result.bam, star_solo_result.bam_index, params.fuscia_mapqual)
	flexiplex_reads = getFusionReadsFlexiplex(
		fusion_target_rows,
		channel.fromPath(params.fastq_r1).collect(),
		channel.fromPath(params.fastq_r2).collect(),
	)
	flexiplex_result = getBarcodesFlexiplex(
                                 flexiplex_reads,
                                 include_list,
                                 flexiplex_demultiplex_options,
				 params.protocol
        )

	arriba_reads = getFusionReadsArriba(
                fusion_target_rows,
		arriba_output,
                channel.fromPath(params.fastq_r1).collect(),
	)
	arriba_result = getBarcodesArriba(
                                 arriba_reads,
                                 include_list,
                                 flexiplex_demultiplex_options,
                                 params.protocol
	)

	// collapse each into a single emission
	fuscia_collected = fuscia_result | collect
	flexiplex_collected = flexiplex_result | collect
	arriba_collected = arriba_result | collect

	// formatting
	fuscia_final = formatFuscia(fuscia_collected, "fuscia_fusion_calls.csv")
	flexiplex_final = formatFlexiplex(flexiplex_collected, "flexiplex_fusion_calls.csv")
	arriba_final = formatArriba(arriba_collected, "arriba_fusion_calls.csv")

	combineFusionCalls(arriba_final,flexiplex_final,fuscia_final)

}
