process buildSTARIndex {
	label 'process_high'
	storeDir "${params.out_dir}/references/"

	input:
	path(ref_fasta)
	path(ref_gtf)
	val(read_length)

	output:
	path("${star_index_dir}"), emit: index

	script:
	sjdb_overhang = read_length.toInteger() - 1
	star_index_dir = "star_index"
	"""
	STAR \
		--runMode genomeGenerate \
		--runThreadN ${task.cpus} \
		--genomeDir ${star_index_dir} \
		--genomeFastaFiles ${ref_fasta} \
		--sjdbGTFfile ${ref_gtf} \
		--sjdbOverhang ${sjdb_overhang}
	"""
}


// Align the cDNA read (R2) with STAR. Barcodes and UMIs come from the flexiplex
// demultiplexing step that runs before alignment, so R1 is not passed to STAR
// at all and the CB/UB tags are written onto the BAM afterwards by
// transferBarcodesToBAM.
// VisiumHD uses different settings to deal with very large files
// (reduced multi-mapping).
process runSTAR {
	label 'process_high'
	publishDir "${params.out_dir}/STAR", mode: 'copy', pattern: "Aligned.sortedByCoord.out.bam*"

	input:
	path(read2)
	path(genome_index)
	val(protocol)

	output:
	path("Aligned.sortedByCoord.out.bam"), emit: bam
	path("Aligned.sortedByCoord.out.bam.bai"), emit: bam_index
	path("Aligned.out.bam"), emit: arriba_bam

	script:
	// STAR takes a comma separated list for multiple files of the same mate
	def read2_files = (read2 instanceof List) ? read2 : [read2]
	def read2_arg = read2_files.join(',')
	def read_files_command = read2_files.every { it.name.endsWith('.gz') } ? "--readFilesCommand zcat " : ""

	STAR_args_common="STAR \
		--runThreadN ${task.cpus} \
		--genomeDir ${genome_index} \
		--genomeLoad NoSharedMemory \
		--readFilesIn ${read2_arg} \
		${read_files_command}--outSAMtype BAM Unsorted SortedByCoordinate \
		--outSAMunmapped Within \
		--outBAMcompression 0 \
		--alignSplicedMateMapLminOverLmate 0.5 \
		--alignSJstitchMismatchNmax 5 -1 5 5 \
		--chimOutType WithinBAM HardClip \
		--outSAMattributes NH HI nM AS \
		--chimScoreJunctionNonGTAG 0 \
		--chimSegmentReadGapMax 3" 

	if(protocol=="10x-3prime-visiumHD"){
	"""
	   ${STAR_args_common} \
		--limitBAMsortRAM 149759137861 \
		--outSAMmultNmax 1 \
		--chimMainSegmentMultNmax 1 \
		--outFilterMultimapNmax 10 \
		--chimMultimapNmax 1 \
		--chimSegmentMin 20 \
		--chimJunctionOverhangMin 20 \
		--chimScoreDropMax 20 \
		--chimScoreSeparation 10

	samtools index Aligned.sortedByCoord.out.bam
	"""
	} else {
	"""
	   ${STAR_args_common} \
		--outFilterMultimapNmax 50 \
		--chimMultimapNmax 50 \
		--chimJunctionOverhangMin 10 \
		--chimSegmentMin 10 \
		--chimScoreDropMax 30 \
		--chimScoreSeparation 1


	samtools index Aligned.sortedByCoord.out.bam
	"""
	}
}

