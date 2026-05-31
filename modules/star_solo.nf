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


// run STAR SOLO use different setting for visium
// to deal with very large files (reduce multi-mapping)
// barcode demultiplexing isn't optimal for visium yet
process runSTARSolo {
	label 'process_high'
	publishDir "${params.out_dir}/STARsolo", mode: 'copy', pattern: "Aligned.sortedByCoord.out.bam*"

	input:
	path(read1)
	path(read2)
	path(genome_index)
	path(include_list)
	val(umi_len)
	val(protocol)

	output:
	path("Aligned.sortedByCoord.out.bam"), emit: bam
	path("Aligned.sortedByCoord.out.bam.bai"), emit: bam_index
	path("Aligned.out.bam"), emit: arriba_bam

	script:

	STAR_args_common="STAR \
		--runThreadN ${task.cpus} \
		--genomeDir ${genome_index} \
		--genomeLoad NoSharedMemory \
		--readFilesIn ${read2} ${read1} \
		--readFilesCommand zcat \
		--outSAMtype BAM Unsorted SortedByCoordinate \
		--outSAMunmapped Within \
		--outBAMcompression 0 \
		--peOverlapNbasesMin 10 \
		--alignSplicedMateMapLminOverLmate 0.5 \
		--alignSJstitchMismatchNmax 5 -1 5 5 \
		--chimOutType WithinBAM HardClip \
		--outSAMattributes NH HI nM AS CB UB \
		--soloUMIdedup NoDedup \
		--chimScoreJunctionNonGTAG 0 \
		--chimSegmentReadGapMax 3" 

	if(protocol=="10x-3prime-visiumHD"){
	"""
	   ${STAR_args_common} \
      		--soloType CB_UMI_Complex \
		--soloCBwhitelist $include_list $include_list \
		--soloCBposition 0_10_0_23 1_-17_1_-4 \
		--soloUMIposition 0_0_0_8 \
		--soloCBmatchWLtype 1MM \
		--soloBarcodeReadLength 0 \
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
		--soloType CB_UMI_Simple \
		--soloCBwhitelist $include_list \
		--soloUMIlen $umi_len \
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

