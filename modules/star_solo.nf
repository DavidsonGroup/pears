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

process runSTARSolo {
	label 'process_high'
	publishDir "${params.out_dir}/STARsolo", mode: 'copy'

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

	script:

	STAR_args_common="STAR \
		--runThreadN ${task.cpus} \
		--genomeDir ${genome_index} \
		--genomeLoad NoSharedMemory \
		--readFilesIn ${read2} ${read1} \
		--readFilesCommand zcat \
		--outSAMtype BAM SortedByCoordinate \
		--outSAMunmapped Within \
		--outBAMcompression 0 \
		--outFilterMultimapNmax 50 \
		--peOverlapNbasesMin 10 \
		--alignSplicedMateMapLminOverLmate 0.5 \
		--alignSJstitchMismatchNmax 5 -1 5 5 \
		--chimSegmentMin 10 \
		--chimOutType WithinBAM HardClip \
		--chimJunctionOverhangMin 10 \
		--chimScoreDropMax 30 \
		--chimScoreJunctionNonGTAG 0 \
		--chimScoreSeparation 1 \
		--chimSegmentReadGapMax 3 \
		--chimMultimapNmax 50 \
		--outSAMattributes NH HI nM AS CB UB \
		--soloUMIdedup NoDedup"

	if(protocol=="10x-3prime-visiumHD"){
	"""
	   ${STAR_args_common} \
      		--soloType CB_UMI_Complex \
		--soloCBwhitelist $include_list $include_list \
		--soloCBposition 0_10_0_23 1_-17_1_-4 \
		--soloUMIposition 0_0_0_8 \
		--soloCBmatchWLtype 1MM \
		--soloBarcodeReadLength 0 \
		--limitBAMsortRAM 149759137861 ;

	samtools index Aligned.sortedByCoord.out.bam
	"""
	} else {
	"""
	   ${STAR_args_common} \
		--soloType CB_UMI_Simple \
		--soloCBwhitelist $include_list \
		--soloUMIlen $umi_len

	samtools index Aligned.sortedByCoord.out.bam
	"""
	}
}

/** process formatBAM {
	label 'process_tiny'
	publishDir "${params.out_dir}/STARsolo"

	output:
	file('*.bam')

	script:
	"""
	samtools view -H Aligned.sortedByCoord.out.bam | sed -E -e 's/SN:([0-9XY])/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - Aligned.sortedByCoord.out.bam > Aligned.sortedByCoord.out_chr.bam

	samtools index Aligned.sortedByCoord.out_chr.bam
	"""
} **/
