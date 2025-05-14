process RUN_STARsolo {
	publishDir "${params.out_dir}/STARsolo", mode: 'copy'
	time = '1d'
	memory = '200'

	input:
	params.genome_index
	params.read1
	params.read2
	params.out_dir
	params.STAR

	output:
	path "Aligned.sortedByCoord.out.bam" // Define outputs
    	path "Aligned.sortedByCoord.out.bam.bai"
    	path "Solo.out/Gene/filtered/barcodes.tsv"
    	path "Solo.out/Gene/filtered/features.tsv"
    	path "Solo.out/Gene/filtered/matrix.mtx"

	script:
	"""
	$params.STAR \
	--runThreadN $task.cpus \
	--genomeDir $params.genome_index \
	--genomeLoad NoSharedMemory \
	--readFilesIn $params.read2 $params.read1 \
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
	--soloType CB_UMI_Simple \
	--soloCBwhitelist $params.barcode_whitelist \
	--soloUMIlen $params.umi_len \
	--soloUMIdedup NoDedup \
	--outSAMattributes NH HI nM AS CB UB \
	--soloBarcodeReadLength 0

	$projectDir/modules/samtools-1.18/samtools index Aligned.sortedByCoord.out.bam

	"""

}

process format_bam {
	publishDir "${params.out_dir}/STARsolo"
	
	input:
	RUN_STARsolo.out

	output:
	file('*.bam') into bam_files, val("aligner done") into aligner_done

	script:
	"""
	$projectDir/modules/samtools-1.18/samtools view -H Aligned.sortedByCoord.out.bam | sed -E -e 's/SN:([0-9XY])/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - Aligned.sortedByCoord.out.bam > Aligned.sortedByCoord.out_chr.bam

	$projectDir/modules/samtools-1.18/samtools index file Aligned.sortedByCoord.out_chr.bam
	"""

}
