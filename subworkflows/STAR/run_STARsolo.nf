process RUN_STARsolo {
	publishDir "${params.out_dir}/STARsolo"
	time = '1d'
	memory = '200'

	input:
	params.reference
	params.reads
	params.out_dir
	params.STAR
	
	output:
	val "STAR done"	

	script:
	"""
	$projectDir/modules/STAR/STAR \
	--runThreadN $params.cpus \
	--genomeDir $params.genome_index --genomeLoad NoSharedMemory \
	--readFilesIn $params.read2 $params.read1 \
	--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outBAMcompression 0\
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
	--chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 \
	--soloType CB_UMI_Simple --soloCBwhitelist $projectDir/modules/barcode_whitelist/$params.protocol_version.txt --soloUMIlen $params.umi_len --outSAMattributes NH HI nM AS CB UB\
	
	"""

}

process format_bam {
	publishDir "${params.out_dir}/STARsolo"
	
	input:
	RUN_STARsolo.out

	script:
	"""
	$projectDir/modules/samtools-1.18/samtools view -H Aligned.sortedByCoord.out.bam | sed -E -e 's/SN:([0-9XY])/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - Aligned.sortedByCoord.out.bam > Aligned.sortedByCoord.out_chr.bam

	$projectDir/modules/samtools-1.18/samtools index file Aligned.sortedByCoord.out_chr.bam
	"""

}
