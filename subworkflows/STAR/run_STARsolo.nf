process RUN_STARsolo {
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
	cd $projectDir/modules/STAR/source/

	STAR \
	--runThreadN $params.cpus \
	--genomeDir $projectDir/STAR/STAR_index_GRCh38_GENCODE38/ --genomeLoad NoSharedMemory \
	--readFilesIn $params.read2 $params.read1 \
	--outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outBAMcompression 0\
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
	--chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 \
	--soloType CB_UMI_Simple --soloCBwhitelist $projectDir/modules/barcode_whitelist/$params.protocol_version.txt --soloUMIlen $params.umi_len --outSAMattributes NH HI nM AS CB UB\

	"""

}