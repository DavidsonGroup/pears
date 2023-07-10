process BUILD_INDEX {
	input:

	output:

	script:
	"""
	STAR 
	--runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir $projectDir/STAR/STAR_index_GRCh38_GENCODE38 \
	--genomeFastaFiles $projectDir/STAR/GRCh38.fa\
	--sjdbGTFfile $projectDir/STAR/GENCODE38.gtf\
	--sjdbOverhang $params.read_length - 1
	"""


}

process RUN_STAR {
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
	cd $params.out_dir

	module load $params.STAR

	STAR \
	--runThreadN 8 \
	--genomeDir /stornext/Bioinf/data/lab_davidson/wu.s/results/5cl/arriba/STAR_index_GRCh38_GENCODE38/ --genomeLoad NoSharedMemory \
	--readFilesIn $params.read2 $params.read1 \
	--outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outBAMcompression 0\
	--outBAMsortingBinsN 200 \
	--limitBAMsortRAM 2500000000 \
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
	--chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 \
	--soloType CB_UMI_Simple --soloCBwhitelist $projectDir/STAR --soloUMIlen 12 --outSAMattributes NH HI nM AS CB UB\

	"""

}
