process RUN_CELLRANGER {
	time = '1d'

	input:
	params.reference
	params.reads
	params.out_dir
	
	script:
	"""
	cd $params.out_dir

	module load cellranger/3.1.0

	cellranger count --id=cellranger_output --transcriptome=$params.reference --fastq=$params.reads 
	
	"""

}
