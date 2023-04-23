project_dir = projectDir


process RUN_CELLRANGER {
	
	input:
	params.reference
	params.reads
	params.out_dir

	output:
 	path '$params.out_dir/cellranger_output/'
	
	script:
	"""
	module load cellranger/3.1.0

	cellranger count --id=cellranger_output --transcriptome=$params.reference --fastq=$params.reads 
	
	"""

}
