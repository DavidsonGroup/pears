project_dir = projectDir

process GEN_MASTERDATA {
	conda './env/pears_env.yml'
	
	input:
	masterdata.out
	params.reference
	params.flexi_searchlen
	params.out_dir

	output:
  	stdout emit: masterdata
	
	script:
	"""
	python $project_dir/scripts/run_fuscia.sh $masterdata.out $params.fuscia_mapqual $params.out_dir
	
	"""

}
