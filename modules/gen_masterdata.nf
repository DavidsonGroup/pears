project_dir = projectDir

process GEN_MASTERDATA {
	conda './env/pears_env.yml'
	
	input:
	params.shr_output
	params.reference
	params.flexi_searchlen
	params.out_dir

	output:
  	stdout emit: masterdata
	
	script:
	"""

	python $project_dir/scripts/gen_masterdata.py $params.shr_output $params.reference $params.flexi_searchlen $params.out_dir 
	
	"""

}
