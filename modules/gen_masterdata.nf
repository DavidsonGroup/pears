project_dir = projectDir

process GEN_MASTERDATA {
	conda './env/pears_env.yml'
	
	input:
	params.shr_output
	params.reference
	params.flexi_searchlen
	params.out_dir
	params.fuscia_up
	params.fuscia_down	

	output:
  	stdout emit: masterdata
	
	script:
	"""

	python $project_dir/submodules/gen_masterdata.py $params.shr_output $params.reference $params.flexi_searchlen $params.out_dir $params.fuscia_up $params.fuscia_down
	
	"""

}
