process GEN_MASTERDATA {
	conda "${projectDir}/env/pears_env.yml"
	
	input:
	params.shr_output
	params.reference
	params.flexi_searchlen
	params.out_dir
	params.fuscia_up
	params.fuscia_down
	
	
	
	script:
	"""

	python $projectDir/submodules/gen_masterdata.py $params.shr_output $params.reference $params.flexi_searchlen $params.out_dir $params.fuscia_up $params.fuscia_down
	
	"""

}
