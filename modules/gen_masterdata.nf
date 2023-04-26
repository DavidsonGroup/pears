process GEN_MASTERDATA {
	conda "${baseDir}/env/pears_env.yml"
	
	output:
	path '$params.out_dir/masterdata.csv'

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
