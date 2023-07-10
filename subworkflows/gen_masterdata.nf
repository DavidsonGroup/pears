process GEN_MASTERDATA {
	conda "${projectDir}/env/pears_env.yml"
	publishDir params.out_dir	

	input:
	params.known_list
	params.reference
	params.flexi_searchlen
	params.out_dir
	params.fuscia_up
	params.fuscia_down

	output:
	path "masterdata.csv"	

	script:
	"""

	python $projectDir/submodules/gen_masterdata.py $params.known_list $params.reference $params.flexi_searchlen . $params.fuscia_up $params.fuscia_down
	
	"""

}
