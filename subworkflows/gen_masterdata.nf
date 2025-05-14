//		conda "${projectDir}/env/pears_env.yml"


process GEN_MASTERDATA {
	publishDir params.out_dir, mode: 'copy'

	input:
	params.known_list
	params.reference
	params.flexiplex_searchlen
	params.out_dir
	params.fuscia_up
	params.fuscia_down

	output:
	path "masterdata.csv"	

	script:
	"""

	python $projectDir/subworkflows/gen_masterdata.py \
	       				$params.known_list \
	       				$params.ref_gene \
					$params.ref_fasta \
					$params.flexiplex_searchlen . $params.fuscia_up $params.fuscia_down
	
	"""

}
