
project_dir = projectDir

process GEN_MASTERDATA {
	conda 'pears'
	
	input:
	params.shr_output
	params.reference
	params.flexi_searchlen
	
	output:
  	file 'masterdata.csv' into ch_output
	
	script:
	"""
	
	python $project_dir/scripts/gen_masterdata.py $params.shr_output $params.reference $params.flexi_searchlen 
	
	"""
	
	

}
