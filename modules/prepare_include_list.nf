process prepareIncludeList {
	label 'process_tiny'

	input:
	path(include_list_file)

	output:
	path("include_list.txt")

	script:
	"""
	if [[ "${include_list_file}" == *.gz ]]; then
		gunzip -c ${include_list_file} > include_list.txt
	else
		cp ${include_list_file} include_list.txt
	fi
	"""
}
