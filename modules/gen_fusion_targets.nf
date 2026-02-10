process genFusionTargets {
	label 'process_low'
	publishDir params.out_dir, mode: 'copy'

	input:
	path(known_list)
	path(ref_gene)
	path(ref_fasta)
	val(flexi_searchlen)
	val(fuscia_up)
	val(fuscia_down)

	output:
	path "fusion_targets.csv"

	script:
	"""
	gen_fusion_targets.py \
		${known_list} \
		${ref_gene} \
		${ref_fasta} \
		${flexi_searchlen} \
		. \
		${fuscia_up} \
		${fuscia_down}
	"""
}
