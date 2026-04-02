process mergeFusionTargetLists {
    tag "merge fusion target lists"

    input:
    val known_list
    val extra_list

    output:
    path "merged_known_fusions.csv"

    script:
    def cmd
    if( known_list && extra_list ) {
        cmd = """
        awk 'NR==1 || FNR>1' ${known_list} ${extra_list} > merged_known_fusions.csv
        """
    }
    else if( known_list ) {
        cmd = """
        cp ${known_list} merged_known_fusions.csv
        """
    }
    else if( extra_list ) {
        cmd = """
        cp ${extra_list} merged_known_fusions.csv
        """
    }
    else {
        cmd = """
        echo "fusion genes,chrom1,base1,strand1,chrom2,base2,strand2" > merged_known_fusions.csv
        """
    }
    cmd
}

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
