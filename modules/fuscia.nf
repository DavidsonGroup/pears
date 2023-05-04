process runFuscia {
	input:
	tuple val (fusion_genes), val (chrom1), val (gene1), val (base1), val (sequence1), val (chrom2), val (gene2), val (base2), val (sequence2)
	val "cellranger done"	 	

	output:
	val "fuscia complete"	

	script:
	"""
	#!/usr/bin/env python3
	import os
	os.makedirs('/stornext/Bioinf/data/lab_davidson/wu.s/nf_pears_test/2-5/fuscia_out', exist_ok = True)
	gr1 = f'$chrom1:{min($gene1, $base1)}-{max($gene1, $base1)}'
	gr2 = f'$chrom2:{min($gene2, $base2)}-{max($gene2, $base2)}'
	fusion_name = f'${fusion_genes}_${task.index}'
	command = f'python $projectDir/submodules/fuscia/discover_chimeric_transcripts.py $params.out_dir/cellranger_output/outs/*.bam {gr1} {gr2} $params.out_dir/fuscia_out/ {fusion_name} $params.fuscia_mapqual'
	os.system(command)
	"""
}

