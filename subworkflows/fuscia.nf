process runFuscia {
	publishDir "${params.out_dir}/fuscia_out", mode: 'copy'

	input:
	tuple val (fusion_genes), val (chrom1), val (gene1), val (base1), val (sequence1), val (chrom2), val (gene2), val (base2), val (sequence2)
    	path bam_file 
	path bam_index

	output:
	path('*')

	script:
	// Use Groovy's conditional logic to handle the min and max comparison
    	def gr1 = "${chrom1}:${Math.min(gene1.toInteger(), base1.toInteger())}-${Math.max(gene1.toInteger(), base1.toInteger())}"
    	def gr2 = "${chrom2}:${Math.min(gene2.toInteger(), base2.toInteger())}-${Math.max(gene2.toInteger(), base2.toInteger())}"
    	def fusion_name = "${fusion_genes}_${task.index}"

    	"""
    	# Define regions and fusion name in shell format
    	gr1=$gr1
    	gr2=$gr2
    	fusion_name=$fusion_name
	
    	# Run the fuscia discovery command
    	python $projectDir/modules/fuscia/discover_chimeric_transcripts.py $bam_file $gr1 $gr2 . $fusion_name $params.fuscia_mapqual
    	"""
}
/**

	"""
	#!/usr/bin/env python3
	import os
	os.makedirs('$params.out_dir/fuscia_out', exist_ok = True)
	gr1 = f'$chrom1:{min($gene1, $base1)}-{max($gene1, $base1)}'
	gr2 = f'$chrom2:{min($gene2, $base2)}-{max($gene2, $base2)}'
	fusion_name = f'${fusion_genes}_${task.index}'
	command = f'python $projectDir/modules/fuscia/discover_chimeric_transcripts.py {bam_file[0]} {gr1} {gr2} . {fusion_name} $params.fuscia_mapqual'
	os.system(command)
	"""


**/