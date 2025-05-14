//FIX -d option set in config file.

process runFlexiplex {
	publishDir "${params.out_dir}/flexiplex_out", mode: 'copy'

	input:
	tuple val (fusion_genes), val (chrom1), val (gene1), val (base1), val (sequence1), val (chrom2), val (gene2), val (base2), val (sequence2)	 
	path barcode_file

	output: 
	path "*"

	script:
	"""
	# Define the fusion name and flexiplex path
    	fusion_name="${fusion_genes}_${task.index}"

    	# Run flexiplex with the specified parameters
    	paste <(gunzip -c ${params.read1}) <(gunzip -c ${params.read2}) | \
	   sed "/^[@+]/! s/^/START/g" | sed "/^[@+]/! s/	//g" | \
	   ${projectDir}/modules/flexiplex/flexiplex -p $task.cpus -n \${fusion_name} \
	   -x ${sequence1}${sequence2} -d grep -f 1 > \${fusion_name}_reads.fastq  

    	${projectDir}/modules/flexiplex/flexiplex -x START \
	   ${params.flexiplex_demultiplex_options} \
	   -k ${barcode_file} -n barcodes_\${fusion_name} \${fusion_name}_reads.fastq 
	   
    """
}

