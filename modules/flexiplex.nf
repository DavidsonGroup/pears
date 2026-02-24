//FIX -d option set in config file.


//path "${fusion_genes}_${task.index}_reads.fastq"
//path barcode_file

process runFlexiplex {
	label 'process_low'
	publishDir "${params.out_dir}/flexiplex_out", mode: 'copy'

	input:
	tuple val(fusion_genes), val(chrom1), val(gene1), val(base1), val(sequence1), val(chrom2), val(gene2), val(base2), val(sequence2)
	path(include_list)
	path(fastq_r1)
	path(fastq_r2)
	val(flexiplex_demultiplex_options)

	output:
	path "barcodes_${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}_reads_barcodes.txt"

	script:
	fusion_name="${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}"
	"""
	# Run flexiplex with the specified parameters
	paste <(gunzip -c ${fastq_r1}) <(gunzip -c ${fastq_r2}) | \
	sed "/^[@+]/! s/^/START/g" | sed "/^[@+]/! s/	//g" | \
	flexiplex -p ${task.cpus} -n ${fusion_name} \
		-x ${sequence1}${sequence2} -d grep -f 1 > ${fusion_name}_reads.fastq

	flexiplex -x START \
		${flexiplex_demultiplex_options} \
		-k ${include_list} -n barcodes_${fusion_name} ${fusion_name}_reads.fastq
    """
}

