//FIX -d option set in config file.

// Search the cDNA read (R2) for reads spanning the fusion junction. Only the
// read IDs are needed downstream - the barcodes are looked up in the
// pipeline-wide demultiplexing table rather than being called again here.
process getFusionReadsFlexiplex {
	label 'process_low'
	publishDir "${params.out_dir}/flexiplex_out", mode: 'copy'

	input:
	tuple val(fusion_genes), val(chrom1), val(gene1), val(base1), val(sequence1), val(chrom2), val(gene2), val(base2), val(sequence2)
	path(fastq_r2)

	output:
	tuple val("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}"),
          path("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}_read_ids.txt"), emit: read_ids
	path("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}_reads.fastq"), emit: reads

	script:
	def fusion_name="${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}"

	"""
	# Grep-like search of the cDNA read for the fusion junction sequence
	for f in ${fastq_r2} ; do
		if [[ "\$f" == *.gz ]] ; then gunzip -c "\$f" ; else cat "\$f" ; fi
	done | \\
	flexiplex -p ${task.cpus} -n ${fusion_name} \\
		-x ${sequence1}${sequence2} -d grep -f 1 > ${fusion_name}_reads.fastq

	awk 'NR%4==1 { sub(/^@/, "", \$0) ; print \$1 }' ${fusion_name}_reads.fastq | \\
		sort -u > ${fusion_name}_read_ids.txt
	"""
}

// Pull the barcodes for the fusion-supporting reads out of the pipeline-wide
// demultiplexing table (see modules/demultiplex.nf). All fusions are done in
// one job: the table is large, and scanning it once beats queueing a short job
// per fusion.
process getBarcodesFlexiplex {
	label 'process_low'
	publishDir "${params.out_dir}/flexiplex_out", mode: 'copy'

	input:
	path(read_ids)
	path(barcode_table)

	output:
	path "barcodes_*_reads_barcodes.txt"

	script:
	"""
	lookup_barcodes.py --barcodes ${barcode_table} ${read_ids}
	"""
}


//	FLEX="/vast/projects/lab_davidson/davidson.n/flexiplex/flexiplex"
