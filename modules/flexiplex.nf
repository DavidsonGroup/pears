//FIX -d option set in config file.

// Search the cDNA read (R2) for reads spanning the fusion junction. Only the
// read IDs are needed downstream - the barcodes are looked up in the
// pipeline-wide demultiplexing table rather than being called again here.
process getFusionReadsFlexiplex {
	label 'process_low'
	publishDir "${params.out_dir}/flexiplex_out", mode: 'copy', pattern: "*_read_ids.txt"

	input:
	tuple val(fusion_genes), val(chrom1), val(gene1), val(base1), val(sequence1), val(chrom2), val(gene2), val(base2), val(sequence2)
	path(fastq_r2)

	output:
	tuple val("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}"),
          path("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}_flexiplex_read_ids.txt"), emit: read_ids
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

	# flexiplex rewrites the read ID: it appends the strand it matched on
	# (_+ or _-), and prepends "<barcode>_<umi>#" when it has called a
	# barcode. Recover the original ID so it can be looked up in the
	# demultiplexing table.
	awk 'NR%4==1' ${fusion_name}_reads.fastq | \\
		sed -E 's/^@//; s/[[:space:]].*\$//; s/^.*#//; s/_[+-]([0-9]+of[0-9]+)?(_C)?\$//' | \\
		sort -u > ${fusion_name}_flexiplex_read_ids.txt
	"""
}

// Look the fusion-supporting reads up in the demultiplexing table (see
// modules/demultiplex.nf) and write the fusion calls. All fusions are done in
// one job, in one pass over the table, and the calls are written straight out
// rather than as a barcode table per fusion that then needs formatting.
process getFusionCallsFlexiplex {
	label 'process_low'
	publishDir "${params.out_dir}", mode: 'copy'

	input:
	path(read_ids)
	path(barcode_table)

	output:
	path("flexiplex_fusion_calls.csv")

	script:
	"""
	lookup_barcodes.py \\
		--barcodes ${barcode_table} \\
		--output flexiplex_fusion_calls.csv \\
		${read_ids}
	"""
}


//	FLEX="/vast/projects/lab_davidson/davidson.n/flexiplex/flexiplex"
