process runArriba {
       	publishDir "${params.out_dir}/arriba_out", mode: 'copy'
	time = '1d'
	memory = '150'

	input:
	path bam_file

	output:
	path "fusions.tsv"

	script:
	"""

	${projectDir}/modules/arriba/arriba \
	    -x $bam_file \
	    -o fusions.tsv \
	    -O fusions.discarded.tsv \
	    -a $params.ref_fasta \
	    -g $params.ref_gene \
	    -f blacklist
	"""
}



process getBarcodes_Arriba {
	echo = true
	publishDir "${params.out_dir}/arriba_out", mode: 'copy'

	input:
        tuple val (fusion_genes), val (chrom1), val (gene1), val (base1), val (sequence1), val (chrom2), val (gene2), val (base2), val (sequence2)
	path fusion_table
	path barcode_file

       	output:
	path "*"

	shell:
	"""
	set +e

	fus=`echo !{fusion_genes} | sed 's/--/\t/g'` ;

	grep "\$fus" !{fusion_table} |\
	   cut -f30 |\
	   sed 's/,/ \\n/g' |\
	   sed 's/^/^@/g' |\
	   grep -f - <(gunzip -c !{params.read1}) -A3 --no-group-separator |\
           sed "/^[@+]/! s/^/START/g" > !{fusion_genes}_${task.index}.fastq ;

	   ${projectDir}/modules/flexiplex/flexiplex -x START \
	   ${params.flexiplex_demultiplex_options} \
	   -k ${barcode_file} -n barcodes_!{fusion_genes}_${task.index} \
	   !{fusion_genes}_${task.index}.fastq

	"""

}

