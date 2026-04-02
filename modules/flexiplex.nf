//FIX -d option set in config file.

process getFusionReadsFlexiplex {
	label 'process_low'
	publishDir "${params.out_dir}/flexiplex_out", mode: 'copy'

	input:
	tuple val(fusion_genes), val(chrom1), val(gene1), val(base1), val(sequence1), val(chrom2), val(gene2), val(base2), val(sequence2)
	path(fastq_r1)
	path(fastq_r2)

	output:
	tuple val("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}"),
          path("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}_reads.fastq")

	script:
	def fusion_name="${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}"

	"""
	# Run flexiplex with the specified parameters
	paste <(gunzip -c ${fastq_r1}) <(gunzip -c ${fastq_r2}) | \
	sed "/^[@+]/! s/^/%%%%%%%/g" | sed "/^[@+]/! s/	//g" | \
	flexiplex -p ${task.cpus} -n ${fusion_name} \
		-x ${sequence1}${sequence2} -d grep -f 1 > ${fusion_name}_reads.fastq
	"""
}

process getBarcodesFlexiplex {
	label 'process_low'
	publishDir "${params.out_dir}/flexiplex_out", mode: 'copy'

	input:
	tuple val(fusion_name), path(reads)
	path(include_list)
	val(flexiplex_demultiplex_options)
	val(protocol)

	output:
	path "barcodes_${fusion_name}_reads_barcodes.txt"

	script:
	FLEX="/vast/projects/lab_davidson/davidson.n/flexiplex/flexiplex"
	
	if (protocol == "10x-3prime-visiumHD") {
	   """
	   cat ${reads} | ${FLEX} \
               -x "%%%%%%%?????????G?????????????" \
               -b "???????????????" \
               -k ${include_list} \
               -e 1 -f 2 -r false | \
    	       ${FLEX} \
               -x "%%%%%%%" \
               -u "?????????" \
               -x "G" \
               -b "?????????????" \
               -k ${include_list} \
               -e 1 -f 2 -n barcodes_${fusion_name} ;
	   """
	} else {

	  """
	  flexiplex -x %%%%%%% \
		  ${flexiplex_demultiplex_options} \
		  -k ${include_list} -n barcodes_${fusion_name} ${reads} ;
          """
	  }
}


