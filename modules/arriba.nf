process runArriba {
	label 'process_high'
	publishDir "${params.out_dir}/arriba_out", mode: 'copy'

	input:
	path(bam_file)
	path(ref_fasta)
	path(ref_gene)

	output:
	path("fusions.tsv")

	script:
	arriba_exclusion_param="-f blacklist" //no filtering
	if(params.arriba_exclusion_file)
	    arriba_exclusion_param="-b ${params.arriba_exclusion_file}"
	"""
	arriba \
		-x ${bam_file} \
		-o fusions.tsv \
		-O fusions.discarded.tsv \
		-a ${ref_fasta} \
		-g ${ref_gene} \
		${arriba_exclusion_param}
	"""
}

process getFusionReadsArriba {
	label 'process_low'
	publishDir "${params.out_dir}/arriba_out", mode: 'copy'

	input:
	tuple val(fusion_genes), val(chrom1), val(gene1), val(base1), val(sequence1), val(chrom2), val(gene2), val(base2), val(sequence2)
	path(fusion_table)
	path(fastq_r1)

	output:
	tuple val("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}"),
          path("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}.fastq")

	script:
	"""
	set +e

	fus=`echo ${fusion_genes} | sed 's/--/\t/g'` ;
	pos=`echo -e "${chrom1}:${base1}\t${chrom2}:${base2}"`
	fusion_name=`echo ${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}`

	grep -e "\$fus" -e "\$pos" ${fusion_table} |\
		cut -f30 |\
		sed 's/,/ \\n/g' |\
		sed 's/^/^@/g' |\
		grep -f - <(gunzip -c ${fastq_r1}) -A3 --no-group-separator |\
			sed "/^[@+]/! s/^/%%%%%%%/g" > "\$fusion_name".fastq ;
	"""
}

process getBarcodesArriba {
	label 'process_low'
	publishDir "${params.out_dir}/arriba_out", mode: 'copy'

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


process get_novel_fusions {
    label 'process_tiny'
    publishDir "${params.out_dir}/arriba_out", mode: 'copy'

    input:
    path("fusions.tsv")

    output:
    path("extra_target.csv")

    script:
    """
    awk -F'\t' -v min_support=${params.min_arriba_support ?: 1} 'BEGIN { OFS="," }
    NR==1 {
        print "fusion genes","chrom1","base1","strand1","chrom2","base2","strand2"
        next
    }
    {
        split(\$5, a, ":")
        split(\$6, b, ":")
        split(\$3, s1, "/")
        split(\$4, s2, "/")

        gene1 = \$1
        gene2 = \$2

        sub(/,.*/, "", gene1)
        sub(/,.*/, "", gene2)

	gsub(/\\(/, "-", gene1); gsub(/\\)/, "", gene1)
	gsub(/\\(/, "-", gene2); gsub(/\\)/, "", gene2)

        support = \$10 + \$11 + \$12
        if (support < min_support) next

        print gene1"--"gene2, a[1], a[2], s1[2], b[1], b[2], s2[2], \$15
    }' fusions.tsv > extra_target.csv
    """
}

