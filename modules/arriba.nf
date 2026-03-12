process runArriba {
	label 'process_medium'
	publishDir "${params.out_dir}/arriba_out", mode: 'copy'

	input:
	path(bam_file)
	path(ref_fasta)
	path(ref_gene)

	output:
	path("fusions.tsv")

	script:
	"""
	arriba \
		-x ${bam_file} \
		-o fusions.tsv \
		-O fusions.discarded.tsv \
		-a ${ref_fasta} \
		-g ${ref_gene} \
		-f blacklist
	"""
}




process getBarcodesArriba {
	label 'process_tiny'
	publishDir "${params.out_dir}/arriba_out", mode: 'copy'

	input:
	tuple val(fusion_genes), val(chrom1), val(gene1), val(base1), val(sequence1), val(chrom2), val(gene2), val(base2), val(sequence2)
	path(fusion_table)
	path(include_list)
	path(fastq_r1)
	val(flexiplex_demultiplex_options)

	output:
	path("barcodes_${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}_reads_barcodes.txt")

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
			sed "/^[@+]/! s/^/START/g" > "\$fusion_name".fastq ;

	flexiplex -x START \
		${flexiplex_demultiplex_options} \
		-k ${include_list} -n barcodes_"\$fusion_name" \
		"\$fusion_name".fastq
	"""
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
        print "fusion genes","chrom1","base1","strand1","chrom2","base2","strand2","classification"
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

