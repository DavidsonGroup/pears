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

// Collect the IDs of the reads arriba assigned to each fusion. The barcodes
// for those reads are looked up in the pipeline-wide demultiplexing table, so
// the reads themselves are no longer pulled out of the FASTQ.
process getFusionReadsArriba {
	label 'process_low'
	publishDir "${params.out_dir}/arriba_out", mode: 'copy'

	input:
	tuple val(fusion_genes), val(chrom1), val(gene1), val(base1), val(sequence1), val(chrom2), val(gene2), val(base2), val(sequence2)
	path(fusion_table)

	output:
	tuple val("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}"),
          path("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}_read_ids.txt")

	script:
	"""
	set +e

	pos=`echo -e "${chrom1}:${base1}\t${chrom2}:${base2}"`
	fusion_name=`echo ${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}`

	grep -e "\$pos" ${fusion_table} |\
		cut -f30 |\
		sed 's/,/\\n/g' |\
		sed 's/[[:space:]]*\$//' |\
		sort -u |\
		grep -v '^\$' > "\$fusion_name"_read_ids.txt

	touch "\$fusion_name"_read_ids.txt
	"""
}

// Pull the barcodes for the arriba fusion reads out of the pipeline-wide
// demultiplexing table (see modules/demultiplex.nf). All fusions are done in
// one job: the table is large, and scanning it once beats queueing a short job
// per fusion.
process getBarcodesArriba {
	label 'process_low'
	publishDir "${params.out_dir}/arriba_out", mode: 'copy'

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


process get_novel_fusions {
    label 'process_tiny'
    publishDir "${params.out_dir}/arriba_out", mode: 'copy'

    input:
    path("fusions.tsv")
    val known_list

    output:
    path("extra_target.csv")

    script:
    """
    echo ${known_list}
    ls -l "${known_list}"
    head ${known_list}

    tail -n +2 ${known_list} | cut -d',' -f1 > allowlist.txt

    awk -F'\t' -v min_support=${params.min_arriba_support ?: 1} '
    BEGIN { OFS="," }

    FNR==NR {
    	    allow[\$1] = 1
    	    next
    }

    FNR==1 {
    	   print "fusion genes","chrom1","base1","strand1","chrom2","base2","strand2","confidence"
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

    	fusion = gene1 "--" gene2
    	support = \$10 + \$11 + \$12

    	if (support < min_support && !(fusion in allow)) next

    	print fusion, a[1], a[2], s1[2], b[1], b[2], s2[2], \$15
	}' allowlist.txt fusions.tsv > extra_target.csv
	"""
}

