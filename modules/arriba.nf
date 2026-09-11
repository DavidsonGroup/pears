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
          path("${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}_arriba_read_ids.txt")

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
		grep -v '^\$' > "\$fusion_name"_arriba_read_ids.txt

	touch "\$fusion_name"_arriba_read_ids.txt
	"""
}

// Look the arriba fusion reads up in the demultiplexing table (see
// modules/demultiplex.nf) and write the fusion calls. All fusions are done in
// one job, in one pass over the table, and the calls are written straight out
// rather than as a barcode table per fusion that then needs formatting.
process getFusionCallsArriba {
	label 'process_low'
	publishDir "${params.out_dir}", mode: 'copy'

	input:
	path(read_ids)
	path(barcode_table)

	output:
	path("arriba_fusion_calls.csv")

	script:
	"""
	lookup_barcodes.py \\
		--barcodes ${barcode_table} \\
		--output arriba_fusion_calls.csv \\
		${read_ids}
	"""
}


// Turn the arriba calls into extra fusion targets for the other tools.
//
// A gene pair on the known list is kept whatever its support, so the
// breakpoints arriba found for it become targets too. With
// --arriba_strict_breakpoints those are dropped instead: the known list
// supplies the coordinates for its own gene pairs, and arriba is only reported
// where its breakpoint matches one of them. Gene pairs that are not on the
// known list always need min_arriba_support to be discovered.
//
// Keep prose out of the awk program below - it is single quoted, so one
// apostrophe in a comment ends the program and awk gets a truncated script.
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

    awk -F'\t' -v min_support=${params.min_arriba_support ?: 1} \\
        -v strict=${params.arriba_strict_breakpoints ? 1 : 0} '
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

    	if (fusion in allow) {
    	   if (strict == 1) next
    	} else if (support < min_support) next

    	print fusion, a[1], a[2], s1[2], b[1], b[2], s2[2], \$15
	}' allowlist.txt fusions.tsv > extra_target.csv
	"""
}

