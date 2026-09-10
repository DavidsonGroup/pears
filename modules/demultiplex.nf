/*
 * -------------------------------------------------
 *  Barcode demultiplexing (run before alignment)
 * -------------------------------------------------
 * Every read is demultiplexed once, with flexiplex, straight off R1. The
 * resulting table (read ID, barcode, edit distances, UMI) is the single source
 * of barcodes for the whole pipeline: it is written onto the STAR BAM as CB/UB
 * tags for fuscia, and joined by read ID for the arriba and flexiplex calls.
 *
 * This replaces STARsolo, which does not assign barcodes and UMIs correctly to
 * chimeric alignment records - exactly the records fusion calling depends on.
 *
 * R1 is a short read holding only the barcode and UMI, so there is no primer or
 * poly-T flank for flexiplex to anchor on. A seven character anchor (%%%%%%%)
 * is prepended to each read (sequence and quality) so that the barcode can be
 * located positionally, which is the same trick the downstream fusion-read
 * demultiplexing used.
 */

// Discover the barcodes present in the data, then intersect with the whitelist.
// Running flexiplex without -k reports every observed barcode with its read
// count; the most abundant params.barcode_discovery_top_n of those are
// intersected with the 10x whitelist to give the final list. Matching against
// this short list rather than the full ~6.8M whitelist is what makes the second
// flexiplex pass tractable on short-read data.
process discoverBarcodes {
	label 'process_high'
	publishDir "${params.out_dir}/demultiplex", mode: 'copy'

	input:
	path(read1)
	path(include_list)
	val(flexiplex_demultiplex_options)

	output:
	path("barcode_list.txt"),               emit: barcode_list
	path("flexiplex_barcodes_counts.txt"),  emit: counts

	script:
	"""
	for f in ${read1} ; do
		if [[ "\$f" == *.gz ]] ; then gunzip -c "\$f" ; else cat "\$f" ; fi
	done | \\
	awk 'NR%4==2 || NR%4==0 { print "%%%%%%%" \$0 ; next } { print }' | \\
	flexiplex -x %%%%%%% \\
		${flexiplex_demultiplex_options} \\
		-p ${task.cpus} -n flexiplex > /dev/null

	head -n ${params.barcode_discovery_top_n} flexiplex_barcodes_counts.txt | \\
		cut -f1 > top_barcodes.txt

	awk 'NR==FNR { wl[\$1] ; next } (\$1 in wl)' ${include_list} top_barcodes.txt \\
		> barcode_list.txt

	n_barcodes=\$(wc -l < barcode_list.txt)
	echo "Kept \$n_barcodes of the top ${params.barcode_discovery_top_n} observed barcodes after intersecting with the whitelist" 1>&2
	if [[ "\$n_barcodes" -eq 0 ]] ; then
		echo "ERROR: no observed barcode matched the whitelist - check --protocol and the read structure" 1>&2
		exit 1
	fi
	"""
}

// Assign a barcode and UMI to every read, against the list from discoverBarcodes
// (or a user supplied --barcode_list).
process demultiplexReads {
	label 'process_high'
	publishDir "${params.out_dir}/demultiplex", mode: 'copy'

	input:
	path(read1)
	path(barcode_list)
	val(flexiplex_demultiplex_options)

	output:
	path("reads_barcodes.txt")

	script:
	"""
	for f in ${read1} ; do
		if [[ "\$f" == *.gz ]] ; then gunzip -c "\$f" ; else cat "\$f" ; fi
	done | \\
	awk 'NR%4==2 || NR%4==0 { print "%%%%%%%" \$0 ; next } { print }' | \\
	flexiplex -x %%%%%%% \\
		${flexiplex_demultiplex_options} \\
		-k ${barcode_list} \\
		-p ${task.cpus} -n flexiplex > /dev/null

	mv flexiplex_reads_barcodes.txt reads_barcodes.txt
	"""
}

// VisiumHD spot barcodes are split in two, either side of the UMI, so flexiplex
// is run twice: the first pass pulls out the second half of the barcode and
// writes it into the read ID, the second pass pulls out the UMI and the first
// half. No discovery step is needed - the slide whitelist is small enough to
// match against directly.
process demultiplexReadsVisiumHD {
	label 'process_high'
	publishDir "${params.out_dir}/demultiplex", mode: 'copy'

	input:
	path(read1)
	path(include_list)

	output:
	path("reads_barcodes.txt")

	script:
	"""
	for f in ${read1} ; do
		if [[ "\$f" == *.gz ]] ; then gunzip -c "\$f" ; else cat "\$f" ; fi
	done | \\
	awk 'NR%4==2 || NR%4==0 { print "%%%%%%%" \$0 ; next } { print }' | \\
	flexiplex \\
		-x "%%%%%%%?????????G?????????????" \\
		-b "???????????????" \\
		-k ${include_list} \\
		-e 1 -f 2 -r false -p ${task.cpus} -n stage1 | \\
	flexiplex \\
		-x "%%%%%%%" \\
		-u "?????????" \\
		-x "G" \\
		-b "?????????????" \\
		-k ${include_list} \\
		-e 1 -f 2 -p ${task.cpus} -n stage2 > /dev/null

	normalise_visium_barcodes.py \\
		--input stage2_reads_barcodes.txt \\
		--output reads_barcodes.txt
	"""
}

// Write the demultiplexed barcodes onto the aligned reads as CB/UB tags, so
// that fuscia sees the same barcodes as the arriba and flexiplex calls.
process transferBarcodesToBAM {
	label 'process_medium'
	publishDir "${params.out_dir}/STAR", mode: 'copy'

	input:
	path(bam_file)
	path(bam_index)
	path(barcode_table)
	path(fusion_targets)

	output:
	path("Aligned.sortedByCoord.tagged.bam"),      emit: bam
	path("Aligned.sortedByCoord.tagged.bam.bai"),  emit: bam_index

	script:
	def region_opts = params.tag_full_bam ? "--all" :
		"--targets ${fusion_targets} --pad ${params.tag_bam_pad}"
	"""
	tag_bam_barcodes.py \\
		--bam ${bam_file} \\
		--barcodes ${barcode_table} \\
		--output Aligned.sortedByCoord.tagged.bam \\
		--threads ${task.cpus} \\
		${region_opts}

	samtools index -@ ${task.cpus} Aligned.sortedByCoord.tagged.bam
	"""
}
