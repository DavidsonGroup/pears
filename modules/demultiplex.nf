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

// Build the barcode list from the reads that are about to be demultiplexed.
// The only barcodes assignable to a read set are the ones observable in it, so
// the list is the observed barcodes that are in the whitelist, plus - unless
// --barcode_list_edit_distance 0 - the whitelist barcodes one error away from
// an observed barcode that is not. Those neighbours are what let a read with a
// sequencing error in its barcode still reach the right cell; without them it
// is corrected to whichever other barcode in the list happens to be close.
// This needs no pass over the library and gives flexiplex a list small enough
// that its fallback comparison against every barcode costs nothing.
process buildBarcodeList {
	label 'process_low'
	publishDir "${params.out_dir}/demultiplex", mode: 'copy'

	input:
	path(reads)
	path(include_list)
	val(barcode_length)

	output:
	path("barcode_list.txt")

	script:
	"""
	demultiplex_reads.py \\
		--whitelist ${include_list} \\
		--barcode-length ${barcode_length} \\
		--edit-distance ${params.barcode_list_edit_distance} \\
		--barcode-list-out barcode_list.txt \\
		${reads}
	"""
}

// Assign the barcodes directly, without handing the list to flexiplex. The
// barcode and UMI are at a fixed offset in R1, and the candidate barcodes for
// each read are already known from building the list, so the assignment is a
// lookup. Where a read's barcode is one error from more than one candidate it
// is resolved by how abundant each candidate is and how likely the implied
// error is given the base quality, rather than being dropped as ambiguous.
// VisiumHD keeps flexiplex, whose two-stage search handles the split barcode.
process demultiplexReadsDirect {
	label 'process_low'
	publishDir "${params.out_dir}/demultiplex", mode: 'copy'

	input:
	path(reads)
	path(include_list)
	val(barcode_length)
	val(umi_length)

	output:
	path("reads_barcodes.txt")

	script:
	"""
	demultiplex_reads.py \\
		--whitelist ${include_list} \\
		--barcode-length ${barcode_length} \\
		--umi-length ${umi_length} \\
		--edit-distance ${params.barcode_list_edit_distance} \\
		--table-out reads_barcodes.txt \\
		--barcode-list-out barcode_list.txt \\
		${reads}
	"""
}

// Gather the read IDs that need a barcode: everything over a fusion target
// region in the BAM (what fuscia will look at) plus every read arriba and
// flexiplex called as fusion supporting.
process collectTargetReadIDs {
	label 'process_low'
	publishDir "${params.out_dir}/demultiplex", mode: 'copy'

	input:
	path(bam_file)
	path(bam_index)
	path(fusion_targets)
	path(fusion_read_ids)

	output:
	path("target_read_ids.txt")

	script:
	"""
	tag_bam_barcodes.py \\
		--bam ${bam_file} \\
		--targets ${fusion_targets} \\
		--pad ${params.tag_bam_pad} \\
		--names-out bam_region_read_ids.txt

	cat bam_region_read_ids.txt ${fusion_read_ids} | \\
		sed 's/^@//' | awk 'NF' | sort -u > target_read_ids.txt

	echo "\$(wc -l < target_read_ids.txt) reads need a barcode" 1>&2
	"""
}

// Cut R1 down to those reads before demultiplexing. Flexiplex spends its time
// on reads that do not match the barcode list exactly - each one costs an edit
// distance against every barcode in the list - so demultiplexing thousands of
// reads instead of the whole library is the difference between seconds and
// hours. Set --demultiplex_all_reads to skip this and demultiplex everything.
process extractTargetReads {
	label 'process_low'
	publishDir "${params.out_dir}/demultiplex", mode: 'copy'

	input:
	path(read1)
	path(target_read_ids)

	output:
	path("target_R1.fastq")

	script:
	"""
	subset_fastq.py \\
		--reads ${target_read_ids} \\
		--output target_R1.fastq \\
		${read1}
	"""
}

// Assign a barcode and UMI to every read given, against the list from
// discoverBarcodes (or a user supplied --barcode_list).
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
