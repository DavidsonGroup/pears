process runFuscia {
	label 'process_tiny'
	publishDir "${params.out_dir}/fuscia_out", mode: 'copy'

	input:
	tuple val(fusion_genes), val(chrom1), val(gene1), val(base1), val(sequence1), val(chrom2), val(gene2), val(base2), val(sequence2)
	path(bam_file)
	path(bam_index)
	val(mapqual)

	output:
	path('*.discovered_discordant_reads.tsv')

	script:
	// Use Groovy's conditional logic to handle the min and max comparison
	def gr1 = "${chrom1}:${Math.min(gene1.toInteger(), base1.toInteger())}-${Math.max(gene1.toInteger(), base1.toInteger())}"
	def gr2 = "${chrom2}:${Math.min(gene2.toInteger(), base2.toInteger())}-${Math.max(gene2.toInteger(), base2.toInteger())}"
	def fusion_name = "${fusion_genes}_${chrom1}_${base1}_${chrom2}_${base2}"

	"""
	# Run the fuscia discovery command
	fuscia_discover_chimeric_transcripts.py ${bam_file} ${gr1} ${gr2} . ${fusion_name} ${mapqual}
	"""
}
