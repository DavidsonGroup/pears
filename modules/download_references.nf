/*
 * Download reference genome and gene annotation based on genome version.
 * Uses storeDir to cache downloads - files are only downloaded once per genome_version.
 */
process downloadReferences {
    label 'process_low'
    storeDir "${params.out_dir}/references/${genome_version.replace('+', '_')}"

    input:
    val genome_version

    output:
    path("*.fa"), emit: fasta
    path("*.gtf"), emit: gtf

    script:
    """
    download_references.py --reference ${genome_version}
    """
}
