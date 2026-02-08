/*
 * Calculate the mode of R2 read lengths from the first N reads.
 * Emits a warning if >10% of reads differ from the mode.
 */
process calculateReadLength {
    label 'process_tiny'

    input:
    path fastq_files

    output:
    stdout

    script:
    """
    calculate_read_length.py ${fastq_files}
    """
}
