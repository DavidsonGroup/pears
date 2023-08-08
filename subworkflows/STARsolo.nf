process{
    publishDir
    input:
    path genomeDir    

    output:
    

    script:
    """
    STAR --runThreadN $params.cpus --genomeDir ${genomeDir}  --genomeLoad NoSharedMemory --readFilesIn $params.reads/*_R1* $params.reads/*_R2* --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outBAMcompression 0 --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --soloType CB_UMI_Simple --soloCBwhitelist $projectDir/modules/STAR/3M-february-2018.txt --soloUMIlen $params.umi_len --outSAMattributes NH HI nM AS CB UB
    """




}
