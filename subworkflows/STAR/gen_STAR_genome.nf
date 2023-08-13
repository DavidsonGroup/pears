process gen_genome{
    publishDir "${projectDir}/modules/STAR/"

    script:
    """
    $projectDir/modules/arriba/download_references.sh $params.genome_version

    $projectDir/modules/STAR/source/STAR --runThreadN $params.cpus \
    --runMode genomeGenerate \
    --genomeDir $projectDir/STAR/STAR_index_GRCh38_GENCODE38 \
    --genomeFastaFiles $projectDir/arriba/GRCh38.fa\
    --sjdbGTFfile $projectDir/arriba/GENCODE38.gtf\
    --sjdbOverhang $params.R2_length - 1
    """

}
