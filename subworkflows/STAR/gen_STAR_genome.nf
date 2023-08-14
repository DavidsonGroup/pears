process gen_genome{
    publishDir "${projectDir}/modules/STAR/"

    script:
    """
<<<<<<< HEAD
    $projectDir/modules/arriba/download_references.sh $params.genome_version

    $projectDir/modules/STAR/source/STAR --runThreadN $params.cpus \
=======
    $projectDir/arriba/download_references.sh $params.genome_version

    STAR --runThreadN $params.cpus \
>>>>>>> a4f53fdf61082f99d2cded4e28076871d04232ef
    --runMode genomeGenerate \
    --genomeDir $projectDir/STAR/STAR_index_GRCh38_GENCODE38 \
    --genomeFastaFiles $projectDir/arriba/GRCh38.fa\
    --sjdbGTFfile $projectDir/arriba/GENCODE38.gtf\
    --sjdbOverhang $params.R2_length - 1
    """

<<<<<<< HEAD
=======

>>>>>>> a4f53fdf61082f99d2cded4e28076871d04232ef
}
