CUR_DIR=$1
REF=$2
READS=$3

cd $CUR_DIR

module load cellranger/3.0.2

cellranger count --id=cellranger_output --transcriptome=$REF \
--fastqs=$READS              
