#needs input (path/to/directory)
#no slash at the end of any
OUTDIR=$1
READS=$2
PEARS=$3

#variables
REF=$OUTDIR/refdata-gex-GRCh38-2020-A
SCRIPTS=$PEARS/scripts
FUSCIA=$PEARS/software/fuscia
FLEXI=$PEARS/software/flexiplex/flexiplex

#reformat data script.sh short-read.csv
python $SCRIPTS/gen_masterdata.sh $OUTDIR/jaffa_results.direct.csv $REF/genes/genes.gtf $REF/
echo 'reformat done'

#run cellranger -count
bash $SCRIPTS/run_cellranger.sh $OUTDIR $REF $READS
echo 'cellranger done'

#run fuscia -detected_chimeric
python $SCRIPTS/run_fuscia.sh $OUTDIR/masterdata.csv 30 $FUSCIA $OUTDIR
echo 'ran fuscia'

#collate fuscia results
python $SCRIPTS/format_fuscia.sh $OUTDIR/fuscia_output/ $OUTDIR
echo 'fuscia reformatted'

#run flexiplex
python $SCRIPTS/run_flexiplex.sh $READS $FLEXI $OUTDIR/masterdata.csv $OUTDIR
echo 'ran flexiplex'

#format flexiplex results
python $SCRIPTS/format_flexiplex.sh $OUTDIR/flexiplex_output $OUTDIR
echo 'flexiplex reformatted'
