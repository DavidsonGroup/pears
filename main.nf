//setting up scrips and software
include { GEN_MASTERDATA } from './subworkflows/gen_masterdata.nf'
include { RUN_STARsolo } from './subworkflows/run_STARsolo.nf'
include { runFuscia } from './subworkflows/fuscia.nf'
include { runFlexiplex } from './subworkflows/flexiplex.nf'
include { runArriba } from './subworkflows/arriba.nf'
include { formatFuscia } from './subworkflows/formatting.nf'
//include { formatFlexiplex } from './subworkflows/formatting.nf'
include { getBarcodes_Arriba } from './subworkflows/arriba.nf'
include {formatFlexiplex as formatFlexiplex1} from './subworkflows/formatting.nf'
include {formatFlexiplex as formatFlexiplex2} from './subworkflows/formatting.nf'

//create channels
//ch_shr_output = params.shr_output ? file(params.shr_output) : file("${params.in_dir}/shr_output.csv")
//ch_reference = params.reference ? params.reference : "${params.in_dir}/reference"

workflow {
	GEN_MASTERDATA()
	masterdata_ch = GEN_MASTERDATA.out
	
        mapped_ch = masterdata_ch \
                | splitCsv(header:true) \
                | map { row -> tuple(row.fusion_genes, row.'chrom1', row.gene1, row.base1, row.sequence1, row.chrom2, row.gene2, row.base2, row.sequence2)}

	STARsolo_result = RUN_STARsolo()

	Fuscia_output_ch = runFuscia(mapped_ch, STARsolo_result[0], STARsolo_result[1]).collect() // Passing only BAM and BAM index
        Flexiplex_output_ch = runFlexiplex(mapped_ch,STARsolo_result[2]).collect()
	Arriba_output_ch = runArriba(STARsolo_result[0])
	ArribaBC_output_ch = getBarcodes_Arriba(mapped_ch,Arriba_output_ch,STARsolo_result[2]).collect()

	formatFuscia(Fuscia_output_ch,"master_fuscia.csv")
	formatFlexiplex1(Flexiplex_output_ch,"flexiplex_out","master_flexiplex.csv")
	formatFlexiplex2(ArribaBC_output_ch,"arriba_out","master_arriba.csv")

}


