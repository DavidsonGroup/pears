//setting up scrips and software
include { GEN_MASTERDATA } from './modules/gen_masterdata.nf'
include { RUN_CELLRANGER } from './modules/run_cellranger.nf'
include { runFuscia } from './modules/fuscia.nf'
include { runFlexiplex } from './modules/flexiplex.nf'
include { FORMATTING } from './modules/formatting.nf'

//create channels
//ch_shr_output = params.shr_output ? file(params.shr_output) : file("${params.in_dir}/shr_output.csv")
//ch_reference = params.reference ? params.reference : "${params.in_dir}/reference"

workflow {

	GEN_MASTERDATA()
	RUN_CELLRANGER()

	masterdata_ch = Channel.fromPath("$params.out_dir/masterdata.csv")
	mapped_ch = masterdata_ch \
	        | splitCsv(header:true) \
	        | map { row -> tuple(row.fusion_genes, row.'chrom1', row.gene1, row.base1, row.sequence1, row.chrom2, row.gene2, row.base2, row.sequence2)} 
		
	mapped_ch | runFuscia
	mapped_ch | runFlexiplex

	FORMATTING()

}
