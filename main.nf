//setting up scrips and software
include { GEN_MASTERDATA } from './subworkflows/gen_masterdata.nf'
include { RUN_STARsolo } from './subworkflows/STAR/run_STARsolo.nf'
include { runFuscia } from './subworkflows/fuscia.nf'
include { runFlexiplex } from './subworkflows/flexiplex.nf'
include { formatFuscia } from './subworkflows/formatting.nf'
include { formatFlexiplex } from './subworkflows/formatting.nf'

//create channels
//ch_shr_output = params.shr_output ? file(params.shr_output) : file("${params.in_dir}/shr_output.csv")
//ch_reference = params.reference ? params.reference : "${params.in_dir}/reference"

workflow {
	if(params.masterdata){GEN_MASTERDATA()}
//	if(params.align){RUN_STARsolo}
	
	if(params.masterdata){masterdata_ch = GEN_MASTERDATA.out}
	else{ masterdata_ch = Channel.fromPath(params.out_dir + '/masterdata.csv')}
        mapped_ch = masterdata_ch \
                | splitCsv(header:true) \
                | map { row -> tuple(row.fusion_genes, row.'chrom1', row.gene1, row.base1, row.sequence1, row.chrom2, row.gene2, row.base2, row.sequence2)}

	if(params.align){
		STARsolo_result = RUN_STARsolo()
		Fuscia_output_ch = runFuscia(mapped_ch, STARsolo_result.out).collect()
	}
	else{ 
		Fuscia_output_ch = runFuscia(mapped_ch, 'aligner done').collect() 
	} 
        Flexiplex_output_ch = runFlexiplex(mapped_ch).collect()
	
	formatFuscia(Fuscia_output_ch)
	formatFlexiplex(Flexiplex_output_ch)

}
