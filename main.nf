//setting up scrips and software
include { GEN_MASTERDATA } from './modules/gen_masterdata.nf'
//include { FUSCIA_WORKFLOW } from '/modules/fuscia_workflow'
//include { FLEXIPLEX_WORKFLOW } from '/modules/flexiplex_workflow'


//create channels
ch_shr_output = params.shr_output ? file(params.shr_output) : file("${params.in_dir}/shr_output.csv")
ch_reference = params.reference ? params.reference : "${params.in_dir}/reference"

workflow PEARS {
	conda './env/environment.yaml'

	GEN_MASTERDATA(params.shr_output, params.reference, params.flexi_searchlen)

}
