//setting up scrips and software
include { GEN_MASTERDATA } from './modules/gen_masterdata.nf'
include { RUN_CELLRANGER } from './modules/run_cellranger.nf'
include { FusciaFlexiplex } from './submodules/FuFl.nf'
//include { runFuscia } from './modules/fuscia.nf'
//include { runFlexiplex } from './modules/flexiplex.nf'
include { FORMATTING } from './modules/formatting.nf'

//create channels
//ch_shr_output = params.shr_output ? file(params.shr_output) : file("${params.in_dir}/shr_output.csv")
//ch_reference = params.reference ? params.reference : "${params.in_dir}/reference"

workflow {
	
	GEN_MASTERDATA()
	RUN_CELLRANGER()
	
	process runModules {
		when:
		file("$params.out_dir/masterdata.csv").exists()

		script:
		"""
		FusciaFlexiplex()
		FORMATTING()
		"""

	}
}
