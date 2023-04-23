//setting up scrips and software
include { GEN_MASTERDATA } from './modules/gen_masterdata.nf'
include { RUN_CELLRANGER } from './modules/run_cellranger.nf'
include { runFuscia } from './modules/fuscia.nf'
include { runFlexiplex } from './modules/flexiplex.nf'


//create channels
//ch_shr_output = params.shr_output ? file(params.shr_output) : file("${params.in_dir}/shr_output.csv")
//ch_reference = params.reference ? params.reference : "${params.in_dir}/reference"

process initFlexiplex{
	script:
	"""
	cd $projectDir/submodules/flexiplex ; make
	"""

}


workflow {

	//GEN_MASTERDATA()
	//RUN_CELLRANGER()
	//runFuscia()
	//if running flexiplex for the first time use: initFlexiplex()
	runFlexiplex()
}
