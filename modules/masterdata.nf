include { gen_masterdata }    from '../submodules/gen_masterdata.py'

workflow GEN_MASTERDATA {
	input:
	params.shr_output
	params.reference
	params.flexi_searchlen

	main:
	masterdata.csv = gen_masterdata(shr_output, reference, flexi_searchlen)

	output:
	file 'masterdata.csv' into ch_output
}
