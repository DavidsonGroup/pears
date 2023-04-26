include { runFuscia } from '../modules/fuscia.nf'
include { runFlexiplex } from '../modules/flexiplex.nf'


workflow FusciaFlexiplex {

	process masterdata_ch{
		when:
		file("$params.out_dir/masterdata.csv").exists()

		script:
		"""
		masterdata_ch = Channel.fromPath(GEN_MASTERDATA.out)
        	mapped_ch = masterdata_ch \
                	| splitCsv(header:true) \
	                | map { row -> tuple(row.fusion_genes, row.'chrom1', row.gene1, row.base1, row.sequence1, row.chrom2, row.gene2, row.base2, row.sequence2)} 
                
        	mapped_ch | runFuscia
		mapped_ch | runFlexiplex
		"""
	}
}
