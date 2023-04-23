include { FUSCIA } from './submodules/fuscia/discover_chimeric_transcripts.py'
project_dir = projectDir

process getVariables {
	input:
	masterdata.out()
	
	output:
	

}

workflow RUN_FUSCIA {
	conda './env/pears_env.yml'
	
	input:
	masterdata.out()
	params.fuscia_mapqual
	params.out_dir

	output:
  	stdout output
	
	script:
	"""
	#!/usr/bin/env python3
	import pandas as pd
	import sys

	def run_fuscia(masterdata, map_qual, fuscia, cur_dir):
	    r = pd.read_csv(masterdata)

	    for index, row in r.iterrows():
		#if map_qual != NULL
		gene1 = f'{row["chrom1"]}:{min(row["gene1"], row["base1"])}-{max(row["gene1"], row["base1"])}'
		gene2 = f'{row["chrom2"]}:{min(row["gene2"], row["base2"])}-{max(row["gene2"], row["base2"])}'
		
	
	masterdata = sys.argv[1]
	map_qual = sys.argv[2]
	fuscia = sys.argv[3]
	cur_dir = sys.argv[4]
	run_fuscia(masterdata, map_qual, fuscia, cur_dir)
	print('fuscia done!')
	
	"""

}


os.system(f'python {FUSCIA} /stornext/Bioinf/data/lab_davidson/wu.s/cellranger/5cl/outs/*.bam {gene1} {gene2} {cur_dir}/fuscia_out {row["fusion genes"]}_{index} {map_qual}')
