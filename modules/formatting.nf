process formatFuscia{

	script:
	"""
	#!/usr/bin/env python3

	import pandas as pd
	import os

	def add_fusion_name(file):
		r = pd.read_table(file)
		if r.empty == False:
			r['fusion'] =  os.path.basename(file).split("_")[0]
			return r

	in_dir = '$params.out_dir/fuscia_out'
	df = pd.DataFrame(columns=['cell_barcode', 'molecular_barcode', 'chrom', 'start', 'end','fusion'])
	for file in os.listdir(in_dir):
		r = add_fusion_name(f'{in_dir}{file}')
		df = pd.concat([df, r], axis = 0, ignore_index = True)
	df['cell_barcode'] = df['cell_barcode'].str.replace('-1', "")
	df.to_csv(f'$params.out_dir/master_fuscia.csv', index=False)

	"""

}

process formatFlexiplex{
	script:
	"""
	#!/usr/bin/env python3
	
	import pandas as pd
	import os
	
	def add_fusion_name(file):
		r = pd.read_table(file)
		if r.empty == False:
			df_temp = pd.DataFrame().assign(cell_barcode = r['Read'], molecular_barcode = r['BarcodeEditDist'])
			df_temp['fusion'] =  os.path.basename(file).split("_")[1]
			return df_temp

	in_dir = f'$params.out_dir/flexiplex_out'
	df = pd.DataFrame(columns = ['cell_barcode', 'molecular_barcode', 'fusion'])
	for file in os.listdir(in_dir):
		if os.path.basename(file)[0:8] == 'barcodes':
			r = add_fusion_name(f'{in_dir}/{file}')
	df = pd.concat([df, r], axis = 0, ignore_index = True)
	df.to_csv(f'$params.out_dir/master_flexiplex.csv', index=False)	

	"""

}


workflow FORMATTING{
	formatFuscia()
	formatFlexiplex()

}
