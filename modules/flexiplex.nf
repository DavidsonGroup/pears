process runFlexiplex {
	input:
	 tuple val (fusion_genes), val (chrom1), val (gene1), val (base1), val (sequence1), val (chrom2), val (gene2), val (base2), val (sequence2)	 

	script:
	"""
	#!/usr/bin/env python3
	import os
	os.makedirs('$params.out_dir/flexiplex_out', exist_ok = True)
	fusion_name = f'${fusion_genes}_${task.index}'
	flexiplex = f'$projectDir/submodules/flexiplex/flexiplex'
	cmd = f'paste -d "&" $params.reads*.fastq | {flexiplex} -n {fusion_name} -l $sequence1 -k $sequence2 -r "" -f 1 -e 1 -u 0 -i false -s false > $params.out_dir/flexiplex_output/{fusion_name}_reads.fastq ; cd $params.out_dir/flexiplex_out/ ; {flexiplex} -l "" -r "&" -e 0 -f 0 -n {fusion_name} {fusion_name}_reads.fastq ; {flexiplex} -l "" -k {fusion_name}_barcodes_counts.txt -b $params.cellbarcode_len -u $params.umi_len -r "&" -e 0 -f 0 -n barcodes_{fusion_name} $params.out_dir/flexiplex_out/{fusion_name}_reads.fastq'
	os.system(cmd)
	
	"""
}
