process runFlexiplex {
    
    input:
    params.reads
    params.out_dir
    params.flexi_searchlen
    params.cellbarcode_len
    params.umi_len

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import sys
    import os
    
    os.makedirs(f'$params.out_dir/flexiplex_output', exist_ok = True)
    df = pd.read_csv('$params.out_dir/masterdata.csv')
    
    for index, row in df.iterrows():
        fusion_name = f'{row["fusion genes"]}_{index}'
        left = row['sequence1']
        right = row['sequence2']
        flexiplex = '$projectDir/submodules/flexiplex/flexiplex'

        os.system(f'paste -d "&" $params.reads/*.fastq | {flexiplex} -n {fusion_name} -l {left} -k {right} -r "" -f 1 -e 1 -u 0 -i false -s false > $params.out_dir/flexiplex_output/{fusion_name}_reads.fastq')
        os.system(f'{flexiplex} -l "" -r "&" -b $params.cellbarcode_len -u {umi_len} -e 0 -f 0 -n {fusion_name} $params.out_dir/flexiplex_output/{fusion_name}_reads.fastq')
        os.system(f'{flexiplex} -l "" -k {fusion_name}_barcodes_counts.txt -b {cellbarcode_len} -u {umi_len} -r "&" -e 0 -f 0 -n barcodes_{fusion_name} $params.out_dir/flexiplex_output/{fusion_name}_reads.fastq')
    """

}
