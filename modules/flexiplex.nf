process runFlexiplex {
    input:
    params.reads
    params.out_dir
    params.flexi_searchlen

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import os

    r = pd.read_csv("$params.out_dir/masterdata.csv")

    for index, row in r.iterrows():
        gene1 = f'{row["chrom1"]}:{min(row["gene1"], row["base1"])}-{max(row["gene1"], row["base1"])}'
        gene2 = f'{row["chrom2"]}:{min(row["gene2"], row["base2"])}-{max(row["gene2"], row["base2"])}'
        command = f'python $projectDir/submodules/fuscia/discover_chimeric_transcripts.py $params.out_dir/cellranger_output/outs/*.bam {gene1} {gene2} $params.out_dir/fuscia_out/ {row["fusion genes"]}_{index} {$params.fuscia_mapqual}'
        os.system(command)
    """

}
