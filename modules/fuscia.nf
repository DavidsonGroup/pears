
process runFuscia {

    input:
    params.fuscia_mapqual
    params.out_dir

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd

    r = pd.read_csv("$params.out_dir/masterdata.csv")

    for index, row in r.iterrows():
        gene1 = f'{row["chrom1"]}:{min(row["gene1"], row["base1"])}-{max(row["gene1"], row["base1"])}'
        gene2 = f'{row["chrom2"]}:{min(row["gene2"], row["base2"])}-{max(row["gene2"], row["base2"])}'
        print(gene1)
	command = f'python ../submodules/fuscia/discover_chimeric_transcripts.py /stornext/Bioinf/data/lab_davidson/wu.s/cellranger/5cl/outs/*.bam {gene1} {gene2} $params.out_dir/fuscia_out/ {row["fusion genes"]}_{index} {map_qual}'
        sh command
    """

}
