include { FUSCIA } from '../modules/fuscia/discover_chimeric_transcripts.py'

process runFuscia {

    input:
    file masterdata
    params.fuscia_mapqual
    params.out_dir

    output:
    file "$params.out_dir/fuscia_out/*"

    script:
    """
    import pandas as pd

    r = pd.read_csv("${masterdata}")

    for index, row in r.iterrows():
        gene1 = f'{row["chrom1"]}:{min(row["gene1"], row["base1"])}-{max(row["gene1"], row["base1"])}'
        gene2 = f'{row["chrom2"]}:{min(row["gene2"], row["base2"])}-{max(row["gene2"], row["base2"])}'
        FUSCIA($params.out_dir/cellranger_output/*.bam,${gene1},${gene2},$params.out_dir/fuscia_out,${row["fusion genes"]}_${index},$params.fuscia_mapqual)
    """
}

runFuscia(masterdata, map_qual, fuscia, out_dir)
