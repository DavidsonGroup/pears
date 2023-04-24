masterdata_ch = Channel.fromPath("${params.out_dir}/masterdata.csv").splitCsv(header: true)
process runFuscia {
    input:
    params.fuscia_mapqual
    params.out_dir

    masterdata_ch.view { row, index -> 
	gene1 = "${row["chrom1"]}:${Math.min(row["base1"].toInteger(), row["gene1"].toInteger())}-${Math.max(row["base1"].toInteger(), row["gene1"].toInteger())}"
    	gene2 = "${row["chrom2"]}:${Math.min(row["base2"].toInteger(), row["gene2"].toInteger())}-${Math.max(row["base2"].toInteger(), row["gene2"].toInteger())}"
	fusion_name = "${row["fusion genes"]}_$index"

	}

    script:
    """
  
    python $projectDir/submodules/fuscia/discover_chimeric_transcripts.py $params.out_dir/cellranger_output/outs/*.bams $gene1 $gene2 $params.out_dir/fuscia_out/ $fusion_name {$params.fuscia_mapqual}'
    """

}
