process run_arriba {
	publishDir "${params.out}/arriba_out"
	time = '1d'
	memory = '150'

	input:
	params.reference
	params.reads
	params.out_dir
	params.STAR
	
	output:
	val "STAR done"	

	script:
	"""
	cd $projectDir/modules/arriba/

	arriba \
	    -x ../reads/Aligned.sortedByCoord.out.bam \
	    -o fusions.tsv -O fusions.discarded.tsv \
	    -a /stornext/Bioinf/data/lab_davidson/wu.s/software/arriba_v2.3.0/GRCh38.fa -g /stornext/Bioinf/data/lab_davidson/wu.s/software/arriba_v2.3.0/GENCODE38.gtf \
	    -b /stornext/Bioinf/data/lab_davidson/wu.s/software/arriba_v2.3.0/database/blacklist_hg38_GRCh38_v2.3.0.tsv.gz
	"""

}

process format_arriba {
	input:
	run_arriba.out

	script:
	"""
	python $projectDir/subworkflows/extract_arriba.py 
	seqtk subseq $params.read2 fusion_ID.txt > out.fq
	seqtk seq -a out.fq > output.fasta
	"""

}
