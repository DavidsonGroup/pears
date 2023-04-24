---
layout: default
---

# What is pears?
### Input
Pears requires:
 - single-cell paired-end reads
 - known list of fusions
 - reference genome

You can either direct pears to a folder (directory) containing the requisites or you can specify different locations.

Reads shoud be named according to cellranger naming convention - [Sample Name]\_S1_L00[Lane Number]\_[Read Type]\_001.fastq.gz. You will need to specify the length of the cell barcode and UMI in your data.

Your known list of fusions should contain:
 - name of gene fusion
 - for both genes: chromosome, base location of breakpoint, strand sense
 - classification **todo check if want to remove?**
You may also need to change the delimiter between genes (i.e. "--" for fusions formatted as such "BCAS4--BCAS3"). The default delimiter is ":".

pears runs `cellranger -count` to get the `.bam` and `.bam.bai` file. However, if you have already have both, you can deselect the cellranger step in **todo add where they can specify this option**.

Similarly, if you would like to use individual sections of pears you may select/deselect modules. See more details. **hyperlink to another section**

```groovy
/* This is the configuration file, here you can adjust parameters and executors. */

params { 
/* pears takes either a directory ("/path/to/in_dir") containing:
 *              - reads: zipped (.fastq.gz) and formatted according to cellranger naming convention
 *                      ([Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz) 
 *                       You will need to specify the length of the cell barcode and unique molecular identifier        
 *              - reference genome: unzipped **TODO: make pears accept zipped format too**
 *              - list of known fusion (shr_ouput): output from JAFFA can be directly piped pears or in the                     
 *                      format of the example.csv file
 * or you can direct pears to each individual input file:
 *              - reads = "/path/to/reads/"
 *              - reference = "/path/to/reference/" 
 *                      - if you have pre-generated .fasta or gene.gtf files, ref_fasta and ref_gene
 *                        can be changed.
 * you will also need to indicate an output directory:
 *              - out_dir = "/path/to/output_directory/"
*/      
        in_dir =
 /*             - reads = "/path/to/reads/"
 *              - reference = "/path/to/reference/"
 *                      - if you have pre-generated .fasta or gene.gtf files, ref_fasta and ref_gene
 *                        can be changed.
 * you will also need to indicate an output directory:
 *              - out_dir = "/path/to/output_directory/"
 */     
        in_dir =
        
        reads = "/stornext/Bioinf/data/lab_davidson/wu.s/datasets/5cl_2R/format/"
        shr_output = "/stornext/Bioinf/data/lab_davidson/wu.s/pears_test/20-4_test/CCLE_formatted.csv"
        
        reference = "/stornext/Bioinf/data/lab_davidson/wu.s/pears_test/refdata-gex-GRCh38-2020-A"
        ref_fasta = "${reference}/fasta/fasta"
        ref_gene = "${reference}/gene/gene.gtf"
        
        out_dir = "/stornext/Bioinf/data/lab_davidson/wu.s/nf_pears_test/20-4"
        
        cellbarcode_len = 16
        umi_len = 10
       
```

#### Adjustable Parameters
The adjustable parameters for pears includes:
 - fuscia_mapqual
 - fuscia_up/fuscia_down
 - flexi_searchlen 

```
//adjustable parameters
/* fuscia_mapqual: minimum mapping quality of reads fuscia will look through
 * fuscia_up/fuscia_down: if there is no gene annotation, you can specify a range (bp) upstream/downstream 
 *                        of the fusion breakpoint
 * flexiplex_searchlen: 2x the length of sequence (bp) flexiplex will use to look for the fusion sequence. 
 *                      Please note, the longer the sequence, the longer flexiplex will take to run.
 *                      See **Davidson et al. (2023) 
*/                      

        fuscia_mapqual = 30
        fuscia_up = 2000 
        fuscia_down = 1000
        flexi_searchlen = 20
        
}
```

### Fuscia
[Link to Fuscia](https://github.com/ding-lab/fuscia)
### Flexiplex
[Link to Flexiplex](https://github.com/DavidsonGroup/flexiplex)
### Output

# Installing pears


# Modules and Requirements
Pears includes three modules: cellranger, fuscia and flexiplex. You may choose to use any combination of the three, however some may have dependencies.

| module       | requirements      |
|:-------------|:------------------|
| cellranger   | reads, reference  |
| fuscia       | reads, list of known fusions (fusion names, chrX:start-end for both genes), `.bam` and `.bam.bai` files   |
| flexiplex    | reads, list of known fusions (fusion names, cell barcode and UMI length, sequence around the fusion breakpoint)  |

Pears also contains additional scripts to help format and streamline between processes.

`gen_masterdata` creates the main file fuscia and flexiplex will drawn from. It will output a dataframe with the columns:
 - fusion name
 - chromosome
 - start of the gene (gene1) to site of breakpoint (base1) or site of breakpoint (base2) to end of gene (gene2)
 - strand sense
 - two sequence +/- n bases around the breakpoint

The header of the file includes: `fusion_name | chrom1 | gene1 | base1 | sequence1 | strand1 | chrom2 | gene2 | base2 | sequence2 | strand2 `

`format_fuscia` and `format_flexiplex` function to collate the `.txt` files flexiplex and fuscia output.

Fuscia will create files that include: `cell_barcode | molecular_barcode | chrom | start | end `.

Flexiplex will create files that include: `Read | Cell_barcode | FlankEditDist | BarcodeEditDist | UMI `. 

