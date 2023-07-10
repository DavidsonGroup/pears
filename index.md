---
layout: default
---

# What is pears?


# Installing pears
1. `git clone https://github.com/DavidsonGroup/pears.git`
2. make sure you have nextflow installed. See [docs](https://www.nextflow.io/docs/latest/getstarted.html).
3. Specify your parameters in the nextflow.config file 
4. run `nextflow run pears -c /path/to/nextflow.config` 


### Input
Pears requires:
 - single-cell paired-end reads
 - known list of fusions
 - reference genome

| input       | requirements      |
|:-------------|:------------------|
| Reads | named according to cellranger: [Sample Name]\_S1_L00[Lane Number]\_[Read Type]\_001.fastq.gz |
| Cell barcode and Unique Molecular Barcode (UMI) | You will need to specify the length of the cell barcode and UMI in your data  |


Your known list of fusions should contain:
 - name of gene fusion
 - for both genes: chromosome, base location of breakpoint, strand sense
You may also need to change the delimiter between genes (i.e. "--" for fusions formatted as such "BCAS4--BCAS3"). The default delimiter is ":".

pears runs STARsolo to align single-cell reads to a reference genome generating alignment and index files. If you have already have both, you can skip the alignment step in the config file. You should place your alignment files in a directory with the path `your_out_dir/STAR/`. 

Similarly, if you would like to use individual sections of pears you may select/deselect tools. See more details. **hyperlink to another section**

#### Adjustable Parameters
The adjustable parameters for pears includes:
 - fuscia_mapqual
 - fuscia_up/fuscia_down
 - flexi_searchlen 

| parameters     | requirements    |
|:-------------|:------------------|
| fuscia_mapqual | the `map_qual` parameter used by fuscia, this will determine the minimum acceptable MAPQ quality of reads fuscia will look for |
| fuscia_up/fuscia_down | if there is no corresponding gene annotation in the reference gene.gtf file, you can specify a range (in bp) upstream/downstream of the fusion breakpoint |
| flexiplex_searchlen | 2x the length of sequence (bp) flexiplex will use to search for the fusion sequence. Please note, the longer the search length the longer flexiplex will take to run. See [Davidson et al.](https://github.com/DavidsonGroup/flexiplex) |

### Arriba
[Link to Arriba](https://github.com/suhrig/arriba)
### Flexiplex
[Link to Flexiplex](https://github.com/DavidsonGroup/flexiplex)
### Fuscia
[Link to Fuscia](https://github.com/ding-lab/fuscia)
### Output

# Modules and Requirements
Pears includes three modules: alignment, fuscia and flexiplex. You may choose to use any combination of the three, however some may have dependencies.

| module       | requirements      |
|:-------------|:------------------|
| STAR   | reads, reference  |
| Arriba | reads, `.bam` and `.bam.bai` files | 
| fuscia       | reads, list of known fusions (fusion names, chrX:start-end for both genes), `.bam` and `.bam.bai` files |
| flexiplex    | reads, list of known fusions (fusion names, cell barcode and UMI length, sequence around the fusion breakpoint)  |

Pears also contains additional scripts to help format and streamline between processes.

`gen_masterdata` creates the main file fuscia and flexiplex will drawn from. It will output a dataframe with the columns:
 - fusion name
 - chromosome
 - start of the gene (gene1) to site of breakpoint (base1) or site of breakpoint (base2) to end of gene (gene2)
 - strand sense
 - two sequence +/- n bases around the breakpoint
The header of the file includes: `fusion_name | chrom1 | gene1 | base1 | sequence1 | strand1 | chrom2 | gene2 | base2 | sequence2 | strand2 `

Fuscia will create files that include: `cell_barcode | molecular_barcode | chrom | start | end `.

Flexiplex will create files that include: `Read | Cell_barcode | FlankEditDist | BarcodeEditDist | UMI `. 

`format_fuscia` and `format_flexiplex` function to collate the `.txt` files flexiplex and fuscia output.

`gen_report` creates the main report file **that is already classified **






