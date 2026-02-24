library(dplyr)
library(tidyverse)

# Rscript --vanilla script.R masterdata.csv fuscia.csv flexiplex.csv arriba.csv arriba_ID.txt CB UMI
args = commandArgs(trailingOnly = TRUE)

masterdata = read.csv(args[1], header = TRUE)
fuscia = read.csv(args[2], header = TRUE)
flexiplex = read.csv(args[3], header = TRUE)
arriba = read.csv(args[4], header =FALSE)
afusion_readID = read.csv(args[5], header = TRUE)
CB = args[6]
UMI = args[7]

#Fuscia 
xfuscia = fuscia %>% 
  group_by(cell_barcode, molecular_barcode, fusion) %>%
  summarize(chrom = paste(chrom, collapse = ", "),
            start = paste(start, collapse = ", "),
            end = paste(end, collapse = ", ")) %>%
  separate(chrom, into = c("chrom1", "chrom2"), sep = ", ") %>%
  separate(start, into = c("start1", "start2"), sep = ", ") %>%
  separate(end, into = c("end1", "end2"), sep = ", " ) %>%
  drop_na()
xfuscia = xfuscia %>% mutate(across(start1:end2, as.integer))

#   2. exclude fusions
mfuscia = subset(xfuscia, !grepl("MT-", fusion))

cfuscia = mfuscia %>% 
  group_by(cell_barcode, molecular_barcode, fusion) %>%
  summarize()
cfuscia$algo = "Fuscia"

#Flexiplex
flexiplex$molecular_barcode = substr(flexiplex$molecular_barcode, 1, CB)

cflexiplex = flexiplex %>%
  group_by(cell_barcode, molecular_barcode, fusion)  %>%
  summarise_all(list(~toString(unique(na.omit(.)))))
cflexiplex$algo = "Flexiplex"

fufl = merge(cflexiplex, cfuscia, all.x = TRUE, all.y = TRUE)

#Arriba
arriba = arriba %>% mutate(ID = gsub(">", "", str_extract(arriba$V1, "^[^ ]+"))) %>% mutate(cell_barcode = substr(V2, 1, CB), UMI = substr(V2, CB+1, CB+1+UMI)) %>% dplyr::select(ID, cell_barcode, UMI)
names(arriba) = c("ID", "cell_barcode", "molecular_barcode")
names(afusion_readID) = c('ID', 'fusion')

arriba_sc = inner_join(afusion_readID, arriba, by="ID")
arriba_sc = arriba_sc %>%
  tidyr::separate_rows(fusion, sep = ",") %>%
  dplyr::select(fusion, cell_barcode, molecular_barcode)
arriba_sc$fusion <- gsub("\\(.*\\)", "", arriba_sc$fusion)

arriba_sc = arriba_sc %>% group_by(cell_barcode, molecular_barcode, fusion) %>%
  distinct() %>% 
  mutate(relevant = ifelse(fusion %in% masterdata$fusion, TRUE, FALSE)) %>% 
  dplyr::filter(relevant == TRUE) %>% 
  dplyr::select(cell_barcode, molecular_barcode, fusion)

arriba_sc$algo = "Arriba"
master = merge(fufl, arriba_sc, all.x = TRUE, all.y = TRUE)

#write.csv(master, file = "./collated_fusions.csv", row.names = FALSE)