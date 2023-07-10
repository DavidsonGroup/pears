library(dplyr)
library(tidyverse)

masterdata = read.csv("/Users/wu.s/Documents/Research/4-5/masterdata.csv", header = TRUE)
fuscia = read.csv("/Users/wu.s/Documents/Research/4-5/master_fuscia.csv", header = TRUE)
flexiplex = read.csv("/Users/wu.s/Documents/Research/4-5/master_flexiplex.csv", header = TRUE)
flexiplex$molecular_barcode = substr(flexiplex$molecular_barcode, 1, 12)
JAFFA = read.csv('/Users/wu.s/Desktop/9cl/JAFFAoutput.csv', header =TRUE)
jfusion_readID = read.csv('/Users/wu.s/Desktop/9cl/JAFFAfusion_readID.txt', sep = "\t", header = FALSE)
ARRIBA = read.csv('/Users/wu.s/Desktop/9cl/Arribaoutput.csv', header =FALSE)
afusion_readID = read.csv('/Users/wu.s/Desktop/9cl/Arribafusions_ID.csv', header = TRUE)

#   1. match UMIs
xfuscia = fuscia %>% 
  group_by(cell_barcode, molecular_barcode, fusion) %>%
  summarize(chrom = paste(chrom, collapse = ", "),
            start = paste(start, collapse = ", "),
            end = paste(end, collapse = ", "))
#dplyr::select(cell_barcode, molecular_barcode, fusion, chrom, start, end)
xfuscia <- xfuscia %>%
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

cflexiplex = flexiplex %>%
  group_by(cell_barcode, molecular_barcode, fusion)  %>%
  summarise_all(list(~toString(unique(na.omit(.)))))
cflexiplex$algo = "Flexiplex"

fufl = merge(cflexiplex, cfuscia, all.x = TRUE, all.y = TRUE)

#JAFFA
JAFFA = JAFFA %>% mutate(ID = gsub(">", "", str_extract(ID, "^[^ ]+"))) %>% mutate(cell_barcode = substr(cell.and.UMI, 1, 16), UMI = substr(cell.and.UMI, 17, 28)) %>% dplyr::select(ID, cell_barcode, UMI)
names(JAFFA) = c("ID", "cell_barcode", "molecular_barcode")
names(jfusion_readID) = c('ID', 'fusion')
JAFFA_sc = inner_join(jfusion_readID, JAFFA, by="ID")

JAFFA_sc = JAFFA_sc %>% dplyr::select(fusion, cell_barcode, molecular_barcode)
JAFFA_sc$fusion = gsub(':', '--', JAFFA_sc$fusion)
JAFFA_sc = JAFFA_sc %>% group_by(cell_barcode, molecular_barcode, fusion) %>%
  summarise_all(list(~toString(unique(na.omit(.))))) %>% 
  mutate(relevant = ifelse(fusion %in% masterdata$fusion, TRUE, FALSE)) %>% 
  dplyr::filter(relevant == TRUE) %>% 
  dplyr::select(cell_barcode, molecular_barcode, fusion)
JAFFA_sc$algo = "JAFFA"
jmaster = merge(fufl, JAFFA_sc, all.x = TRUE, all.y = TRUE)
jmaster %>% dplyr::filter(algo == "JAFFA")

#Arriba
ARRIBA = ARRIBA %>% mutate(ID = gsub(">", "", str_extract(ARRIBA$V1, "^[^ ]+"))) %>% mutate(cell_barcode = substr(V2, 1, 16), UMI = substr(V2, 17, 28)) %>% dplyr::select(ID, cell_barcode, UMI)
names(ARRIBA) = c("ID", "cell_barcode", "molecular_barcode")
names(afusion_readID) = c('ID', 'fusion')

ARRIBA_sc = inner_join(afusion_readID, ARRIBA, by="ID")
ARRIBA_sc = ARRIBA_sc %>%
  tidyr::separate_rows(fusion, sep = ",") %>%
  dplyr::select(fusion, cell_barcode, molecular_barcode)
ARRIBA_sc$fusion <- gsub("\\(.*\\)", "", ARRIBA_sc$fusion)

ARRIBA_sc = ARRIBA_sc %>% group_by(cell_barcode, molecular_barcode, fusion) %>%
  distinct() %>% 
  mutate(relevant = ifelse(fusion %in% jmaster$fusion, TRUE, FALSE)) %>% 
  dplyr::filter(relevant == TRUE) %>% 
  dplyr::select(cell_barcode, molecular_barcode, fusion)

ARRIBA_sc$algo = "Arriba"
ARRIBA_sc %>% arrange(cell_barcode)

master = merge(jmaster, ARRIBA_sc, all.x = TRUE, all.y = TRUE)
master = master %>% mutate(fuscia = ifelse(algo == "Fuscia", TRUE, FALSE), flexiplex = ifelse(algo == "Flexiplex", TRUE, FALSE), JAFFA = ifelse(algo == "JAFFA", TRUE, FALSE), ARRIBA = ifelse(algo == "Arriba", TRUE, FALSE))

#write.csv(master, file = "/Users/wu.s/Desktop/9cl/mcf7freq_cluster.csv", row.names = FALSE)
