#libraries
library(dplyr)

#take in master_fuscia

#filtering steps 

#   1. threshold
      #one to multiple reads
fusrange = fuscia %>% group_by(fusion, chrom) %>% summarize(across(start:end, ~ mean(.)))
      #within 10,000 bp +- of average
threshold <- 10000
      # filter rows in df1 that have a value within threshold range of df2
merged_fusrange = left_join(fuscia, fusrange, by = "fusion")

merged_fusrange$between_start <- ifelse(abs(merged_fusrange$start.x - merged_fusrange$start.y) <= threshold, TRUE, FALSE) 
merged_fusrange$between_end <- ifelse(abs(merged_fusrange$end.x - merged_fusrange$end.y) <= threshold, TRUE, FALSE) 
merged_filtered = merged_fusrange %>% dplyr::filter(between_start != FALSE & between_end != FALSE)
mfuscia = merged_filtered %>% dplyr::select(cell_barcode, molecular_barcode, chrom.x, start.x, end.x, fusion)

#   2. exclude fusions
mfuscia = subset(mfuscia, !grepl("MT-", fusion))

#   3. match UMIs
cfuscia = mfuscia %>% 
  group_by(cell_barcode, molecular_barcode, fusion) %>%
  summarise_all(list(~toString(unique(na.omit(.))))) %>%
  dplyr::select(cell_barcode, molecular_barcode, fusion)
cfuscia$algo = "Fuscia"

# take in master_flexiplex
cflexiplex = flexiplex %>%
  group_by(cell_barcode, molecular_barcode, fusion)  %>%
  summarise_all(list(~toString(unique(na.omit(.)))))
cflexiplex$algo = "Flexiplex"

fufl = merge(cflexiplex, cfuscia, all.x = TRUE, all.y = TRUE)
fufl

