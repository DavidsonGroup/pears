# library
library(dplyr)
library(tidyr)

# parameters
obs_freq = 1
cluster_threshold = 0.8

# classification system
fusioncluster = dplyr::left_join(fufl, umap) %>% tidyr::drop_na()

freq_cluster = fusioncluster %>% group_by(fusion) %>% count(clusters)
freq_fusion = fusioncluster %>% group_by(fusion) %>% count(fusion)
clustering = dplyr::left_join(freq_cluster, freq_fusion, by = "fusion")
clustering$n.z = clustering$n.x/clustering$n.y

fusioncluster = left_join(fusioncluster, clustering, by = c("fusion", "clusters"))
fusioncluster = fusioncluster %>% group_by(cell_barcode, molecular_barcode, fusion) %>% mutate(both = all(c("Fuscia", "Flexiplex") %in% algo))

output <- fusioncluster %>% 
  mutate(
    confidence = ifelse( n.y <= obs_freq, "Low Confidence", ifelse(
      both & n.z > cluster_threshold, "High Confidence", ifelse( both | n.z > cluster_threshold, "Medium Confidence", "Low Confidence")
    ))
  ) %>% 
  mutate(confidence = factor(confidence, levels = c("High Confidence", "Medium Confidence", "Low Confidence"))) %>%
  arrange(confidence) %>% dplyr::select(!c(n.x, both))
names(output) = c("cell_barcode", "molecular_barcode", "fusion", "algorithm", "cluster", "frequency", "on_cluster", "confidence")
output

