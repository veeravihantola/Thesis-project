####### Combining liver and muscle data #########

rm(list = ls())  # clear environment

library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db) 
library(pathview)
library(tidyverse)
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)
library(enrichplot)
library(fgsea)
library(scales)
library(readxl)
library(RColorBrewer)
library(stringr)


#### changing the default plot colours
# Extract a 5-color vector from the Road/Purple palette
rdpu_colors <- rev(brewer.pal(5, "RdPu"))

# Set the session option with the vector
options(enrichplot.colours = rdpu_colors)


## Read in differential expression results
liver_res12 <- read_xlsx("/Users/s3040370/Downloads/Liver_proteins_3m_for_R.xlsx")

## Create a tibble
#liver_res12_tb <- liver_res12 %>%
#  as_tibble()

## clusterProfiler does not work as easily using gene names, so we will turn gene names into Ensembl IDs using 
## clusterProfiler::bitr and merge the IDs back with the DE results
##keytypes(org.Mm.eg.db)

ids_liver <- bitr(liver_res12$Genes, 
            fromType = "SYMBOL", 
            toType = c("ENSEMBL", "ENTREZID"), 
            OrgDb = "org.Mm.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_liver <- which(duplicated(ids_liver$SYMBOL) == FALSE)


ids_liver <- ids_liver[non_duplicates_liver, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_liver <- merge(x=liver_res12, y=ids_liver, by.x="Genes", by.y="SYMBOL")  
## Why are there now more obs. than in the filtered one?

## Extract significant results
sig_liver <- subset(merged_gene_ids_liver, p < 0.05)


################# Then let's do all that to muscle as well #######################

## Read in differential expression results
muscle_res12 <- read_xlsx("/Users/s3040370/Downloads/Muscle_3m_protein_for_R.xlsx")

## Create a tibble
#muscle_res12_tb <- muscle_res12 %>%
#  as_tibble()

## clusterProfiler does not work as easily using gene names, so we will turn gene names into Ensembl IDs using 
## clusterProfiler::bitr and merge the IDs back with the DE results
##keytypes(org.Mm.eg.db)

ids_muscle <- bitr(muscle_res12$Genes, 
            fromType = "SYMBOL", 
            toType = c("ENSEMBL", "ENTREZID"), 
            OrgDb = "org.Mm.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_muscle <- which(duplicated(ids_muscle$SYMBOL) == FALSE)

ids_muscle <- ids_muscle[non_duplicates_muscle, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_muscle <- merge(x=muscle_res12, y=ids_muscle, by.x="Genes", by.y="SYMBOL")  
## Why are there now more obs. than in the filtered one?

## Extract significant results
sig_muscle <- subset(merged_gene_ids_muscle, p < 0.05)


########### And same for spinal cord ############

## Read in differential expression results
spinalcord_res12 <- read_xlsx("/Users/s3040370/Downloads/Spinalcord_12m_for_R.xlsx")


## Create a tibble
#spinalcord_res12_tb <- liver_res12 %>%
#  as_tibble()

## clusterProfiler does not work as easily using gene names, so we will turn gene names into Ensembl IDs using 
## clusterProfiler::bitr and merge the IDs back with the DE results
##keytypes(org.Mm.eg.db)

ids_sc <- bitr(spinalcord_res12$Genes, 
            fromType = "SYMBOL", 
            toType = c("ENSEMBL", "ENTREZID"), 
            OrgDb = "org.Mm.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_sc <- which(duplicated(ids_sc$SYMBOL) == FALSE)

ids_sc <- ids_sc[non_duplicates_sc, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_sc <- merge(x=spinalcord_res12, y=ids_sc, by.x="Genes", by.y="SYMBOL")  
## Why are there now more obs. than in the filtered one?

## Extract significant results
sig_spinalcord <- subset(merged_gene_ids_sc, merged_gene_ids_sc$`p` < 0.05)

########### For brain ###########

## Read in differential expression results
brain_res12 <- read_xlsx("/Users/s3040370/Downloads/Brain_proteins_3m_for_R.xlsx")

## clusterProfiler does not work as easily using gene names, so we will turn gene names into Ensembl IDs using 
## clusterProfiler::bitr and merge the IDs back with the DE results
##keytypes(org.Mm.eg.db)

ids_brain <- bitr(brain_res12$Genes, 
                  fromType = "SYMBOL", 
                  toType = c("ENSEMBL", "ENTREZID"), 
                  OrgDb = "org.Mm.eg.db")

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates_brain <- which(duplicated(ids_brain$SYMBOL) == FALSE)


ids_brain <- ids_brain[non_duplicates_brain, ] 

## Merge the Ensembl IDs with the results     
merged_gene_ids_brain <- merge(x=brain_res12, y=ids_brain, by.x="Genes", by.y="SYMBOL")  
## Why are there now more obs. than in the filtered one?

## Extract significant results
sig_brain <- subset(merged_gene_ids_brain, p < 0.05)

###################### Now combine the two ###########################

# obs this only works for two data sets!

shared_genes <- intersect(sig_liver$Genes, sig_muscle$Genes)
only_liver <- setdiff(sig_liver$Genes, sig_muscle$Genes)
only_muscle <- setdiff(sig_muscle$Genes, sig_liver$Genes)

########## trying it for all the 4 at once #######

# Put all your data frames into a list
listing <- list(sig_liver$Genes, sig_muscle$Genes, sig_brain$Genes, sig_spinalcord$Genes)

# Find the common rows across all data frames
shared_genes <- Reduce(intersect, listing)

print(shared_genes)


## Making them into lists for ggvenn visualisation

gene_list <- list(
  Liver = sig_liver$Genes,
  Muscle = sig_muscle$Genes)

ggvenn(gene_list, fill_color = heat_hcl(2))

############### Doing an enrichment analysis for the shared genes ##############

#### no point of doing if there are very few shared genes

shared_genes_ensembl <- bitr(shared_genes,
                             fromType = "SYMBOL", 
                             toType = "ENSEMBL", 
                             OrgDb = "org.Mm.eg.db")

shared_genes_ensembl <- as.character(shared_genes_ensembl$ENSEMBL)

enrich_shared <- enrichGO(gene = shared_genes_ensembl, 
                keyType = "ENSEMBL", 
                OrgDb = org.Mm.eg.db, 
                ont = "ALL", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

shared_enrich_result <- data.frame(enrich_shared@result)

############## And the non-shared ones to see if there are still similar pathways ##############

#### this is better when doing all four

sig_muscle_genes <- as.character(sig_muscle$ENSEMBL)

ego_muscle <- enrichGO(gene = sig_muscle_genes, 
                          keyType = "ENSEMBL", 
                          OrgDb = org.Mm.eg.db, 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05, 
                          readable = TRUE)

results_ego_muscle <- data.frame(ego_muscle@result)



sig_liver_genes <- as.character(sig_liver$ENSEMBL)

ego_liver <- enrichGO(gene = sig_liver_genes, 
                      keyType = "ENSEMBL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = TRUE)

results_ego_liver <- data.frame(ego_liver@result)



sig_brain_genes <- as.character(sig_brain$ENSEMBL)

ego_brain <- enrichGO(gene = sig_brain_genes, 
                      keyType = "ENSEMBL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = TRUE)

results_ego_brain <- data.frame(ego_brain@result)



sig_sc_genes <- as.character(sig_spinalcord$ENSEMBL)

ego_sc <- enrichGO(gene = sig_sc_genes, 
                      keyType = "ENSEMBL", 
                      OrgDb = org.Mm.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      qvalueCutoff = 0.05, 
                      readable = TRUE)

results_ego_sc <- data.frame(ego_sc@result)

################ Enrichmap analysis ###############

map_muscle <- pairwise_termsim(ego_muscle)

emapplot_muscle <- emapplot(map_muscle, showCategory=40,
                            node_label_size = 4)
emapplot_muscle + ggtitle("Muscle 12m") +
  theme(
    plot.title = element_text(size = 20)
  )


map_liver <- pairwise_termsim(ego_liver)

emapplot_liver <- emapplot(map_liver, showCategory=35,
                           node_label_size = 4)
emapplot_liver + ggtitle("Liver 12m") +
  theme(
    plot.title = element_text(size = 20)
  )


map_brain <- pairwise_termsim(ego_brain)

emapplot_brain <- emapplot(map_brain, showCategory=35,
                           node_label_size = 4)
emapplot_brain + ggtitle("Brain 12m") +
  theme(
    plot.title = element_text(size = 20)
  )


map_sc <- pairwise_termsim(ego_sc)

emapplot_sc <- emapplot(map_sc, 
                        showCategory=38,
                        node_label_size = 4)
emapplot_sc + ggtitle("Spinal cord 12m") +
  theme(
    plot.title = element_text(size = 20)
  )

####### Trying to reduce the number of similar pathway annotations ##############

ego_muscle_simplified <- clusterProfiler::simplify(
  ego_muscle,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
  )

muscle_simplified_df <- ego_muscle_simplified@result

ego_liver_simplified <- clusterProfiler::simplify(
  ego_liver,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

liver_simplified_df <- ego_liver_simplified@result

ego_brain_simplified <- clusterProfiler::simplify(
  ego_brain,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

brain_simplified_df <- ego_brain_simplified@result


ego_sc_simplified <- clusterProfiler::simplify(
  ego_sc,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)

sc_simplified_df <- ego_sc_simplified@result




map_muscle_simplified <- pairwise_termsim(ego_muscle_simplified)

emapplot_muscle <- emapplot(map_muscle_simplified, showCategory=50)
emapplot_muscle + ggtitle("Muscle 12m") +
  theme(
    plot.title = element_text(size = 20)
  )

############# Comparing shared pathways between organs ###############

shared_pathways <- intersect(ego_liver@result$Description, ego_muscle@result$Description)

## or for more than two

# Put all your data frames into a list
pathway_listing <- list(ego_liver_simplified@result$Description, ego_muscle_simplified@result$Description, ego_brain_simplified@result$Description, ego_sc_simplified@result$Description)

# Find the common rows across all data frames
shared_pathways <- Reduce(intersect, pathway_listing)



#### Combining them into one table and adding a column "organ" ######

shared_muscle <- subset(
  results_ego_muscle,
  Description %in% shared_pathways
  )

shared_muscle$organ <- "muscle"


shared_liver <- subset(
  results_ego_liver,
  Description %in% shared_pathways
  )

shared_liver$organ <- "liver"


shared_brain <- subset(
  results_ego_brain,
  Description %in% shared_pathways
)

shared_brain$organ <- "brain"


shared_sc <- subset(
  results_ego_sc,
  Description %in% shared_pathways
)

shared_sc$organ <- "spinalcord"


shared_pathways_df <- bind_rows(shared_liver, shared_muscle, shared_brain, shared_sc)


#### building a table with only top 20 enrichment ratio of BP terms from every organ

top_muscle <- subset(
  muscle_simplified_df, muscle_simplified_df$ONTOLOGY == "BP") %>%
  slice_max(order_by = FoldEnrichment, n = 20)

top_muscle$organ <- "muscle"


top_liver <- subset(
  liver_simplified_df, liver_simplified_df$ONTOLOGY == "BP") %>%
  slice_max(order_by = FoldEnrichment, n = 20)
  

top_liver$organ <- "liver"


top_brain <- subset(
  brain_simplified_df, brain_simplified_df$ONTOLOGY == "BP") %>%
  slice_max(order_by = FoldEnrichment, n = 20)

top_brain$organ <- "brain"


top_sc <- subset(
  sc_simplified_df, sc_simplified_df$ONTOLOGY == "BP") %>%
  slice_max(order_by = FoldEnrichment, n = 20)

top_sc$organ <- "spinalcord"


top_pathways <- bind_rows(top_liver, top_muscle, top_brain, top_sc)




################ Building a bubbleplot ##################


## with all of the terms if there were not that many

ggplot(shared_pathways_df, aes(x = organ,
                                                                    y = Description,
                                                                    showCategory = 20,
                                                                    size = FoldEnrichment,
                                                                    fill = p.adjust)) +
  geom_point(shape = 21, colour = "black") +
  scale_size(range = c(2,10),
             name = expression("Enrichment ratio")) +
  scale_fill_gradient(low = "hotpink",
                      high = "yellow",
                      name = "p.adjust",
                      trans = pseudo_log_trans(sigma = 0.001)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1),
    axis.title = element_blank(),
    panel.grid.minor.y = element_line(colour = "grey50"),
    panel.grid.minor.x = element_line(colour = "grey50"))

#### for changing the scaling, you can also call

library(scales)

# and add 

pseudo_log_trans(sigma = 0.001) #or

trans = "sqrt"

# and try to customize the p-value points with

breaks = c()

#### if you want to use ready palettes ###

ggplot(top_pathways, aes(x = organ,
                               y = Description,
                               showCategory = 20,
                               size = FoldEnrichment,
                               fill = p.adjust)) +
  geom_point(shape = 21, colour = "black") +
  scale_size(range = c(2,10),
             name = expression("Enrichment ratio")) +
scale_fill_distiller(
  palette = "RdPu",
  name = "p.adjust",
  trans = pseudo_log_trans(sigma = 0.001),
  direction = -1     # Change to -1 to reverse the gradient direction
) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_blank(),
    panel.grid.minor.y = element_line(colour = "grey50"),
    panel.grid.minor.x = element_line(colour = "grey50"))

## nice palettes: RdPu, YlGnBu, 

### Displaying pathways with highest enrichment scores ###
top_enrich <- subset(shared_pathways_df, FoldEnrichment>5) %>% dplyr::pull(Description)

ggplot(subset(shared_pathways_df, Description %in% top_enrich), aes(x = organ,
               y = Description,
               showCategory = 20,
               size = FoldEnrichment,
               fill = p.adjust)) +
  geom_point(shape = 21, colour = "black") +
  scale_size(range = c(2,10),
             name = expression("Enrichment ratio")) +
  scale_fill_gradient(low = "darkgreen",
                      high = "white",
                      name = "p.adjust") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1),
    axis.title = element_blank(),
    panel.grid.minor.y = element_line(colour = "grey50"),
    panel.grid.minor.x = element_line(colour = "grey50"))


### Displaying pathways with highest p-values ###

top_padjust <- subset(shared_pathways_df, p.adjust<0.001) %>% dplyr::pull(Description)

ggplot(subset(shared_pathways_df, Description %in% top_padjust), aes(x = organ,
                                                                    y = Description,
                                                                    size = FoldEnrichment,
                                                                    fill = p.adjust)) +
  geom_point(shape = 21, colour = "black") +
  scale_size(range = c(2,10),
             name = expression("Enrichment ratio")) +
  scale_fill_gradient(low = "brown",
                      high = "white",
                      name = "p.adjust") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1),
    axis.title = element_blank(),
    panel.grid.minor.y = element_line(colour = "grey50"),
    panel.grid.minor.x = element_line(colour = "grey50"))




############## Separating up- and downregulated proteins into their own categories for the ORA analysis ############

### creating a function for the visuals of dotplots

custom_formatter <- function(x) {
  # 1. Capitalise the very first letter of the GO term
  x_formatted <- str_to_sentence(x)
  
  # Standardized case-insensitive matching ensures nothing is missed
  x_formatted <- str_replace_all(x_formatted, "(?i)mrna", "mRNA")
  x_formatted <- str_replace_all(x_formatted, "(?i)trna", "tRNA")
  x_formatted <- str_replace_all(x_formatted, "(?i)gtpase", "GTPase")
  x_formatted <- str_replace_all(x_formatted, "(?i)rrna", "rRNA")
  x_formatted <- str_replace_all(x_formatted, "(?i)rna", "RNA")
  x_formatted <- str_replace_all(x_formatted, "(?i)atp", "ATP")
  x_formatted <- str_replace_all(x_formatted, "(?i)nadph", "NADPH")
  x_formatted <- str_replace_all(x_formatted, "(?i)ampa", "AMPA")
  x_formatted <- str_replace_all(x_formatted, "(?i)g0", "G0")
  x_formatted <- str_replace_all(x_formatted, "(?i)g1", "G1")
  x_formatted <- str_replace_all(x_formatted, "(?i)nf-kappab", "NF-kappaB")
  x_formatted <- str_replace_all(x_formatted, "(?i)ii", "II")
  
  # 2. Fix specific acronyms using word boundaries (\\b) so they match exactly
  x_formatted <- str_replace_all(x_formatted, "\\bdna\\b", "DNA")
  x_formatted <- str_replace_all(x_formatted, "\\brna\\b", "RNA")
  x_formatted <- str_replace_all(x_formatted, "\\batp\\b", "ATP")
  # Add any other acronyms here following the same pattern
  
  # 3. Maintain text wrapping at 30 characters
  str_wrap(x_formatted, width = 30)
}


#### analysis for the liver ###

sig_up_liver <- subset(sig_liver, log.2.FC > 0)

sig_up_liver <- as.character(sig_up_liver$ENSEMBL)

ego_up_liver <- enrichGO(gene = sig_up_liver, 
                   keyType = "ENSEMBL", 
                   OrgDb = org.Mm.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

up_simplified_liver <- clusterProfiler::simplify(ego_up_liver,
                                           cutoff = 0.8,          # Similarity threshold (0 to 1)
                                           by = "p.adjust",       # Metric used to select representative term
                                           select_fun = min)

sig_down_liver <- subset(sig_liver, sig_liver$`log.2.FC`< 0)

sig_down_liver <- as.character(sig_down_liver$ENSEMBL)

ego_down_liver <- enrichGO(gene = sig_down_liver, 
                     keyType = "ENSEMBL", 
                     OrgDb = org.Mm.eg.db, 
                     ont = "BP", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

down_simplified_liver <- clusterProfiler::simplify(ego_down_liver,
                                             cutoff = 0.55,          # Similarity threshold (0 to 1)
                                             by = "p.adjust",       # Metric used to select representative term
                                             select_fun = min)

dp_up_liver <- dotplot(ego_up_liver, showCategory = 15, title = "Liver", label_format = custom_formatter)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dp_down_liver <- dotplot(ego_down_liver, showCategory = 15, title = "Downregulated in Liver", label_format = custom_formatter)
dp_up_simp_liver <- dotplot(up_simplified_liver, showCategory = 15, title = "Liver", label_format = custom_formatter)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dp_down_simp_liver <- dotplot(down_simplified_liver, showCategory = 15, title = "Liver", label_format = custom_formatter)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

dp_up_liver
dp_down_liver
dp_up_simp_liver
dp_down_simp_liver


up_summary_liver <- data.frame(ego_up_liver)
up_simp_summary_liver <- data.frame(up_simplified_liver)
down_simp_summary_liver <- data.frame(down_simplified_liver)


#### same for muscle ###

sig_up_muscle <- subset(sig_muscle, sig_muscle$`log.2.FC` > 0)

sig_up_muscle <- as.character(sig_up_muscle$ENSEMBL)

ego_up_muscle <- enrichGO(gene = sig_up_muscle, 
                         keyType = "ENSEMBL", 
                         OrgDb = org.Mm.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)

up_simplified_muscle <- clusterProfiler::simplify(ego_up_muscle,
                                                 cutoff = 0.6,          # Similarity threshold (0 to 1)
                                                 by = "p.adjust",       # Metric used to select representative term
                                                 select_fun = min)

sig_down_muscle <- subset(sig_muscle, sig_muscle$`log.2.FC`< 0)

sig_down_muscle <- as.character(sig_down_muscle$ENSEMBL)

ego_down_muscle <- enrichGO(gene = sig_down_muscle, 
                           keyType = "ENSEMBL", 
                           OrgDb = org.Mm.eg.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = TRUE)

down_simplified_muscle <- clusterProfiler::simplify(ego_down_muscle,
                                             cutoff = 0.6,          # Similarity threshold (0 to 1)
                                             by = "p.adjust",       # Metric used to select representative term
                                             select_fun = min)

dp_up_muscle <- dotplot(ego_up_muscle, showCategory = 15, title = "Upregulated in Muscle", label_format = custom_formatter) +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    
dp_down_muscle <- dotplot(ego_down_muscle, showCategory = 15, title = "Downregulated in Muscle", label_format = custom_formatter)
dp_up_simp_muscle <- dotplot(up_simplified_muscle, showCategory = 15, title = "Muscle", label_format = custom_formatter)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dp_down_simp_muscle <- dotplot(down_simplified_muscle, showCategory = 15, title = "Muscle", label_format = custom_formatter)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

dp_up_muscle
dp_down_muscle
dp_up_simp_muscle
dp_down_simp_muscle

up_simp_summary_muscle <- data.frame(up_simplified_muscle)
down_simp_summary_muscle <- data.frame(down_simplified_muscle)
up_summary_muscle <- data.frame(ego_up_muscle)
down_summary_muscle <- data.frame(ego_down_muscle)



#### same for brain ###

sig_up_brain <- subset(sig_brain, sig_brain$`log.2.FC` > 0)

sig_up_brain <- as.character(sig_up_brain$ENSEMBL)

ego_up_brain <- enrichGO(gene = sig_up_brain, 
                          keyType = "ENSEMBL", 
                          OrgDb = org.Mm.eg.db, 
                          ont = "BP", 
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05, 
                          readable = TRUE)

up_simplified_brain <- clusterProfiler::simplify(ego_up_brain,
                                                  cutoff = 0.6,          # Similarity threshold (0 to 1)
                                                  by = "p.adjust",       # Metric used to select representative term
                                                  select_fun = min)

sig_down_brain <- subset(sig_brain, sig_brain$`log.2.FC`< 0)

sig_down_brain <- as.character(sig_down_brain$ENSEMBL)

ego_down_brain <- enrichGO(gene = sig_down_brain, 
                            keyType = "ENSEMBL", 
                            OrgDb = org.Mm.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = TRUE)

down_simplified_brain <- clusterProfiler::simplify(ego_down_brain,
                                                    cutoff = 0.5,          # Similarity threshold (0 to 1)
                                                    by = "p.adjust",       # Metric used to select representative term
                                                    select_fun = min)

dp_up_brain <- dotplot(ego_up_brain, showCategory = 15, title = "Upregulated in Brain", label_format = custom_formatter)
dp_down_brain <- dotplot(ego_down_brain, showCategory = 15, title = "Downregulated in Brain", label_format = custom_formatter)
dp_up_simp_brain <- dotplot(up_simplified_brain, showCategory = 15, title = "Brain", label_format = custom_formatter)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dp_down_simp_brain <- dotplot(down_simplified_brain, showCategory = 15, title = "Brain", label_format = custom_formatter)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

dp_up_brain
dp_down_brain
dp_up_simp_brain
dp_down_simp_brain


up_summary_brain <- data.frame(ego_up_brain)
down_summary_brain <- data.frame(ego_down_brain)
up_simp_summary_brain <- data.frame(up_simplified_brain)
down_simp_summary_brain <- data.frame(down_simplified_brain)


#### same for spinal cord ###

sig_up_sc <- subset(sig_spinalcord, sig_spinalcord$`log.2.FC` > 0)

sig_up_sc <- as.character(sig_up_sc$ENSEMBL)

ego_up_sc <- enrichGO(gene = sig_up_sc, 
                         keyType = "ENSEMBL", 
                         OrgDb = org.Mm.eg.db, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         qvalueCutoff = 0.05, 
                         readable = TRUE)

up_simplified_sc <- clusterProfiler::simplify(ego_up_sc,
                                                 cutoff = 0.52,          # Similarity threshold (0 to 1)
                                                 by = "p.adjust",       # Metric used to select representative term
                                                 select_fun = min)

sig_down_sc <- subset(sig_spinalcord, sig_spinalcord$`log.2.FC`< 0)

sig_down_sc <- as.character(sig_down_sc$ENSEMBL)

ego_down_sc <- enrichGO(gene = sig_down_sc, 
                           keyType = "ENSEMBL", 
                           OrgDb = org.Mm.eg.db, 
                           ont = "BP", 
                           pAdjustMethod = "BH", 
                           qvalueCutoff = 0.05, 
                           readable = TRUE)

down_simplified_sc <- clusterProfiler::simplify(ego_down_sc,
                                                   cutoff = 0.6,          # Similarity threshold (0 to 1)
                                                   by = "p.adjust",       # Metric used to select representative term
                                                   select_fun = min)


# creating dotplots
dp_up_sc <- dotplot(ego_up_sc, showCategory = 15, title = "Upregulated in Spinal cord", label_format = custom_formatter)
dp_down_sc <- dotplot(ego_down_sc, showCategory = 15, title = "Spinal cord", label_format = custom_formatter)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dp_up_simp_sc <- dotplot(up_simplified_sc, showCategory = 15, title = "Spinal cord", label_format = custom_formatter)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
dp_down_simp_sc <- dotplot(down_simplified_sc, showCategory = 15, title = "Spinal cord", label_format = custom_formatter)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

# having a look at the dotplots
dp_up_sc
dp_down_sc
dp_up_simp_sc
dp_down_simp_sc


# saving results into a data frame
up_summary_sc <- data.frame(ego_up_sc)
down_summary_sc <- data.frame(ego_down_sc)
up_simp_summary_sc <- data.frame(up_simplified_sc)
down_simp_summary_sc <- data.frame(down_simplified_sc)



################ Emapplots for up and downregulated separately ###############

map_up_muscle <- pairwise_termsim(ego_up_muscle)

emapplot_up_muscle <- emapplot(map_up_muscle, showCategory=40,
                            node_label_size = 4)
emapplot_up_muscle + ggtitle("Upregulated in Muscle 3m") +
  theme(
    plot.title = element_text(size = 20)
  )

map_down_muscle <- pairwise_termsim(ego_down_muscle)

emapplot_down_muscle <- emapplot(map_down_muscle, showCategory=35,
                               node_label_size = 4)
emapplot_down_muscle + ggtitle("Downregulated in Muscle 3m") +
  theme(
    plot.title = element_text(size = 20)
  )




map_up_liver <- pairwise_termsim(ego_up_liver)

emapplot_up_liver <- emapplot(map_up_liver, showCategory=40,
                           node_label_size = 4)
emapplot_up_liver + ggtitle("Upregulated in Liver 3m") +
  theme(
    plot.title = element_text(size = 20)
  )

map_down_liver <- pairwise_termsim(ego_down_liver)

emapplot_down_liver <- emapplot(map_down_liver, showCategory=40,
                              node_label_size = 4)
emapplot_down_liver + ggtitle("Downregulated in Liver 3m") +
  theme(
    plot.title = element_text(size = 20)
  )




map_up_brain <- pairwise_termsim(ego_up_brain)

emapplot_up_brain <- emapplot(map_up_brain, showCategory=40,
                           node_label_size = 4)
emapplot_up_brain + ggtitle("Upregulated in Brain 3m") +
  theme(
    plot.title = element_text(size = 20)
  )

map_down_brain <- pairwise_termsim(ego_down_brain)

emapplot_down_brain <- emapplot(map_down_brain, showCategory=35,
                              node_label_size = 4)
emapplot_down_brain + ggtitle("Downregulated in Brain 3m") +
  theme(
    plot.title = element_text(size = 20)
  )



map_up_sc <- pairwise_termsim(ego_up_sc)

emapplot_up_sc <- emapplot(map_up_sc, 
                        showCategory=70,
                        node_label_size = 4)
emapplot_up_sc + ggtitle("Upregulated in Spinal cord 12m") +
  theme(
    plot.title = element_text(size = 20)
  )

map_down_sc <- pairwise_termsim(ego_down_sc)

emapplot_down_sc <- emapplot(map_down_sc, 
                           showCategory=40,
                           node_label_size = 4)
emapplot_down_sc + ggtitle("Downregulated in Spinal cord 3m") +
  theme(
    plot.title = element_text(size = 20)
  )


### checking how many shared pathways there are in up- and downregulated groups

# Put all your data frames into a list
pathways_up <- list(up_simplified_liver@result$Description, up_simplified_muscle@result$Description, up_simplified_brain@result$Description, up_simplified_sc@result$Description)

# Find the common rows across all data frames
shared_ups <- Reduce(intersect, pathways_up)


# Put all your data frames into a list
pathways_down <- list(ego_down_liver@result$Description, ego_down_muscle@result$Description, ego_down_brain@result$Description, ego_down_sc@result$Description)

# Find the common rows across all data frames
shared_down <- Reduce(intersect, pathways_down)


write.csv(up_simp_summary_brain, "/Users/s3040370/Documents/Upregulated_brain.csv")
write.csv(up_simp_summary_muscle, "/Users/s3040370/Documents/Upregulated_muscle.csv")
write.csv(up_simp_summary_sc, "/Users/s3040370/Documents/Upregulated_spinalcord.csv")
write.csv(up_simp_summary_liver, "/Users/s3040370/Documents/Upregulated_liver.csv")

write.csv(down_simp_summary_brain, "/Users/s3040370/Documents/Downregulated_brain.csv")
write.csv(down_simp_summary_muscle, "/Users/s3040370/Documents/Downregulated_muscle.csv")
write.csv(down_simp_summary_sc, "/Users/s3040370/Documents/Downregulated_spinalcord.csv")
write.csv(down_simp_summary_liver, "/Users/s3040370/Documents/Downregulated_liver.csv")

#### visualising with a bubble plot

ggplot(top_30_up, aes(x = organ,
                         y = Description,
                         showCategory = 5,
                         size = FoldEnrichment,
                         fill = p.adjust)) +
  geom_point(shape = 21, colour = "black") +
  scale_size(range = c(2,10),
             name = expression("Enrichment ratio")) +
  scale_fill_distiller(
    palette = "RdPu",
    name = "p.adjust",
    trans = pseudo_log_trans(sigma = 0.001),
    direction = -1     # Change to -1 to reverse the gradient direction
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1,
                               size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_blank(),
    panel.grid.minor.y = element_line(colour = "grey50"),
    panel.grid.minor.x = element_line(colour = "grey50"))




top_30_up <- slice_min(shared_up_df, order_by = p.adjust, n = 30)



#### patchwork plots ####

## Combine the graphs using patchwork
library(patchwork)
combined_up <- dp_up_simp_brain + dp_up_simp_sc + dp_up_simp_muscle + dp_up_simp_liver

combined_down <- dp_down_simp_brain + dp_down_sc + dp_down_simp_muscle + dp_down_simp_liver


# Add the global title and style it
combined_up + 
  plot_annotation(
    title = "Upregulated pathways at 12m
    ",
    theme = theme(
      # Adjust size, bolding, and center the text (hjust = 0.5)
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )

combined_down + 
  plot_annotation(
    title = "Downregulated pathways at 12m
    ",
    theme = theme(
      # Adjust size, bolding, and center the text (hjust = 0.5)
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  )




dp_down <- dp_down + scale_y_discrete(position = "right")

(dp_up | dp_down) +
  plot_layout(guides = "collect",
              axis_titles = "collect_x") &
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.2, "cm"),
    plot.margin = margin(2, 2, 2, 2)
  )
