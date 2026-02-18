library(DESeq2)
library(tidyverse)
library(here)

get.gos <- function(NAME, obj = dds, species = "hs"){
  if(species == "hs"){
    library(org.Hs.eg.db)
    DB <- org.Hs.eg.db
  }else {
    library(org.Mm.eg.db)
    DB <- org.Mm.eg.db
  }
  
  library(org.Hs.eg.db)
  res <- metadata(obj)$de_results[[NAME]]
  
  up_genes <- res %>% as.data.frame() %>% 
    rownames_to_column("gene") %>% 
    dplyr::filter(padj < 0.05, log2FoldChange > 1) %>% pull(gene)
  dn_genes <- res %>% as.data.frame() %>% 
    rownames_to_column("gene") %>% 
    dplyr::filter(padj < 0.05, log2FoldChange < -1) %>% pull(gene)
  
  universe <- res %>% as.data.frame() %>% 
    rownames_to_column("gene") %>% 
    dplyr::filter(!is.na(padj)) %>% pull(gene)
  
  up_go <- clusterProfiler::enrichGO(up_genes, DB, keyType = "SYMBOL", ont = "BP", universe = universe) %>% 
    as.data.frame()
  dn_go <- clusterProfiler::enrichGO(dn_genes, DB, keyType = "SYMBOL", ont = "BP", universe = universe) %>% 
    as.data.frame()
  
  if(!is.list(metadata(obj)$fe_results)){
    metadata(obj)$fe_results <- list()
  }
  if(!is.list(metadata(obj)$fe_results[[NAME]])){
    metadata(obj)$fe_results[[NAME]] <- list()
  }
  
  metadata(obj)$fe_results[[NAME]][["up_go"]] <- up_go
  metadata(obj)$fe_results[[NAME]][["dn_go"]] <- dn_go
  return(obj)
}
## EXAMPLE
# dds_2 <- get.gos("KO effect in untreated")



get.gsea <- function(NAME, obj = dds, type = "HALLMARK", conditions, species = "hs"){
  if(species == "hs"){
    SPECIES <- "HS"
  }else {
    SPECIES <- "MM"
  }
  
  if(type == "HALLMARK"){
    gene_sets <- msigdbr::msigdbr(collection = "H", db_species = SPECIES)
  }else if(type == "REACTOME"){
    gene_sets <- msigdbr::msigdbr(collection = "C2", subcollection = "CP:REACTOME", db_species = SPECIES)
  }
  
  gene_sets <- gene_sets %>% 
    dplyr::select(gs_name, gene_symbol) %>% 
    distinct()
  
  rankings <- counts(dds, normalized = T) %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") %>% 
    pivot_longer(-gene) %>% 
    left_join(as.data.frame(colData(dds)), join_by(name == sample)) %>% 
    dplyr::filter(condition %in% conditions) %>% 
    dplyr::select(gene, condition, value) %>% 
    group_by(condition, gene) %>% 
    summarise(sd = sd(value, na.rm = T), mean_expr = mean(value, na.rm = T)) %>% 
    mutate(sd = ifelse(sd >= 0.2 * abs(mean_expr), sd,  0.2 * abs(mean_expr)),
           mean_expr = ifelse(mean_expr == 0, 1, mean_expr))%>% 
    mutate(condition = ifelse(condition == conditions[2], "group2", "group1")) %>% 
    pivot_wider(names_from = condition, values_from = mean_expr:sd) %>% 
    dplyr::filter(if_all(sd_group1:sd_group2,\(x) x > 0),
                  if_any(mean_expr_group1:mean_expr_group2,\(x) x > 0)) %>%
    mutate(ranking = (mean_expr_group1 - mean_expr_group2) / (sd_group1 + sd_group2),
           gene = fct_reorder(gene, ranking, .desc = T)) %>% 
    filter(!is.na(ranking))%>% 
    arrange(desc(ranking)) %>% 
    pull(ranking, name = gene)
  
  gsea <- clusterProfiler::GSEA(rankings, TERM2GENE = gene_sets)
  
  if(!is.list(metadata(obj)$fe_results)){
    metadata(obj)$fe_results <- list()
  }
  if(!is.list(metadata(obj)$fe_results[[NAME]])){
    metadata(obj)$fe_results[[NAME]] <- list()
  }
  metadata(obj)$fe_results[[NAME]][[paste0("gsea_", type)]] <- as.data.frame(gsea)
  
  
  return(obj)
}

## Example
# 
# dds_3 <- get.gsea(NAME = "KO effect in untreated", obj = dds_2, type = "HALLMARK", conditions = c("ko_untreated", "wt_untreated"))
# dds_3 <- get.gsea(NAME = "KO effect in untreated", obj = dds_3, type = "REACTOME", conditions = c("ko_untreated", "wt_untreated"))
