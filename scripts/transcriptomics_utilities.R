here::i_am("scripts/transcriptomics_utilities.R")
library(DESeq2)
library(tidyverse)
library(here)

prefilter.dds <- function(dds, min_expression = 10, smallest_group = 3) {
  keep <- rowSums(counts(dds) >= min_expression) >= smallest_group
  dds <- dds[keep, ]
  return(dds)
}

filter.biotype <- function(dds, filter_cat = c("protein_coding")) {
  keep <- rowData(dds)$gene_biotype %in% filter_cat
  dds <- dds[keep, ]
  return(dds)
}

add.ids <- function(se, add = c("SYMBOL", "ENTREZID", "GENETYPE")){
  library(org.Hs.eg.db)
  
  for (ID_TYPE in add) {
    se <- addIds(se, column = ID_TYPE)
  }
  names(rowData(se)) <- str_to_lower(names(rowData(se)))
  return(se)
}


filter.and.activate <- function(dds){
  keep <- !is.na(rowData(dds)$symbol)
  dds <- dds[keep,]
  rownames(dds) <- rowData(dds)$symbol
  return(dds)
}





merge.names <- function(obj) {
  res <- obj %>%
    as.data.frame() %>%
    {
      if (str_length(rownames(obj)[1]) > 1) {
        rownames_to_column(., "gene_id")
      } else {
        .
      }
    } %>%
    left_join(rowData(dds) %>% as.data.frame() %>% dplyr::select(gene_id, symbol, entrezid)) %>%
    dplyr::select(gene_id, symbol, entrezid, everything())
  return(res)
}


plot.volcano <- function(RES, FC_CO = log2(1.5), PVAL_CO = 0.05, TITLE = "", FC = log2FoldChange, PV = padj, INT = F) {
  res_update <- RES %>%
    merge.names() %>%
    dplyr::mutate(
      fc = {{ FC }},
      pv = {{ PV }},
      fc = round(fc, 2),
      signif = case_when(
        is.na(fc) | is.na(pv) ~ "no",
        abs(fc) > FC_CO & pv < PVAL_CO ~ "both",
        abs(fc) > FC_CO ~ "fold change",
        pv < PVAL_CO ~ "p-value",
        T ~ "no"
      ),
      signif = fct_relevel(signif, "no", "fold change", "p-value")
    ) %>%
    dplyr::filter(!is.na(pv))



  p <- res_update %>%
    ggplot(aes(fc, -log10(pvalue), color = signif, alpha = signif, text = paste("log2 FC:", fc, "\nGene:", symbol))) +
    geom_point() +
    annotate("text",
      y = 0,
      x = c(max(res_update %>% filter(!is.na(padj)) %>% pull(fc), na.rm = T), min(res_update %>% filter(!is.na(padj)) %>% pull(fc), na.rm = T)),
      label = c(paste("higher in\n", GROUP[1]), paste("higher in\n", GROUP[2])),
      alpha = 0.8
    ) +
    theme_minimal() +
    scale_alpha_manual(values = c(0.1, 0.4, 0.4, 1)) +
    scale_color_manual(values = c("darkgrey", "orange", "darkblue", "red")) +
    labs(title = TITLE, x = "log2 FC", y = "-log10 p-value", color = "Significant by:", alpha = "Significant by:")

  ggsave(
    plot = p,
    here("results", COMPARISON, "plots", paste0("volcano_plot_", COMPARISON, ".pdf")),
    width = 14, height = 8
  )
  ggsave(
    plot = p,
    here("results", COMPARISON, "plots", paste0("volcano_plot_", COMPARISON, ".png")),
    width = 14, height = 8, bg = "white"
  )

  if (INT) {
    library(plotly)
    p <- ggplotly(p, tooltip = c("y", "size", "text"), width = 1600, height = 900) %>%
      layout(
        margin = list(l = 30, r = 8, b = 100, t = 50),
        annotations = list(
          x = 1, y = -0.2,
          text = paste(dim(res_update)[1], "genes present"),
          xref = "paper", yref = "paper", showarrow = F,
          xanchor = "right", yanchor = "auto", xshift = 0, yshift = 0,
          font = list(size = 10)
        )
      )
  }
  return(p)
}



plot.lots.of.categories <- function(){
  ggplot(dat, aes(x = cat, y = num))+
    geom_boxplot()+
    scale_x_discrete(guide = guide_axis(n.dodge = 2))
}