here::i_am("scripts/transcriptomics_utilities.R")
library(DESeq2)
library(tidyverse)
library(here)

prefilter.dds <- function(dds, min_expression = 10, smallest_group = 3){
  keep <- rowSums(counts(dds) >= min_expression) >= smallest_group
  dds <- dds[keep,]
  return(dds)
}

filter.biotype <- function(dds, filter_cat = c("protein_coding")){
  keep <- rowData(dds)$gene_biotype %in% filter_cat
  dds <- dds[keep,]
  return(dds)
}

merge.names <- function(obj){
  res <- obj %>% 
    as.data.frame() %>% 
    { if(str_length(rownames(obj)[1]) > 1)
      
      rownames_to_column(., "gene_id") else
        .
    } %>% 
    left_join(rowData(dds) %>% as.data.frame() %>% dplyr::select(gene_id, symbol, entrezid)) %>%  
    dplyr::select(gene_id, symbol, entrezid, everything())
  return(res)
}


plot.volcano <- function(RES, FC_CO = log2(1.5), PVAL_CO = 0.05, TITLE = "", FC = log2FoldChange, PV = padj, INT = F){
  res_update <- RES %>% 
    merge.names() %>% 
    dplyr::mutate(fc =  {{FC}},
                  pv = {{PV}},
                  fc = round(fc, 2),
                  signif = case_when(
                    is.na(fc) | is.na(pv) ~ "no",
                    abs(fc) > FC_CO & pv < PVAL_CO ~ "both",
                    abs(fc) > FC_CO ~ "fold change",
                    pv < PVAL_CO ~ "p-value",
                    T ~ "no"),
                  signif = fct_relevel(signif, "no", "fold change", "p-value")) %>% 
    dplyr::filter(!is.na(pv))
  
  
  
  p <- res_update %>% 
    ggplot(aes(fc, -log10(pvalue), color = signif, alpha = signif, text = paste("log2 FC:", fc, "\nGene:", symbol)))+
    geom_point()+
    annotate("text", y = 0, 
             x = c(max(res_update %>% filter(!is.na(padj)) %>% pull(fc), na.rm = T), min(res_update %>% filter(!is.na(padj)) %>% pull(fc), na.rm = T)), 
             label = c(paste("higher in\n", GROUP[1]), paste("higher in\n", GROUP[2])),
             alpha = 0.8)+
    theme_minimal()+
    scale_alpha_manual(values = c(0.1, 0.4, 0.4, 1))+
    scale_color_manual(values = c("darkgrey", "orange", "darkblue", "red"))+
    labs(title = TITLE, x = "log2 FC", y = "-log10 p-value", color = "Significant by:", alpha = "Significant by:")
  
  ggsave(plot = p, 
         here("results", COMPARISON, "plots", paste0("volcano_plot_", COMPARISON, ".pdf")), 
         width = 14, height = 8)
  ggsave(plot = p, 
         here("results", COMPARISON, "plots", paste0("volcano_plot_", COMPARISON, ".png")), 
         width = 14, height = 8, bg = "white")
  
  if(INT){
    library(plotly)
    p <- ggplotly(p, tooltip = c("y", "size", "text"), width = 1600, height = 900) %>%
      layout(
        margin = list(l = 30, r = 8, b = 100, t = 50),
        annotations = list(
          x = 1, y = -0.2,
          text = paste(dim(res_update)[1], "genes present"),
          xref = 'paper', yref = 'paper', showarrow = F,
          xanchor = 'right', yanchor = 'auto', xshift = 0, yshift = 0,
          font = list(size = 10)
        )
      )
  }
  return(p)
}




hm_js <- function(df,          
                  clust_samples = T, # should the samples (columns)  be clustered
                  clust_params = T, # should the parameters (rows) be clustered
                  param_order = NULL, # if the parameters are not clustered, you can supply a custom order
                  linkage = "complete", #what linkage method for clustering
                  
                  show_dend_params = F, # show a dendrogram for parameter clustering
                  show_dend_samples = F,  # show a dendrogram for sample clustering
                  
                  show_sample_names = T, #show the sample names 
                  
                  normalise_params = T,  # normalisation across parameters (z-score)
                  normalise_samples = F, # normalisation across samples
                  norm_method = "zscore", # which normmilastion to use, zscore or max 
                  
                  excluded_vars = c(), # numeric variables included as columns which are not to be included in the heatmap(i.e. batch 1,2,3)
                  id_col = "", # column identifying each sample, defaults to first colum if empty
                  
                  color_code = c(high = "#FF1c00", low = "darkblue"), # color code for the heatmap
                  custom_threshold = NULL, #custom threshold for z-score 
                  outlier.removal = T, # should the zscore be cleared of outliers (outliers are set to the percentile detailed below)
                  outlier.threshold = 0.95, # want percentile of z-scores should be trimmed
                  
                  add_annotation = F, # should color bar annotation be included?
                  anno_col = "", # which columns contains annotations, defaults to everything but the identifier or the not excluded numeric columns
                  annotation_colors = list(), # list of color vectors, expects one for every anno col
                  
                  .plot = T, # should a plot be generated or only the data frame be returned 
                  return_list = F # if a plot is generated, should it be returned as a plot or as a list of subplots 
){
  require(tidyverse, quietly = T)
  require(patchwork, quietly = T)
  require(tidymodels, quietly = T)
  require(tidyclust, quietly = T)
  
  if(show_dend_params | show_dend_samples){
    library(factoextra)
  }
  df_input <- df
  
  if(normalise_params) { #normalise params
    df_input <- df_input %>% 
      {if(norm_method == "zscore") 
        mutate(., across(.cols = where(is.numeric) & !excluded_vars, scale)) 
        else .} %>% 
      {if(norm_method == "max") 
        mutate(., across(.cols = where(is.numeric) & !excluded_vars, 
                         function(x){(x-min(x))/(max(x)-min(x))})) 
        else .}}
  
  if(normalise_samples){ #normalize samples
    df_input <- df_input %>% 
      pivot_longer(cols = where(is.numeric) & !excluded_vars,
                   names_to = "params", 
                   values_to = "val") %>% 
      pivot_wider(names_from = ifelse(nchar(id_col) > 0, id_col, names(df)[1]),
                  values_from = "val") %>% 
      {if(norm_method == "zscore") 
        mutate(., across(.cols = where(is.numeric) & !excluded_vars, scale)) 
        else .} %>% 
      {if(norm_method == "max") 
        mutate(., across(.cols = where(is.numeric) & !excluded_vars, 
                         function(x){(x-min(x))/(max(x)-min(x))})) 
        else .} %>% 
      pivot_longer(cols = where(is.numeric) & !excluded_vars,
                   names_to = ifelse(nchar(id_col) > 0, id_col, names(df)[1]),
                   values_to = "val") %>% 
      pivot_wider(names_from = "params", values_from = "val")
  }
  
  message("normalisation done")
  if(clust_samples){
    samples_rec <- recipe(~., data = df_input) %>%
      {if(length(excluded_vars) >0) 
        step_rm(., excluded_vars) 
        else .} %>% 
      step_rm(-c(all_numeric_predictors())) 
    
    
    samples_wf <- workflow() %>%
      add_recipe(samples_rec) %>%
      add_model(hier_clust(linkage_method = linkage))
    
    samples_fit <- samples_wf %>%
      fit(data = df_input)  
    
    samples_fit <- samples_fit %>% 
      extract_fit_engine() 
    
    order_samples <-samples_fit$order
    
    
    if(show_dend_samples){
      dend_samples <- samples_fit %>%
        fviz_dend()
    }
    message("clustering samples done")}
  
  if(clust_params){
    df_input_mod <- df_input %>% 
      dplyr::select(ifelse(nchar(id_col) > 0, id_col, names(df_input)[1]),where(is.numeric), -excluded_vars) %>% 
      pivot_longer(cols = where(is.numeric) & !excluded_vars,
                   names_to = "params", 
                   values_to = "val") %>% 
      pivot_wider(names_from = ifelse(nchar(id_col) > 0, id_col, names(df)[1]),
                  values_from = "val")
    
    
    params_rec <- recipe(~., data = df_input_mod) %>%
      {if(length(excluded_vars) >0) step_rm(., excluded_vars) else .} %>% 
      step_rm(-c(all_numeric_predictors())) 
    
    params_wf <- workflow() %>%
      add_recipe(params_rec) %>%
      add_model(hier_clust(linkage_method = linkage))
    
    params_fit <- params_wf %>%
      fit(data = df_input_mod)  
    
    
    
    params_fit <- params_fit %>% 
      extract_fit_engine() 
    
    order_params <- params_fit$order
    
    
    if(show_dend_params){
      dend_params <- params_fit %>%
        fviz_dend() %>% attr("dendrogram")
    }
    message("clustering params done")
  }
  
  if(is.numeric(custom_threshold) & norm_method == "zscore"){
    df_plot <- df_input %>% 
      pivot_longer(cols = where(is.numeric) & !excluded_vars,
                   names_to = "params", 
                   values_to = "val")
    
    df_plot <- df_plot %>% 
      mutate(val = case_when(
        val > custom_threshold ~ custom_threshold,
        val < -custom_threshold ~ -custom_threshold,
        T ~ val)) 
    
    fill_labels =  round(seq(-custom_threshold, custom_threshold, length.out = 5),2) %>% as.numeric()
    fill_names = c(paste0("< -", custom_threshold), 
                   paste0("- ", custom_threshold/2 %>% round(2)), 
                   "0",
                   paste0(custom_threshold/2 %>% round(2)),
                   paste0("> ", custom_threshold))
  }else if(outlier.removal & norm_method == "zscore"){
    df_plot <- df_input %>% 
      pivot_longer(cols = where(is.numeric) & !excluded_vars,
                   names_to = "params", 
                   values_to = "val")
    
    thres <- df_plot$val %>% quantile(probs = outlier.threshold) %>% round(2)
    df_plot <- df_plot %>% 
      mutate(val = case_when(
        val > thres ~ thres,
        val < -thres ~ -thres,
        T ~ val)) 
    
    fill_labels =  round(seq(-thres, thres, length.out = 5),2) %>% as.numeric()
    fill_names = c(paste0("< -", thres), 
                   paste0("- ", thres/2 %>% round(2)), 
                   "0",
                   paste0(thres/2 %>% round(2)),
                   paste0("> ", thres))
  }else{
    df_plot <- df_input %>% 
      pivot_longer(cols = where(is.numeric) & !excluded_vars,
                   names_to = "params", 
                   values_to = "val") 
  }
  
  
  df_plot <- df_plot %>% 
    {if(nchar(id_col) > 0) dplyr::rename(.,SAMPLE = id_col)  else dplyr::rename(.,SAMPLE = 1)} %>% 
    mutate(SAMPLE = as_factor(SAMPLE),
           params = as_factor(params)) %>% 
    {if (clust_params) mutate(., params = fct_relevel(params, levels(params)[order_params])) else . }%>%
    {if (clust_samples) mutate(., SAMPLE = fct_relevel(SAMPLE, levels(SAMPLE)[order_samples])) else . } %>% 
    {if (!clust_params & !is.null(param_order)) mutate(., params = fct_relevel(params, param_order)) else . }
  
  if(.plot){
    message("plotting now")
    
    
    hm <- df_plot %>% 
      ggplot(aes(x = SAMPLE, 
                 y = params,
                 fill = val))+
      geom_tile()+
      {if(outlier.removal | (is.numeric(custom_threshold) & norm_method == "zscore"))
        scale_fill_gradient2(high = color_code[1], low = color_code[2],
                             breaks = fill_labels,
                             labels = fill_names) 
        else
          scale_fill_gradient2(high = color_code[1], low = color_code[2])}+
      scale_x_discrete(position = "top")+
      # theme_void()+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "bottom")+
      {if(length(unique(df_plot$params)) < 100) 
        theme(axis.text.y = element_text())}+
      {if(!show_sample_names)
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) 
        else if(length(unique(df_plot$SAMPLE)) > 50) 
          theme(axis.text.x = element_text(hjust = 1, angle = 90))}+
      {if(norm_method == "zscore") labs(fill = "z-score norm.\nmRNA expression")}+
      {if(norm_method == "max") labs(fill = "Maximum scaled\nmRNA expression")}
    
    if(add_annotation){
      message("adding annotation")
      
      hm <- hm +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
      plot_list <- list(hm)
      
      
      if(length(anno_col) > 1){
        if(!all(anno_col %in% colnames(df_plot))){
          stop("not all supplied columns are in the data frame")
        }
        anno_cols <- anno_col
      }else{
        anno_cols <- df_plot %>% 
          dplyr::select(-c(SAMPLE, params, val)) %>% names()
      }
      
      for(col in 1:length(anno_cols)){
        anno <- df_plot %>% 
          ggplot(aes(y = 1, x = SAMPLE, fill = .data[[anno_cols[col]]]))+
          geom_tile()+
          scale_y_discrete(breaks = seq(from = 0, to = 1, by = 0.25), labels = rep("", 5))+
          {if(length(annotation_colors) >= col && is.character(annotation_colors[[col]])) 
            scale_fill_manual(values = annotation_colors[[col]])}+
          theme_void()+
          theme(legend.position = "bottom")
        plot_list <- prepend(plot_list, list(anno))
      }
      
      if(return_list) {
        return(plot_list)
      }
      
      message("wrapping plots up")
      plot <- wrap_plots(plot_list, 
                         heights = c(if(clust_samples & show_dend_samples){0.3}, rep(0.05, length(anno_cols)), 1), 
                         ncol = 1, guides = "collect") & theme(legend.position = "bottom")
      return(plot)
    }else{
      return(hm) 
    }
    
  }else{
    return(df_plot)
  }}


