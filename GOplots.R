go_enrich <- function(gene_ids, types, org_db){
  eid <- bitr(gene_ids, fromType = "SYMBOL", toType = types, OrgDb = org_db)
  ego_BP <- enrichGO(gene = eid$ENTREZID, #as.character(unique(inp_sub$EnsemblID)),
                     OrgDb         = org_db,
                     #keytype       = 'ENTREZID',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     #pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05)
  return(ego_BP)
}

kegg_enrich <- function(gene_ids, types, org_db, orgname){
  eid <- bitr(gene_ids, fromType = "SYMBOL", toType = types, OrgDb = org_db)
  ego_kegg <- enrichKEGG(gene = eid$ENTREZID, #as.character(unique(inp_sub$EnsemblID)),
                         org=orgname,
                         #keytype       = 'ENTREZID',
                         #ont           = "BP",
                         pAdjustMethod = "BH",
                         #pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05)
  return(ego_kegg)
}


for(clus in unique(ami_clusters$cluster)){
  cat("Working for AMI GOplot ClusNo:", as.character(clus),"/", as.character(length(unique(ami_clusters$cluster))),"\n")
  if(as.integer(clus)%in%ami_human_clusterids){
    org_db = "org.Hs.eg.db"
    types <- c('ENTREZID','ENSEMBL')
  } else{
    org_db = "org.Ss.eg.db"
    types <- c('ENTREZID')
  }
  plt_name = paste0("AMI_GOplot_clus_", as.character(clus),".png")
  file_name = file.path("data", "GOplots", plt_name)
  # png(file_name)
  gene_list <- unique(ami_clusters$geneid[which(ami_clusters$cluster==as.numeric(clus))])
  ego <- go_enrich(gene_list, types, org_db)
  dplot <- dotplot(ego, showCategory=15)
  ggsave(plot = dplot, filename = file_name)
  # dev.off()
  # Sys.sleep(5)
}
  

for(clus in unique(cmi_clusters$cluster)){
  cat("Working for CMI GOplot ClusNo:", as.character(clus),"/", as.character(length(unique(cmi_clusters$cluster))),"\n")
  if(as.integer(clus)%in%cmi_human_clusterids){
    org_db = "org.Hs.eg.db"
    types <- c('ENTREZID','ENSEMBL')
  } else{
    org_db = "org.Ss.eg.db"
    types <- c('ENTREZID')
  }
  gene_list <- unique(cmi_clusters$geneid[which(cmi_clusters$cluster==as.numeric(clus))])
  ego <- go_enrich(gene_list, types, org_db)
  plt_name = paste0("CMI_GOplot_clus_", as.character(clus),".png")
  file_name = file.path("data", "GOplots", plt_name)
  dplot <- dotplot(ego, showCategory=15)
  ggsave(plot = dplot, filename = file_name)
}

#### KEGG plot #########

for(clus in unique(ami_clusters$cluster)){
  cat("Working for AMI KEGGplot ClusNo:", as.character(clus),"/", as.character(length(unique(ami_clusters$cluster))),"\n")
  if(as.integer(clus)%in%ami_human_clusterids){
    org_db = "org.Hs.eg.db"
    types <- c('ENTREZID','ENSEMBL')
    org_name = 'hsa'
  } else{
    org_db = "org.Ss.eg.db"
    types <- c('ENTREZID')
    org_name <- 'ssc'
  }
  gene_list <- unique(ami_clusters$geneid[which(ami_clusters$cluster==as.numeric(clus))])
  egokegg <- kegg_enrich(gene_list, types, org_db, org_name)
  plt_name = paste0("AMI_KEGGplot_clus_", as.character(clus),".png")
  file_name = file.path("data", "KEGGplots", plt_name)
  dplot <- dotplot(egokegg, showCategory=15)
  ggsave(plot = dplot, filename = file_name)
}

for(clus in unique(cmi_clusters$cluster)){
  cat("Working for CMI KEGGplot ClusNo:", as.character(clus),"/", as.character(length(unique(cmi_clusters$cluster))),"\n")
  if(as.integer(clus)%in%cmi_human_clusterids){
    org_db = "org.Hs.eg.db"
    types <- c('ENTREZID','ENSEMBL')
    org_name = 'hsa'
  } else{
    org_db = "org.Ss.eg.db"
    types <- c('ENTREZID')
    org_name <- 'ssc'
  }
  gene_list <- unique(cmi_clusters$geneid[which(cmi_clusters$cluster==as.numeric(clus))])
  egokegg <- kegg_enrich(gene_list, types, org_db, org_name)
  plt_name = paste0("CMI_KEGGplot_clus_", as.character(clus),".png")
  file_name = file.path("data", "KEGGplots", plt_name)
  # png(file_name)
  dplot <- dotplot(egokegg, showCategory=15)
  ggsave(plot = dplot, filename = file_name)
}
